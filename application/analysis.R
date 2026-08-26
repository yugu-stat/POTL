# Run from the application directory. The script is configured below for
# 100 replicates, 100 CPU workers, and evaluation at the target-study endpoint:
#   Rscript analysis.R
# TRANSCOX_CONDA_ENV can optionally identify a non-default Python environment.

# Avoid oversubscribing BLAS/OpenMP/TensorFlow inside replicate workers. Any
# values explicitly set by the caller are respected.
if(!nzchar(Sys.getenv("OMP_NUM_THREADS"))) Sys.setenv(OMP_NUM_THREADS = "1")
if(!nzchar(Sys.getenv("TF_NUM_INTRAOP_THREADS"))) {
  Sys.setenv(TF_NUM_INTRAOP_THREADS = "1")
}
if(!nzchar(Sys.getenv("TF_NUM_INTEROP_THREADS"))) {
  Sys.setenv(TF_NUM_INTEROP_THREADS = "1")
}

Rcpp::sourceCpp("EM.cpp")

library(reticulate)
transcox_conda_env = Sys.getenv(
  "TRANSCOX_CONDA_ENV",
  unset = "/home/ygu/anaconda3/envs/TransCoxEnvi"
)
if(dir.exists(transcox_conda_env)) {
  use_condaenv(transcox_conda_env)
}
source_python(system.file("python", "TransCoxFunction.py", package = "TransCox"))

library(intsurv)
library(survival)
library(dplyr)
library(parallel)
library(TransCox)
library(Hmisc)
library(rms)

nrep = 100L
ncore = 100L

# read target and source data
data = read.table("processed_target_data.txt", header = T)
datas = read.table("processed_source_data.txt", header = T)

# standardize continuous covariates
data[,c(3,4,10:19)] = scale(data[,c(3,4,10:19)], center = F)
datas[,c(3,4,10:19)] = scale(datas[,c(3,4,10:19)], center = F)

# sort data by time  
data = data %>%
  filter(time>0) %>%
  arrange(time) 
datas = datas %>%
  filter(time>0) %>%
  arrange(time) 

# Use the full target-study endpoint for covariate intervals, model fitting,
# and evaluation in every split and method.
target_study_tau = max(data$time)

# change the last few subjects' status to 0 to achieve better convergence 
ns = nrow(datas)
datas$status[c(ns-1,ns)] = 0
data$status[nrow(data)] = 0

# proportion of training data
ptrain = 0.7

# plot KM curve
alldata = data.frame("time" = c(data$time, datas$time),
                     "status" = c(data$status, datas$status),
                     "group" = c(rep(0, nrow(data)),
                                 rep(1, nrow(datas)))) %>%
  mutate(group = factor(group))
fit_km = survfit(Surv(time, status) ~ group, data = alldata)
pdf(file = "km.pdf", width = 8, height = 5.5)
plot(fit_km, lty = c("solid", "dashed"),
     xlab = "Time (years)", ylab = "Survival Probability",
     main = NULL,
     lwd = 2)
legend("topright", legend = c("TCGA-BRCA", "METABRIC"), lty = c("solid", "dashed"),
       lwd = 2)
dev.off()

make_censor_km = function(d) {
  fit = survfit(Surv(d$time, 1-d$status) ~ 1)
  time = as.numeric(fit$time)
  censor_survival = as.numeric(fit$surv)
  
  # Carry the last positive censoring estimate through the tail, matching the
  # C++ metric convention used for evaluation at the target-study endpoint.
  keep = censor_survival > 0
  list(time = time[keep], surv = censor_survival[keep])
}

km_value = function(km, q, left_limit = FALSE) {
  q = as.numeric(q)
  index = findInterval(q, km$time)
  if(left_limit && length(km$time)) {
    exact = index > 0L
    exact[exact] = abs(km$time[index[exact]]-q[exact]) <=
      64*.Machine$double.eps*pmax(1, abs(q[exact]))
    index[exact] = index[exact]-1L
  }
  c(1, km$surv)[index+1L]
}

make_cpp_censor_matrix = function(censor_km) {
  rbind(c(0, 1.0), cbind(censor_km$time, censor_km$surv))
}

survival_predictions = function(beta, baseline, newx, grid, r) {
  beta = as.numeric(beta)
  baseline = as.matrix(baseline)
  newx = as.matrix(newx)
  if(ncol(baseline) != 2L || nrow(baseline) == 0L ||
     any(!is.finite(baseline))) {
    stop("the fitted baseline hazard must be a finite two-column matrix")
  }
  if(length(beta) != ncol(newx) || any(!is.finite(beta)) ||
     any(!is.finite(newx))) {
    stop("beta and the validation covariates are incompatible or non-finite")
  }
  if(length(grid) < 2L || any(!is.finite(grid)) || grid[1] != 0 ||
     is.unsorted(grid, strictly = TRUE) || !is.finite(r) || r < 0) {
    stop("the prediction grid or transformation parameter is invalid")
  }
  if(any(baseline[,2] < -sqrt(.Machine$double.eps))) {
    stop("the fitted baseline contains a negative hazard increment")
  }
  baseline[,2] = pmax(baseline[,2], 0)
  baseline = baseline[order(baseline[,1]),,drop = FALSE]
  cumulative_hazard = c(0, cumsum(baseline[,2]))
  baseline_at_grid = cumulative_hazard[
    findInterval(grid, baseline[,1])+1L
  ]
  eta = drop(newx %*% beta)
  conditional_hazard = outer(exp(pmin(eta, 700)), baseline_at_grid)
  if(abs(r) < 1e-8) {
    pred = exp(-conditional_hazard)
  } else {
    pred = exp(-log1p(r*conditional_hazard)/r)
  }
  if(any(!is.finite(pred)) || any(pred < 0) || any(pred > 1)) {
    stop("the fitted model produced invalid survival probabilities")
  }
  pred
}

# Standard Graf-style IPCW Brier score.  An event observed by u is weighted by
# 1/G(Y-), while a subject still observed beyond u is weighted by 1/G(u).
# metric.hpp incorrectly uses 1/G(Y) for both terms, so its intBS is discarded.
corrected_ipcw_intBS = function(res, datav, censor_km, r, ngrid = 101L) {
  grid = seq(0, target_study_tau, length.out = ngrid)
  pred = survival_predictions(
    res$beta, res$lambda, as.matrix(datav[,-c(1,2)]), grid, r
  )
  gy_left = km_value(censor_km, datav$time, left_limit = TRUE)
  needed_event = datav$status == 1L
  if(any(!is.finite(gy_left[needed_event])) ||
     any(gy_left[needed_event] <= 0)) {
    stop("the censoring KM has no support for an event by the target-study endpoint")
  }
  brier = vapply(seq_along(grid), function(j) {
    u = grid[j]
    gu = km_value(censor_km, u)
    if(length(gu) != 1L || !is.finite(gu) || gu <= 0) {
      return(NA_real_)
    }
    contribution = numeric(nrow(datav))
    event_index = which(datav$time <= u & datav$status == 1L)
    alive_index = which(datav$time > u)
    contribution[event_index] = pred[event_index,j]^2/
      gy_left[event_index]
    contribution[alive_index] = (1-pred[alive_index,j])^2/gu
    mean(contribution)
  }, numeric(1L))
  if(any(!is.finite(brier))) {
    stop("the corrected IPCW Brier curve is non-finite")
  }
  sum(diff(grid)*(head(brier, -1L)+tail(brier, -1L))/2)/target_study_tau
}

# Uno concordance at the target-study endpoint. A comparable pair is anchored
# by an observed event i and a subject j observed strictly later. The pair
# receives the training-sample censoring weight 1/G(Y_i-)^2.
corrected_ipcw_cindex = function(time, status, risk_score, censor_km) {
  time = as.numeric(time)
  status = as.numeric(status)
  risk_score = as.numeric(risk_score)
  n = length(time)
  if(n < 2L || length(status) != n || length(risk_score) != n ||
     any(!is.finite(time)) || any(!is.finite(status)) ||
     any(!is.finite(risk_score)) || any(!status %in% c(0, 1))) {
    stop("the data for the IPCW C-index are invalid")
  }
  
  # Events at the largest observed time cannot form a comparable pair and do
  # not require censoring support.
  event_index = which(status == 1 & time < max(time))
  if(!length(event_index)) {
    stop("the validation data have no comparable event by the target-study endpoint")
  }
  event_g_left = km_value(
    censor_km, time[event_index], left_limit = TRUE
  )
  if(any(!is.finite(event_g_left)) || any(event_g_left <= 0)) {
    stop("the censoring KM has no support for a C-index event")
  }
  event_weight = 1/event_g_left^2
  if(any(!is.finite(event_weight))) {
    stop("the IPCW C-index weights are non-finite")
  }
  
  numerator = 0
  denominator = 0
  for(k in seq_along(event_index)) {
    i = event_index[k]
    later = which(time > time[i])
    if(!length(later)) next
    
    score_difference = risk_score[i]-risk_score[later]
    concordance_credit = as.numeric(score_difference > 0) +
      0.5*as.numeric(score_difference == 0)
    numerator = numerator + event_weight[k]*sum(concordance_credit)
    denominator = denominator + event_weight[k]*length(later)
  }
  if(!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    stop("the validation data have no finite comparable C-index pairs")
  }
  
  value = numerator/denominator
  if(!is.finite(value) || value < 0 || value > 1) {
    stop("the IPCW C-index is outside [0, 1]")
  }
  value
}

sum_metric = function(res, datav, censor_km, r) {
  res$intBS = corrected_ipcw_intBS(
    res, datav, censor_km, r
  )
  
  score = drop(as.matrix(datav[,-c(1,2)]) %*% as.numeric(res$beta))
  Cindex = corrected_ipcw_cindex(
    time = datav$time,
    status = datav$status,
    risk_score = score,
    censor_km = censor_km
  )
  
  res = c("Cindex" = Cindex,
          "intBS" = res$intBS,
          "RMST" = res$RMST)
  
  return(res)
}

trans_fit = function(datat, datav, rawXt, rawXv, S, weight, r_opt, xi,
                     kmcensor, censor_km, N) {
  # transfer learning algorithm
  weight = weight/sum(weight)
  res = TransFit(datat$time, datat$status, rawXt, S, weight, 
                 datav$time, datav$status, rawXv, kmcensor,
                 r_opt, N, xi, target_study_tau, verbose = 0)
  
  res = sum_metric(res, datav, censor_km, r_opt)
  
  return(res)
}

evaluate_external_fit = function(beta, lambda, time, datav, rawXv,
                                 kmcensor, censor_km, r) {
  beta = as.numeric(beta)
  lambda = as.numeric(lambda)
  time = as.numeric(time)
  if(length(beta) != ncol(datav)-2L || any(!is.finite(beta)) ||
     length(lambda) != length(time) || any(!is.finite(lambda)) ||
     any(!is.finite(time))) {
    stop("external baseline times and increments are invalid")
  }
  if(any(lambda < -sqrt(.Machine$double.eps))) {
    stop("external fit has a negative baseline-hazard increment")
  }
  lambda = pmax(lambda, 0)
  if(any(lambda > 0 & time <= 0)) {
    stop("external fit has a nonpositive baseline-hazard time")
  }
  keep = lambda > 0 & time <= target_study_tau
  lambda = lambda[keep]
  time = time[keep]
  if(!length(time)) {
    stop("external fit has no nonzero baseline increment by the target-study endpoint")
  }
  order_index = order(time)
  lambda = lambda[order_index]
  time = time[order_index]
  res = Metric(beta, lambda, time, datav$time, datav$status,
               rawXv, kmcensor, r, target_study_tau)
  sum_metric(res, datav, censor_km, r)
}

# choose the value of xi by 5-fold cross validation
cross_valid = function(datat, rawXt, S, weight, r_opt, xi, kmcensor, N, tau) {
  nt = nrow(datat)
  nfold = 5
  nrepeat = 5
  foldsize = round(nt/nfold)
  # assume that the data is already sorted by Y
  ave_loss = 0.0
  for(b in 1:nrepeat) {
    # create 5 folds
    folds = list()
    remain_sid = 1:nt
    for(fid in 1:nfold) {
      if(fid==nfold) {
        sid = remain_sid
      } else {
        sid = sample(remain_sid, size = foldsize)
      }
      folds[[fid]] = sort(sid)
      remain_sid = setdiff(remain_sid, sid)
    }
    b_loss = 0.0
    for(fid in 1:nfold) {
      sid = folds[[fid]]
      test = datat[sid,]
      test_rawX = rawXt[sid,]
      test_rawX[,1] = 0:(nrow(test_rawX)-1)
      train = datat[-sid,]
      train_rawX = rawXt[-sid,]
      train_rawX[,1] = 0:(nrow(train_rawX)-1)
      train_S = S[-sid] 
      train_weight = weight[-sid]
      train_weight = train_weight/sum(train_weight)
      
      # fit transformation model based on training data
      train_res = TransFit(train$time, train$status, train_rawX, train_S, train_weight, 
                           test$time, test$status, test_rawX, kmcensor, 
                           r_opt, N, xi, tau, verbose = 0, pll = 1, maxit = 100)
      
      # criterion: log profile likelihood
      b_loss = (b_loss*(fid-1)-train_res$logPL)/fid
      
      # # criterion: C-index
      # score = as.matrix(test[,c("age", "nlymph", "stage34", "IDC", "ER", "PR", "HER2",
      #                 paste0("PC", 1:10))]) %*% train_res$beta
      # if(sum(test[,"status"])){
      #   cindex = cIndex(time = test[,"time"],
      #                   event = test[,"status"],
      #                   risk_score = score)
      #   b_loss = (b_loss*(fid-1)-cindex["index"])/fid
      # }
      
      # # criterion: IBS
      # b_loss = (b_loss*(fid-1)+train_res$intBS)/fid
    }
    ave_loss = (ave_loss*(b-1)+b_loss)/b
  } 
  return(ave_loss)
}

eval_pred = function(repnum, data, datas, ptrain) {
  # sample size
  n = nrow(data)
  ns = nrow(datas)
  nt = round(n*ptrain)
  
  # source-study end time
  taus = max(datas$time)
  
  # Gau-Lag approx
  N = 20
  # candidate tuning parameter xi
  xi_cand = c(0.0, 2^seq(-15,15,1))
  # candidate value of r
  r_cand = rs_cand = seq(0, 1.5, 0.05)
  
  # divide data into training and testing sets
  set.seed(repnum)
  train_ind = sort(sample(1:n, round(n*ptrain)))
  # training data
  datat = data[train_ind,]
  taut = max(datat$time)
  rawXt = cbind(0:(nt-1),
                rep(0,nt),
                rep(target_study_tau,nt),
                datat[,-c(1,2)])
  rawXt = matrix(unlist(rawXt), nrow = nt)
  # validation data
  datav = data[-train_ind,]
  rawXv = cbind(0:(n-nt-1),
                rep(0,n-nt),
                rep(target_study_tau,n-nt),
                datav[,-c(1,2)])
  rawXv = matrix(unlist(rawXv), nrow = n-nt)
  # source data
  rawXs = cbind(0:(ns-1),
                rep(0,ns),
                rep(taus,ns),
                datas[,-c(1,2)])
  rawXs = matrix(unlist(rawXs), nrow = ns)
  
  # # select the optimal rs
  # aic_rs = sapply(rs_cand, function(rs_try) {
  #   cat("rs = ", rs_try, "\n")
  #   fit_source = sourceFit(datas$time, datas$status, rawXs,
  #                          datat$time, rawXt, rs_try, N, verbose = 0)
  #   -2*fit_source$logL
  # })
  # rs_opt = rs_cand[which.min(aic_rs)]
  rs_opt = 0.1
  
  # # select the optimal r
  # aic = sapply(r_cand, function(r_try) {
  #   fit_target = sourceFit(datat$time, datat$status, rawXt,
  #                          datat$time, rawXt, r_try, N, verbose = 0)
  #   -2*fit_target$logL
  # })
  # r_opt = r_cand[which.min(aic)]
  r_opt = 1.5
  
  # obtain source prediction 
  fit_source = sourceFit(datas$time, datas$status, rawXs, 
                         datat$time, rawXt, rs_opt, N,
                         verbose = 0, compVar = 1)
  S = fit_source$predSY
  
  # no weights
  weight = rep(1.0, length(S))
  
  # Estimate the censoring distribution from this target-training split only.
  # The list is used by corrected_ipcw_intBS; the matrix is retained for the
  # C++ RMST calculation and profile-likelihood CV interface.
  censor_km = make_censor_km(datat)
  kmcensor = make_cpp_censor_matrix(censor_km)
  
  ### Proposed method
  # tune parameter xi
  cv_loss = tryCatch({
    sapply(xi_cand, function(x) cross_valid(datat, rawXt, S, weight, r_opt, x, kmcensor, N, taut))
  }, error = function(e) {
    stop(sprintf("rep %d: CV for xi failed: %s",
                 repnum, conditionMessage(e)), call. = FALSE)
  })
  xi = xi_cand[which.min(cv_loss)]
  
  # fit model with optimal xi
  res = tryCatch({
    trans_fit(datat, datav, rawXt, rawXv, S, weight, r_opt, xi,
              kmcensor, censor_km, N)
  }, error = function(e) {
    stop(sprintf("rep %d: POTL fit/evaluation failed: %s",
                 repnum, conditionMessage(e)), call. = FALSE)
  })
  
  ### Target-only method
  res_tar = tryCatch({
    trans_fit(datat, datav, rawXt, rawXv, S, weight, r_opt, xi = 0.0,
              kmcensor, censor_km, N)
  }, error = function(e) {
    stop(sprintf("rep %d: target-only fit/evaluation failed: %s",
                 repnum, conditionMessage(e)), call. = FALSE)
  })
  
  ### TransCox
  tmpdata = data.frame(datat) %>%
    mutate(status = status+1) 
  tmpdatas = data.frame(datas) %>%
    mutate(status = status+1) 
  LRres <- SelLR_By_BIC(primData = tmpdata,
                        auxData = tmpdatas,
                        cov = c("age", "nlymph", "stage34", "IDC", "ER", "PR", "HER2",
                                paste0("PC", 1:10)),
                        statusvar = "status", lambda1 = 0.1, lambda2 = 0.1,
                        learning_rate_vec = 10^(seq(-5,0,1)),
                        nsteps_vec = as.integer(c(100, 200)))
  # select the best tuning parameter using BIC
  if(!is.na(LRres$best_lr) & !is.na(LRres$best_nsteps)) {
    best_lr = LRres$best_lr
    best_nsteps = as.integer(LRres$best_nsteps)
  } else {
    best_lr = 0.001
    best_nsteps = 100L
  }
  
  BICres <- SelParam_By_BIC(primData = tmpdata,
                            auxData = tmpdatas,
                            cov = c("age", "nlymph", "stage34", "IDC", "ER", "PR", "HER2",
                                    paste0("PC", 1:10)),
                            statusvar = "status",
                            lambda1_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)),
                            lambda2_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)),
                            learning_rate = best_lr, nsteps = best_nsteps)
  Cout <- GetAuxSurv(tmpdatas, cov = c("age", "nlymph", "stage34", "IDC", "ER", "PR", "HER2",
                                       paste0("PC", 1:10)))
  Pout <- GetPrimaryParam(tmpdata, q = Cout$q, estR = Cout$estR)
  Tres <- runTransCox_one(Pout, l1 = BICres$best_la1, l2 = BICres$best_la2,
                          learning_rate = best_lr, nsteps = best_nsteps,
                          cov = c("age", "nlymph", "stage34", "IDC", "ER", "PR", "HER2",
                                  paste0("PC", 1:10)))
  beta_est = Tres$new_beta
  lambda_est = Tres$new_IntH
  uniqt = Tres$time
  res_li23 = evaluate_external_fit(
    beta_est, lambda_est, uniqt, datav, rawXv, kmcensor, censor_km,
    r = 0.0
  )
  
  ## CoxTL (Lu et al.)
  source("coxTL.R")
  coxtl_cov = c("age", "nlymph", "stage34", "IDC", "ER", "PR", "HER2",
                paste0("PC", 1:10))
  p = length(coxtl_cov)
  data_t = datat[,c(coxtl_cov, "time", "status")]
  data_s = datas[,c(coxtl_cov, "time", "status")]
  data_t$R = factor(0, levels = c(0, 1))
  data_s$R = factor(1, levels = c(0, 1))
  coxtl_weights = CoxTL_density_weights(data_t, data_s, p)
  cox_tl<-run_CoxTL(data_t, data_s, weights = coxtl_weights, di = p)
  beta_est = cox_tl$coefficients
  cumlam = CoxTL_target_basehaz(cox_tl)
  uniqt = cumlam$time
  lambda_est = cumlam$hazard-c(0.0, cumlam$hazard[-nrow(cumlam)])
  res_coxtl = evaluate_external_fit(
    beta_est, lambda_est, uniqt, datav, rawXv, kmcensor, censor_km,
    r = 0.0
  )
  
  ### Pooled method
  datac = rbind(datat, datas)
  rawXc = rbind(rawXt, rawXs)
  # sort source data by Y
  ind = order(datac$time)
  datac = datac[ind,]
  rawXc = rawXc[ind,]
  rawXc[,1] = 0:(nt+ns-1)
  rawXc[,3] = taus
  # # select the optimal r for combined data
  # aic_com = sapply(r_cand, function(r_try) {
  #   fit_com = sourceFit(datac$time, datac$status, rawXc,
  #                       datat$time, rawXt, r_try, N, verbose = 0)
  #   -2*fit_com$logL
  # })
  # rc_opt = r_cand[which.min(aic_com)]
  rc_opt = 0.1
  res_com = trans_fit(datac, datav, rawXc, rawXv, rep(1,nt+ns), rep(1,nt+ns),
                      rc_opt, xi = 0.0, kmcensor, censor_km, N)
  
  allres = rbind(res, res_tar, res_li23, res_coxtl, res_com)
  allres = cbind(1:5, allres)
  attr(allres, "selected_xi") = xi
  
  return(allres)
}

run_replicate = function(i) {
  tryCatch({
    metrics = eval_pred(i, data, datas, ptrain)
    selected_xi = attr(metrics, "selected_xi", exact = TRUE)
    if(!is.matrix(metrics) || !identical(dim(metrics), c(5L, 4L)) ||
       any(!is.finite(metrics)) ||
       !identical(as.numeric(metrics[,1]), as.numeric(seq_len(5L))) ||
       length(selected_xi) != 1L || !is.finite(selected_xi)) {
      stop("replicate returned malformed/non-finite metrics or tuning result")
    }
    list(ok = TRUE, rep = i, metrics = metrics,
         selected_xi = selected_xi,
         error = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, rep = i, metrics = NULL,
         error = conditionMessage(e))
  })
}

replicate_runs = mclapply(seq_len(nrep), run_replicate, mc.cores = ncore)

valid_run = vapply(replicate_runs, function(x) {
  is.list(x) && isTRUE(x$ok) && is.matrix(x$metrics) &&
    identical(dim(x$metrics), c(5L, 4L)) && all(is.finite(x$metrics)) &&
    length(x$selected_xi) == 1L && is.finite(x$selected_xi)
}, logical(1))

if(any(!valid_run)) {
  failure_summary = data.frame(
    rep = which(!valid_run),
    error = vapply(replicate_runs[!valid_run], function(x) {
      if(is.list(x) && length(x$error) == 1L && !is.na(x$error)) {
        x$error
      } else if(inherits(x, "try-error") && length(x) >= 1L) {
        as.character(x)[1]
      } else {
        "malformed worker result"
      }
    }, character(1))
  )
  failure_file = "allres_failures.RData"
  save(replicate_runs, failure_summary, target_study_tau,
       nrep, ncore, file = failure_file)
  preview_n = min(5L, nrow(failure_summary))
  preview = paste(
    sprintf("rep %d: %s", failure_summary$rep[seq_len(preview_n)],
            failure_summary$error[seq_len(preview_n)]),
    collapse = " | "
  )
  stop(sprintf(
    "%d replicate(s) failed; diagnostics saved to %s. %s",
    nrow(failure_summary), failure_file, preview
  ), call. = FALSE)
}

# Aggregate only validated, successful replicate results.
allres = as.data.frame(do.call(rbind, lapply(replicate_runs, `[[`, "metrics")))
replicate_diagnostics = data.frame(
  rep = vapply(replicate_runs, `[[`, integer(1), "rep"),
  selected_xi = vapply(replicate_runs, `[[`, numeric(1), "selected_xi")
)
result_file = "allres.RData"
save(allres, replicate_diagnostics, target_study_tau, nrep, ncore,
     file = result_file)
nmeth = 5
me = sapply(1:nmeth, function(i) apply(allres[allres[,1]==i,-1], 2, median, na.rm=T))
se = sapply(1:nmeth, function(i) apply(allres[allres[,1]==i,-1], 2, mad, na.rm=T))
out = sapply(1:nmeth, function(i) sprintf("%.3f (%.3f)", me[,i], se[,i]))

rownames(out) = c("Cindex", "intBS", "RMST")
colnames(out) = c("POTL", "Target-only", "TransCox", "CoxTL", "Pooled")
print(out, quote = F)

library(xtable)
print(xtable(out))

# plot prediction errors
suppressPackageStartupMessages(library(ggplot2))
load("allres.RData")

method_levels = c("POTL", "Target-only", "TransCox", "CoxTL", "Pooled")
metric_variables = c("Cindex", "intBS", "RMST")
metric_labels = c(
  "C-index",
  "IBS",
  "RMST"
)
method_id = as.integer(allres[[1L]])
metric_lower_bound = c(Cindex = 0, intBS = 0, RMST = 0)
metric_upper_bound = c(Cindex = 1, intBS = 1, RMST = Inf)

summary_rows = lapply(seq_along(metric_variables), function(metric_index) {
  variable = metric_variables[metric_index]
  do.call(rbind, lapply(seq_along(method_levels), function(id) {
    values = allres[[variable]][method_id == id]
    center = stats::median(values)
    spread = stats::mad(values, center = center, constant = 1.4826)
    data.frame(
      method = method_levels[id],
      metric = metric_labels[metric_index],
      median = center,
      mad = spread,
      lower = max(metric_lower_bound[[variable]], center-spread),
      upper = min(metric_upper_bound[[variable]], center+spread),
      stringsAsFactors = FALSE
    )
  }))
})
plot_data = do.call(rbind, summary_rows)
plot_data$method = factor(plot_data$method, levels = rev(method_levels))
plot_data$metric = factor(plot_data$metric, levels = metric_labels)
plot_data$median_label = ifelse(
  plot_data$metric == "RMST",
  sprintf("%.2f", plot_data$median),
  sprintf("%.3f", plot_data$median)
)

method_colors = c(
  "POTL" = "blue",
  "Target-only" = "black",
  "TransCox" = "red",
  "CoxTL" = "green",
  "Pooled" = "purple"
)
method_shapes = c(
  "POTL" = 16,
  "Target-only" = 17,
  "TransCox" = 15,
  "CoxTL" = 18,
  "Pooled" = 8
)

# Keep the C-index panel anchored at chance performance. For IBS and RMST,
# round down just below the smallest error-bar endpoint so the intervals fill
# their panels without clipping any results.
ibs_axis_start = max(
  0,
  floor(100 * min(plot_data$lower[plot_data$metric == "IBS"])) / 100
)
rmst_axis_start = max(
  0,
  floor(2 * min(plot_data$lower[plot_data$metric == "RMST"])) / 2
)
axis_anchors = data.frame(
  method = factor(rep("Pooled", 3L), levels = rev(method_levels)),
  metric = factor(metric_labels, levels = metric_labels),
  value = c(0.5, ibs_axis_start, rmst_axis_start)
)

prediction_error_plot = ggplot(
  plot_data,
  aes(y = method, x = median, color = method, shape = method)
) +
  geom_blank(
    data = axis_anchors,
    aes(y = method, x = value),
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y", width = 0.18, linewidth = 0.5
  ) +
  geom_point(size = 2.5, stroke = 0.45) +
  geom_text(
    aes(label = median_label),
    nudge_y = 0.25, color = "grey20", size = 2.2,
    show.legend = FALSE
  ) +
  facet_grid(cols = vars(metric), scales = "free_x") +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_shape_manual(values = method_shapes, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(add = c(0.4, 0.75))
  ) +
  labs(x = "Prediction error", y = NULL) +
  theme_bw(base_size = 9, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(linewidth = 0.45),
    strip.background = element_rect(fill = "grey93", linewidth = 0.45),
    strip.text = element_text(size = 8.5, lineheight = 1.05),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.margin = margin(5, 6, 4, 5, unit = "pt")
  )

pdf_file = file.path(".", "pred_error.pdf")
ggsave(
  pdf_file, prediction_error_plot, device = "pdf",
  width = 7.25, height = 3.25, units = "in", useDingbats = FALSE
)

message("Created: ", pdf_file)
