Rcpp::sourceCpp("EM_diffcovar.cpp")
source("coxTL.R")
source("individual_survival_bias.R")

library(reticulate)
use_condaenv("/home/ygu/anaconda3/envs/TransCoxEnvi")
source_python(system.file("python", "TransCoxFunction.py", package = "TransCox"))

library(intsurv)
library(survival)
library(nleqslv)
library(parallel)
library(TransCox)
library(dplyr)
library(Hmisc)
library(rms)
library(MASS)
library(splines)

# target sample size
n = 100
# validation sample size
nv = 10000
# source sample size
ns = 1000
# study end time
tau = 2
taus = 5

# Grid and fixed validation cohort used only for individual survival-curve bias
bias_times = seq(tau/4, tau, by = tau/4)
bias_validation_seed = 20260821L
simulation_metric_names = c("L2distS", "Stau", "Cindex", "intBS", "RMST")
method_names = c("POTL", "Target-only", "TransCox", "CoxTL", "Pooled")

# parameters
beta1 = 0.5
beta2 = -0.5
beta3 = 0.2
beta4 = -0.1
Lambda = function(t) {
  Lambda_t = log(1+0.5*t)
  return(Lambda_t)
}
r = 0.0
rs = 0.0
# transformation function
G = function(x, r) {
  if(abs(r)<1e-8) {
    G_x = x
  } else {
    G_x = log(1+r*x)/r
  }
  return(G_x)
}

# varying coefficient for time-covariate interaction
ia = function(t) {
  t
}

# Gau-Lag approx
N = 20
# candidate tuning parameter xi
xi_cand = c(0.0, 2^seq(-15,15,1))

# Candidate numbers of non-intercept cubic B-spline basis columns for each
# time-varying coefficient. Selection uses source-only AIC.
# The fitted model lets both coefficients vary; in the current SC4/SC5 data
# generators, the true X2 coefficient remains constant.
source_vc_degree = 3L
source_vc_df_set = 4:5
source_vc_diagnostic_names = c(
  "n_events", "selected_df", "selected_basis_df", "selected_warning_count",
  paste0("AIC_df", source_vc_df_set)
)

# inverse transformation function
G_inv = function(g, r) {
  if(abs(r)<1e-8) {
    G_inv_g = g
  } else {
    G_inv_g = (exp(g*r)-1)/r
  }
  return(G_inv_g)
}

log_sum_exp = function(x) {
  if(length(x)==0L) return(-Inf)
  xmax = max(x)
  if(!is.finite(xmax)) return(xmax)
  xmax+log(sum(exp(x-xmax)))
}

source_vc_basis_spec = function(event_times, spline_df, degree, boundary) {
  n_inner = spline_df-degree
  if(n_inner<0L) {
    stop("spline_df must be at least the spline degree")
  }
  knots = if(n_inner==0L) {
    numeric(0)
  } else {
    unique(as.numeric(stats::quantile(
      event_times,
      probs = seq_len(n_inner)/(n_inner+1),
      names = FALSE
    )))
  }
  list(knots = knots, degree = degree, boundary = boundary)
}

source_vc_basis = function(t, spec) {
  unclass(splines::bs(
    as.numeric(t),
    knots = spec$knots,
    degree = spec$degree,
    intercept = FALSE,
    Boundary.knots = spec$boundary,
    warn.outside = FALSE
  ))
}

fit_source_vc_candidate = function(train, event_times, spline_df, degree,
                                   boundary) {
  spec = source_vc_basis_spec(event_times, spline_df, degree, boundary)
  basis_at = function(t) source_vc_basis(t, spec)
  tt_basis = function(x, t, ...) {
    z = x*basis_at(t)
    colnames(z) = paste0("b", seq_len(ncol(z)))
    z
  }

  warning_messages = character(0)
  fit = withCallingHandlers(
    survival::coxph(
      survival::Surv(time, status) ~ X1+X2+tt(X1)+tt(X2),
      data = train,
      ties = "breslow",
      tt = list(tt_basis, tt_basis),
      singular.ok = FALSE,
      model = FALSE,
      x = FALSE,
      y = FALSE,
      # Generated event times are continuous. Disabling near-tie adjustment
      # keeps coxph risk sets identical to the raw times used below for the
      # manual Breslow baseline-hazard calculation.
      control = survival::coxph.control(iter.max = 50, timefix = FALSE)
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if(length(warning_messages)>0L) {
    stop(
      "varying-coefficient Cox fit issued warning(s): ",
      paste(unique(warning_messages), collapse = " | ")
    )
  }
  coefficients = stats::coef(fit)
  if(length(coefficients)==0L || any(!is.finite(coefficients))) {
    stop("varying-coefficient Cox fit returned non-finite coefficients")
  }
  if(length(fit$loglik)<2L || !is.finite(fit$loglik[2])) {
    stop("varying-coefficient Cox fit returned a non-finite log likelihood")
  }
  basis_df = ncol(basis_at(event_times))
  expected_coef = 2L+2L*basis_df
  if(length(coefficients)!=expected_coef) {
    stop("unexpected varying-coefficient Cox model dimension")
  }
  aic = -2*fit$loglik[2]+2*length(coefficients)

  # Return only the small objects needed for prediction. The coxph fit uses an
  # expanded risk-set representation and should be released before fitting the
  # next candidate.
  list(
    coefficients = coefficients,
    spec = spec,
    spline_df = spline_df,
    basis_df = basis_df,
    aic = aic,
    warning_messages = warning_messages
  )
}

fit_source_varying_cox = function(datas, target_data,
                                  df_set = source_vc_df_set,
                                  degree = source_vc_degree) {
  required_columns = c("X1", "X2", "Y", "Delta")
  if(!all(required_columns %in% colnames(datas)) ||
     !all(required_columns %in% colnames(target_data))) {
    stop("source and target data must contain X1, X2, Y, and Delta")
  }
  train = data.frame(
    X1 = as.numeric(datas[,"X1"]),
    X2 = as.numeric(datas[,"X2"]),
    time = as.numeric(datas[,"Y"]),
    status = as.integer(datas[,"Delta"])
  )
  target = data.frame(
    X1 = as.numeric(target_data[,"X1"]),
    X2 = as.numeric(target_data[,"X2"]),
    time = as.numeric(target_data[,"Y"])
  )
  if(any(!is.finite(as.matrix(train))) || any(!is.finite(as.matrix(target)))) {
    stop("varying-coefficient Cox inputs must be finite")
  }
  if(any(train$time<0) || any(target$time<0) ||
     !all(train$status %in% c(0L, 1L))) {
    stop("invalid survival times or source event indicators")
  }
  df_set = sort(unique(as.integer(df_set)))
  if(length(df_set)==0L || any(df_set<degree)) {
    stop("df_set must contain integers at least as large as degree")
  }
  event_times = sort(unique(train$time[train$status==1L]))
  n_events = sum(train$status==1L)
  if(length(event_times)<2L || n_events<=2L+2L*min(df_set)) {
    stop("too few source events for the varying-coefficient Cox model")
  }
  boundary = c(0, max(train$time))
  if(boundary[2]<=boundary[1]) {
    stop("source follow-up has no positive time range")
  }

  candidates = vector("list", length(df_set))
  candidate_errors = rep(NA_character_, length(df_set))
  aic_values = rep(Inf, length(df_set))
  for(j in seq_along(df_set)) {
    candidate = tryCatch(
      fit_source_vc_candidate(
        train, event_times, df_set[j], degree, boundary
      ),
      error = function(e) {
        candidate_errors[j] <<- conditionMessage(e)
        NULL
      }
    )
    # Single-bracket assignment preserves the candidate position when a fit
    # fails and candidate is NULL.
    candidates[j] = list(candidate)
    if(!is.null(candidate)) {
      aic_values[j] = candidate$aic
    }
  }
  if(!any(is.finite(aic_values))) {
    stop(
      "all varying-coefficient Cox candidates failed: ",
      paste(candidate_errors, collapse = " | ")
    )
  }
  selected_index = which.min(aic_values)
  selected = candidates[[selected_index]]
  coefficients = selected$coefficients
  event_basis = source_vc_basis(event_times, selected$spec)
  gamma1 = coefficients[startsWith(names(coefficients), "tt(X1)")]
  gamma2 = coefficients[startsWith(names(coefficients), "tt(X2)")]
  if(length(gamma1)!=ncol(event_basis) || length(gamma2)!=ncol(event_basis)) {
    stop("could not identify both time-varying coefficient bases")
  }
  beta1_time = unname(coefficients["X1"]+drop(event_basis%*%gamma1))
  beta2_time = unname(coefficients["X2"]+drop(event_basis%*%gamma2))
  if(any(!is.finite(beta1_time)) || any(!is.finite(beta2_time))) {
    stop("estimated source coefficient curves are non-finite")
  }

  event_count = tabulate(
    match(train$time[train$status==1L], event_times),
    nbins = length(event_times)
  )
  log_basehaz_increment = vapply(seq_along(event_times), function(j) {
    risk = train$time>=event_times[j]
    eta = beta1_time[j]*train$X1[risk]+beta2_time[j]*train$X2[risk]
    log(event_count[j])-log_sum_exp(eta)
  }, numeric(1))
  if(any(!is.finite(log_basehaz_increment))) {
    stop("estimated source baseline-hazard increments are non-finite")
  }

  S = vapply(seq_len(nrow(target)), function(i) {
    keep = event_times<=target$time[i]
    if(!any(keep)) return(1.0)
    log_hazard_terms = log_basehaz_increment[keep]+
      beta1_time[keep]*target$X1[i]+beta2_time[keep]*target$X2[i]
    log_cumulative_hazard = log_sum_exp(log_hazard_terms)
    if(log_cumulative_hazard==-Inf) return(1.0)
    if(!is.finite(log_cumulative_hazard) ||
       log_cumulative_hazard>log(.Machine$double.xmax)) return(0.0)
    exp(-exp(log_cumulative_hazard))
  }, numeric(1))
  if(any(!is.finite(S)) || any(S<0 | S>1)) {
    stop("source survival predictions are invalid")
  }

  aic_diagnostics = stats::setNames(
    aic_values, paste0("AIC_df", df_set)
  )
  diagnostics = c(
    n_events = n_events,
    selected_df = selected$spline_df,
    selected_basis_df = selected$basis_df,
    selected_warning_count = length(selected$warning_messages),
    aic_diagnostics
  )
  diagnostics = diagnostics[source_vc_diagnostic_names]
  list(S = S, diagnostics = diagnostics)
}

# simulate data from transformation model with r
sim_data = function(i, censor = TRUE) {
  X1 = ifelse(runif(1)>0.5, 1, 0)
  X2 = runif(1)
  X34 = mvrnorm(n = 1, 
                mu = rep(0, 2), 
                Sigma = matrix(c(1,2^(-0.5),2^(-0.5),1), nrow = 2, ncol = 2))
  X3 = X34[1]
  X4 = X34[2]
  Ti = (exp(G_inv(-log(runif(1)), r)*exp(-beta1*X1-beta2*X2-beta3*X3-beta4*X4))-1)/0.5
  
  # generate censoring time from Unif(1.5,4), 
  # such that censoring rate is around 50%
  if(censor) {
    Ci = runif(1, min = 1.5, max = 4)
    Ci = min(Ci, tau)
    Y = min(Ti, Ci)
    Delta = I(Ti <= Ci)
  } else {
    Y = Ti
    Delta = 1
  }
  
  data = c(i, X1, X2, X3, X4, Y, Delta)
  return(data)
}

make_bias_validation = function(seed = bias_validation_seed) {
  set.seed(seed)
  data = t(sapply(1:nv, sim_data, censor = FALSE))
  colnames(data) = c("ID", "X1", "X2", "X3", "X4", "Y", "Delta")
  rawX = cbind(0:(nv-1),
               rep(0,nv),
               rep(tau,nv),
               data[,c("X1", "X2", "X3", "X4")])
  list(data = data, rawX = rawX)
}

sim_source_data = function(i, sc, r, nu) {
  X1 = ifelse(runif(1)>0.5, 1, 0)
  X2 = runif(1)
  
  if(sc==1) {
    Ti = (exp(G_inv(-log(runif(1)), r)*exp(-beta1*X1-beta2*X2))-1)/nu
  } 
  
  if(sc==2) {
    # Lambda0(t) = nu*t
    Ti = G_inv(-log(runif(1)), r)*exp(-beta1*X1-beta2*X2)/nu
  }
  
  if(sc==3) {
    # different beta
    # Lambda0(t) = nu*t
    Ti = G_inv(-log(runif(1)), r)*exp(-0.7*X1+0.7*X2)/nu
  } 
  
  if(sc==4) {
    # Cox model with time-covariate interaction after tau
    lam_t = function(t) {
      lam0 = 0.5/(1+0.5*t)
      lamt = lam0*exp(beta1*X1+beta2*X2-0.3*I(t>tau)*ia(t)*X1)
      return(lamt)
    }
    cumlam_t = function(t) {
      G(integrate(lam_t, lower=0, upper=t)$value, 0.0)
    }
    tmp = -log(runif(1))
    Ti = nleqslv(0, function(t) cumlam_t(t) - tmp)$x
  }
  
  if(sc==5) {
    # Cox model with time-covariate interaction throughout
    lam_t = function(t) {
      lam0 = 0.5/(1+0.5*t)
      lamt = lam0*exp(0.5*X1-0.5*X2-0.3*ia(t)*X1)
      return(lamt)
    }
    cumlam_t = function(t) {
      G(integrate(lam_t, lower=0, upper=t)$value, 0.0)
    }
    tmp = -log(runif(1))
    # Ti = uniroot(function(t) cumlam_t(t) - tmp, lower=0, upper=9999)$root
    if(cumlam_t(taus)<tmp) {
      Ti = taus+0.1
    } else {
      Ti = nleqslv(0, function(t) cumlam_t(t) - tmp)$x
    }
  }
  
  # generate censoring time from Unif(3.5,7)
  Ci = runif(1, min = 3.5, max = 7)
  Ci = min(Ci, taus)
  Y = min(Ti, Ci)
  Delta = I(Ti <= Ci)
  data = c(i, X1, X2, Y, Delta)
  return(data)
}

sum_metric = function(res, datav, rawXv_bias, pred_r = r) {
  individual_survival_bias = IndividualSurvBias(
    res$beta, res$lambda[,2], res$lambda[,1],
    rawXv_bias, pred_r, r, bias_times
  )
  
  # C-index
  if(length(res$beta)>2) {
    score = datav[,c("X1", "X2", "X3", "X4")] %*% res$beta
  } else {
    score = datav[,c("X1", "X2")] %*% res$beta
  }
  Cindex = cIndex(time = datav[,"Y"],
                  event = datav[,"Delta"],
                  risk_score = score)
  Cindex = Cindex["index"]
  attributes(Cindex) = NULL
  
  metrics = c("L2distS" = res$L2distS,
              "Stau" = res$Stau,
              "Cindex" = Cindex,
              "intBS" = res$intBS,
              "RMST" = res$RMST)
  
  list(metrics = metrics,
       individual_survival_bias = individual_survival_bias)
}

trans_fit = function(data, datav, rawX, rawXv, rawXv_bias,
                     S, weight, r_opt, xi, kmcensor) {
  # transfer learning algorithm
  weight = weight/sum(weight)
  res = TransFit(data[,"Y"], data[,"Delta"], rawX, S, weight, 
                 datav[,"Y"], datav[,"Delta"], rawXv, kmcensor,
                 r_opt, N, xi, tau)
  
  res = sum_metric(res, datav, rawXv_bias, pred_r = r_opt)
  
  return(res)
}

# choose the value of xi by 5-fold cross validation
cross_valid = function(datat, rawXt, S, weight, r_opt, xi, kmcensor) {
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
      train_res = TransFit(train[,"Y"], train[,"Delta"], train_rawX, train_S, train_weight, 
                           test[,"Y"], test[,"Delta"], test_rawX, kmcensor, 
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

run_all_meth = function(rep, sc, bias_validation) {
  set.seed(rep)
  data = t(sapply(1:n, sim_data))
  colnames(data) = c("ID", "X1", "X2", "X3", "X4", "Y", "Delta")
  datav = t(sapply(1:nv, sim_data, censor = FALSE))
  colnames(datav) = c("ID", "X1", "X2", "X3", "X4", "Y", "Delta")
  
  # sort data by Y
  ind = order(data[,"Y"])
  data = data[ind,]
  rawX = cbind(0:(n-1),
               rep(0,n),
               rep(tau,n),
               data[,c("X1", "X2", "X3", "X4")])
  
  # subset of data by removing X3 and X4
  sub_data = data[,-c(4,5)]
  sub_rawX = rawX[,1:(ncol(rawX)-2)]
  
  # KM estimator for censoring distribution
  # no censoring for validation data
  kmcensor = matrix(nrow = 0, ncol = 2)
  
  # sort validation data by Y
  ind = order(datav[,"Y"])
  datav = datav[ind,]
  rawXv = cbind(0:(nv-1),
                rep(0,nv),
                rep(tau,nv),
                datav[,c("X1", "X2", "X3", "X4")])
  
  # simulate source data
  if(sc %in% c(1,4:5)) {
    nus = 0.5
  } else {
    # sc2/3: source model = target model with different parameters
    nus = 0.4
  }
  
  datas = t(sapply(1:ns, sim_source_data, sc=sc, r=rs, nu=nus))
  colnames(datas) = c("ID", "X1", "X2", "Y", "Delta")
  # sort source data by Y
  ind = order(datas[,"Y"])
  datas = datas[ind,]
  rawXs = cbind(0:(ns-1),
                rep(0,ns),
                rep(taus,ns),
                datas[,c("X1", "X2")])
  
  # fit source model
  if(sc %in% 1:3) {
    # obtain source prediction and weights
    fit_source = sourceFit(datas[,"Y"], datas[,"Delta"], rawXs, 
                           data[,"Y"], rawX, rs, N)
    S = fit_source$predSY
  } else if(sc %in% 4:5) {
    # Estimate both source covariate effects as smooth functions of time. The
    # spline complexity is selected using source-only partial-likelihood AIC.
    source_vc_fit = fit_source_varying_cox(datas, data)
    S = source_vc_fit$S
    message(
      "source varying-coefficient spline df = ",
      source_vc_fit$diagnostics["selected_df"],
      " based on source-only AIC."
    )
  } else {
    stop("source prediction model is not defined for this scenario")
  }
  
  # no weight for penalty
  weight = rep(1.0, n)
  
  # tune parameter xi
  cv_loss = tryCatch({
    sapply(xi_cand, function(x) cross_valid(data, rawX, S, weight, r, x, kmcensor))
  }, error = function(e) {
    message("rep ", rep, ": error in CV for selecting xi")
    e
  })
  xi = xi_cand[which.min(cv_loss)]
  message("xi = ", xi, " based on 5-fold cross validation.")
 
  # fit model with optimal xi
  res = trans_fit(data, datav, rawX, rawXv, bias_validation$rawX,
                  S, weight, r, xi, kmcensor)
  
  ##############################################################################
  # compare existing methods
  ## Cox model with target data only
  res_tar = trans_fit(data, datav, rawX, rawXv, bias_validation$rawX,
                      S, weight, r, xi = 0.0, kmcensor)
  
  ## Cox model with both target and source data
  datac = rbind(sub_data, datas)
  rawXc = rbind(sub_rawX, rawXs)
  # sort source data by Y
  ind = order(datac[,"Y"])
  datac = datac[ind,]
  rawXc = rawXc[ind,]
  rawXc[,1] = 0:(n+ns-1)
  rawXc[,3] = 5
  res_com = trans_fit(datac, datav, rawXc, rawXv, bias_validation$rawX,
                      rep(1,n+ns), rep(1,n+ns), r, xi = 0.0, kmcensor)
  
  ## Li et al (2023, JASA)
  # select the best learning rate and number of steps using BIC
  tmpdata = data.frame(data[,c(2,3,6,7)]) %>%
    mutate(time = Y,
           status = Delta+1) %>%
    dplyr::select(-Y, -Delta)
  tmpdatas = data.frame(datas[,2:5]) %>%
    mutate(time = Y,
           status = Delta+1) %>%
    dplyr::select(-Y, -Delta)
  LRres <- SelLR_By_BIC(primData = tmpdata,
                        auxData = tmpdatas,
                        cov = c("X1", "X2"),
                        statusvar = "status", lambda1 = 0.1, lambda2 = 0.1,
                        learning_rate_vec = 10^(seq(-5,0,1)),
                        nsteps_vec = c(100, 200))
  # select the best tuning parameter using BIC
  BICres <- SelParam_By_BIC(primData = tmpdata,
                            auxData = tmpdatas,
                            cov = c("X1", "X2"),
                            statusvar = "status",
                            lambda1_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)),
                            lambda2_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)),
                            learning_rate = LRres$best_lr, nsteps = LRres$best_nsteps)
  Cout <- GetAuxSurv(tmpdatas, cov = c("X1", "X2"))
  Pout <- GetPrimaryParam(tmpdata, q = Cout$q, estR = Cout$estR)
  Tres <- runTransCox_one(Pout, l1 = BICres$best_la1, l2 = BICres$best_la2,
                          learning_rate = LRres$best_lr, nsteps = LRres$best_nsteps,
                          cov = c("X1", "X2"))
  beta_est = Tres$new_beta
  lambda_est = Tres$new_IntH
  uniqt = Tres$time
  res_li23 = Metric(beta_est, lambda_est, uniqt, datav[,"Y"], datav[,"Delta"],
                    rawXv, kmcensor, r, tau)
  res_li23 = sum_metric(res_li23, datav, bias_validation$rawX)
  
  ## CoxTL (Lu et al.)
  p = 2
  data_t = CoxTL_prepare_data(
    data[,-c(1,4,5),drop = FALSE], p, population = 0,
    covariates = c("X1", "X2")
  )
  data_s = CoxTL_prepare_data(
    datas[,-1,drop = FALSE], p, population = 1,
    covariates = c("X1", "X2")
  )
  coxtl_weights = CoxTL_density_weights(data_t, data_s, p)
  cox_tl<-run_CoxTL(data_t, data_s, weights = coxtl_weights, di = p)
  beta_est = cox_tl$coefficients
  cumlam = CoxTL_target_basehaz(cox_tl)
  uniqt = cumlam$time
  lambda_est = cumlam$hazard-c(0.0, cumlam$hazard[-nrow(cumlam)])
  ind = which(lambda_est!=0 & uniqt<tau)
  uniqt = uniqt[ind]
  lambda_est = lambda_est[ind]
  res_coxtl = Metric(beta_est, lambda_est, uniqt, datav[,"Y"], datav[,"Delta"],
                     rawXv, kmcensor, r, tau)
  res_coxtl = sum_metric(res_coxtl, datav, bias_validation$rawX)
  ##############################################################################
  
  # 1: POTL, 2: Target-only, 3: TransCox, 4: CoxTL, 5: Pooled
  method_results = list(res = res, res_tar = res_tar, res_li23 = res_li23,
                        res_coxtl = res_coxtl, res_com = res_com)
  combined_result = combine_method_results(
    method_results, method_names, bias_validation$data[,"ID"], bias_times
  )
  if(sc %in% 4:5) {
    combined_result$source_vc_diagnostics = source_vc_fit$diagnostics
  }
  combined_result
}

# simulation in single scenario
single_sc_sim = function(sc) {
  nrep = 200
  ncore = 200
  mc_preschedule = !(sc %in% 4:5)
  
  bias_validation = make_bias_validation()
  sim_results = mclapply(
    1:nrep, run_all_meth, sc = sc,
    bias_validation = bias_validation, mc.cores = ncore,
    mc.preschedule = mc_preschedule
  )
  if(sc %in% 4:5) {
    source_vc_diagnostics = t(vapply(sim_results, function(x) {
      values = stats::setNames(
        rep(NA_real_, length(source_vc_diagnostic_names)),
        source_vc_diagnostic_names
      )
      if(is.list(x) && !is.null(x$source_vc_diagnostics)) {
        available = intersect(names(x$source_vc_diagnostics), names(values))
        values[available] = as.numeric(x$source_vc_diagnostics[available])
      }
      values
    }, stats::setNames(
      rep(NA_real_, length(source_vc_diagnostic_names)),
      source_vc_diagnostic_names
    )))
    source_vc_diagnostics = data.frame(
      replication = seq_len(nrep),
      scenario = sc,
      source_vc_diagnostics,
      check.names = FALSE
    )
  }
  expected_bias_dim = c(nv, length(bias_times), length(method_names))
  valid_rep = vapply(
    sim_results, valid_simulation_result, logical(1),
    expected_bias_dim = expected_bias_dim
  )
  if(!any(valid_rep)) stop("All simulation replications failed.")
  allres = lapply(sim_results[valid_rep], function(x) x$metrics)
  if(sc %in% 4:5) {
    save(allres, source_vc_diagnostics,
         file = sprintf("simres_sc%d.RData", sc))
    utils::write.csv(
      source_vc_diagnostics,
      file = sprintf("source_vc_diagnostics_sc%d.csv", sc),
      row.names = FALSE
    )
  } else {
    save(allres, file = sprintf("simres_sc%d.RData", sc))
  }
  individual_bias = summarize_individual_survival_bias(
    sim_results, bias_validation$data, bias_times, method_names,
    scenario = sc
  )
  save(individual_bias,
       file = sprintf("individual_survival_bias_single_source_diffcovar_sc%d.RData", sc))
  plot_individual_survival_bias(
    individual_bias,
    file = sprintf("individual_survival_bias_single_source_diffcovar_sc%d.pdf", sc)
  )
  rm(sim_results)
  allres <- as.data.frame(do.call(rbind, allres))
  allres = apply(allres, c(1,2), as.numeric) 
  allres = na.omit(allres)
  
  # aggregate results over all replicates
  nmeth = 5
  me = sapply(1:nmeth, function(i) apply(allres[allres[,1]==i,-1], 2, median, na.rm=T))
  # me[1:2,] = me[1:2,]-tval
  se = sapply(1:nmeth, function(i) apply(allres[allres[,1]==i,-1], 2, mad, na.rm=T))
  out = sapply(1:nmeth, function(i) sprintf("%.3f (%.3f)", me[,i], se[,i]))
  
  metric = c("L2D", "Dtau", "C-index", "IBS", "RMST")
  scenario = c(sc, rep("", 4))
  out = cbind(scenario, metric, out)
  colnames(out) = c("SC", "Metric", "POTL", "Target-only", "TransCox", "CoxTL", "Pooled")
  print(out, quote = F)
  return(out)
}

# formal simulation over all scenarios
res_by_sc = lapply(1:5, single_sc_sim)
res = do.call(rbind, res_by_sc)
write.table(res, file = "res.txt", quote = F, row.names = F, col.names = T)
library(xtable)
print(xtable(res), include.rownames = F)
