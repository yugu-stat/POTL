# Predict survival probabilities for two representative future patients using
# the same data processing, models, tuning criterion, and target-study endpoint
# as analysis.R. Run this script from the application directory:
#   Rscript prediction.R

# Avoid oversubscribing BLAS/OpenMP/TensorFlow. Values explicitly set by the
# caller are respected.
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

library(survival)
library(dplyr)
library(TransCox)
library(ggplot2)
library(ggpubr)

covariate_names = c(
  "age", "nlymph", "stage34", "IDC", "ER", "PR", "HER2",
  paste0("PC", 1:10)
)
method_levels = c("POTL", "Target-only", "TransCox", "CoxTL", "Pooled")

# Match the palette in simulations/plot_main_simulation_figures.R.
method_colors = c(
  "POTL" = "blue",
  "Target-only" = "black",
  "TransCox" = "red",
  "CoxTL" = "green",
  "Pooled" = "purple"
)

# Read and preprocess exactly the files used by analysis.R.
data = read.table("processed_target_data.txt", header = TRUE)
datas = read.table("processed_source_data.txt", header = TRUE)

data[,c(3,4,10:19)] = scale(data[,c(3,4,10:19)], center = FALSE)
datas[,c(3,4,10:19)] = scale(datas[,c(3,4,10:19)], center = FALSE)

data = data %>%
  filter(time > 0) %>%
  arrange(time)
datas = datas %>%
  filter(time > 0) %>%
  arrange(time)

target_study_tau = max(data$time)
source_study_tau = max(datas$time)

# Match the endpoint status changes used for numerical convergence.
ns = nrow(datas)
datas$status[c(ns-1, ns)] = 0
data$status[nrow(data)] = 0

make_censor_km = function(d) {
  fit = survfit(Surv(d$time, 1-d$status) ~ 1)
  time = as.numeric(fit$time)
  censor_survival = as.numeric(fit$surv)

  # Match analysis.R by carrying the last positive censoring estimate through
  # the tail for evaluation at the target-study endpoint.
  keep = censor_survival > 0
  list(time = time[keep], surv = censor_survival[keep])
}

make_cpp_censor_matrix = function(censor_km) {
  rbind(c(0, 1.0), cbind(censor_km$time, censor_km$surv))
}

# Construct a right-continuous step survival curve for the two representative
# patients. The baseline's second column must contain hazard increments.
survival_curve = function(newX, beta, baseline, r, method_name) {
  newX = as.matrix(newX)
  beta = as.numeric(beta)
  baseline = as.matrix(baseline)

  if(nrow(newX) != 2L || length(beta) != ncol(newX) ||
     any(!is.finite(newX)) || any(!is.finite(beta))) {
    stop(method_name, ": beta and representative-patient covariates are invalid")
  }
  if(ncol(baseline) != 2L || any(!is.finite(baseline)) ||
     !is.finite(r) || r < 0) {
    stop(method_name, ": baseline or transformation parameter is invalid")
  }
  if(any(baseline[,2] < -sqrt(.Machine$double.eps))) {
    stop(method_name, ": fitted baseline has a negative hazard increment")
  }

  baseline[,2] = pmax(baseline[,2], 0)
  if(any(baseline[,2] > 0 & baseline[,1] <= 0)) {
    stop(method_name, ": fitted baseline has a nonpositive hazard time")
  }
  keep = baseline[,2] > 0 & baseline[,1] > 0 &
    baseline[,1] <= target_study_tau
  baseline = baseline[keep,,drop = FALSE]
  if(nrow(baseline)) {
    baseline = baseline[order(baseline[,1]),,drop = FALSE]
    curve_time = c(0, baseline[,1])
    baseline_hazard = c(0, cumsum(baseline[,2]))
  } else {
    curve_time = 0
    baseline_hazard = 0
  }

  eta = drop(newX %*% beta)
  conditional_hazard = outer(exp(pmin(eta, 700)), baseline_hazard)
  if(abs(r) < 1e-8) {
    prediction = exp(-conditional_hazard)
  } else {
    prediction = exp(-log1p(r*conditional_hazard)/r)
  }
  if(any(!is.finite(prediction)) || any(prediction < 0) ||
     any(prediction > 1)) {
    stop(method_name, ": fitted model produced invalid survival probabilities")
  }

  result = data.frame(
    time = curve_time,
    surv0 = prediction[1,],
    surv1 = prediction[2,]
  )
  if(tail(result$time, 1L) < target_study_tau) {
    endpoint = tail(result, 1L)
    endpoint$time = target_study_tau
    result = rbind(result, endpoint)
  }
  result
}

# Fit POTL, target-only, or pooled transformation models. TransFit requires a
# validation argument to compute metrics, but those artificial metrics are not
# used or reported here; only the training-fit beta and baseline are retained.
fit_transformation_curve = function(fit_data, fit_rawX, newX, S, weight,
                                    r_opt, xi, kmcensor, N, method_name) {
  if(length(S) != nrow(fit_data) || any(!is.finite(S)) ||
     any(S < 0) || any(S > 1) || length(weight) != nrow(fit_data) ||
     any(!is.finite(weight)) || any(weight < 0) || sum(weight) <= 0) {
    stop(method_name, ": source predictions or weights are invalid")
  }
  weight = weight/sum(weight)
  first_event = which(fit_data$status == 1L)[1L]
  first_censored = which(fit_data$status == 0L)[1L]
  if(is.na(first_event) || is.na(first_censored)) {
    stop(method_name, ": fit data must contain an event and a censored case")
  }
  validation_index = sort(unique(c(
    seq_len(min(6L, nrow(fit_data))), first_event, first_censored
  )))
  validation_data = fit_data[validation_index,,drop = FALSE]
  validation_rawX = fit_rawX[validation_index,,drop = FALSE]
  validation_rawX[,1] = 0:(nrow(validation_rawX)-1L)

  fit = TransFit(
    fit_data$time, fit_data$status, fit_rawX, S, weight,
    validation_data$time, validation_data$status, validation_rawX,
    kmcensor, r_opt, N, xi, target_study_tau, verbose = 0
  )
  list(
    fit = fit,
    curve = survival_curve(
      newX, fit$beta, fit$lambda, r_opt, method_name
    )
  )
}

# This is intentionally identical to analysis.R: xi is chosen by the
# five-fold, five-repeat negative profile-likelihood criterion.
cross_valid = function(datat, rawXt, S, weight, r_opt, xi, kmcensor, N) {
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
                           r_opt, N, xi, target_study_tau,
                           verbose = 0, pll = 1, maxit = 100)

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

add_method = function(curve, method_name) {
  curve$method = method_name
  curve
}

# Full target/source design matrices.
n = nrow(data)
ns = nrow(datas)
N = 20
xi_cand = c(0.0, 2^seq(-15,15,1))

rawX = cbind(
  0:(n-1L), rep(0,n), rep(target_study_tau,n),
  data[,-c(1,2)]
)
rawX = matrix(unlist(rawX), nrow = n)
rawXs = cbind(
  0:(ns-1L), rep(0,ns), rep(source_study_tau,ns),
  datas[,-c(1,2)]
)
rawXs = matrix(unlist(rawXs), nrow = ns)

# Two representative patients differ only in disease stage.
median_covariates = apply(data[,covariate_names,drop = FALSE], 2, median)
newX = rbind(median_covariates, median_covariates)
newX[1,"stage34"] = 0
newX[2,"stage34"] = 1
storage.mode(newX) = "double"
rownames(newX) = c("Early Stage", "Advanced Stage")

rs_opt = 0.1
r_opt = 1.5
rc_opt = 0.1

# Source survival predictions used by POTL.
fit_source = sourceFit(
  datas$time, datas$status, rawXs,
  data$time, rawX, rs_opt, N, verbose = 0, compVar = 1
)
S = fit_source$predSY
weight = rep(1.0, length(S))

# Match analysis.R's positive-tail censoring convention.
censor_km = make_censor_km(data)
kmcensor = make_cpp_censor_matrix(censor_km)

### POTL
# Set a seed so the full-data profile-likelihood CV result and figure are
# reproducible. The CV algorithm itself is unchanged from analysis.R.
set.seed(1L)
cv_loss = tryCatch({
  sapply(
    xi_cand,
    function(x) cross_valid(
      data, rawX, S, weight, r_opt, x, kmcensor, N
    )
  )
}, error = function(e) {
  stop("CV for xi failed: ", conditionMessage(e), call. = FALSE)
})
if(length(cv_loss) != length(xi_cand) || !any(is.finite(cv_loss))) {
  stop("CV for xi returned no finite candidate")
}
xi = xi_cand[which.min(ifelse(is.finite(cv_loss), cv_loss, Inf))]
message("xi = ", xi, " based on 5-fold profile-likelihood cross-validation.")

potl = fit_transformation_curve(
  data, rawX, newX, S, weight, r_opt, xi, kmcensor, N,
  "POTL"
)
res = add_method(potl$curve, "POTL")

### Target-only
target_only = fit_transformation_curve(
  data, rawX, newX, S, weight, r_opt, xi = 0.0, kmcensor, N,
  "Target-only"
)
res_tar = add_method(target_only$curve, "Target-only")

### TransCox
tmpdata = data.frame(data) %>% mutate(status = status+1)
tmpdatas = data.frame(datas) %>% mutate(status = status+1)
LRres = SelLR_By_BIC(
  primData = tmpdata,
  auxData = tmpdatas,
  cov = covariate_names,
  statusvar = "status", lambda1 = 0.1, lambda2 = 0.1,
  learning_rate_vec = 10^(seq(-5,0,1)),
  nsteps_vec = as.integer(c(100,200))
)
if(!is.na(LRres$best_lr) && !is.na(LRres$best_nsteps)) {
  best_lr = LRres$best_lr
  best_nsteps = as.integer(LRres$best_nsteps)
} else {
  best_lr = 0.001
  best_nsteps = 100L
}

BICres = SelParam_By_BIC(
  primData = tmpdata,
  auxData = tmpdatas,
  cov = covariate_names,
  statusvar = "status",
  lambda1_vec = c(0.1, 0.5, seq(1,10,by = 0.5)),
  lambda2_vec = c(0.1, 0.5, seq(1,10,by = 0.5)),
  learning_rate = best_lr, nsteps = best_nsteps
)
Cout = GetAuxSurv(tmpdatas, cov = covariate_names)
Pout = GetPrimaryParam(tmpdata, q = Cout$q, estR = Cout$estR)
Tres = runTransCox_one(
  Pout, l1 = BICres$best_la1, l2 = BICres$best_la2,
  learning_rate = best_lr, nsteps = best_nsteps,
  cov = covariate_names
)
res_li23 = survival_curve(
  newX, Tres$new_beta, cbind(Tres$time, Tres$new_IntH),
  r = 0.0, method_name = "TransCox"
)
res_li23 = add_method(res_li23, "TransCox")

### CoxTL (Lu et al.)
source("coxTL.R")
p = length(covariate_names)
data_t = data[,c(covariate_names, "time", "status")]
data_s = datas[,c(covariate_names, "time", "status")]
data_t$R = factor(0, levels = c(0,1))
data_s$R = factor(1, levels = c(0,1))
coxtl_weights = CoxTL_density_weights(data_t, data_s, p)
cox_tl = run_CoxTL(data_t, data_s, weights = coxtl_weights, di = p)
cumlam = CoxTL_target_basehaz(cox_tl)
lambda_est = cumlam$hazard-c(0.0, cumlam$hazard[-nrow(cumlam)])
res_coxtl = survival_curve(
  newX, cox_tl$coefficients, cbind(cumlam$time, lambda_est),
  r = 0.0, method_name = "CoxTL"
)
res_coxtl = add_method(res_coxtl, "CoxTL")

### Pooled
datac = rbind(data, datas)
rawXc = rbind(rawX, rawXs)
pooled_order = order(datac$time)
datac = datac[pooled_order,]
rawXc = rawXc[pooled_order,]
rawXc[,1] = 0:(n+ns-1L)
rawXc[,3] = source_study_tau

pooled = fit_transformation_curve(
  datac, rawXc, newX, rep(1,n+ns), rep(1,n+ns),
  rc_opt, xi = 0.0, kmcensor, N, "Pooled"
)
res_com = add_method(pooled$curve, "Pooled")

allres = rbind(res, res_tar, res_li23, res_coxtl, res_com)
allres$method = factor(allres$method, levels = method_levels)

if(any(!is.finite(as.matrix(allres[,c("time", "surv0", "surv1")]))) ||
   any(allres$time < 0) || any(allres$time > target_study_tau) ||
   any(allres$surv0 < 0 | allres$surv0 > 1) ||
   any(allres$surv1 < 0 | allres$surv1 > 1)) {
  stop("assembled survival curves failed validation")
}

# Save the curve data and tuning result so the expensive fit is reproducible
# without rerunning all methods.
save(
  allres, newX, cv_loss, xi, target_study_tau,
  file = "surv_pred.RData"
)

# Figure generation ---------------------------------------------------------
make_survival_panel = function(survival_column, panel_title) {
  ggplot(
    data = allres,
    aes(x = time, y = .data[[survival_column]], color = method)
  )+
    geom_step(linewidth = 0.7)+
    scale_x_continuous(name="Time (years)", 
                       limits=c(0, target_study_tau),
                       breaks = seq(0,target_study_tau,5),
                       labels = seq(0,target_study_tau,5))+
    scale_y_continuous(name="Survival Probability", 
                       limits = c(floor(min(allres$surv1)/0.1)*0.1, 1), 
                       breaks = seq(floor(min(allres$surv1)/0.1)*0.1,1,0.1))+
    scale_color_manual(
      name = NULL, values = method_colors,
      breaks = method_levels, drop = FALSE
    )+
    labs(title = panel_title)+
    theme_bw()+
    theme(axis.title.x = element_text(margin = margin(t=6)),
          axis.title.y = element_text(margin = margin(r=6)),
          legend.background = element_rect(fill="transparent",colour=NA),
          legend.key = element_blank(),
          plot.margin = margin(l=10, r=10, t=5, b=5))
}

fig0 = make_survival_panel("surv0", "Early Stage")
fig1 = make_survival_panel("surv1", "Advanced Stage")
fig = ggarrange(
  plotlist = list(fig0,fig1), ncol = 2,
  legend = "right", common.legend = TRUE
)
ggsave(
  filename = "surv_pred.pdf", plot = fig, device = "pdf",
  width = 10, height = 5, units = "in", useDingbats = FALSE
)
message("Saved survival curves to surv_pred.pdf and surv_pred.RData")
