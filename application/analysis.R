Rcpp::sourceCpp("EM.cpp")

library(reticulate)
use_condaenv("/home/ygu/anaconda3/envs/TransCoxEnvi")
source_python(system.file("python", "TransCoxFunction.py", package = "TransCox"))

library(intsurv)
library(survival)
library(dplyr)
library(parallel)
library(TransCox)
library(Hmisc)
library(rms)

nrep = 20
ncore = detectCores()-1
is_parallel = 1

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

sum_metric = function(res, datav, rawXv) {
  # C-index
  score = as.matrix(datav[,-c(1,2)]) %*% res$beta
  Cindex = cIndex(time = datav$time,
                  event = datav$status,
                  risk_score = score)
  Cindex = Cindex["index"]
  attributes(Cindex) = NULL
  
  res = c("Cindex" = Cindex,
          "intBS" = res$intBS,
          "RMST" = res$RMST)
  
  return(res)
}

trans_fit = function(datat, datav, rawXt, rawXv, S, weight, r_opt, xi, kmcensor, N, tau) {
  # transfer learning algorithm
  weight = weight/sum(weight)
  res = TransFit(datat$time, datat$status, rawXt, S, weight, 
                 datav$time, datav$status, rawXv, kmcensor,
                 r_opt, N, xi, tau)
  
  res = sum_metric(res, datav, rawXv)
  
  return(res)
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
  
  # study end time
  tau = max(data$time)
  taus = max(datas$time)
  
  # Gau-Lag approx
  N = 20
  # candidate tuning parameter xi
  xi_cand = c(0.0, 2^seq(-15,15,1))
  # candidate value of r
  r_cand = rs_cand = seq(0, 1.5, 0.05)
  
  # divide data into training and testing sets
  set.seed(20+repnum)
  train_ind = sort(sample(1:n, round(n*ptrain)))
  # training data
  datat = data[train_ind,]
  taut = max(datat$time)
  rawXt = cbind(0:(nt-1),
                rep(0,nt),
                rep(tau,nt),
                datat[,-c(1,2)])
  rawXt = matrix(unlist(rawXt), nrow = nt)
  # validation data
  datav = data[-train_ind,]
  rawXv = cbind(0:(n-nt-1),
                rep(0,n-nt),
                rep(tau,n-nt),
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
  message("optimal r for source study based on AIC: ", rs_opt)
  
  # # select the optimal r
  # aic = sapply(r_cand, function(r_try) {
  #   fit_target = sourceFit(datat$time, datat$status, rawXt,
  #                          datat$time, rawXt, r_try, N, verbose = 0)
  #   -2*fit_target$logL
  # })
  # r_opt = r_cand[which.min(aic)]
  r_opt = 1.5
  message("optimal r for target study based on AIC: ", r_opt)
  
  # obtain source prediction 
  fit_source = sourceFit(datas$time, datas$status, rawXs, 
                         datat$time, rawXt, rs_opt, N, compVar = 1)
  S = fit_source$predSY
  
  # no weights
  weight = rep(1.0, length(S))
  
  # obtain KM estimator for censoring distribution
  km_fit = survfit(Surv(datat$time, 1-datat$status) ~ 1)
  kmcensor = cbind(km_fit$time, km_fit$surv)
  kmcensor = rbind(c(0, 1.0), kmcensor)
  kmcensor = kmcensor[-nrow(kmcensor),]
  
  ### Proposed method
  # tune parameter xi
  cv_loss = tryCatch({
    sapply(xi_cand, function(x) cross_valid(datat, rawXt, S, weight, r_opt, x, kmcensor, N, taut))
  }, error = function(e) {
    message("rep ", repnum, ": error in CV for selecting xi")
    e
  })
  xi = xi_cand[which.min(cv_loss)]
  message("xi = ", xi, " based on 5-fold cross validation.")
  
  # fit model with optimal xi
  res = tryCatch({
    trans_fit(datat, datav, rawXt, rawXv, S, weight, r_opt, xi, kmcensor, N, taut)
  }, error = function(e) {
    message("rep ", repnum, ": error in TransFit with best xi")
    e
  })
  
  ### Target-only method
  res_tar = tryCatch({
    trans_fit(datat, datav, rawXt, rawXv, S, weight, r_opt, xi = 0.0, kmcensor, N, taut)
  }, error = function(e) {
    message("rep ", repnum, ": error in target_fit")
    e
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
                        nsteps_vec = c(100, 200))
  # select the best tuning parameter using BIC
  if(!is.na(LRres$best_lr) & !is.na(LRres$best_nsteps)) {
    best_lr = LRres$best_lr
    best_nsteps = LRres$best_nsteps
  } else {
    best_lr = 0.001
    best_nsteps = 100
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
  res_li23 = Metric(beta_est, lambda_est, uniqt, datav$time, datav$status,
                    rawXv, kmcensor, 0.0, taut)
  res_li23 = sum_metric(res_li23, datav, rawXv)
  
  ## CoxTL (Lu et al.)
  source("coxTL.R")
  p = 17
  data_t = data
  data_s = datas
  cox_tl<-run_CoxTL(data_t,data_s, p)
  beta_est = cox_tl$coefficients
  cumlam = basehaz(cox_tl, centered = F)
  uniqt = cumlam$time
  lambda_est = cumlam$hazard-c(0.0, cumlam$hazard[-nrow(cumlam)])
  ind = which(lambda_est!=0 & uniqt<tau)
  uniqt = uniqt[ind]
  lambda_est = lambda_est[ind]
  res_coxtl = Metric(beta_est, lambda_est, uniqt, datav$time, datav$status,
                     rawXv, kmcensor, 0.0, taut)
  res_coxtl = sum_metric(res_coxtl, datav, rawXv)
  
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
  message("optimal r for combined data based on AIC: ", rc_opt)
  res_com = trans_fit(datac, datav, rawXc, rawXv, rep(1,nt+ns), rep(1,nt+ns),
                      rc_opt, xi = 0.0, kmcensor, N, taut)
  
  allres = rbind(res, res_tar, res_li23, res_coxtl, res_com)
  allres = cbind(1:5, allres)
    
  return(allres)
}

if(is_parallel) {
  allres = mclapply(1:nrep, function(i) eval_pred(i, data, datas, ptrain), mc.cores = ncore)
} else {
  allres = sapply(1:nrep, function(i) eval_pred(i, data, datas, ptrain))
}

# aggregate results from all replicates
allres <- as.data.frame(do.call(rbind, allres))
nmeth = 5
me = sapply(1:nmeth, function(i) apply(allres[allres[,1]==i,-1], 2, median, na.rm=T))
se = sapply(1:nmeth, function(i) apply(allres[allres[,1]==i,-1], 2, mad, na.rm=T))
out = sapply(1:nmeth, function(i) sprintf("%.3f (%.3f)", me[,i], se[,i]))

rownames(out) = c("Cindex", "intBS", "RMST")
colnames(out) = c("POTL", "Target-only", "TransCox", "CoxTL", "Pooled")
print(out, quote = F)

library(xtable)
print(xtable(out))

