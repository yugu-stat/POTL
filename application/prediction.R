# predict survival probability for future patients using various methods

Rcpp::sourceCpp("EM_invvar.cpp")

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

# read target and source data
data = read.table("target_data_pca.txt", header = T)
datas = read.table("source_data_pca.txt", header = T)

# standardize continuous covariates
data[,c(3,4,10:19)] = scale(data[,c(3,4,10:19)], center = F)
datas[,c(3,4,10:19)] = scale(datas[,c(3,4,10:19)], center = F)
# subject 155 has median age and nlymph

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

# transformation function
G = function(x, r) {
  if(abs(r)<1e-8) {
    G_x = x
  } else {
    G_x = log(1+r*x)/r
  }
  return(G_x)
} 

# compute survival function estimator based on beta and lambda estimates 
surv_prob = function(newX, beta, lambda, r) {
  # lambda: unique time | lambda est
  cumlam = c(0, cumsum(lambda[,2]))
  uniqt = c(0, lambda[,1])
  ebx = exp(newX %*% beta)
  St = data.frame("time" = uniqt,
                  "surv0" = exp(-G(cumlam*ebx[1,1], r)),
                  "surv1" = exp(-G(cumlam*ebx[2,1], r)))
  return(St)
}

trans_fit = function(data, rawX, newX, S, weight, r_opt, xi, kmcensor, N, tau) {
  # transfer learning algorithm
  weight = weight/sum(weight)
  # create artificial validation set, never used
  datav = head(data)
  rawXv = head(rawX)
  res = TransFit(data$time, data$status, rawX, S, weight, 
                 datav$time, datav$status, rawXv, kmcensor,
                 r_opt, N, xi, tau)
  St = surv_prob(newX, res$beta, res$lambda, r_opt)
  
  return(St)
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

# predict survival probabilities
# sample size
n = nrow(data)
ns = nrow(datas)

# study end time
tau = max(data$time)
taus = max(datas$time)

# Gau-Lag approx
N = 20
# candidate tuning parameter xi
xi_cand = c(0.0, 2^seq(-15,15,1))
# candidate value of r
r_cand = rs_cand = seq(0, 1.5, 0.05)

# target and source covariates
rawX = cbind(0:(n-1),
              rep(0,n),
              rep(tau,n),
              data[,-c(1,2)])
rawX = matrix(unlist(rawX), nrow = n)
rawXs = cbind(0:(ns-1),
              rep(0,ns),
              rep(taus,ns),
              datas[,-c(1,2)])
rawXs = matrix(unlist(rawXs), nrow = ns)

# new patient's covariates
medX = apply(data[,-c(1,2)], 2, median)
newX = matrix(medX, nrow = 2, ncol = length(medX), byrow = T)
newX[2,3] = 1

# # select the optimal rs
# aic_rs = sapply(rs_cand, function(rs_try) {
#   cat("rs = ", rs_try, "\n")
#   fit_source = sourceFit(datas$time, datas$status, rawXs,
#                          datat$time, rawXt, rs_try, N, verbose = 1)
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

# obtain source prediction and weights
fit_source = sourceFit(datas$time, datas$status, rawXs, 
                       data$time, rawX, rs_opt, N, compVar = 1)
S = fit_source$predSY
weight = rep(1.0, length(S))

# obtain KM estimator for censoring distribution
km_fit = survfit(Surv(data$time, 1-data$status) ~ 1)
kmcensor = cbind(km_fit$time, km_fit$surv)
kmcensor = rbind(c(0, 1.0), kmcensor)
kmcensor = kmcensor[-nrow(kmcensor),]

### Proposed method
# tune parameter xi
cv_loss = tryCatch({
  sapply(xi_cand, function(x) cross_valid(data, rawX, S, weight, r_opt, x, kmcensor, N, tau))
}, error = function(e) {
  message("Error in CV for selecting xi")
  e
})
#print(cv_loss)
xi = xi_cand[which.min(cv_loss)]
# xi = 100
message("xi = ", xi, " based on 5-fold cross validation.")
# xi = 100.0
# fit model with optimal xi
res = tryCatch({
  trans_fit(data, rawX, newX, S, weight, r_opt, xi, kmcensor, N, tau)
}, error = function(e) {
  message("Error in TransFit with best xi")
  e
})
res$method = 1
if(max(res$time)<tau) {
  res = rbind(res, tail(res, n=1))
  res$time[nrow(res)] = tau
}

### Target-only method
res_tar = tryCatch({
  trans_fit(data, rawX, newX, S, weight, r_opt, xi = 0.0, kmcensor, N, tau)
}, error = function(e) {
  message("Error in target_fit")
  e
})
res_tar$method = 2
if(max(res_tar$time)<tau) {
  res_tar = rbind(res_tar, tail(res_tar, n=1))
  res_tar$time[nrow(res_tar)] = tau
}

### TransCox
tmpdata = data.frame(data) %>%
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
lambda_est = cbind(Tres$time, Tres$new_IntH)
res_li23 = surv_prob(newX, beta_est, lambda_est, 0.0)
res_li23$method = 3
if(max(res_li23$time)<tau) {
  res_li23 = rbind(res_li23, tail(res_li23, n=1))
  res_li23$time[nrow(res_li23)] = tau
}

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
lambda_est = cbind(uniqt[ind], lambda_est[ind])
res_coxtl = surv_prob(newX, beta_est, lambda_est, 0.0)
res_coxtl$method = 4
if(max(res_coxtl$time)<tau) {
  res_coxtl = rbind(res_coxtl, tail(res_coxtl, n=1))
  res_coxtl$time[nrow(res_coxtl)] = tau
}

### Combined method
datac = rbind(data, datas)
rawXc = rbind(rawX, rawXs)
# sort source data by Y
ind = order(datac$time)
datac = datac[ind,]
rawXc = rawXc[ind,]
rawXc[,1] = 0:(n+ns-1)
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
res_com = trans_fit(datac, rawXc, newX, rep(1,n+ns), rep(1,n+ns),
                    rc_opt, xi = 0.0, kmcensor, N, tau)
res_com$method = 5
if(max(res_com$time)<tau) {
  res_com = rbind(res_com, tail(res_com, n=1))
  res_com$time[nrow(res_com)] = tau
}

allres = rbind(res, res_tar, res_li23, res_coxtl, res_com)

# plot survival curves
library(ggplot2)
library(ggpubr)

allres = allres %>%
  filter(time<=20) %>%
  mutate(method = factor(method))

fig0 = ggplot(data = allres)+
  geom_step(aes(x=time, y=surv0, color=method))+
  scale_x_continuous(name="Time (years)", 
                     limits=c(0, 20),
                     breaks = seq(0,20,5),
                     labels = seq(0,20,5))+
  scale_y_continuous(name="Survival Probability", 
                     limits = c(floor(min(allres$surv1)/0.1)*0.1, 1), 
                     breaks = seq(floor(min(allres$surv1)/0.1)*0.1,1,0.1))+
  scale_color_manual(name = NULL,
                     values = c("red", "blue", "green", "magenta", "cyan"),
                     labels = c("POTL", "Target-only", "TransCox", "CoxTL", "Pooled"))+
  labs(title = "Early Stage")+
  theme_bw()+
  theme(axis.title.x = element_text(margin = margin(t=6)),
        axis.title.y = element_text(margin = margin(r=6)),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.key = element_blank(),
        plot.margin = margin(l=10, r=10, t=5, b=5))

fig1 = ggplot(data = allres)+
  geom_step(aes(x=time, y=surv1, color=method))+
  scale_x_continuous(name="Time (years)", 
                     limits=c(0, 20),
                     breaks = seq(0,20,5),
                     labels = seq(0,20,5))+
  scale_y_continuous(name="Survival Probability", 
                     limits = c(floor(min(allres$surv1)/0.1)*0.1, 1), 
                     breaks = seq(floor(min(allres$surv1)/0.1)*0.1,1,0.1))+
  scale_color_manual(name = NULL,
                     values = c("red", "blue", "green", "magenta", "cyan"),
                     labels = c("POTL", "Target-only", "TransCox", "CoxTL", "Pooled"))+
  labs(title = "Advanced Stage")+
  theme_bw()+
  theme(axis.title.x = element_text(margin = margin(t=6)),
        axis.title.y = element_text(margin = margin(r=6)),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.key = element_blank(),
        plot.margin = margin(l=10, r=10, t=5, b=5))

fig = ggarrange(plotlist = list(fig0, fig1), 
                ncol = 2, 
                legend = "right", 
                common.legend = T)
ggsave(fig, file = "surv_pred.pdf",
       width = 10, height = 5, units = "in")
