Rcpp::sourceCpp("EM.cpp")
source("coxTL.R")

library(intsurv)
library(survival)
library(nleqslv)
library(parallel)
library(dplyr)

# target sample size
n = 100
# validation sample size
nv = 10000
# study end time
tau = 2
# source sample size
ns = n*c(2,3,4,5,6)
# source study end time
taus = c(4,4,5,6,6)
nsc = 5

# pooling method: cv or ns
pool = "cv"

# parameters
beta1 = 0.5
beta2 = -0.5
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
# Gau-Lag approx
N = 20
# candidate tuning parameter xi
xi_cand = c(0.0, 2^seq(-15,15,1))
# candidate xi_k for each source
xi_cand_K = 100*seq(1.0, 2.0, 0.5)
# candidate value of r
r_cand = seq(0, 1.5, 0.05)

# inverse transformation function
G_inv = function(g, r) {
  if(abs(r)<1e-8) {
    G_inv_g = g
  } else {
    G_inv_g = (exp(g*r)-1)/r
  }
  return(G_inv_g)
}

# varying coefficient for time-covariate interaction
ia = function(t) {
  t
}

# simulate data from transformation model with r
sim_data = function(i, censor = TRUE) {
  X1 = ifelse(runif(1)>0.5, 1, 0)
  X2 = runif(1)
  Ti = (exp(G_inv(-log(runif(1)), r)*exp(-beta1*X1-beta2*X2))-1)/0.5
  # # generate censoring time from Exp(0.1)
  # Ci = -log(runif(1))/0.1
  # generate censoring time from Unif(1.5,4), 
  # such that censoring rate is around 50%
  if(censor) {
    Ci = runif(1, min = 1.5, max = 4)
    Ci = min(Ci, tau)
    # Ci = tau
    Y = min(Ti, Ci)
    Delta = I(Ti <= Ci)
  } else {
    Y = Ti
    Delta = 1
  }
  
  data = c(i, X1, X2, Y, Delta)
  return(data)
}

sim_source_data = function(i,sc) {
  X1 = ifelse(runif(1)>0.5, 1, 0)
  X2 = runif(1)
  # # change X2 to beta(1,2) for source data
  # X2 = rbeta(1, 1, 2)
  
  if(sc==1) {
    rs = r
    nus = 0.5
    Ti = (exp(G_inv(-log(runif(1)), rs)*exp(-beta1*X1-beta2*X2))-1)/nus
  } 
  
  if(sc==2) {
    rs = r
    nus = 0.4
    # Lambda0(t) = nu*t
    Ti = G_inv(-log(runif(1)), rs)*exp(-beta1*X1-beta2*X2)/nus
  }
  
  if(sc==3) {
    rs = r
    nus = 0.4
    # different beta values
    # Lambda0(t) = nu*t
    Ti = G_inv(-log(runif(1)), rs)*exp(-0.7*X1+0.7*X2)/nus
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
    if(cumlam_t(taus[sc])<tmp) {
      Ti = taus[sc]+0.1
    } else {
      Ti = nleqslv(0, function(t) cumlam_t(t) - tmp)$x
    }
  }
  
  # generate censoring time from Unif(3.5,7)
  Ci = runif(1, min = 3.5, max = 7)
  Ci = min(Ci, taus[sc])
  # Ci = taus
  Y = min(Ti, Ci)
  Delta = I(Ti <= Ci)
  data = c(i, X1, X2, Y, Delta)
  return(data)
}

get_source_predictor = function(sc, data, rawX, kmcensor) {
  # simulate source data
  tmp_n = ns[sc]
  tmp_tau = taus[sc]
  datas = t(sapply(1:tmp_n, sim_source_data, sc=sc))
  colnames(datas) = c("ID", "X1", "X2", "Y", "Delta")
  # sort source data by Y
  ind = order(datas[,"Y"])
  datas = datas[ind,]
  rawXs = cbind(0:(tmp_n-1),
                rep(0,tmp_n),
                rep(tmp_tau,tmp_n),
                datas[,c("X1", "X2")])
  
  if(sc %in% 1:3) {
    # obtain source prediction and weights
    fit_source = sourceFit(datas[,"Y"], datas[,"Delta"], rawXs, 
                           data[,"Y"], rawX, rs, N)
    S = fit_source$predSY
  } else if(sc==4) {
    # use event times + tau as split points
    train = as.data.frame(datas)
    colnames(train) = c("ID", "X1", "X2", "time", "status")
    cuts <- sort(unique(train$time[train$status == 1]))
    cuts <- sort(unique(c(cuts, tau)))
    cuts <- cuts[cuts < max(train$time)]
    
    train_long <- survSplit(
      Surv(time, status) ~ .,
      data  = train,
      cut   = cuts,
      start = "tstart",
      end   = "tstop",
      event = "status"
    )
    
    # explicit time-dependent covariate
    train_long$Z <- with(train_long, X1 * (tstop > tau) * ia(tstop))
    
    fit_cp <- coxph(
      Surv(tstart, tstop, status) ~ X1 + X2 + Z,
      data = train_long,
      ties = "breslow",
      x = TRUE
    )
    
    bh <- basehaz(fit_cp, centered = FALSE)
    bh$dh <- c(bh$hazard[1], diff(bh$hazard))
    pred_surv_one <- function(x1, x2, y, fit, bh, tau) {
      cf <- coef(fit)
      
      b1 <- cf["X1"]
      b2 <- cf["X2"]
      g  <- cf["Z"]
      
      idx <- bh$time <= y
      if (!any(idx)) return(1)
      
      tt <- bh$time[idx]
      dh <- bh$dh[idx]
      
      eta_t <- b1 * x1 + b2 * x2 + g * x1 * (tt > tau) * ia(tt)
      Hhat <- sum(exp(eta_t) * dh)
      
      exp(-Hhat)
    }
    newdat = as.data.frame(data)
    S <- mapply(
      pred_surv_one,
      x1 = newdat$X1,
      x2 = newdat$X2,
      y  = newdat$Y,
      MoreArgs = list(fit = fit_cp, bh = bh, tau = tau)
    )
  } else {
    train = as.data.frame(datas)
    colnames(train) = c("ID", "X1", "X2", "time", "status")
    event_times <- sort(unique(train$time[train$status == 1]))
    cuts <- event_times[event_times < max(train$time)]
    
    train_long <- survSplit(
      Surv(time, status) ~ .,
      data  = train,
      cut   = cuts,
      start = "tstart",
      end   = "tstop",
      event = "status_split"
    )
    
    train_long$Z <- with(train_long, X1 * ia(tstop))
    
    fit_cp <- coxph(
      Surv(tstart, tstop, status_split) ~ X1 + X2 + Z,
      data = train_long,
      ties = "breslow",
      x = TRUE
    )
    
    bh <- basehaz(fit_cp, centered = FALSE)
    bh$dh <- c(bh$hazard[1], diff(bh$hazard))
    
    pred_surv_one <- function(x1, x2, y, fit, bh) {
      cf <- coef(fit)
      b1 <- cf["X1"]
      b2 <- cf["X2"]
      g  <- cf["Z"]
      
      idx <- bh$time <= y
      if (!any(idx)) return(1)
      
      tt <- bh$time[idx]
      dh <- bh$dh[idx]
      
      eta_t <- b1 * x1 + b2 * x2 + g * x1 * ia(tt)
      Hhat <- sum(exp(eta_t) * dh)
      
      exp(-Hhat)
    }
    newdat = as.data.frame(data)
    S = mapply(
      pred_surv_one,
      x1 = newdat$X1,
      x2 = newdat$X2,
      y  = newdat$Y,
      MoreArgs = list(fit = fit_cp, bh = bh)
    )
  }
  
  return(list("S" = S,
              "datas" = datas,
              "rawXs" = rawXs))
}

sum_metric = function(res, datav, rawXv) {
  # C-index
  score = datav[,c("X1", "X2")] %*% res$beta
  Cindex = cIndex(time = datav[,"Y"],
                  event = datav[,"Delta"],
                  risk_score = score)
  Cindex = Cindex["index"]
  attributes(Cindex) = NULL
  
  res = c("L2distS" = res$L2distS,
          "Stau" = res$Stau,
          "Cindex" = Cindex,
          "intBS" = res$intBS,
          "RMST" = res$RMST)
  
  return(res)
}

trans_fit = function(data, datav, rawX, rawXv, S, weight, r_opt, xi, kmcensor) {
  # transfer learning algorithm
  weight = weight/sum(weight)
  res = TransFit(data[,"Y"], data[,"Delta"], rawX, S, weight, 
                 datav[,"Y"], datav[,"Delta"], rawXv, kmcensor,
                 r_opt, N, xi, tau)
  
  res = sum_metric(res, datav, rawXv)
  
  return(res)
}

# choose the value of xi by 5-fold cross validation
cross_valid = function(datat, rawXt, all_S, r_opt, all_xi, kmcensor, single = 1) {
  
  if(single) {
    # for single or pooled source predictor
    S = all_S
    xi = all_xi
  } else {
    # for aggregating multiple source predictors
    c = all_xi/sum(all_xi)
    S = all_S %*% c
    xi = sum(all_xi)
  }
  
  # no weighting 
  weight = rep(1.0, length(S))
  
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

run_all_meth = function(rep) {
  set.seed(rep)
  data = t(sapply(1:n, sim_data))
  colnames(data) = c("ID", "X1", "X2", "Y", "Delta")
  datav = t(sapply(1:nv, sim_data, censor = FALSE))
  colnames(datav) = c("ID", "X1", "X2", "Y", "Delta")
  
  # sort data by Y
  ind = order(data[,"Y"])
  data = data[ind,]
  rawX = cbind(0:(n-1),
               rep(0,n),
               rep(tau,n),
               data[,c("X1", "X2")])
  
  # KM estimator for censoring distribution
  # no censoring for validation data
  kmcensor = matrix(nrow = 0, ncol = 2)
  
  # sort validation data by Y
  ind = order(datav[,"Y"])
  datav = datav[ind,]
  rawXv = cbind(0:(nv-1),
                rep(0,nv),
                rep(tau,nv),
                datav[,c("X1", "X2")])
  
  # obtain source prediction and weights
  fit_source = lapply(1:nsc, get_source_predictor, data = data, rawX = rawX, kmcensor = kmcensor)
  all_S = lapply(fit_source, function(x) x$S)
  all_S = do.call(cbind, all_S)
  
  if(pool=="cv") {
    # tune xi_k for each source study altogether
    all_xi_cand <- t(expand.grid(xi1 = xi_cand_K,
                                 xi2 = xi_cand_K,
                                 xi3 = xi_cand_K,
                                 xi4 = xi_cand_K,
                                 xi5 = xi_cand_K))
    # each column of all_xi_cand is a set of xi_k (k=1,...,5)
    cv_loss = tryCatch({
      sapply(1:ncol(all_xi_cand), function(i) cross_valid(data, rawX, all_S, r, all_xi_cand[,i], kmcensor, single = 0))
    }, error = function(e) {
      message("rep ", rep, ": error in CV for selecting xi")
      e
    })
    all_xi = all_xi_cand[,which.min(cv_loss)]
    message("all xi based on 5-fold cross validation:")
    print(all_xi)
    
    # optimal source weights
    c = all_xi/sum(all_xi)
  } else if(pool=="ns") {
    # weight each source by sample size
    c = ns/sum(ns)
  }
  
  S = all_S %*% c
  weight = rep(1.0, length(S))
  
  # tune xi with optimal source predictor
  cv_loss = tryCatch({
    sapply(xi_cand, function(x) cross_valid(data, rawX, S, r, x, kmcensor))
  }, error = function(e) {
    message("rep ", rep, ": error in CV for selecting xi")
    e
  })
  xi = xi_cand[which.min(cv_loss)]
  message("xi = ", xi, " based on 5-fold cross validation.")
  
  # fit model with optimal xi
  res = trans_fit(data, datav, rawX, rawXv, S, weight, r, xi, kmcensor)
  
  ##############################################################################
  # compare existing methods
  ## Cox model with target data only
  res_tar = trans_fit(data, datav, rawX, rawXv, S, weight, r, xi = 0.0, kmcensor)
  
  # ## Cox model with all target and source data
  # all_datas = lapply(fit_source, function(x) x$datas)
  # datas = do.call(rbind, all_datas)
  # all_rawXs = lapply(fit_source, function(x) x$rawXs)
  # rawXs = do.call(rbind, all_rawXs)
  # 
  # datac = rbind(data, datas)
  # rawXc = rbind(rawX, rawXs)
  # # sort source data by Y
  # ind = order(datac[,"Y"])
  # datac = datac[ind,]
  # rawXc = rawXc[ind,]
  # rawXc[,1] = 0:(n+sum(ns)-1)
  # rawXc[,3] = max(taus)
  # 
  # ## select the optimal r for combined data
  # aic_com = sapply(r_cand, function(r_try) {
  #   fit_com = sourceFit(datac[,"Y"], datac[,"Delta"], rawXc,
  #                       data[,"Y"], rawX, r_try, N, verbose = 0)
  #   -2*fit_com$logL
  # })
  # rc_opt = r_cand[which.min(aic_com)]
  # message("optimal rc based on AIC: ", rc_opt)
  # # rc_opt = r
  # res_com = trans_fit(datac, datav, rawXc, rawXv, rep(1,n+sum(ns)), rep(1,n+sum(ns)),
  #                     rc_opt, xi = 0.0, kmcensor)
  # 
  # # CoxTL 
  # p = 2
  # data_t = data[,-1]
  # data_s = lapply(all_datas, function(x) {y = x[,-1]; colnames(y) = c("X1", "X2", "time", "status"); y})
  # colnames(data_t) = c("X1", "X2", "time", "status")
  # res_coxtl = rep(NA, 5)
  # tryCatch({
  #   cox_tl<-run_CoxTL_ms(data_t, data_s, p)
  #   beta_est = cox_tl$coefficients
  #   cumlam = basehaz(cox_tl, centered = F)
  #   uniqt = cumlam$time
  #   lambda_est = cumlam$hazard-c(0.0, cumlam$hazard[-nrow(cumlam)])
  #   ind = which(lambda_est!=0 & uniqt<tau)
  #   uniqt = uniqt[ind]
  #   lambda_est = lambda_est[ind]
  #   res_coxtl = Metric(beta_est, lambda_est, uniqt, datav[,"Y"], datav[,"Delta"],
  #                      rawXv, kmcensor, r, tau)
  #   res_coxtl = sum_metric(res_coxtl, datav, rawXv)
  # }, error = function(e) {
  #   # message("CoxTL: failed Cox regression due to extreme weights")
  #   message("CoxTL failed: ", conditionMessage(e))
  # })
  ##############################################################################
  res_coxtl = res_com = rep(NA, 5)
  sumres = rbind(res, res_tar, res_coxtl, res_com)
  # 1: POTL, 2: target, 3: CoxTL, 4: pooled
  sumres = cbind(1:4, sumres)
  
  return(sumres)
}

# simulation in single scenario
nrep = 200
# ncore = detectCores()/2
ncore = 200
is_parallel = 1

if(is_parallel) {
  allres = mclapply(1:nrep, run_all_meth, mc.cores = ncore)
} else {
  allres = sapply(1:nrep, run_all_meth)
}

save(allres, file = sprintf("res_multi_source_%s.RData", pool))
allres <- as.data.frame(do.call(rbind, allres))
allres = apply(allres, c(1,2), as.numeric) 
allres = na.omit(allres)

# aggregate results over all replicates
nmeth = 4
me = sapply(1:nmeth, function(i) apply(allres[allres[,1]==i,-1], 2, median, na.rm=T))
# me[1:2,] = me[1:2,]-tval
se = sapply(1:nmeth, function(i) apply(allres[allres[,1]==i,-1], 2, mad, na.rm=T))
out = sapply(1:nmeth, function(i) sprintf("%.3f (%.3f)", me[,i], se[,i]))

metric = c("L2D", "Dtau", "C-index", "IBS", "RMST")
out = cbind(metric, out)
colnames(out) = c("Metric", "POTL", "Target-only", "CoxTL", "Pooled")
print(out, quote = F)
write.table(out, file = sprintf("res_%s.txt", pool), 
            quote = F, row.names = F, col.names = T)
library(xtable)
print(xtable(out), include.rownames = F)
