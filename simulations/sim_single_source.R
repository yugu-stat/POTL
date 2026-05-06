Rcpp::sourceCpp("EM.cpp")
source("coxTL.R")

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

# target sample size
n = 100
# validation sample size
nv = 10000
# source sample size
ns = 1000
# study end time
tau = 2
taus = 5

# whether there's covariate shift in source data
shift_s = FALSE
# whether there's covariate shift in validation data
shift_v = FALSE

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

# varying coefficient for time-covariate interaction
ia = function(t) {
  t
}

# Gau-Lag approx
N = 20
# candidate tuning parameter xi
# xi_cand = seq(0, 150, 5)
xi_cand = c(0.0, 2^seq(-15,15,1))
# candidate value of r
r_cand = rs_cand = seq(0, 1.5, 0.05)

# inverse transformation function
G_inv = function(g, r) {
  if(abs(r)<1e-8) {
    G_inv_g = g
  } else {
    G_inv_g = (exp(g*r)-1)/r
  }
  return(G_inv_g)
}

# simulate data from transformation model with r
sim_data = function(i, censor = TRUE, shift = FALSE) {
  if(shift) {
    # covariate distribution shift for validation data 
    X1 = ifelse(runif(1)>0.5, 1, 0)
    X2 = rbeta(1, 1, 2)
  } else {
    # target sample
    X1 = ifelse(runif(1)>0.5, 1, 0)
    X2 = runif(1)
  }
  
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

sim_source_data = function(i, sc, r, nu, shift = FALSE) {
  X1 = ifelse(runif(1)>0.5, 1, 0)
  if(!shift) {
    X2 = runif(1)
  } else {
    X2 = rbeta(1, 1, 2) 
  }
  
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

run_all_meth = function(rep, sc) {
  set.seed(rep)
  data = t(sapply(1:n, sim_data))
  colnames(data) = c("ID", "X1", "X2", "Y", "Delta")
  datav = t(sapply(1:nv, sim_data, censor = FALSE, shift = shift_v))
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
  
  # simulate source data
  if(sc %in% c(1,4:5)) {
    nus = 0.5
  } else {
    # sc2/3: source model = target model with different parameters
    nus = 0.4
  }
  
  datas = t(sapply(1:ns, sim_source_data, sc=sc, r=rs, nu=nus, shift=shift_s))
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
  # xi = 100.0
  # fit model with optimal xi
  res = trans_fit(data, datav, rawX, rawXv, S, weight, r, xi, kmcensor)
  
  ##############################################################################
  # compare existing methods
  ## Cox model with target data only
  res_tar = trans_fit(data, datav, rawX, rawXv, S, weight, r, xi = 0.0, kmcensor)
  
  ## Cox model with both target and source data
  datac = rbind(data, datas)
  rawXc = rbind(rawX, rawXs)
  # sort source data by Y
  ind = order(datac[,"Y"])
  datac = datac[ind,]
  rawXc = rawXc[ind,]
  rawXc[,1] = 0:(n+ns-1)
  rawXc[,3] = 5
  res_com = trans_fit(datac, datav, rawXc, rawXv, rep(1,n+ns), rep(1,n+ns),
                      r, xi = 0.0, kmcensor)
  
  ## Li et al (2023, JASA)
  # select the best learning rate and number of steps using BIC
  tmpdata = data.frame(data[,2:5]) %>%
    mutate(time = Y,
           status = Delta+1) %>%
    select(-Y, -Delta)
  tmpdatas = data.frame(datas[,2:5]) %>%
    mutate(time = Y,
           status = Delta+1) %>%
    select(-Y, -Delta)
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
  res_li23 = sum_metric(res_li23, datav, rawXv)
  
  ## CoxTL (Lu et al.)
  p = 2
  data_t = data[,-1]
  data_s = datas[,-1]
  colnames(data_t) = colnames(data_s) = c("X1", "X2", "time", "status")
  res_coxtl = rep(NA, 5)
  tryCatch({
    cox_tl<-run_CoxTL(data_t,data_s, p)
    beta_est = cox_tl$coefficients
    cumlam = basehaz(cox_tl, centered = F)
    uniqt = cumlam$time
    lambda_est = cumlam$hazard-c(0.0, cumlam$hazard[-nrow(cumlam)])
    ind = which(lambda_est!=0 & uniqt<tau)
    uniqt = uniqt[ind]
    lambda_est = lambda_est[ind]
    res_coxtl = Metric(beta_est, lambda_est, uniqt, datav[,"Y"], datav[,"Delta"],
                       rawXv, kmcensor, r, tau)
    res_coxtl = sum_metric(res_coxtl, datav, rawXv)
  }, error = function(e) {
    message("CoxTL: failed Cox regression due to extreme weights")
  })
  
  ##############################################################################
  
  sumres = rbind(res, res_tar, res_li23, res_coxtl, res_com)
  # 1: POTL, 2: Target-only, 3: TransCox, 4: CoxTL, 5: Pooled
  sumres = cbind(1:5, sumres)
  
  return(sumres)
}

# simulation in single scenario
single_sc_sim = function(sc) {
  nrep = 200
  ncore = 200
  
  allres = mclapply(1:nrep, run_all_meth, sc=sc, mc.cores = ncore)
  save(allres, file = sprintf("simres_sc%d.RData", sc))
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
