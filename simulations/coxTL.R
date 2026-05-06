# 9/30/2025: tune both lambda and nu using candidate set 10^(seq(-5,5,1))

# optimal density
density_opt<-function(Xs,Xt,lam){
  mm<-function(theta_0){
    theta_0<-as.matrix(theta_0,ncol=1)
    ff<-mean(Xt%*%theta_0)-mean(exp(Xs%*%theta_0))+lam*norm(theta_0,type = "2")^2
    return(ff)
  }
  return(mm)
}

## fit CoxTL proposed in Lu et al. (2025+)
run_CoxTL<-function(data_t,data_s,p,folds_num=5,
                    lam_set = 10^(seq(-4,1,1)),
                    nu_set=10^(seq(-4,1,1))){
  data_t = as.data.frame(data_t)
  data_s = as.data.frame(data_s)
  num_t<-dim(data_t)[1]
  num_s<-dim(data_s)[1]

  cv_record = matrix(NA, nrow = length(lam_set), ncol = length(nu_set))
  for(lam_idx in 1:length(lam_set)) {
    lambda = lam_set[lam_idx]
    for(nu_idx in 1:length(nu_set)) {
      nu = nu_set[nu_idx]
      # cross validation
      split_set <- sample(rep(1:folds_num, length.out = num_t))
      cindex = rep(NA,folds_num)
      for(k in 1:folds_num){
        data_t1<-data_t[split_set==k,] # validation
        data_t2<-data_t[split_set!=k,] # train
        data_ts<-rbind(data_s,data_t2)
        
        # obtain weights
        theta_ini<-rep(0,p+1) # set initial value
        Xs_new<-as.matrix(cbind(rep(1,num_s),data_s[,1:p])) # add intercept term
        Xt_new<-as.matrix(cbind(rep(1,nrow(data_t2)),data_t2[,1:p])) # add intercept term
        theta = rep(NA, p+1)
        tryCatch({
          w_opt <- nlm(density_opt(Xs_new,Xt_new,lambda), theta_ini, iterlim = 500)
          theta<-w_opt$estimate
        }, error = function(e) {
          message("lambda = ", lambda, ": error in nlm")
        })
        weights<-exp(Xs_new%*%matrix(theta,ncol=1))
        
        # fit weighted Cox regression
        tryCatch({
          cox_train<- rms::cph(survival::Surv(time, status) ~ X1+X2, 
                               data = data_ts,weights=c(nu*weights,rep(1,dim(data_t2)[1])),
                               method = "breslow", x=TRUE,y=TRUE,surv=TRUE)
          cindex[k]<-1-rcorr.cens(exp(as.matrix(data_t1[,1:p])%*%coef(cox_train)),Surv(data_t1$time,data_t1$status))['C Index']
        }, error = function(e) {
          message("lambda = ", lambda, ": results in extreme weights and failed Cox regression")
        })
      }
      cv_record[lam_idx, nu_idx] = mean(cindex)
    }
  }

  ind = arrayInd(which.max(cv_record), .dim=dim(cv_record))
  lam_best<-lam_set[ind[1]]
  nu_best = nu_set[ind[2]]
  message("optimal lambda: ", lam_best)
  message("optimal nu: ", nu_best)
  
  # obtain weights
  theta_ini<-rep(0,p+1) # set initial value
  Xs_new<-as.matrix(cbind(rep(1,num_s),data_s[,1:p])) # add intercept term
  Xt_new<-as.matrix(cbind(rep(1,num_t),data_t[,1:p])) # add intercept term
  w_opt <- nlm(density_opt(Xs_new,Xt_new,lam_best), theta_ini, iterlim = 500)
  theta<-w_opt$estimate
  weights<-exp(Xs_new%*%matrix(theta,ncol=1))
  cox_final<- rms::cph(survival::Surv(time, status) ~ X1+X2, 
                       data=rbind(data_s,data_t),
                       weights=c(nu_best*weights,rep(1,num_t)),
                       method = "breslow",x=TRUE,y=TRUE,surv=TRUE)
  return(cox_final)
}

## fit CoxTL-multi-source proposed in Lu et al. (2025+)
run_CoxTL_ss<-function(sc, data_t,data_s_list,p,folds_num=5,
                       lam_set = 10^(seq(-4,1,1)),
                       nu_set=10^(seq(-4,1,1))){
  data_t = as.data.frame(data_t)
  data_s = as.data.frame(data_s_list[[sc]])
  num_t<-dim(data_t)[1]
  num_s<-dim(data_s)[1]
  
  cv_record = matrix(NA, nrow = length(lam_set), ncol = length(nu_set))
  for(lam_idx in 1:length(lam_set)) {
    lambda = lam_set[lam_idx]
    for(nu_idx in 1:length(nu_set)) {
      nu = nu_set[nu_idx]
      # cross validation
      split_set <- sample(rep(1:folds_num, length.out = num_t))
      cindex = rep(NA,folds_num)
      for(k in 1:folds_num){
        data_t1<-data_t[split_set==k,] # validation
        data_t2<-data_t[split_set!=k,] # train
        data_ts<-rbind(data_s,data_t2)
        
        # obtain weights
        theta_ini<-rep(0,p+1) # set initial value
        Xs_new<-as.matrix(cbind(rep(1,num_s),data_s[,1:p])) # add intercept term
        Xt_new<-as.matrix(cbind(rep(1,nrow(data_t2)),data_t2[,1:p])) # add intercept term
        theta = rep(NA, p+1)
        tryCatch({
          w_opt <- nlm(density_opt(Xs_new,Xt_new,lambda), theta_ini, iterlim = 500)
          theta<-w_opt$estimate
        }, error = function(e) {
          message("lambda = ", lambda, ": error in nlm")
        })
        weights<-exp(Xs_new%*%matrix(theta,ncol=1))
        
        # fit weighted Cox regression
        tryCatch({
          cox_train<- rms::cph(survival::Surv(time, status) ~ X1+X2, 
                               data = data_ts,weights=c(nu*weights,rep(1,dim(data_t2)[1])),
                               method = "breslow", x=TRUE,y=TRUE,surv=TRUE)
          cindex[k]<-1-rcorr.cens(exp(as.matrix(data_t1[,1:p])%*%coef(cox_train)),Surv(data_t1$time,data_t1$status))['C Index']
        }, error = function(e) {
          message("lambda = ", lambda, ": results in extreme weights and failed Cox regression")
        })
      }
      cv_record[lam_idx, nu_idx] = mean(cindex)
    }
  }
  
  ind = arrayInd(which.max(cv_record), .dim=dim(cv_record))
  lam_best<-lam_set[ind[1]]
  nu_best = nu_set[ind[2]]
  message("optimal lambda: ", lam_best)
  message("optimal nu: ", nu_best)
  
  # obtain weights
  if(!is.na(lam_best)) {
    theta_ini<-rep(0,p+1) # set initial value
    Xs_new<-as.matrix(cbind(rep(1,num_s),data_s[,1:p])) # add intercept term
    Xt_new<-as.matrix(cbind(rep(1,num_t),data_t[,1:p])) # add intercept term
    w_opt <- nlm(density_opt(Xs_new,Xt_new,lam_best), theta_ini, iterlim = 500)
    theta<-w_opt$estimate
    weights<-exp(Xs_new%*%matrix(theta,ncol=1))
  } else {
    weights = rep(1.0, ns[sc])
  }
  
  if(is.na(nu_best)) {
    # ratio in sample size
    nu_best = ns[sc]/n
  }

  return(list(weights = weights,
              nu = nu_best))
}

run_CoxTL_ms<-function(data_t,data_s_list,p,folds_num=5,
                    lam_set = 10^(seq(-4,1,1)),
                    nu_set=10^(seq(-4,1,1))){
  all_sc_res = lapply(1:nsc, run_CoxTL_ss, data_t = data_t, data_s_list = data_s_list, p = p)
  all_weights = lapply(all_sc_res, function(x) x$weights)
  all_nu = lapply(all_sc_res, function(x) x$nu)
  
  data_t = as.data.frame(data_t)
  data_s = as.data.frame(do.call(rbind, data_s_list)) 
  num_t<-dim(data_t)[1]
  num_s<-dim(data_s)[1]
  
  for(sc in 1:nsc) {
    # tune nu for source sc while fixing nu at initial values for all other sources 
    cv_record = rep(NA, nrow = length(nu_set))
    for(nu_idx in 1:length(nu_set)) {
      nu_sc = nu_set[nu_idx]
      all_nu[[sc]] = nu_sc
      
      # cross validation
      split_set <- sample(rep(1:folds_num, length.out = num_t))
      cindex = rep(NA,folds_num)
      for(k in 1:folds_num){
        data_t1<-data_t[split_set==k,] # validation
        data_t2<-data_t[split_set!=k,] # train
        data_ts<-rbind(data_s,data_t2)
        cox_weights = unlist(lapply(1:nsc, function(x) all_nu[[x]]*all_weights[[x]]))
        
        # remove subjects with 0 weights
        rm_ind = which(cox_weights==0.0)
        if(length(rm_ind)) {
          data_ts = data_ts[-rm_ind,]
          cox_weights = cox_weights[-rm_ind]
        }
        
        # fit weighted Cox regression
        tryCatch({
          cox_train<- rms::cph(survival::Surv(time, status) ~ X1+X2, 
                               data = data_ts,weights=c(cox_weights,rep(1,dim(data_t2)[1])),
                               method = "breslow", x=TRUE,y=TRUE,surv=TRUE)
          cindex[k]<-1-rcorr.cens(exp(as.matrix(data_t1[,1:p])%*%coef(cox_train)),Surv(data_t1$time,data_t1$status))['C Index']
        }, error = function(e) {
          message("nu = ", nu_sc, ": weighted Cox regression failed")
        })
      }
      cv_record[nu_idx] = mean(cindex, na.rm = T)
    }
    ind = which.max(cv_record)
    if(length(ind)) {
      nu_sc = nu_set[ind]
      all_nu[[sc]] = nu_sc
    } else {
      # ratio in sample size
      all_nu[[sc]] = ns[sc]/n
    }
  }
  
  message("best nu: ")
  print(unlist(all_nu))
  cox_weights = unlist(lapply(1:nsc, function(x) all_nu[[x]]*all_weights[[x]]))
 
  # remove subjects with 0 weights
  rm_ind = which(cox_weights==0.0)
  if(length(rm_ind)) {
    data_s = data_s[-rm_ind,]
    cox_weights = cox_weights[-rm_ind]
  }
  cox_final<- rms::cph(survival::Surv(time, status) ~ X1+X2, 
                       data=rbind(data_s,data_t),
                       weights=c(cox_weights,rep(1,num_t)),
                       method = "breslow",x=TRUE,y=TRUE,surv=TRUE)
  return(cox_final)
}

# ## fit CoxTL proposed in Lu et al. (2025+)
# run_CoxTL<-function(data_t,data_s,p,weights=NULL,folds_num=5,lam_set=c(0.001,0.005,0.01,0.05,0.1,1:10)){
#   data_t = as.data.frame(data_t)
#   data_s = as.data.frame(data_s)
#   num_t<-dim(data_t)[1]
#   num_s<-dim(data_s)[1]
#   # size ratio = 0.1
#   size_ratio<-num_t/num_s
#   if (is.null(weights)){
#     weights=rep(1,num_s)
#   }
#   cv_record = rep(0, length(lam_set))
#   for(lam_idx in 1:length(lam_set)) {
#     lambda<-size_ratio*lam_set[lam_idx]
#     for(b in 1:folds_num) {
#       split_set<-sample(1:folds_num,num_t,replace=TRUE)
#       cv_lam = rep(NA, folds_num)
#       for(k in 1:folds_num){
#         data_t1<-data_t[split_set==k,] # validation
#         data_t2<-data_t[split_set!=k,] # train
#         data_ts<-rbind(data_s,data_t2)
#         cox_train<- rms::cph(survival::Surv(time, status) ~ X1+X2,
#                              data = data_ts,weights=c(lambda*weights,rep(1,dim(data_t2)[1])),x=TRUE,y=TRUE,surv=TRUE)
#         cv_lam[k]<-1-rcorr.cens(exp(as.matrix(data_t1[,1:p])%*%coef(cox_train)),Surv(data_t1$time,data_t1$status))['C Index']
#       }
#       cv_record[lam_idx] = (cv_record[lam_idx]*(b-1)+mean(cv_lam))/b
#     }
#   }
#   lam_best<-size_ratio*lam_set[which.max(cv_record)]
#   message("optimal lambda: ", lam_best)
#   cox_final<- rms::cph(survival::Surv(time, status) ~ X1+X2,
#                        data=rbind(data_s,data_t),
#                        weights=c(lam_best*weights,rep(1,num_t)),x=TRUE,y=TRUE,surv=TRUE)
#   return(cox_final)
# }
