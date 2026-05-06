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
      split_set<-sample(1:folds_num,num_t,replace=TRUE)
      cindex = rep(NA,folds_num)
      for(k in 1:folds_num){
        data_t1<-data_t[split_set==k,] # validation
        data_t2<-data_t[split_set!=k,] # train
        data_ts<-rbind(data_s,data_t2)
        
        # # obtain weights
        # theta_ini<-rep(0,p+1) # set initial value
        # Xs_new<-as.matrix(cbind(rep(1,num_s),data_s[,2+(1:p)])) # add intercept term
        # Xt_new<-as.matrix(cbind(rep(1,nrow(data_t2)),data_t2[,2+(1:p)])) # add intercept term
        # theta = rep(NA, p+1)
        # tryCatch({
        #   w_opt <- nlm(density_opt(Xs_new,Xt_new,lambda), theta_ini, iterlim = 500)
        #   theta<-w_opt$estimate
        # }, error = function(e) {
        #   message("lambda = ", lambda, ": error in nlm")
        # })
        # weights<-exp(Xs_new%*%matrix(theta,ncol=1))
        
        # lambda<=1 always causes errors due to extreme weights; use equal weights instead
        weights = rep(1,num_s)
        
        # fit weighted Cox regression
        tryCatch({
          cox_train<- rms::cph(survival::Surv(time, status) ~ age+nlymph+stage34+IDC+ER+PR+HER2+
                                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
                               data = data_ts,weights=c(nu*weights,rep(1,dim(data_t2)[1])),
                               method = "breslow", x=TRUE,y=TRUE,surv=TRUE)
          cindex[k]<-1-rcorr.cens(exp(as.matrix(data_t1[,2+(1:p)])%*%coef(cox_train)),Surv(data_t1$time,data_t1$status))['C Index']
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
  # theta_ini<-rep(0,p+1) # set initial value
  # Xs_new<-as.matrix(cbind(rep(1,num_s),data_s[,2+(1:p)])) # add intercept term
  # Xt_new<-as.matrix(cbind(rep(1,num_t),data_t[,2+(1:p)])) # add intercept term
  # w_opt <- nlm(density_opt(Xs_new,Xt_new,lam_best), theta_ini, iterlim = 500)
  # theta<-w_opt$estimate
  # weights<-exp(Xs_new%*%matrix(theta,ncol=1))
  
  # lambda<=1 always causes errors due to extreme weights; use equal weights instead
  weights = rep(1,num_s)
  cox_final<- rms::cph(survival::Surv(time, status) ~ age+nlymph+stage34+IDC+ER+PR+HER2+
                         PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
                       data=rbind(data_s,data_t),
                       weights=c(nu_best*weights,rep(1,num_t)),
                       method = "breslow",x=TRUE,y=TRUE,surv=TRUE)
  return(cox_final)
}
