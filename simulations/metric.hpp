#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double trueS(double u, // time
             const vec& x,
             double r) {
  double Lambda_u, true_surv;
  rowvec beta = {0.5, -0.5};
  
  Lambda_u = log(1+0.5*u);
  Lambda_u *= exp(as_scalar(beta*x));
  true_surv = StX(Lambda_u, r);
  return true_surv;
}

double diffS2(double u,
              const rowvec& pred,
              const vec& t,
              const vec& x,
              double r,
              int L) {
  // compute the square of Shat(u)-S(u), 
  double diff2, Shat;
  
  Shat = 1.0;
  for(int l=L-1; l>=0; --l) {
    if(u >= t(l)) {
      Shat = pred(l);
      break;
    }
  }
  diff2 = pow(Shat-trueS(u, x, r), 2.0);
  return diff2;
}

double diffS2_hat(double u,
                  const rowvec& pred1,
                  const rowvec& pred2,
                  const vec& t1,
                  const vec& t2) {
  // compute {S1(u|X)-S2(u|X)}^2
  double diff, Shat1, Shat2;
  int L1, L2, l1, l2, flag = 0;
  
  L1 = t1.n_elem;
  L2 = t2.n_elem;
  
  Shat1 = Shat2 = 1.0;
  for(l1=L1-1, l2=L2-1; l1>=0, l2>=0; --l1, --l2) {
    if(u >= t1(l1)) {
      Shat1 = pred1(l1);
      flag++;
    }
    if(u >= t2(l2)) {
      Shat2 = pred2(l2);
      flag++;
    }
    if(flag==2) {
      break;
    }
  }
  diff = pow(Shat1-Shat2, 2.0);
  return diff;
}

double trap(double a, 
            double b, 
            const rowvec& pred,
            const vec& t,
            const vec& x,
            double r,
            int L,
            int m = 100) {
  // m is number of grids. Higher value means more accuracy
  double h = (b-a)/(m-1);
  // Evaluate endpoints
  double value = 0.5*(diffS2(a,pred,t,x,r,L) + diffS2(b,pred,t,x,r,L));
  // Now the midpoints
  for(int k=2; k < m; k++){
    value += diffS2(a + h*(k-1),pred,t,x,r,L);
  }
  value*=h;
  return value;
}

double trap_hat(double a,
                double b,
                const rowvec& pred1,
                const rowvec& pred2,
                const vec& t1,
                const vec& t2, 
                int m = 100) {
  // m is number of grids. Higher value means more accuracy
  double h = (b-a)/(m-1);
  // Evaluate endpoints
  double value = 0.5*(diffS2_hat(a,pred1,pred2,t1,t2) + diffS2_hat(b,pred1,pred2,t1,t2));
  // Now the midpoints
  for(int k=2; k < m; k++){
    value += diffS2_hat(a + h*(k-1),pred1,pred2,t1,t2);
  }
  value*=h;
  return value;
}

double Stau(const vec& t,
            const mat& pred, // output from PredSurv
            const field<vec>& Xv,
            double r,
            int nv,
            int L,
            double endt) {
  // compute difference between S(tau) and predicated S(tau)
  double error;
  
  error = 0.0;
  for(int i=0; i<nv; ++i) {
    error += sqrt(diffS2(endt, pred.row(i), t, Xv(i,0), r, L));
  }
  error/=nv;
  return error;
}

double supdistS(const vec& t,
                const mat& pred, // output from PredSurv
                const field<vec>& Xv,
                double r,
                int nv,
                int L,
                double endt) {
  // compute the supremum distance between the predicted and true survival function
  // over [0,tau]
  double error, dist1, dist2, supdist, curS;
  
  error = 0.0;
  for(int i=0; i<nv; ++i) {
    curS = trueS(t(0), Xv(i,0), r);
    supdist = fabs(1-curS);
    for(int l=0; l<L-1; ++l) {
      // compute dist1 = |Shat(t_l)-S(t_l)|
      dist1 = fabs(pred(i,l)-curS);
      if(dist1 > supdist) {
        supdist = dist1;
      }
      // compute dist2 = |Shat(t_l)-S(t_{l+1})|
      curS = trueS(t(l+1), Xv(i,l+1), r);
      dist2 = fabs(pred(i,l)-curS);
      if(dist2 > supdist) {
        supdist = dist2;
      }
    } // end l
    error += supdist; 
  } // end i 
  error/=nv;
  return error;
}

double L2distS(const vec& t,
               const mat& pred, // output from PredSurv
               const field<vec>& Xv,
               double r,
               int nv,
               int L,
               double endt) {
  // compute L2 distance between the predicted and true survival function
  double error, int_diffS2;
  
  error = 0.0;
  for(int i=0; i<nv; ++i) {
    // integral of diffS2 from 0 to tau
    int_diffS2 = trap(0.0, endt, pred.row(i), t, Xv(i,0), r, L);
    error += sqrt(int_diffS2);
  }
  error/=nv;
  return error;
}

double L2distS_2r(const vec& t1,
                   const vec& t2,
                   const mat& pred1,
                   const mat& pred2,
                   double endt,
                   double rho = 1.0) {
  // compute the 2rth power of L2 distance between Shat1 and Shat2 based on generated data
  double diff, tmp;
  int np;
  
  np = pred1.n_rows;
  if(np != pred2.n_rows) {
    Rcout << "two predmat do not contain the same number of subjects." << endl;
  }
  
  diff = 0.0;
  for(int i=0; i<np; ++i) {
    tmp = trap_hat(0.0, endt, pred1.row(i), pred2.row(i), t1, t2);
    diff += pow(tmp, rho);
  }
  diff /= np;
  return diff;
}

double Brier(double u, // time point
             const vec& t,
             const vec& Yv, // validation data
             const ivec& Deltav,
             const mat& pred, // output from PredSurv
             const mat& kmcensor,
             int nv) {
  // compute Brier score at time point u
  // kmcenor is the KM estimator for the censoring distribution.
  // kmcensor is a matrix with two columns, the first column
  // represents time, the second column represents KM estimator,
  // the first row is (0,1)
  double brier, Ghat;
  int Lc, lc; // number of distinct censoring times
  uvec ind;
  vec Shat;
  
  Lc = kmcensor.n_rows;
  if(t(0)<=u) {
    ind = find(t<=u);
    Shat = pred.col(ind(ind.n_elem-1));
  } else {
    Shat = ones<vec>(nv);
  }
  
  brier = 0.0;
  lc = 0;
  for(int i=0; i<nv; ++i) {
    if(Lc>0) {
      while(lc<Lc && Yv(i)>=kmcensor(lc, 0)) {
        lc++;
      }
      lc--;
      Ghat = kmcensor(lc, 1);
    } else {
      Ghat = 1.0;
    }

    if(Yv(i)<=u && Deltav(i)==1) {
      brier += pow(Shat(i), 2.0)/Ghat;
    }
    if(Yv(i)>u) {
      brier += pow(1.0-Shat(i), 2.0)/Ghat;
    }
  } // end i
  brier /= nv;
  return brier;
}

double intBrier(const vec& t,
                const vec& Yv, // validation data
                const ivec& Deltav,
                const mat& pred, // output from PredSurv
                const mat& kmcensor,
                int nv,
                double endt,
                int m = 100) {
  // compute integrated brier score over [0, endt] using trap method
  double h, value;
  
  h = endt/(m-1);
  value = 0.5*(Brier(0.0,t,Yv,Deltav,pred,kmcensor,nv) + Brier(endt,t,Yv,Deltav,pred,kmcensor,nv));
  // Now the midpoints
  for(int k=2; k < m; k++){
    value += Brier(h*(k-1),t,Yv,Deltav,pred,kmcensor,nv);
  }
  value *= h;
  value /=  endt;
  return value;
}

vec medT(const vec& t,
         const mat& pred, // output from PredSurv
         int nv,
         int L,
         double endt) {
  // compute median survival time
  vec med;
  
  med.set_size(nv);
  med.fill(endt);
  for(int i=0; i<nv; ++i) {
    for(int l=0; l<L; ++l) {
      if(pred(i,l)<=0.5) {
        med(i) = t(l);
        break;
      }
    }
  }
  return med;
}

vec RMST(const vec& t,
         const mat& pred, // output from PredSurv
         int nv,
         int L,
         double endt) {
  // compute restricted mean survival time
  vec rmst, dt, newt;
  mat newpred;
  
  // add time zero
  newt.zeros(L+2);
  newt(span(1,L)) = t.head(L);
  newt(L+1) = endt;
  newpred = join_horiz(ones<vec>(nv), pred);
  
  rmst.zeros(nv);
  dt = newt.tail(L+1)-newt.head(L+1);
  for(int i=0; i<nv; ++i) {
    for(int l=0; l<L+1; ++l) {
      rmst(i) += newpred(i,l)*dt(l);
    }
  }
  return rmst;
}

double diffT(const vec& That,
             const vec& Yv, // validation data
             const ivec& Deltav,
             const mat& kmcensor,
             int nv,
             double endt) {
  double diff, den, Ghat;
  int Lc, lc;
  
  Lc = kmcensor.n_rows;
  diff = den = 0.0;
  lc = 0;
  for(int i=0; i<nv; ++i) {
    if(Lc>0) {
      while(lc<Lc && Yv(i)>=kmcensor(lc, 0)) {
        lc++;
      }
      lc--;
      Ghat = kmcensor(lc, 1);
    } else {
      Ghat = 1.0;
    }
    
    if(Yv(i)<=endt && Deltav(i)==1) {
      diff += fabs(Yv(i)-That(i))/Ghat;
      den += 1.0/Ghat;
    }
  } // end i
  diff /= den;
  return diff;
}
