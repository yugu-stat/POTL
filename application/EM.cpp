#include <RcppArmadillo.h>
#include "gaulag.hpp"
#include "internal.hpp"
#include "metric.hpp"
#include <ctime>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// 5/2/2024 changes: 
// 1. weight penalty by inverse of variance of source prediction
// 2. compute variance estimator for S(Y) based on source data

// data should be sorted based on Y

// conditional expectation
struct conde{
  mat wbar;
  vec z;
};

conde Estep(const vec& beta,
            const vec& lambda,
            const vec& t, // unique time points
            const vec& Y,
            const ivec& Delta,
            const field<vec>& X, // covariates (subject, time)
            const vec& S, // predicted survival prob
            const vec& weight, // weight of S
            gaulag gl,
            double r, // variance of Gamma dist
            int N,
            int n,
            int L,
            double xi) {
  // E-step: compute conditional expectations of Wbar and z
  conde e;
  vec ru, ewz, numw;
  mat likeli;
  int l;
  double sum0, ebx, numz, den, tmp;
  
  e.wbar.zeros(n, L);
  e.z.ones(n);
  likeli.ones(n, N);
  ru = r*gl.node;
  
  if(fabs(r)>1e-8) {
    for(int i=0; i<n; ++i) {
      sum0 = den = numz = 0.0;
      ewz.zeros(L);
      numw.zeros(L);
      l = 0;
      while(l<L && t(l)<=Y(i)) {
        ebx = exp(as_scalar(beta.t()*X(i,l)));
        sum0 += lambda(l)*ebx;
        ewz(l) = (1.0-S(i))*lambda(l)*ebx;
        l++;
      }
      l--;
      if(l>=0) {
        for(int g=0; g<N; ++g) {
          if(Delta(i)==1) {
            likeli(i,g) *= ru(g)*lambda(l)*ebx;
          }
          likeli(i,g) *= exp(-(1.0+xi*weight(i)*S(i))*ru(g)*sum0);
          likeli(i,g) *= pow(1.0-exp(-ru(g)*sum0), xi*weight(i)*(1-S(i)));
          tmp = pow(ru(g), 1.0/r-1.0)*likeli(i,g);
          den += gl.weight(g)*tmp;
          numz += gl.weight(g)*ru(g)*tmp;
          if(fabs(sum0)>1e-8) {
            numw(span(0,l)) += gl.weight(g)*ewz(span(0,l))*ru(g)/(1.0-exp(-ru(g)*sum0))*tmp;
          }
        } // end g
        if(fabs(den) > 1e-8) {
          e.z(i) = numz/den;
          e.wbar(span(i), span(0,l)) = numw(span(0,l)).t()/den;
        }
      }
    } // end i
  } else {
    // r=0: Cox model: only need to compute EW
    for(int i=0; i<n; ++i) {
      sum0 = 0.0;
      l = 0;
      ewz.zeros(L);
      while(l<L && t(l)<=Y(i)) {
        ebx = exp(as_scalar(beta.t()*X(i,l)));
        sum0 += lambda(l)*ebx;
        ewz(l) = (1.0-S(i))*lambda(l)*ebx;
        l++;
      } 
      l--;
      if(l>=0 && fabs(sum0)>1e-8) {
        e.wbar(span(i), span(0,l)) = ewz(span(0,l)).t()/(1.0-exp(-sum0));
      }
    } // end i
  }
  
  return e;
} 

vec profile_Estep(const vec& beta,
                  const vec& lambda,
                  const vec& t, // unique time points
                  const vec& Y,
                  const ivec& Delta,
                  const field<vec>& X, // covariates (subject, time)
                  const vec& S, // predicted survival prob
                  gaulag gl,
                  double r, // variance of Gamma dist
                  int N,
                  int n,
                  int L) {
  // E-step for computing the profile likelihood. 
  // No penalty term, thus no need to compute EW.
  vec ru, ez;
  mat likeli;
  int l;
  double sum0, ebx, numz, den, tmp;
  
  ez.ones(n);
  likeli.ones(n, N);
  ru = r*gl.node;
  
  if(fabs(r)>1e-8) {
    for(int i=0; i<n; ++i) {
      sum0 = den = numz = 0.0;
      l = 0;
      while(l<L && t(l)<=Y(i)) {
        ebx = exp(as_scalar(beta.t()*X(i,l)));
        sum0 += lambda(l)*ebx;
        l++;
      }
      l--;
      
      if(l>=0) {
        for(int g=0; g<N; ++g) {
          if(Delta(i)==1) {
            likeli(i,g) *= ru(g)*lambda(l)*ebx;
          }
          likeli(i,g) *= exp(-ru(g)*sum0);
          tmp = pow(ru(g), 1.0/r-1.0)*likeli(i,g);
          den += gl.weight(g)*tmp;
          numz += gl.weight(g)*ru(g)*tmp;
        } // end g
        if(fabs(den) > 1e-8) {
          ez(i) = numz/den;
        }
      }
    } // end i
  }
  return ez;
}

bool Mstep(vec& beta,
           vec& lambda,
           const vec& t, // unique time points
           const vec& Y,
           const ivec& Delta,
           const field<vec>& X, // covariates (subject, time)
           const vec& weight, // weight of source prediction
           int n,
           int p,
           int L,
           conde e,
           double xi) {
  // M-step: update beta and lambda
  int i, nfail;
  double ebx;
  vec score, sum0, sumw;
  mat sum1, sumwx, info, info_inv;
  cube sum2;
  bool isinv;
  
  score.zeros(p);
  info.zeros(p,p);
  lambda.zeros(L);
  i = n-1;
  sum0.zeros(L);
  sum1.zeros(p, L);
  sum2.zeros(p, p, L);
  sumw.zeros(L);
  sumwx.zeros(p, L);
  for(int l=L-1; l>=0; --l) {
    nfail = 0;
    while(i>=0 && Y(i)>=t(l)) {
      for(int ll=0; ll<=l; ++ll) {
        ebx = exp(as_scalar(beta.t()*X(i,ll)));
        sum0(ll) += (1.0+xi*weight(i))*e.z(i)*ebx;
        sum1.col(ll) += (1.0+xi*weight(i))*e.z(i)*ebx*X(i,ll);
        sum2.slice(ll) += (1.0+xi*weight(i))*e.z(i)*ebx*X(i,ll)*X(i,ll).t();
        sumw(ll) += xi*weight(i)*e.wbar(i,ll);
        sumwx.col(ll) += xi*weight(i)*e.wbar(i,ll)*X(i,ll);
      }
      if(fabs(Y(i)-t(l)) < 1e-8 && Delta(i)==1) {
        nfail++;
        score += X(i,l);
      }
      i--;
    } // end i
    
    // update lambda
    score += sumwx.col(l);
    if(fabs(sum0(l)) > 1e-8) {
      lambda(l) = nfail + sumw(l);
      lambda(l) /= sum0(l);
      score -= nfail*sum1.col(l)/sum0(l);
      score -= sumw(l)*sum1.col(l)/sum0(l);
      info += nfail*(sum2.slice(l)/sum0(l)-sum1.col(l)*sum1.col(l).t()/pow(sum0(l), 2));
      info += sumw(l)*(sum2.slice(l)/sum0(l)-sum1.col(l)*sum1.col(l).t()/pow(sum0(l),2));
    }
  } // end l
  
  // update beta
  isinv = pinv(info_inv, info);
  if(isinv) {
    beta += info_inv*score;
  } 
  
  return isinv;
  
}

void profile_Mstep(vec& beta,
                   vec& lambda,
                   const vec& t, // unique time points
                   const vec& Y,
                   const ivec& Delta,
                   const field<vec>& X, // covariates (subject, time)
                   int n,
                   int L,
                   vec ez) {
  // M-step for computing the profile likelihood
  // keep beta fixed.
  int i;
  double nfail, ebx;
  vec sum0;
  
  lambda.zeros(L);
  i = n-1;
  sum0.zeros(L);
  for(int l=L-1; l>=0; --l) {
    nfail = 0.0;
    while(i>=0 && Y(i)>=t(l)) {
      for(int ll=0; ll<=l; ++ll) {
        ebx = exp(as_scalar(beta.t()*X(i,ll)));
        sum0(ll) += ez(i)*ebx;
      }
      if(fabs(Y(i)-t(l)) < 1e-8 && Delta(i)==1) {
        nfail = nfail + 1.0;
      }
      i--;
    } // end i
    
    // update lambda
    if(fabs(sum0(l)) > 1e-8) {
      lambda(l) = nfail/sum0(l);
    }
  } // end l
}

void EM(vec& beta,
        vec& lambda,
        const vec& t, // unique time points
        const vec& Y,
        const ivec& Delta,
        const field<vec>& X, // covariates (subject, time)
        const vec& S, // predicted survival prob
        const vec& weight, // weight of S
        gaulag gl,
        double r, // variance of Gamma dist
        int N,
        int n,
        int p,
        int L,
        double xi,
        double tol,
        int maxit,
        bool profile, // whether to profile lambda out
        bool verbose) {
  // EM algorithm
  double diff;
  vec beta0, lambda0, ez;
  int iter;
  conde e;
  bool isinv;
  
  iter = 0;
  diff = 1.0;
  while(iter<maxit && diff>tol) {
    iter++;
    beta0 = beta;
    lambda0 = lambda;
    if(!profile) {
      e = Estep(beta, lambda, t, Y, Delta, X, S, weight, gl, r, N, n, L, xi);
      isinv = Mstep(beta, lambda, t, Y, Delta, X, weight, n, p, L, e, xi);
      if(isinv) {
        diff = norm(beta-beta0, 2) + norm(lambda-lambda0, "inf");
        if(verbose && iter % 100 == 1) {
          time_t now = time(0);
          char* dt = ctime(&now);
          Rcout << dt << "Iteration " << iter << ": difference = " << diff << endl;
        }
      } else {
        break;
      }
    } else {
      // EM for profile likelihood
      ez = profile_Estep(beta, lambda, t, Y, Delta, X, S, gl, r, N, n, L);
      profile_Mstep(beta, lambda, t, Y, Delta, X, n, L, ez);
      diff = norm(lambda-lambda0, "inf");
    }
  }
  
  if(!profile) {
    if(isinv) {
      if(diff>tol) {
        Rcout << "Nonconvergence, diff = " << diff << "\n xi = " << xi << endl;
      } else if(verbose) {
        Rcout << "EM algorithm converged after " << iter << " iterations." << endl;
      }
    } else {
      Rcout << "EM algorithm failed due to singular information matrix!" << endl;
    }
  }
}

double logLikeli(const vec& beta,
                 const vec& lambda,
                 const vec& t, // unique time points
                 const vec& Y,
                 const ivec& Delta,
                 const field<vec>& X, // covariates (subject, time)
                 double r, // variance of Gamma dist
                 int n,
                 int L) {
  // compute the likelihood function without penalty
  double loglikeli, sum0, ebx;
  int l;
  vec Li;
  
  Li.ones(n);
  
  for(int i=0; i<n; ++i) {
    sum0 = 0.0;
    l = 0;
    while(l<L && t(l)<=Y(i)) {
      ebx = exp(as_scalar(beta.t()*X(i,l)));
      sum0 += lambda(l)*ebx;
      l++;
    } // end l
    l--;
    Li(i) *= StX(sum0, r);
    if(l>=0 && fabs(Y(i)-t(l)) < 1e-8 && Delta(i)==1) {
      Li(i) *= lambda(l)*ebx;
      Li(i) *= Gprime(sum0, r);
    }
  } // end i
  
  loglikeli = sum(log(Li(find(Li))));
  return loglikeli;
}

mat compCov(const vec& beta,
            const vec& lambda,
            conde e,
            const vec& t,
            const vec& Y,
            const ivec& Delta,
            const field<vec>& X, // covariates (subject, time)
            int n,
            int L,
            int p) {
  // compute covariance estimator of beta and lambda
  mat cov, info;
  vec parlc, sum1;
  double ebx;
  
  info.zeros(p+L, p+L);
  for(int i=0; i<n; ++i) {
    parlc.zeros(p+L);
    sum1.zeros(p);
    for(int l=0; l<L; ++l) {
      if(t(l)>Y(i)) {
        break;
      }
      ebx = exp(as_scalar(beta.t()*X(i,l)));
      sum1 += lambda(l)*e.z(i)*ebx*X(i,l);
      parlc(p+l) -= e.z(i)*ebx;
      if(fabs(Y(i)-t(l)) < 1e-8 && Delta(i)==1) {
        parlc.head(p) += X(i,l);
        parlc(p+l) += 1.0/lambda(l);
      }
    } // end l
    parlc.head(p) -= sum1;
    info += parlc * parlc.t();
  } // end i
  if(!pinv(cov, info)) {
    cov.fill(datum::nan);
    Rcout << "variance cannot be estimated due to singular information matrix!" << endl;
  }
  return cov;
}

mat PredSurv_Y(const vec& beta,
               const vec& lambda,
               const mat& cov, // covariance matrix of beta and lambda
               const vec& t,
               const vec& Yv, // validation data
               const field<vec>& Xv,
               double r,
               int nv,
               int L,
               int p,
               bool compVar) {
  // predict survival probability for validation data
  vec pred, var, sum1, parJ;
  double ebx, sum0, J, varJ;
  mat pred_res;
  
  pred.set_size(nv);
  var.set_size(nv);
  for(int i=0; i<nv; ++i) {
    sum0 = 0.0;
    sum1.zeros(p);
    parJ.zeros(p+L);
    for(int l=0; l<L; ++l) {
      if(t(l)>Yv(i)) {
        break;
      }
      ebx = exp(as_scalar(beta.t()*Xv(i,l)));
      sum0 += lambda(l)*ebx;
      sum1 += lambda(l)*ebx*Xv(i,l);
      parJ(p+l) = (fabs(r)>1e-8 ? r : 1.0)*ebx;
    }
    pred(i) = StX(sum0, r);
    if(compVar) {
      if(fabs(r)>1e-8) {
        parJ.head(p) += r*sum1;
        J = 1.0+r*sum0;
        varJ = as_scalar(parJ.t()*cov*parJ);
        var(i) = pow(1.0/r*pow(J, -1.0/r-1.0), 2.0)*varJ;
      } else {
        // J is defined as S
        parJ.head(p) = sum1;
        parJ *= -pred(i);
        var(i) = as_scalar(parJ.t()*cov*parJ);
      }
    } 
  } // end i
  pred_res = join_horiz(pred, var);
  return pred_res;
}

mat PredSurv(const vec& beta,
             const vec& lambda,
             const field<vec>& Xv,
             double r,
             int nv,
             int L) {
  // predict covariate-specific survival function for validation data
  // output is a n-by-L matrix, each row represents a subject,
  // each column represents a distinct time point.
  mat pred;
  double ebx, sum0;
  
  pred.set_size(nv, L);
  
  for(int i=0; i<nv; ++i) {
    sum0 = 0.0;
    for(int l=0; l<L; ++l) {
      ebx = exp(as_scalar(beta.t()*Xv(i,l)));
      sum0 += lambda(l)*ebx;
      pred(i,l) = StX(sum0, r);
    } // end l
  } // end i 
  return pred;
}

// [[Rcpp::export]]
List TransFit(const vec& Y,
              const ivec& Delta,
              const mat& rawX, // raw covariates, each subject may have multiple rows
              const vec& S, // predicted survival prob
              const vec& weight, // weight of S
              const vec& Yv, // validation data
              const ivec& Deltav,
              const mat& rawXv, 
              const mat& kmcensor, // KM estimator for censoring distribution
              double r, // variance of Gamma dist
              int N,
              double xi,
              double endt,
              double tol = 1e-4,
              int maxit = 2000,
              int verbose = 1,
              int pll = 0) {
  // Fit the transformation model
  int n, nv, p, L, Lv;
  vec t, tv, beta, lambda, plambda, med, rmst, true_rmst;
  uvec ind;
  field<vec> X, Xv, Xvtv;
  gaulag gl;
  List res;
  mat pred;
  
  n = Y.n_elem;
  nv = Yv.n_elem;
  t = unique(Y);
  L = t.n_elem;
  
  // convert rawX to X
  p = getX(X, rawX, t, n, L);
  p = getX(Xv, rawXv, t, nv, L);
  
  // initialize beta and lambda
  beta.zeros(p);
  lambda.set_size(L);
  lambda.fill(1.0/L);
  
  // generate GL weights and nodes
  gl = GLcompute(N);
  
  // EM algorithm
  EM(beta, lambda, t, Y, Delta, X, S, weight, gl, r, N, n, p, L, xi, tol, maxit, FALSE, verbose);

  // output beta and lambda estimates
  res["beta"] = beta;
  res["lambda"] = join_horiz(t, lambda);
  
  // difference in S(tau)
  ind = find(t<=endt);
  Lv = ind.n_elem;
  pred = PredSurv(beta, lambda, Xv, r, nv, Lv);
  // res["Stau"] = Stau(t, pred, Xv, r, nv, Lv, endt);
  
  // supremum distance between Shat and true S
  // res["supdistS"] = supdistS(t, pred, Xv, r, nv, Lv, endt);
  
  // L2 distance between Shat and true S
  // res["L2distS"] = L2distS(t, pred, Xv, r, nv, Lv, endt);
  
  // integrated Brier score
  res["intBS"] = intBrier(t, Yv, Deltav, pred, kmcensor, nv, endt);
  
  // median survival time
  med = medT(t, pred, nv, Lv, endt);
  res["medT"] = diffT(med, Yv, Deltav, kmcensor, nv, endt);
  
  // restricted mean survival time
  rmst = RMST(t, pred, nv, Lv, endt);
  true_rmst = min(Yv, ones<vec>(nv)*endt);
  res["RMST"] = diffT(rmst, true_rmst, Deltav, kmcensor, nv, endt);
  
  // likelihood for training data
  res["logL"] = logLikeli(beta, lambda, t, Y, Delta, X, r, n, L);
  
  // profile likelihood for testing data
  if(pll) {
    tv = unique(Yv);
    Lv = tv.n_elem;
    p = getX(Xvtv, rawXv, tv, nv, Lv);
    plambda.zeros(Lv);
    EM(beta, plambda, tv, Yv, Deltav, Xvtv, S, weight, gl, r, N, nv, p, Lv, xi, tol, maxit, TRUE, FALSE);
    res["logPL"] = logLikeli(beta, plambda, tv, Yv, Deltav, Xvtv, r, nv, Lv);
  }
  
  return res;
}

// [[Rcpp::export]]
List sourceFit(const vec& Ys,
               const ivec& Deltas,
               const mat& rawXs, // raw covariates, each subject may have multiple rows
               const vec& Y, // target data
               const mat& rawX,
               double r, // variance of Gamma dist
               int N,
               double tol = 1e-4,
               int maxit = 2000,
               int verbose = 1,
               int pll = 0,
               int compVar = 0,
               int predmat = 0) {
  // Fit the transformation model
  int ns, n, p, Ls;
  vec ts, beta, lambda, S, weight;
  uvec ind;
  field<vec> Xs, X;
  gaulag gl;
  List res;
  mat cov, pred_res;
  conde e;
  
  ns = Ys.n_elem;
  n = Y.n_elem;
  ts = unique(Ys(find(Deltas==1)));
  Ls = ts.n_elem;
  
  // convert rawX to X
  p = getX(Xs, rawXs, ts, ns, Ls);
  p = getX(X, rawX, ts, n, Ls);
  
  // initialize beta and lambda
  beta.zeros(p);
  lambda.set_size(Ls);
  lambda.fill(1.0/Ls);
  
  // create artificial S and weight (since they are used)
  S.ones(ns);
  weight.ones(ns);
  
  // generate GL weights and nodes
  gl = GLcompute(N);
  
  // EM algorithm
  EM(beta, lambda, ts, Ys, Deltas, Xs, S, weight, gl, r, N, ns, p, Ls, 0.0, tol, maxit, FALSE, verbose);
  
  // output beta and lambda estimates
  res["beta"] = beta;
  res["lambda"] = join_horiz(ts, lambda);
  
  // output source prediction and variance 
  e = Estep(beta, lambda, ts, Ys, Deltas, Xs, S, weight, gl, r, N, ns, Ls, 0.0);
  cov = compCov(beta, lambda, e, ts, Ys, Deltas, Xs, ns, Ls, p);
  pred_res = PredSurv_Y(beta, lambda, cov, ts, Y, X, r, n, Ls, p, compVar);
  res["predSY"] = pred_res.col(0);
  if(compVar) {
    res["varSY"] = pred_res.col(1);
  }
  
  if(predmat) {
    res["predmat"] = PredSurv(beta, lambda, X, r, n, Ls);
    res["uniqt"] = ts;
  }
  
  // likelihood for training data
  res["logL"] = logLikeli(beta, lambda, ts, Ys, Deltas, Xs, r, ns, Ls);
  
  return res;
}

// [[Rcpp::export]]
double xi_opt(int n,
              int ns,
              const vec& t1,
              const vec& t2,
              const mat& pred1,
              const mat& pred2,
              double endt,
              double q = 1.0,
              double c = 1000.0) {
  double xi, diff;
  
  xi = c*pow(ns, q/2.0)/pow(n, 1.0/2.0);
  diff = L2distS_2r(t1, t2, pred1, pred2, endt) + pow(ns, -q);
  if(diff/c > 1.0/n) {
    xi = 0.0;
  }
  return xi;
}

// [[Rcpp::export]]
double xi_prac(int n,
               int ns,
               const vec& t1,
               const vec& t2,
               const mat& pred1,
               const mat& pred2,
               double endt,
               double q = 1.0,
               double rho = 1.0,
               double c = 5.0) {
  double xi, diff;
  
  diff = L2distS_2r(t1, t2, pred1, pred2, endt);
  Rcout << "distS2: " << diff << endl;
  diff += pow(ns, -q);
  if(diff<=c) {
    xi = c*L2distS_2r(t1, t2, pred1, pred2, endt, -rho);
  } else {
    xi = 0.0;
  }

  return xi;
}

// [[Rcpp::export]]
List Metric(const vec& beta,
            const vec& lambda,
            const vec& t,
            const vec& Yv, // validation data
            const ivec& Deltav,
            const mat& rawXv, 
            const mat& kmcensor, // KM estimator for censoring distribution
            double r, // variance of Gamma dist
            double endt) {
  // compute different metrics based on beta and lambda estimates from input
  List res;
  mat pred;
  int p, nv, L;
  field<vec> Xv;
  vec med, rmst, true_rmst;
  
  nv = Yv.n_elem;
  L = t.n_elem;
  
  // convert rawX to X
  p = getX(Xv, rawXv, t, nv, L);
  
  res["beta"] = beta;
  res["lambda"] = join_horiz(t, lambda);
  
  // difference in S(tau)
  pred = PredSurv(beta, lambda, Xv, r, nv, L);
  // res["Stau"] = Stau(t, pred, Xv, r, nv, L, endt);
  
  // supremum distance between Shat and true S
  // res["supdistS"] = supdistS(t, pred, Xv, r, nv, L, endt);
  
  // L2 distance between Shat and true S
  // res["L2distS"] = L2distS(t, pred, Xv, r, nv, L, endt);
  
  // integrated Brier score
  res["intBS"] = intBrier(t, Yv, Deltav, pred, kmcensor, nv, endt);
  
  // median survival time
  med = medT(t, pred, nv, L, endt);
  res["medT"] = diffT(med, Yv, Deltav, kmcensor, nv, endt);
  
  // restricted mean survival time
  rmst = RMST(t, pred, nv, L, endt);
  true_rmst = min(Yv, ones<vec>(nv)*endt);
  res["RMST"] = diffT(rmst, true_rmst, Deltav, kmcensor, nv, endt);
  
  return res;
}
