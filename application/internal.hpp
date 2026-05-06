#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Gauss-Lag quadrature
struct gaulag{
  vec node;
  vec weight;
};

gaulag GLcompute(int N) {
  // output Gauss-Lag nodes and weights
  double x[N], w[N];
  gaulag gl;
  gl.node.set_size(N);
  gl.weight.set_size(N);
  
  // compute nodes and weights
  laguerre_compute(N, x, w);
  
  // copy nodes and weights
  for(int i=0; i<N; ++i) {
    gl.node(i) = x[i];
    gl.weight(i) = w[i];
  }
  
  return gl;
}

int getX(field<vec>& X, 
         const mat& rawX, 
         const vec& t, 
         int n,
         int L) {
  /* convert raw X to a n-by-L field with each element of p dimensions. 
   rawX: id | start time | end time | X1 | X2 | ... */
  int p, i, l;
  rowvec Xil;
  
  X.set_size(n,L);
  p = rawX.n_cols-3;
  
  i = 0; // i: row index
  for(int curi=0; curi<n; ++curi) {
    // curi: current subject id
    l = 0;
    while(i < rawX.n_rows && rawX(i,0) == curi) {
      while(l<L && rawX(i,1)<t(l) && t(l)<=rawX(i,2)) {
        Xil = rawX(i, span(3,p+2));
        X(curi, l) = Xil.t();
        l++;
      }
      i++;
    }
  }
  
  return p;
}

double StX(double x, // Lambda(t|X)
           double r) {
  // return S(t|X) 
  double St, Gx;
  
  if(fabs(r)<1e-8) {
    Gx = x;
  } else {
    Gx = log(1+r*x)/r;
  }
  
  St = exp(-Gx);
  return St;
}

double Gprime(double x, 
              double r) {
  // compute the derivative of G(x)
  double deriv;
  if(fabs(r)<1e-8) {
    deriv = 1.0;
  } else {
    deriv = 1.0/(1.0+r*x);
  }
  return deriv;
}
