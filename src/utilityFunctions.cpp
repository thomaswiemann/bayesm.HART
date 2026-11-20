#include "bayesm.HART.h"

unireg runiregG(vec const& y, mat const& X, mat const& XpX, vec const& Xpy, double sigmasq, mat const& A, 
              vec const& Abetabar, double nu, double ssq) {

// Keunwoo Kim 09/16/2014

// Purpose: 
//  perform one Gibbs iteration for Univ Regression Model
//  only does one iteration so can be used in rhierLinearModel

// Model:
//  y = Xbeta + e  e ~N(0,sigmasq)
//  y is n x 1
//  X is n x k
//  beta is k x 1 vector of coefficients

// Prior:  
//  beta ~ N(betabar,A^-1)
//  sigmasq ~ (nu*ssq)/chisq_nu

  unireg out_struct;
  
  int n = y.size();
  int k = XpX.n_cols;
  
  //first draw beta | sigmasq
  mat IR = solve(trimatu(chol(XpX/sigmasq+A)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  vec btilde = (IR*trans(IR)) * (Xpy/sigmasq + Abetabar);
  vec beta = btilde + IR*vec(rnorm(k));
  
  //now draw sigmasq | beta
  double s = sum(square(y-X*beta));
  sigmasq = (s + nu*ssq)/rchisq(1,nu+n)[0]; //rchisq returns a vectorized object, so using [0] allows for the conversion to double
  
  out_struct.beta = beta;
  out_struct.sigmasq = sigmasq;  

  return (out_struct);
}

//Used in rnegbinRW and rhierNegbinRw-------------------------------------------------------------------------------------
double llnegbin(vec const& y, vec const& lambda, double alpha, bool constant){

// Keunwoo Kim 11/02/2014

// Computes the log-likelihood

  int i;
  int nobs = y.size();  
  vec prob = alpha/(alpha+lambda);    
  vec logp(nobs);
  if (constant){
    // normalized log-likelihood
    for (i=0; i<nobs; i++){
      logp[i] = R::dnbinom(y[i], alpha, prob[i], 1);
    }    
  }else{
    // unnormalized log-likelihood
    logp = sum(alpha*log(prob) + y % log(1-prob));
  }
  return (sum(logp));
}

double lpostbeta(double alpha, vec const& beta, mat const& X, vec const& y, vec const& betabar, mat const& rooti){

// Keunwoo Kim 11/02/2014  (bayesm)
// rooti is bayesm's root with rooti * trans(rooti) = Sigma^{-1} (rmixGibbs,
// cov::rootpi, lndMvn convention), not chol(Sigma^{-1}).

  vec lambda = exp(X*beta);
  double ll = llnegbin(y, lambda, alpha, FALSE);

  // unnormalized log prior  N(beta | betabar, Sigma)  with rooti rooti^T = Sigma^{-1}
  vec z = trans(rooti)*(beta-betabar);
  double lprior = - 0.5*sum(z%z);

  return (ll+lprior);
}
