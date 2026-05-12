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

double lpostbeta(double alpha, vec const& beta, mat const& X, vec const& y, vec const& betabar, mat const& rootA){

// Keunwoo Kim 11/02/2014

// Computes log posterior for beta | alpha

  vec lambda = exp(X*beta);
  double ll = llnegbin(y, lambda, alpha, FALSE);

  // unormalized prior
  vec z = rootA*(beta-betabar);
  double lprior = - 0.5*sum(z%z);
  
  return (ll+lprior);
}

double lpostalpha(double alpha, vec const& beta, mat const& X, vec const& y, double a, double b){

// Keunwoo Kim 11/02/2014

// Computes log posterior for alpha | beta

  vec lambda = exp(X*beta);
  double ll = llnegbin(y, lambda, alpha, TRUE);
  // unormalized prior
  double lprior = (a-1)*log(alpha) - b*alpha;  
  
  return (ll+lprior);
}

//Used in rhierMnlDP and rhierMnlRwMixture------------------------------------------------------------------------
mnlMetropOnceOut mnlMetropOnce(vec const& y, mat const& X, vec const& oldbeta, 
                                                 double oldll,double s, mat const& incroot, 
                                                 vec const& betabar, mat const& rootpi){ 
// Wayne Taylor 10/01/2014

// function to execute rw metropolis for the MNL
// y is n vector with element = 1,...,j indicating which alt chosen
// X is nj x k matrix of xvalues for each of j alt on each of n occasions
// RW increments are N(0,s^2*t(inc.root)%*%inc.root)
// prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
//  inc.root, rootpi are upper triangular
//  this means that we are using the UL decomp of Sigma^-1 for prior 
// oldbeta is the current


mnlMetropOnceOut metropout_struct;

double unif;
vec betadraw, alphaminv(2);

int stay = 0;
vec betac = oldbeta + s*trans(incroot)*as<vec>(rnorm(X.n_cols));
double cll = llmnl(betac,y,X);
double clpost = cll+lndMvn(betac,betabar,rootpi);
double ldiff = clpost-oldll-lndMvn(oldbeta,betabar,rootpi);
//alphaminv << 1 << exp(ldiff);
alphaminv = { 1, exp(ldiff) };
double alpha = min(alphaminv);

     if(alpha < 1) {
       unif = runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double
      } else { 
        unif=0;}
     if (unif <= alpha) {
       betadraw = betac;
       oldll = cll;
      } else {
        betadraw = oldbeta;
        stay = 1;
      }

metropout_struct.betadraw = betadraw;
metropout_struct.stay = stay;  
metropout_struct.oldll = oldll;

return (metropout_struct);
}