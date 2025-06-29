#include "bayesm.HART.h"

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