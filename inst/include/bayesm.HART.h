#ifndef __BAYESM_HART_H__
#define __BAYESM_HART_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include "bart.h"    
#include "heterbart.h"   
#include "rn.h"      
#include "common.h"

#include "logging.h"    

using namespace arma;
using namespace Rcpp;

//CUSTOM STRUCTS--------------------------------------------------------------------------------------------------
//Used in rhierLinearMixture, rhierLinearModel, rhierMnlDP, rhierMnlRwMixture, rhierNegbinRw, and rsurGibbs
struct moments{
  vec y;
  mat X;
  mat XpX;
  vec Xpy;
  mat hess;
};

//Used in rhierLinearMixture, rhierLinearModel, rhierMnlRWMixture, and utilityFunctions.cpp
struct unireg{
    vec beta;
    double sigmasq;
  };

//Used in rhierMnlDP, rhierMnlRwMixture, and utilityFunctions.cpp
struct mnlMetropOnceOut{
  vec betadraw;
  int stay;
  double oldll;
};  

//EXPOSED FUNCTIONS-----------------------------------------------------------------------------------------------
List rwishart(double nu, mat const& V);

List rmultireg(mat const& Y, mat const& X, mat const& Bbar, mat const& A, double nu, mat const& V);

vec rdirichlet(vec const& alpha);

double llmnl(vec const& beta, vec const& y, mat const& X);

double lndMvn(vec const& x, vec const& mu, mat const& rooti);

double lndIWishart(double nu, mat const& V, mat const& IW);

vec breg(vec const& y, mat const& X, vec const& betabar, mat const& A);

List rmixGibbs( mat const& y,  mat const& Bbar, mat const& A, double nu, mat const& V,  vec const& a, vec const& p,  vec const& z);
  //rmixGibbs contains the following support functions, which are called ONLY THROUGH rmixGibbs: drawCompsFromLabels, drawLabelsFromComps, and drawPFromLabels
List rmixGibbs_BART(mat const& y, mat const& Bbar, mat const& A, double nu, mat const& V, vec const& a, vec const& p, vec const& z);

//SUPPORT FUNCTIONS (contained in utilityFunctions.cpp and trunNorm.cpp)-----------------------------------------------------------

//Used in rhierLinearModel, rhierLinearMixture and rhierMnlRWMixture
mat drawDelta(mat const& x,mat const& y,vec const& z,List const& comps,vec const& deltabar,mat const& Ad);

void drawBART(mat const& x, mat const& y, std::vector<bart>& models, vec const& ind, List const& comps, arn gen);
void update_stdoldbetas(mat const& oldbetas, std::vector<double*>& pstd_oldbetas_cols, vec const& ind, List const& comps);
std::vector<bart> initializeBART(double* pZt, size_t nlgt, size_t nz, std::vector<double*> const& pstd_oldbetas_cols, List const& bart_params, arn gen);

unireg runiregG(vec const& y, mat const& X, mat const& XpX, vec const& Xpy, double sigmasq, mat const& A, vec const& Abetabar, double nu, double ssq);

//Used in rhierMnlDP and rhierMnlRwMixture
mnlMetropOnceOut mnlMetropOnce(vec const& y, mat const& X, vec const& oldbeta, double oldll,double s, mat const& incroot, vec const& betabar, mat const& rootpi);

//FUNCTION TIMING (contained in functionTiming.cpp)---------------------------------------------------------------
void startMcmcTimer();
void infoMcmcTimer(int rep, int R);
void endMcmcTimer();

//DEBUG FUNCTIONS
inline void checkMatrix(mat const& M, const char* name) {
	cout << "Matrix " << name << " dimensions: " << M.n_rows << "x" << M.n_cols << endl;
	cout << "Matrix " << name << " condition number: " << cond(M) << endl;
	cout << "Matrix " << name << " determinant: " << det(M) << endl;
	cout << "First few elements: " << M(0, 0) << " " << M(0, 1) << endl;
}

#endif