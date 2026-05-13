#ifndef __BAYESM_HART_H__
#define __BAYESM_HART_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include "bart.h"    
#include "heterbart.h"   
#include "varbart.h"
#include "cov_helpers.h"
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

// Heteroscedastic standardization: tilde_theta_i = D(Z_i)^{-1/2} L(Z_i) (theta_i - mu).
// Diagonal case: pass an empty phi_models.
void update_stdoldbetas_het(
    mat const& oldbetas,
    std::vector<double*>& pstd_oldbetas_cols,
    vec const& mu,
    std::vector<varbart> const& d_models,
    std::vector<std::vector<heterbart>> const& phi_models = std::vector<std::vector<heterbart>>());

// Construct D variance-tree ensembles, one per dimension.  var_params carries
// num_trees, power, base, nu (baseline), lambda (baseline), plus optional
// DART hyperparameters.  Pratola per-tree calibration is applied internally.
std::vector<varbart> initializeVarBART(
    double* pZt, size_t nlgt, size_t nz,
    std::vector<double*> const& peta_sq_cols,
    List const& var_params, arn gen);

// Construct D(D-1)/2 phi-tree ensembles (one per off-diagonal Cholesky entry).
// phi_models[j][k] is the heterbart for phi_{jk}(.), defined for j > k.
// phi_models[j].size() == j; phi_models[0] is empty.
//
// py_cols[j][k] points to a length-nlgt buffer that will hold the per-iteration
// pseudo-response Y_i = r_i / c_i^{(k)} (caller-owned, mutated in place each
// MCMC iteration).
//
// phi_params: num_trees, power, base, tau, [numcut, sparse, ...].  Sets
// pi.nmin = 2 and pi.ess_min = ess_min (default 5.0) per the
// varying-coefficient transform; see discussions/2026-05-12-phibart-vs-heterbart.md.
std::vector<std::vector<heterbart>> initializePhiBART(
    double* pZt, size_t nlgt, size_t nz,
    std::vector<std::vector<double*>> const& py_cols,
    List const& phi_params, arn gen);

// Conjugate normal draw of the single posterior mean mu under heteroscedastic
// Sigma(Z_i):  mu | {theta_i}, Sigma(.) ~ N(m_post, V_post) with
//   V_post^{-1} = Amu + sum_i Sigma(Z_i)^{-1}
//   m_post     = V_post (Amu * mubar0 + sum_i Sigma(Z_i)^{-1} theta_i).
// mubar0 is the prior mean vector (length p).  Diagonal-only path is handled
// transparently by passing an empty phi_models.
vec drawMuHeterCov(
    mat const& oldbetas,
    vec const& mubar0,
    mat const& Amu,
    std::vector<varbart> const& d_models,
    std::vector<std::vector<heterbart>> const& phi_models,
    arn& gen);

unireg runiregG(vec const& y, mat const& X, mat const& XpX, vec const& Xpy, double sigmasq, mat const& A, vec const& Abetabar, double nu, double ssq);

//Used in rhierMnlDP and rhierMnlRwMixture
mnlMetropOnceOut mnlMetropOnce(vec const& y, mat const& X, vec const& oldbeta, double oldll,double s, mat const& incroot, vec const& betabar, mat const& rootpi);

//Used in rnegbinRW and rhierNegbinRw
double llnegbin(vec const& y, vec const& lambda, double alpha, bool constant);
// `rooti` is the bayesm-`rooti` upper-triangular matrix with rooti * rooti^T = Sigma^{-1}.
// (Same convention as lndMvn / mnlMetropOnce_con.  See utilityFunctions.cpp for details.)
double lpostbeta(double alpha, vec const& beta, mat const& X, vec const& y, vec const& betabar, mat const& rooti);
double lpostalpha(double alpha, vec const& beta, mat const& X, vec const& y, double a, double b);

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