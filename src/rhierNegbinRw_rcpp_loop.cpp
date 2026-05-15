#include "bayesm.HART.h"
#include "global_covariance_backend.h"
#include "hetercov_backend.h"
#include "negbin_theta_updater.h"
#include "unified_sampler_driver.h"
#include <memory>

//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
double llnegbinpooled(std::vector<moments> const& regdata_vector,
                      mat const& oldbetas,
                      double alpha){

// Wayne Taylor 12/01/2014

// "Unlists" the regdata and calculates the negative binomial loglikelihood using individual-level betas

  int nreg = regdata_vector.size();
  double ll = 0.0;

  for(int reg = 0; reg<nreg; reg++){
  vec lambda = exp(regdata_vector[reg].X*trans(oldbetas(reg,span::all)));
  ll = ll + llnegbin(regdata_vector[reg].y,lambda,alpha,TRUE);
  }

  return(ll);
}

// [[Rcpp::export]]
List rhierNegbinRw_rcpp_loop(List const& regdata, List const& hessdata, mat const& Z,
                             mat oldbetas, vec const& deltabar, mat const& Ad,
                             mat const& mubar, mat const& Amu,
                             double nu, mat const& V, double a, double b,
                             int R, int keep, double sbeta, double alphacroot, int nprint,
                             bool drawdelta, mat olddelta, vec const& a_mix,
                             vec oldprob, vec ind,
                             double alpha, bool fixalpha,
                             bool useBART, List const& bart_params,
                             bool useHeterCov,
                             List const& var_params,
                             List const& phi_params){

// Wayne Taylor 12/01/2014
// Modified by Thomas Wiemann 2025 -- HART/BART support, mixture-of-normals,
//   plus heter-cov branch (modified-Cholesky tree ensembles for Sigma(Z_i)).
//
// Three execution modes (mutually exclusive):
//   * Heter-cov: hetercov_block ensembles for d_j(Z_i) [+ phi_jk(Z_i)].
//   * Mix-BART : sum-of-trees Delta(Z) on top of mixture-of-normals u_i.
//   * Linear   : original bayesm linear-Delta + rmixGibbs path.

  // Heter-cov path validation.
  if (useHeterCov && !useBART)   stop("useHeterCov requires useBART = TRUE");
  if (useHeterCov && !drawdelta) stop("useHeterCov requires drawdelta = TRUE");

  int nreg = regdata.size();
  int nvar = V.n_cols;

  // convert regdata and hessdata Lists to std::vector of struct
  std::vector<moments> regdata_vector;
  moments regdatai_struct;
  List regdatai, hessi;

  for (int reg = 0; reg<nreg; reg++){
    regdatai = regdata[reg];
    hessi = hessdata[reg];

    regdatai_struct.y = as<vec>(regdatai["y"]);
    regdatai_struct.X = as<mat>(regdatai["X"]);
    regdatai_struct.hess = as<mat>(hessi["hess"]);
    regdata_vector.push_back(regdatai_struct);
  }

  // allocate driver-owned storage
  cube betadraw = zeros<cube>(nreg, nvar, R/keep);
  vec loglike = zeros<vec>(R/keep);

  // ====================================================================
  // Create backend + updater adapters, then run the unified MCMC loop.
  // ====================================================================
  std::unique_ptr<CovarianceBackend> backend;
  if (useHeterCov) {
    backend = std::make_unique<HeterCovBackend>(
        oldbetas, Z, mubar, Amu, R, keep,
        bart_params, var_params, phi_params);
  } else {
    backend = std::make_unique<GlobalCovarianceBackend>(
        oldbetas, Z, mubar, Amu, nu, V, a_mix,
        deltabar, Ad, drawdelta, useBART,
        olddelta, oldprob, ind, R, keep, bart_params);
  }

  NegbinThetaUpdaterAdapter updater(oldbetas, regdata_vector,
                                     sbeta, alpha, alphacroot,
                                     a, b, fixalpha, R, keep);

  run_unified_sampler(R, keep, nprint, nreg,
                      *backend, updater, oldbetas, betadraw, loglike);

  // ====================================================================
  // Output assembly
  // ====================================================================
  BackendPackPayload bp = backend->pack();
  UpdaterPackPayload up = updater.pack();
  AcceptanceMetrics am = updater.acceptance_metrics();

  List out = bp.fields;
  out["betadraw"]     = betadraw;
  out["loglike"]      = loglike;
  out["acceptrbeta"]  = am.acceptrbeta;
  out["acceptralpha"] = am.acceptralpha;

  // Merge updater fields (alphadraw)
  CharacterVector upnames = up.fields.names();
  for (int i = 0; i < up.fields.size(); i++) {
    out[as<std::string>(upnames[i])] = up.fields[i];
  }

  return out;
}
