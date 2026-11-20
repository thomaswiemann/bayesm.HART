#include "bayesm.HART.h"
#include "hart_backend.h"
#include "linear_theta_updater.h"
#include "unified_sampler_driver.h"
#include <memory>

//[[Rcpp::export]]
List rhierLinearMixture_rcpp_loop(List const& regdata, mat const& Z,
                                  vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,
                                  double nu, mat const& V, double nu_e, vec const& ssq,
                                  int R, int keep, int nprint, bool drawdelta,
                                  mat olddelta, vec const& a, vec oldprob, vec ind, vec tau,
                                  std::string const& hart_model, std::string const& hart_cov,
                                  List const& bart_params,
                                  List const& var_params,
                                  List const& phi_params,
                                  mat const& Beta_init){

// Wayne Taylor 10/02/2014
// Modified by Thomas Wiemann 2025 -- HART support.
// Covariance backend by hart_model (make_backend in hart_backend.h).

  if (hart_model != "none" && !drawdelta)
    stop("hart_model = '%s' requires Z.", hart_model.c_str());

  int nreg = regdata.size();
  int nvar = V.n_cols;

  // convert List to std::vector of type "moments"
  std::vector<moments> regdata_vector;
  moments regdatai_struct;
  List regdatai;

  for (int reg = 0; reg<nreg; reg++){
    regdatai = regdata[reg];

    regdatai_struct.y = as<vec>(regdatai["y"]);
    regdatai_struct.X = as<mat>(regdatai["X"]);
    regdatai_struct.XpX = as<mat>(regdatai["XpX"]);
    regdatai_struct.Xpy = as<vec>(regdatai["Xpy"]);
    regdata_vector.push_back(regdatai_struct);
  }

  // allocate driver-owned storage
  mat oldbetas = zeros<mat>(nreg, nvar);
  cube betadraw(nreg, nvar, R/keep);
  vec loglike(R/keep);

  // Warm start at the per-unit OLS slopes under a HART model.
  if (hart_model != "none" && Beta_init.n_rows == (uword)nreg &&
      Beta_init.n_cols == (uword)nvar) {
    oldbetas = Beta_init;
  }

  // ====================================================================
  // Create backend + updater adapters, then run the unified MCMC loop.
  // ====================================================================
  std::unique_ptr<CovarianceBackend> backend = make_backend(
      hart_model, hart_cov, oldbetas, Z, mubar, Amu, nu, V, a,
      deltabar, Ad, drawdelta, olddelta, oldprob, ind, R, keep,
      bart_params, var_params, phi_params);

  LinearThetaUpdaterAdapter updater(oldbetas, regdata_vector, tau, nu_e, ssq, R, keep);

  run_unified_sampler(R, keep, nprint, nreg,
                      *backend, updater, oldbetas, betadraw, loglike);

  // ====================================================================
  // Output assembly: backend.pack() + updater.pack() + betadraw/loglike
  // ====================================================================
  BackendPackPayload bp = backend->pack();
  UpdaterPackPayload up = updater.pack();

  List out = bp.fields;
  out["betadraw"] = betadraw;
  out["loglike"]  = loglike;

  // Merge updater fields (taudraw)
  CharacterVector upnames = up.fields.names();
  for (int i = 0; i < up.fields.size(); i++) {
    out[as<std::string>(upnames[i])] = up.fields[i];
  }

  return out;
}
