#ifndef HETERCOV_BLOCK_H
#define HETERCOV_BLOCK_H

#include <RcppArmadillo.h>
#include <vector>
#include <sstream>
#include "bayesm.HART.h"

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// HeterCovState  --  modified-Cholesky tree-ensemble Sigma(Z_i) sampler state.
//
// Owns the mean / diagonal-variance / off-diagonal-phi tree ensembles, all
// associated working buffers, output cubes, and the per-iteration mu and
// delta_Z vectors.  Designed to be the single source of truth for heter-cov
// MCMC across rhierMnlRwMixture, rhierNegbinRw, and rhierLinearMixture.
//
// The model-specific Step E (per-unit theta_i posterior draw) is the ONLY
// piece that lives outside this struct; it accesses {var_models, phi_models,
// mu_post, delta_Z} read-only via cov::cov_eval / cov::rootpi.
// ---------------------------------------------------------------------------
struct HeterCovState {
    // ---- Shape parameters --------------------------------------------------
    int  nlgt, nz, nvar, R, keep;
    bool useFullCov;          // true iff phi_params is non-empty (Phase 2)
    size_t n_pairs;           // = nvar*(nvar-1)/2 if useFullCov else 0

    // DART configuration: mean / variance / phi trees independently configured.
    int  burn, var_burn, phi_burn;
    bool sparse, var_sparse, phi_sparse;

    // ---- Tree ensembles (owned) -------------------------------------------
    std::vector<bart>                          bart_models;  // mean trees
    std::vector<varbart>                       var_models;   // diagonal d_j
    std::vector<std::vector<heterbart>>        phi_models;   // phi_jk (jagged)

    // ---- Tree-string serialization streams (owned) ------------------------
    std::vector<std::stringstream>             treess;       // mean
    std::vector<std::stringstream>             var_treess;   // var
    std::vector<std::vector<std::stringstream>> phi_treess;  // phi (jagged)
    bool store_trees;

    // ---- Output cubes / matrices (owned) ----------------------------------
    arma::cube  varcount,    varprob;       // mean trees
    arma::cube  var_varcount, var_varprob;  // var  trees
    arma::mat   mu_draw;                    // R/keep x nvar

    // ---- Working buffers (owned) ------------------------------------------
    arma::mat              eta_sq_buf;      // nlgt x nvar (var-tree response)
    std::vector<double*>   peta_sq_cols;    // ptrs into eta_sq_buf cols

    arma::mat              Y_phi;           // nlgt x n_pairs (phi-tree response)
    arma::vec              sigma_phi_buf;   // nlgt (per-iter phi sigmas)
    std::vector<std::vector<double*>> py_phi_cols;  // ptrs into Y_phi (jagged)

    // ---- Per-iter state consumed by Step E --------------------------------
    arma::vec   mu_post;     // current mu draw (length nvar)
    arma::mat   delta_Z;     // nlgt x nvar; delta_Z(i,.) = Sigma(Z_i)^{1/2} delta_i

    // ---- Optional keep-time cache payload ----------------------------------
    bool        has_unique_map;
    arma::uvec  unique_first_idx;       // n_unique, 0-based row indices
    arma::cube  delta_z_unique_draws;   // n_unique x nvar x (R/keep)
    arma::cube  sigma_z_unique_flat;    // n_unique x (nvar*nvar) x (R/keep)
};

// ---------------------------------------------------------------------------
// hetercov_init
//   Pre-loop allocation.  Sets up shape parameters, DART config, tree-stream
//   headers, output cubes, working buffers, and per-iter state (zeroed).
//
//   Parameters:
//     bart_params, var_params, phi_params -- parsed Lists from the R wrapper
//     (see R/heter_cov_priors.R).  bart_params drives mean-tree configuration
//     even though the trees themselves are initialized inside the loop on
//     rep == 0.
// ---------------------------------------------------------------------------
void hetercov_init(HeterCovState& s,
                   int nlgt, int nz, int nvar, int R, int keep,
                   Rcpp::List const& bart_params,
                   Rcpp::List const& var_params,
                   Rcpp::List const& phi_params);

// ---------------------------------------------------------------------------
// hetercov_draw_iter
//   Per-iteration MCMC step.  Internally dispatches on rep == 0 to perform
//   first-iter initialization (baseline mu, baseline eta_sq, var/phi/mean
//   tree initialization).  For rep >= 1 executes the Wiemann (2025) Sec. 3.1
//   Gibbs ordering:
//
//     Step A : drawMuHeterCov                   -> s.mu_post
//     Step C : variance trees d_j               -> var_models[j]
//     Step D': phi_jk trees (Phase 2 only)      -> phi_models[j][k]
//     Step B : mean trees delta_j               -> bart_models[j]
//     end    : recompute s.delta_Z = Sigma(Z_i)^{1/2} delta_i
//
//   OWNERSHIP CONTRACT -- pstd_oldbetas_cols:
//     Caller MUST allocate a writable working buffer (e.g.
//     `mat std_oldbetas = oldbetas;`) and populate
//     `pstd_oldbetas_cols[j] = std_oldbetas.colptr(j)` BEFORE the first call.
//     This function writes standardized betas into that storage via
//     `update_stdoldbetas_het` and reads them back through bart::setdata
//     pointers established in `initializeBART` on rep == 0.
// ---------------------------------------------------------------------------
void hetercov_draw_iter(HeterCovState& s,
                        arma::mat const& oldbetas,
                        arma::vec const& mubar0,
                        arma::mat const& Amu,
                        double* pZt,
                        std::vector<double*>& pstd_oldbetas_cols,
                        Rcpp::List const& bart_params,
                        Rcpp::List const& var_params,
                        Rcpp::List const& phi_params,
                        int rep, arn& gen);

// ---------------------------------------------------------------------------
// hetercov_store
//   Per-keep serialization: appends each tree's structure to the per-dimension
//   tree-string stream and slices varcount/varprob into the (mkeep-1) layer.
//   Writes the current mu_post into mu_draw row mkeep-1.
// ---------------------------------------------------------------------------
void hetercov_store(HeterCovState& s, int mkeep);

// ---------------------------------------------------------------------------
// hetercov_pack
//   Post-loop: returns a List with all heter-cov outputs:
//     bart_models, var_models, phi_models (NULL if diagonal-only),
//     mu_draw, varcount, varprob, var_varcount, var_varprob.
//
//   The caller merges its own model-specific outputs (betadraw, loglike,
//   SignRes, taudraw, alphadraw, ...) into the returned List before returning
//   to R.
// ---------------------------------------------------------------------------
Rcpp::List hetercov_pack(HeterCovState& s);

#endif
