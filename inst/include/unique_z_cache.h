#ifndef UNIQUE_Z_CACHE_H
#define UNIQUE_Z_CACHE_H

#include <RcppArmadillo.h>

// Keep-time cache of per-unit quantities at the unique rows of Z.
//
// The R wrapper adds `z_index` (1-based map from units to unique Z rows) and
// `n_unique` to the tree-parameter List (see R/utils.R::.build_unique_z_map).
// At each keep step the backend stores the coefficient-scale mean offset
// Delta(Z_u) for the first unit of every unique row and, for covariate-
// dependent covariance, Sigma(Z_u) flattened column-major.  R reassembles the
// draws into the `hart_cache` attribute (R/utils.R::.assemble_hart_cache).
//
// Inactive (no `z_index` supplied) is a valid state: store_* and pack() are
// then no-ops, so callers need no branching.
struct UniqueZCache {
    bool       active = false;
    arma::uvec first_idx;         // n_unique, 0-based row of first appearance
    arma::cube delta_draws;       // n_unique x nvar x nkeep
    arma::cube sigma_draws_flat;  // n_unique x (nvar*nvar) x nkeep (optional)

    // `who` names the caller in error messages.
    void init(Rcpp::List const& params, int n, int nvar, int nkeep,
              bool with_sigma, const char* who);

    void store_delta(int mkeep, arma::mat const& delta_Z);
    void store_sigma(int mkeep, arma::uword u, arma::mat const& Sigma);
    void pack(Rcpp::List& out) const;

    arma::uword n_unique() const { return first_idx.n_elem; }
};

#endif
