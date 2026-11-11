#ifndef UNIQUE_Z_CACHE_H
#define UNIQUE_Z_CACHE_H

#include <RcppArmadillo.h>

// Keep-time cache of per-unit quantities at the unique rows of Z.  The R
// wrapper supplies `z_index` (1-based unit -> unique row) and `n_unique` in
// the tree-parameter List (.build_unique_z_map); each keep stores Delta(Z_u)
// for the first unit of every unique row and, with covariate-dependent
// covariance, Sigma(Z_u) flattened column-major (.assemble_hart_cache reads
// them back).  Without `z_index` the cache is inactive and store_*/pack() are
// no-ops.
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
