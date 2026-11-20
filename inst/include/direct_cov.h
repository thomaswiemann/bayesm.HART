#ifndef DIRECT_COV_H
#define DIRECT_COV_H

#include <RcppArmadillo.h>
#include "rn.h"
#include "unique_z_cache.h"

// Covariance model of the direct HART backend:  theta_i = mu + f(Z_i) + eps_i
// with eps_i ~ N(0, Sigma_i).  The backend owns the mean blocks and defers
// Sigma to one of these (appendix_sampler, app:hs_sweep block 3).
//
// Lifetime:  init() once with the backend's dimensions and Z buffer, start()
// once with the initial residuals (no RNG), then draw() every sweep with
// E = theta - mu - f.  precisions() and root(i) stay valid until the next
// draw().
class DirectCov {
public:
    virtual ~DirectCov() = default;

    virtual void init(int n, int D, int nz, double* pZt, int nkeep) = 0;
    virtual void start(arma::mat const& E) = 0;
    // DART switch at rep == burn + 1, per ensemble.
    virtual void start_dart(int rep) = 0;
    virtual void draw(arma::mat const& E, rn& gen) = 0;

    virtual std::vector<arma::mat> const& precisions() const = 0;
    // Upper-tri R_i with R_i R_i^T = Sigma_i^{-1}.
    virtual arma::mat const& root(int i) const = 0;

    virtual void store(int mkeep, arma::vec const& mu, bool store_trees,
                       UniqueZCache& zc) = 0;
    virtual void pack(Rcpp::List& out, bool store_trees) = 0;
    virtual Rcpp::List checkpoint() = 0;
    virtual void restore(Rcpp::List const& state) = 0;
    // Sigma varies with Z: drives capabilities() and the Sigma(Z) cache.
    virtual bool covariate_sigma() const = 0;
};

#endif
