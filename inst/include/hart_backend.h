#ifndef HART_BACKEND_H
#define HART_BACKEND_H

#include <memory>
#include "bayesm.HART.h"
#include "covariance_backend.h"

// Shared state for the tree-based (HART) covariance backends.
//
// Holds what every HART backend needs regardless of preference model: the
// unit coefficients the updater writes, the transposed characteristics Z^T
// (unit i is column i, the layout bart::setdata expects), the problem
// dimensions, the keep-time tree flag, the RNG, and the per-unit prior-mean
// cache the unified driver reads through prior_mean().
//
// Subclasses implement draw_iteration(), precision_root(), store(), pack()
// and capabilities().  They must refresh prior_mean_ at the end of every
// draw_iteration().
//
// Implementations:
//   HeterCovBackend        - stage 0: pre-plan modified-Cholesky kernel
//   NormalizedHartBackend  - stage 1: theta_i = mu + S{f*(Z_i) + eps_i}
//   DirectHartBackend      - stage 4: theta_i = mu + f(Z_i) + S_i eps_i
class HartBackend : public CovarianceBackend {
protected:
    arma::mat const&       theta_;       // n x D unit coefficients (driver-owned)
    arma::mat              Zt_;          // nz x n
    int                    n_, nz_, D_, nkeep_;
    bool                   store_trees_;
    std::vector<arma::vec> prior_mean_;  // n, refreshed each iteration
    arn                    gen_;

    HartBackend(arma::mat const& theta, arma::mat const& Z, int R, int keep,
                bool store_trees)
        : theta_(theta), Zt_(Z.t()),
          n_(theta.n_rows), nz_(Z.n_cols), D_(theta.n_cols), nkeep_(R / keep),
          store_trees_(store_trees),
          prior_mean_(theta.n_rows, arma::vec(theta.n_cols, arma::fill::zeros)) {}

    // Column-major Z^T buffer handed to bart::setdata; stable for the
    // lifetime of the backend.
    double* pZt() { return Zt_.memptr(); }

public:
    arma::vec const& prior_mean(int i) const override { return prior_mean_[i]; }
};

// Select the covariance backend for one sampler run (shared by the MNL,
// Linear, and NegBin loop files):
//   useHeterCov            -> HeterCovBackend        (direct, pre-plan kernel)
//   useBART && !useHeterCov -> NormalizedHartBackend  (normalized HART)
//   otherwise              -> GlobalCovarianceBackend (bayesm paths)
std::unique_ptr<CovarianceBackend> make_backend(
    bool useBART, bool useHeterCov,
    arma::mat& oldbetas, arma::mat const& Z,
    arma::mat const& mubar, arma::mat const& Amu,
    double nu, arma::mat const& V, arma::vec const& a,
    arma::vec const& deltabar, arma::mat const& Ad, bool drawdelta,
    arma::mat const& olddelta, arma::vec const& oldprob, arma::vec const& ind,
    int R, int keep,
    Rcpp::List const& bart_params, Rcpp::List const& var_params,
    Rcpp::List const& phi_params);

#endif
