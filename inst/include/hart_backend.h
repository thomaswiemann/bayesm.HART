#ifndef HART_BACKEND_H
#define HART_BACKEND_H

#include <memory>
#include "bayesm.HART.h"
#include "covariance_backend.h"

// Shared state of the tree-based (HART) covariance backends: the unit
// coefficients the updater writes, Z^T in the column layout bart::setdata
// expects, the dimensions, the keep-time tree flag, the RNG, and the per-unit
// prior means the driver reads through prior_mean().  Subclasses implement
// draw_iteration(), precision_root(), store(), pack(), capabilities(), and
// refresh prior_mean_ at the end of every draw_iteration().
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

    // Optional logical entry of a parameter List.
    static bool list_flag(Rcpp::List const& p, const char* name, bool def) {
        return p.containsElementNamed(name) ? Rcpp::as<bool>(p[name]) : def;
    }

public:
    arma::vec const& prior_mean(int i) const override { return prior_mean_[i]; }
};

// Covariance backend for one sampler run (shared by the MNL, Linear, and
// NegBin loop files):
//   useHeterCov             -> HeterCovBackend         (direct)
//   useBART && !useHeterCov -> NormalizedHartBackend   (normalized)
//   otherwise               -> GlobalCovarianceBackend (bayesm paths)
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
