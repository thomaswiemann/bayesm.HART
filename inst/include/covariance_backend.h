#ifndef COVARIANCE_BACKEND_H
#define COVARIANCE_BACKEND_H

#include <RcppArmadillo.h>

// Capability flags reported by each backend.
struct BackendCapabilities {
    bool has_mixture_draws;      // global: true, heter: false
    bool has_covariate_sigma;    // global: false, heter: true
    bool has_delta_z;            // true when mean-tree or linear delta active
};

// Typed payload returned by pack().
struct BackendPackPayload {
    Rcpp::List fields;
    BackendCapabilities caps;
};

// Abstract interface for covariance-regime sampling.
//
// Contract:
//   precision_root(i) returns upper-tri R_i with R_i * R_i^T = Sigma_i^{-1}.
//   prior_mean(i) returns m_i such that theta_i | . ~ N(m_i, Sigma_i).
//   Returned references remain valid until the next draw_iteration() call.
//
// Implementations:
//   GlobalCovarianceBackend   - mixture-of-normals + optional BART/linear delta
//   CovarDepCovarianceBackend - modified-Cholesky tree ensembles for Sigma(Z_i)
class CovarianceBackend {
public:
    virtual ~CovarianceBackend() = default;

    // Run all covariance sub-blocks for one MCMC iteration.
    // Must be called before precision_root / prior_mean are read.
    virtual void draw_iteration(int rep) = 0;

    // Upper-tri R_i with R_i R_i^T = Sigma_i^{-1}.
    virtual arma::mat const& precision_root(int i) const = 0;

    // Prior mean m_i for unit i.
    virtual arma::vec const& prior_mean(int i) const = 0;

    // Per-keep serialization.
    virtual void store(int mkeep) = 0;

    // Post-loop: return all backend outputs.
    virtual BackendPackPayload pack() = 0;

    // Report capabilities (static after init).
    virtual BackendCapabilities capabilities() const = 0;
};

#endif
