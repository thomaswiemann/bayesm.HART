#ifndef THETA_UPDATER_H
#define THETA_UPDATER_H

#include <RcppArmadillo.h>

struct UpdaterPackPayload {
    Rcpp::List fields;
};

struct AcceptanceMetrics {
    double acceptrbeta  = -1.0;   // -1 means not applicable
    double acceptralpha = -1.0;   // NegBin only
};

// Abstract interface for likelihood-specific per-unit theta updates.
//
// Contract:
//   update_unit(i, rootpi_i, prior_mean_i) draws theta_i from its
//   full conditional, using rootpi_i and prior_mean_i as the prior
//   kernel parameters.  The updater must not recompute rootpi or
//   prior_mean independently.
//
// Implementations:
//   MnlThetaUpdater    - RW-MH with SignRes support
//   LinearThetaUpdater - conjugate Gaussian via runiregG
//   NegbinThetaUpdater - RW-MH for beta_i + alpha MH block
class ThetaUpdater {
public:
    virtual ~ThetaUpdater() = default;

    // Draw theta_i from full conditional.
    virtual void update_unit(int i, arma::mat const& rootpi_i,
                             arma::vec const& prior_mean_i) = 0;

    // Called once per iteration after all update_unit calls, before store.
    // Default no-op. Override for per-iteration blocks (e.g. NegBin alpha MH).
    virtual void post_iteration() {}

    // Per-keep serialization (taudraw for linear, alphadraw for negbin).
    virtual void store(int mkeep) = 0;

    // Post-loop: return updater-specific outputs.
    virtual UpdaterPackPayload pack() = 0;

    // Acceptance rates (aggregated).
    virtual AcceptanceMetrics acceptance_metrics() const = 0;

    // Current log-likelihood (sum over units), called at keep steps.
    virtual double current_loglike() const = 0;
};

#endif
