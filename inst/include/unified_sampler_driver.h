#ifndef UNIFIED_SAMPLER_DRIVER_H
#define UNIFIED_SAMPLER_DRIVER_H

#include <RcppArmadillo.h>
#include "covariance_backend.h"
#include "theta_updater.h"

// Run the unified MCMC loop.
//
// The driver owns: rep loop, betadraw/loglike storage, keep/store cadence.
// The driver does NOT own: covariance algebra (backend), likelihood logic (updater).
//
// oldbetas is passed by mutable reference: updater writes, backend reads.
// betadraw and loglike are filled at keep steps.
//
// After this returns, the caller handles output assembly (Phase 1: existing
// code; Phase 4+: SamplerOutputBuilder).
void run_unified_sampler(
    int R, int keep, int nprint, int nreg,
    CovarianceBackend& backend,
    ThetaUpdater& updater,
    arma::mat& oldbetas,
    arma::cube& betadraw,
    arma::vec& loglike);

#endif
