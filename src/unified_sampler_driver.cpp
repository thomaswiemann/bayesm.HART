#include "bayesm.HART.h"
#include "unified_sampler_driver.h"

void run_unified_sampler(
    int R, int keep, int nprint, int nreg,
    CovarianceBackend& backend,
    ThetaUpdater& updater,
    arma::mat& oldbetas,
    arma::cube& betadraw,
    arma::vec& loglike) {

    if (nprint > 0) startMcmcTimer();

    for (int rep = 0; rep < R; rep++) {

        // Step A-D: backend draws covariance/mean-tree blocks (atomic).
        backend.draw_iteration(rep);

        // Step E: per-unit theta update using backend-provided prior kernel.
        for (int i = 0; i < nreg; i++) {
            arma::mat rootpi_i    = backend.precision_root(i);
            arma::vec prior_mean_i = backend.prior_mean(i);
            updater.update_unit(i, rootpi_i, prior_mean_i);
        }

        // Per-iteration hook (e.g. NegBin alpha MH). No-op by default.
        updater.post_iteration();

        if (nprint > 0 && (rep + 1) % nprint == 0) infoMcmcTimer(rep, R);

        if ((rep + 1) % keep == 0) {
            int mkeep = (rep + 1) / keep;
            betadraw.slice(mkeep - 1) = oldbetas;
            loglike[mkeep - 1] = updater.current_loglike();
            backend.store(mkeep);
            updater.store(mkeep);
        }
    }

    if (nprint > 0) endMcmcTimer();
}
