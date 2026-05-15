#ifndef LINEAR_THETA_UPDATER_H
#define LINEAR_THETA_UPDATER_H

#include "bayesm.HART.h"
#include "theta_updater.h"

// Linear theta updater: conjugate Gaussian via runiregG + tau draw.
class LinearThetaUpdaterAdapter : public ThetaUpdater {
public:
    arma::mat& oldbetas_;
    std::vector<moments> const& regdata_;
    arma::vec tau_;
    arma::mat taudraw_;
    double nu_e_;
    arma::vec const& ssq_;
    int nreg_, nvar_, R_, keep_;

    LinearThetaUpdaterAdapter(
        arma::mat& oldbetas,
        std::vector<moments> const& regdata,
        arma::vec tau, double nu_e, arma::vec const& ssq,
        int R, int keep);

    void update_unit(int i, arma::mat const& rootpi_i,
                     arma::vec const& prior_mean_i) override;
    void store(int mkeep) override;
    UpdaterPackPayload pack() override;
    AcceptanceMetrics acceptance_metrics() const override;
    double current_loglike() const override;
};

#endif
