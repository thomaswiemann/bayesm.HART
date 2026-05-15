#ifndef MNL_THETA_UPDATER_H
#define MNL_THETA_UPDATER_H

#include "bayesm.HART.h"
#include "theta_updater.h"

// Forward declare (defined in rhierMnlRwMixture_rcpp_loop.cpp)
double llmnl_con(arma::vec const& betastar, arma::vec const& y,
                 arma::mat const& X, arma::vec const& SignRes);
mnlMetropOnceOut mnlMetropOnce_con(arma::vec const& y, arma::mat const& X,
    arma::vec const& oldbeta, double oldll, double s, arma::mat const& incroot,
    arma::vec const& betabar, arma::mat const& rootpi, arma::vec const& SignRes);

// MNL theta updater: RW-MH with SignRes support + incroot proposal tuning.
class MnlThetaUpdaterAdapter : public ThetaUpdater {
public:
    arma::mat& oldbetas_;
    std::vector<moments> const& lgtdata_;
    arma::vec const& SignRes_;
    arma::vec oldll_;
    int nacceptbeta_;
    double s_;
    int nlgt_, nvar_, R_;

    MnlThetaUpdaterAdapter(
        arma::mat& oldbetas,
        std::vector<moments> const& lgtdata,
        arma::vec const& SignRes,
        double s, int R);

    void update_unit(int i, arma::mat const& rootpi_i,
                     arma::vec const& prior_mean_i) override;
    void store(int mkeep) override;
    UpdaterPackPayload pack() override;
    AcceptanceMetrics acceptance_metrics() const override;
    double current_loglike() const override;
};

#endif
