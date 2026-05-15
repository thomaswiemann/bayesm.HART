#ifndef NEGBIN_THETA_UPDATER_H
#define NEGBIN_THETA_UPDATER_H

#include "bayesm.HART.h"
#include "theta_updater.h"

// Forward declare (defined in rhierNegbinRw_rcpp_loop.cpp)
double llnegbinpooled(std::vector<moments> const& regdata_vector,
                      arma::mat const& oldbetas,
                      double alpha);

// NegBin theta updater: RW-MH for beta_i + alpha MH block.
class NegbinThetaUpdaterAdapter : public ThetaUpdater {
public:
    arma::mat& oldbetas_;
    std::vector<moments> const& regdata_;
    arma::vec oldlpostbeta_, clpostbeta_;
    int nacceptbeta_, nacceptalpha_;
    double sbeta_, alpha_, alphacroot_;
    double a_prior_, b_prior_;
    bool fixalpha_;
    arma::vec alphadraw_;
    int nreg_, nvar_, R_, keep_;

    NegbinThetaUpdaterAdapter(
        arma::mat& oldbetas,
        std::vector<moments> const& regdata,
        double sbeta, double alpha, double alphacroot,
        double a_prior, double b_prior, bool fixalpha,
        int R, int keep);

    void update_unit(int i, arma::mat const& rootpi_i,
                     arma::vec const& prior_mean_i) override;
    void store(int mkeep) override;
    UpdaterPackPayload pack() override;
    AcceptanceMetrics acceptance_metrics() const override;
    double current_loglike() const override;
    void post_iteration() override { draw_alpha(); }

    // Alpha MH step (called via post_iteration after per-unit loop)
    void draw_alpha();
};

#endif
