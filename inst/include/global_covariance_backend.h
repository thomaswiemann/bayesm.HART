#ifndef GLOBAL_COVARIANCE_BACKEND_H
#define GLOBAL_COVARIANCE_BACKEND_H

#include "bayesm.HART.h"
#include "covariance_backend.h"

// bayesm paths (Prior$hart_model = "none"): mixture-of-normals prior on
// theta_i with an optional linear Delta(Z) = Z Delta, drawn with rmixGibbs
// and drawDelta.  Shared by the MNL, Linear, and NegBin loop files.
class GlobalCovarianceBackend : public CovarianceBackend {
public:
    GlobalCovarianceBackend(
        arma::mat& oldbetas, arma::mat const& Z,
        arma::mat const& mubar, arma::mat const& Amu,
        double nu, arma::mat const& V, arma::vec const& a,
        arma::vec const& deltabar, arma::mat const& Ad,
        bool drawdelta,
        arma::mat olddelta, arma::vec oldprob, arma::vec ind,
        int R, int keep);

    void draw_iteration(int rep) override;
    arma::mat const& precision_root(int i) const override;
    arma::vec const& prior_mean(int i) const override;
    void store(int mkeep) override;
    BackendPackPayload pack() override;
    BackendCapabilities capabilities() const override;

protected:
    arma::mat& oldbetas_;
    Rcpp::List oldcomp_;
    arma::vec  oldprob_, ind_;
    arma::mat  olddelta_, delta_Z_;

    arma::mat  probdraw_;
    Rcpp::List compdraw_;
    arma::mat  Deltadraw_;
    std::vector<arma::vec> comp_mu_cache_;
    std::vector<arma::mat> comp_root_cache_;
    std::vector<arma::vec> unit_prior_mean_cache_;

    arma::mat const& Z_;
    arma::mat const& mubar_;
    arma::mat const& Amu_;
    arma::mat const& V_;
    arma::vec const& a_;
    arma::vec const& deltabar_;
    arma::mat const& Ad_;
    double nu_;
    bool drawdelta_;
    int nreg_, nz_, nvar_, R_, keep_;

    void refresh_component_cache();
};

#endif
