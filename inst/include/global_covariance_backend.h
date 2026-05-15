#ifndef GLOBAL_COVARIANCE_BACKEND_H
#define GLOBAL_COVARIANCE_BACKEND_H

#include "bayesm.HART.h"
#include "mixbart_block.h"
#include "covariance_backend.h"

// Model-agnostic global covariance backend.
// Wraps mixbart_block + rmixGibbs + drawDelta behind the CovarianceBackend
// interface.  Shared by MNL, Linear, and NegBin loop files.
class GlobalCovarianceBackend : public CovarianceBackend {
public:
    arma::mat& oldbetas_;
    Rcpp::List oldcomp_;
    arma::vec  oldprob_, ind_;
    arma::mat  olddelta_, delta_Z_;

    MixBartState mbs_;
    arma::mat std_oldbetas_;
    std::vector<double*> pstd_oldbetas_cols_;
    arma::mat Zt_;
    double* pZt_;
    arn gen_;

    arma::mat probdraw_;
    Rcpp::List compdraw_;
    arma::mat Deltadraw_;
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
    Rcpp::List const& bart_params_;
    double nu_;
    bool drawdelta_, useBART_;
    int nreg_, nz_, nvar_, R_, keep_;

    GlobalCovarianceBackend(
        arma::mat& oldbetas, arma::mat const& Z,
        arma::mat const& mubar, arma::mat const& Amu,
        double nu, arma::mat const& V, arma::vec const& a,
        arma::vec const& deltabar, arma::mat const& Ad,
        bool drawdelta, bool useBART,
        arma::mat olddelta, arma::vec oldprob, arma::vec ind,
        int R, int keep, Rcpp::List const& bart_params);

    void draw_iteration(int rep) override;
    arma::mat const& precision_root(int i) const override;
    arma::vec const& prior_mean(int i) const override;
    void store(int mkeep) override;
    BackendPackPayload pack() override;
    BackendCapabilities capabilities() const override;
    void refresh_component_cache();
};

#endif
