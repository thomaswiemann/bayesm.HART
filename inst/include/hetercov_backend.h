#ifndef HETERCOV_BACKEND_H
#define HETERCOV_BACKEND_H

#include "bayesm.HART.h"
#include "hetercov_block.h"
#include "covariance_backend.h"

// Model-agnostic covariate-dependent covariance backend.
// Wraps hetercov_block behind the CovarianceBackend interface.
// Shared by MNL, Linear, and NegBin loop files.
class HeterCovBackend : public CovarianceBackend {
public:
    arma::mat& oldbetas_;
    HeterCovState hcs_;
    arma::mat std_oldbetas_;
    std::vector<double*> pstd_oldbetas_cols_;
    arma::mat Zt_;
    double* pZt_;
    arn gen_;
    std::vector<arma::mat> rootpi_cache_;
    std::vector<arma::vec> prior_mean_cache_;

    arma::mat const& mubar_;
    arma::mat const& Amu_;
    Rcpp::List const& bart_params_;
    Rcpp::List const& var_params_;
    Rcpp::List const& phi_params_;
    int nreg_, nz_, nvar_;

    HeterCovBackend(
        arma::mat& oldbetas, arma::mat const& Z,
        arma::mat const& mubar, arma::mat const& Amu,
        int R, int keep,
        Rcpp::List const& bart_params,
        Rcpp::List const& var_params,
        Rcpp::List const& phi_params);

    void draw_iteration(int rep) override;
    arma::mat const& precision_root(int i) const override;
    arma::vec const& prior_mean(int i) const override;
    void store(int mkeep) override;
    BackendPackPayload pack() override;
    BackendCapabilities capabilities() const override;
};

#endif
