#ifndef HETERCOV_BACKEND_H
#define HETERCOV_BACKEND_H

#include "hart_backend.h"
#include "hetercov_block.h"

// Covariate-dependent covariance backend on the pre-plan modified-Cholesky
// kernel (hetercov_block).  Selected by Prior$hart_model = "direct" until
// DirectHartBackend replaces it in stage 4.  Shared by MNL, Linear, and
// NegBin loop files.
class HeterCovBackend : public HartBackend {
public:
    HeterCovState hcs_;
    arma::mat std_oldbetas_;
    std::vector<double*> pstd_oldbetas_cols_;
    std::vector<arma::mat> rootpi_cache_;

    arma::mat const& mubar_;
    arma::mat const& Amu_;
    Rcpp::List const& bart_params_;
    Rcpp::List const& var_params_;
    Rcpp::List const& phi_params_;

    HeterCovBackend(
        arma::mat& oldbetas, arma::mat const& Z,
        arma::mat const& mubar, arma::mat const& Amu,
        int R, int keep,
        Rcpp::List const& bart_params,
        Rcpp::List const& var_params,
        Rcpp::List const& phi_params);

    void draw_iteration(int rep) override;
    arma::mat const& precision_root(int i) const override;
    void store(int mkeep) override;
    BackendPackPayload pack() override;
    BackendCapabilities capabilities() const override;
};

#endif
