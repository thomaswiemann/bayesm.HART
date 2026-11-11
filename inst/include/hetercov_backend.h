#ifndef HETERCOV_BACKEND_H
#define HETERCOV_BACKEND_H

#include "hart_backend.h"
#include "hetercov_block.h"

// Direct HART (Prior$hart_model = "direct") on the modified-Cholesky kernel
// in hetercov_block.  Shared by the MNL, Linear, and NegBin loop files.
class HeterCovBackend : public HartBackend {
public:
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

protected:
    HeterCovState hcs_;
    arma::mat std_oldbetas_;
    std::vector<double*> pstd_oldbetas_cols_;
    std::vector<arma::mat> rootpi_cache_;

    arma::mat const& mubar_;
    arma::mat const& Amu_;
    Rcpp::List const& bart_params_;
    Rcpp::List const& var_params_;
    Rcpp::List const& phi_params_;
};

#endif
