#include "hart_backend.h"
#include "global_covariance_backend.h"
#include "normalized_hart_backend.h"
#include "hetercov_backend.h"

std::unique_ptr<CovarianceBackend> make_backend(
    bool useBART, bool useHeterCov,
    arma::mat& oldbetas, arma::mat const& Z,
    arma::mat const& mubar, arma::mat const& Amu,
    double nu, arma::mat const& V, arma::vec const& a,
    arma::vec const& deltabar, arma::mat const& Ad, bool drawdelta,
    arma::mat const& olddelta, arma::vec const& oldprob, arma::vec const& ind,
    int R, int keep,
    Rcpp::List const& bart_params, Rcpp::List const& var_params,
    Rcpp::List const& phi_params)
{
    if (useHeterCov)
        return std::make_unique<HeterCovBackend>(
            oldbetas, Z, mubar, Amu, R, keep,
            bart_params, var_params, phi_params);
    if (useBART)
        return std::make_unique<NormalizedHartBackend>(
            oldbetas, Z, arma::vec(mubar.row(0).t()), Amu, nu, V, R, keep,
            bart_params);
    return std::make_unique<GlobalCovarianceBackend>(
        oldbetas, Z, mubar, Amu, nu, V, a, deltabar, Ad, drawdelta,
        olddelta, oldprob, ind, R, keep);
}
