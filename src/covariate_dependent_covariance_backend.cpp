#include "hetercov_backend.h"

using namespace arma;
using namespace Rcpp;

HeterCovBackend::HeterCovBackend(
    mat& oldbetas, mat const& Z,
    mat const& mubar, mat const& Amu,
    int R, int keep,
    List const& bart_params,
    List const& var_params,
    List const& phi_params)
    : HartBackend(oldbetas, Z, R, keep, list_flag(bart_params, "store_trees", true)),
      mubar_(mubar), Amu_(Amu),
      bart_params_(bart_params), var_params_(var_params),
      phi_params_(phi_params)
{
    hetercov_init(hcs_, n_, nz_, D_, R, keep,
                  bart_params, var_params, phi_params);
    std_oldbetas_ = theta_;
    pstd_oldbetas_cols_.resize(D_);
    for (int j = 0; j < D_; j++)
        pstd_oldbetas_cols_[j] = std_oldbetas_.colptr(j);
    rootpi_cache_.assign(n_, mat(D_, D_, fill::zeros));
}

void HeterCovBackend::draw_iteration(int rep) {
    hetercov_draw_iter(hcs_, theta_, mubar_.row(0).t(), Amu_,
                       pZt(), pstd_oldbetas_cols_,
                       bart_params_, var_params_, phi_params_,
                       rep, gen_);

    // Per-unit kernels read by the driver.
    for (int i = 0; i < n_; i++) {
        cov::cov_eval ev{hcs_.var_models, hcs_.phi_models, (size_t)i};
        rootpi_cache_[i] = cov::rootpi(ev);
        prior_mean_[i]   = hcs_.mu_post + hcs_.delta_Z.row(i).t();
    }
}

mat const& HeterCovBackend::precision_root(int i) const {
    return rootpi_cache_[i];
}

void HeterCovBackend::store(int mkeep) {
    hetercov_store(hcs_, mkeep);
}

BackendPackPayload HeterCovBackend::pack() {
    BackendPackPayload bp;
    bp.fields = hetercov_pack(hcs_);
    bp.caps = capabilities();
    return bp;
}

BackendCapabilities HeterCovBackend::capabilities() const {
    return {false, true, true};
}
