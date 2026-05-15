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
    : oldbetas_(oldbetas), mubar_(mubar), Amu_(Amu),
      bart_params_(bart_params), var_params_(var_params),
      phi_params_(phi_params),
      nreg_(oldbetas.n_rows), nz_(Z.n_cols), nvar_(oldbetas.n_cols)
{
    hetercov_init(hcs_, nreg_, nz_, nvar_, R, keep,
                  bart_params, var_params, phi_params);
    std_oldbetas_ = oldbetas_;
    pstd_oldbetas_cols_.resize(nvar_);
    for (int j = 0; j < nvar_; j++)
        pstd_oldbetas_cols_[j] = std_oldbetas_.colptr(j);
    Zt_ = Z.t();
    pZt_ = Zt_.memptr();
    rootpi_cache_.resize(nreg_);
    prior_mean_cache_.resize(nreg_);
    for (int i = 0; i < nreg_; i++) {
        rootpi_cache_[i].set_size(nvar_, nvar_);
        prior_mean_cache_[i].set_size(nvar_);
    }
}

void HeterCovBackend::draw_iteration(int rep) {
    hetercov_draw_iter(hcs_, oldbetas_, mubar_.row(0).t(), Amu_,
                       pZt_, pstd_oldbetas_cols_,
                       bart_params_, var_params_, phi_params_,
                       rep, gen_);

    // Cache per-unit kernels once per iteration for reference-based access
    // from the unified sampler driver.
    for (int i = 0; i < nreg_; i++) {
        cov::cov_eval ev{hcs_.var_models, hcs_.phi_models, (size_t)i};
        rootpi_cache_[i] = cov::rootpi(ev);
        prior_mean_cache_[i] = hcs_.mu_post + hcs_.delta_Z.row(i).t();
    }
}

mat const& HeterCovBackend::precision_root(int i) const {
    return rootpi_cache_[i];
}

vec const& HeterCovBackend::prior_mean(int i) const {
    return prior_mean_cache_[i];
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
