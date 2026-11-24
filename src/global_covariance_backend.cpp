#include "global_covariance_backend.h"

using namespace arma;
using namespace Rcpp;

void GlobalCovarianceBackend::refresh_component_cache() {
    int ncomp = oldcomp_.size();
    comp_mu_cache_.resize(ncomp);
    comp_root_cache_.resize(ncomp);
    for (int k = 0; k < ncomp; k++) {
        List comp = oldcomp_[k];
        comp_mu_cache_[k] = as<vec>(comp[0]);
        comp_root_cache_[k] = as<mat>(comp[1]);
    }
}

GlobalCovarianceBackend::GlobalCovarianceBackend(
    mat& oldbetas, mat const& Z,
    mat const& mubar, mat const& Amu,
    double nu, mat const& V, vec const& a,
    vec const& deltabar, mat const& Ad,
    bool drawdelta,
    mat olddelta, vec oldprob, vec ind,
    int R, int keep)
    : oldbetas_(oldbetas), oldprob_(oldprob), ind_(ind),
      olddelta_(olddelta), Z_(Z), mubar_(mubar), Amu_(Amu),
      V_(V), a_(a), deltabar_(deltabar), Ad_(Ad),
      nu_(nu), drawdelta_(drawdelta),
      nreg_(oldbetas.n_rows), nz_(Z.n_cols), nvar_(V.n_cols),
      R_(R), keep_(keep)
{
    delta_Z_.zeros(nreg_, nvar_);
    probdraw_.zeros(R / keep, oldprob_.size());
    compdraw_ = List(R / keep);
    Deltadraw_ = drawdelta ? zeros<mat>(R / keep, nz_ * nvar_) : zeros<mat>(1, 1);
    unit_prior_mean_cache_.resize(nreg_);
    for (int i = 0; i < nreg_; i++) unit_prior_mean_cache_[i].set_size(nvar_);
}

void GlobalCovarianceBackend::draw_iteration(int /*rep*/) {
    if (drawdelta_) {
        olddelta_.reshape(nvar_, nz_);
        List mgout = rmixGibbs(oldbetas_ - delta_Z_, mubar_, Amu_,
                               nu_, V_, a_, oldprob_, ind_);
        oldcomp_ = mgout["comps"];
        oldprob_ = as<vec>(mgout["p"]);
        ind_     = as<vec>(mgout["z"]);
        olddelta_ = drawDelta(Z_, oldbetas_, ind_, oldcomp_, deltabar_, Ad_);
        olddelta_.reshape(nvar_, nz_);
        delta_Z_ = Z_ * trans(olddelta_);
    } else {
        List mgout = rmixGibbs(oldbetas_, mubar_, Amu_,
                               nu_, V_, a_, oldprob_, ind_);
        oldcomp_ = mgout["comps"];
        oldprob_ = as<vec>(mgout["p"]);
        ind_     = as<vec>(mgout["z"]);
    }
    refresh_component_cache();

    // Per-unit prior means read by the driver.
    for (int i = 0; i < nreg_; i++) {
        int k = static_cast<int>(ind_[i]) - 1;
        unit_prior_mean_cache_[i] = comp_mu_cache_[k];
        if (drawdelta_) unit_prior_mean_cache_[i] += delta_Z_.row(i).t();
    }
}

mat const& GlobalCovarianceBackend::precision_root(int i) const {
    int k = static_cast<int>(ind_[i]) - 1;
    return comp_root_cache_[k];
}

vec const& GlobalCovarianceBackend::prior_mean(int i) const {
    return unit_prior_mean_cache_[i];
}

void GlobalCovarianceBackend::store(int mkeep) {
    probdraw_(mkeep - 1, span::all) = trans(oldprob_);
    compdraw_[mkeep - 1] = oldcomp_;
    if (drawdelta_) Deltadraw_(mkeep - 1, span::all) = trans(vectorise(olddelta_));
}

BackendPackPayload GlobalCovarianceBackend::pack() {
    List nmix = List::create(Named("probdraw") = probdraw_,
                             Named("zdraw") = R_NilValue,
                             Named("compdraw") = compdraw_);
    BackendPackPayload bp;
    bp.fields = drawdelta_
        ? List::create(Named("Deltadraw") = Deltadraw_, Named("nmix") = nmix)
        : List::create(Named("nmix") = nmix);
    bp.caps = capabilities();
    return bp;
}

BackendCapabilities GlobalCovarianceBackend::capabilities() const {
    return {true, false, drawdelta_};
}

List GlobalCovarianceBackend::checkpoint() {
    return List::create(Named("delta") = olddelta_,
                        Named("prob") = NumericVector(oldprob_.begin(),
                                                       oldprob_.end()),
                        Named("ind") = NumericVector(ind_.begin(), ind_.end()));
}

void GlobalCovarianceBackend::restore(List const& state) {
    olddelta_ = as<mat>(state["delta"]);
    oldprob_ = as<vec>(state["prob"]);
    ind_ = as<vec>(state["ind"]);
    if ((int)ind_.n_elem != nreg_ || oldprob_.n_elem != a_.n_elem)
        stop("GlobalCovarianceBackend: state dimension mismatch.");
    if (drawdelta_) {
        if ((int)olddelta_.n_rows != nvar_ || (int)olddelta_.n_cols != nz_)
            stop("GlobalCovarianceBackend: delta dimension mismatch.");
        delta_Z_ = Z_ * olddelta_.t();
    }
}
