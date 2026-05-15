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
    bool drawdelta, bool useBART,
    mat olddelta, vec oldprob, vec ind,
    int R, int keep, List const& bart_params)
    : oldbetas_(oldbetas), oldprob_(oldprob), ind_(ind),
      olddelta_(olddelta), Z_(Z), mubar_(mubar), Amu_(Amu),
      V_(V), a_(a), deltabar_(deltabar), Ad_(Ad),
      bart_params_(bart_params), nu_(nu),
      drawdelta_(drawdelta), useBART_(useBART),
      nreg_(oldbetas.n_rows), nz_(Z.n_cols), nvar_(V.n_cols),
      R_(R), keep_(keep)
{
    delta_Z_.zeros(nreg_, nvar_);
    probdraw_.zeros(R / keep, oldprob_.size());
    compdraw_ = List(R / keep);
    Deltadraw_ = (drawdelta && !useBART) ?
        zeros<mat>(R / keep, nz_ * nvar_) : zeros<mat>(1, 1);

    if (useBART && drawdelta) {
        mixbart_init(mbs_, nreg_, nz_, nvar_, R, keep, bart_params, oldprob_, ind_);
        std_oldbetas_ = oldbetas_;
        pstd_oldbetas_cols_.resize(nvar_);
        for (int j = 0; j < nvar_; j++)
            pstd_oldbetas_cols_[j] = std_oldbetas_.colptr(j);
        Zt_ = Z_.t();
        pZt_ = Zt_.memptr();
    }
    unit_prior_mean_cache_.resize(nreg_);
    for (int i = 0; i < nreg_; i++) {
        unit_prior_mean_cache_[i].set_size(nvar_);
    }
}

void GlobalCovarianceBackend::draw_iteration(int rep) {
    if (useBART_ && drawdelta_) {
        mixbart_draw_iter(mbs_, oldbetas_, mubar_, Amu_, nu_, V_, a_,
                          pZt_, pstd_oldbetas_cols_, bart_params_, rep, gen_);
        oldcomp_ = mbs_.oldcomp;
        oldprob_ = mbs_.oldprob;
        ind_     = mbs_.ind;
        delta_Z_ = mbs_.delta_Z;
    } else if (drawdelta_) {
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

    // Build per-unit prior means once per iteration so the driver can consume
    // const references without per-unit allocation/copy overhead.
    for (int i = 0; i < nreg_; i++) {
        int k = static_cast<int>(ind_[i]) - 1;
        unit_prior_mean_cache_[i] = comp_mu_cache_[k];
        if (drawdelta_) {
            unit_prior_mean_cache_[i] += delta_Z_.row(i).t();
        }
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
    if (useBART_ && drawdelta_) {
        mixbart_store(mbs_, mkeep);
    } else if (drawdelta_) {
        Deltadraw_(mkeep - 1, span::all) = trans(vectorise(olddelta_));
    }
}

BackendPackPayload GlobalCovarianceBackend::pack() {
    List nmix = List::create(Named("probdraw") = probdraw_,
                              Named("zdraw") = R_NilValue,
                              Named("compdraw") = compdraw_);
    BackendPackPayload bp;
    if (useBART_ && drawdelta_) {
        bp.fields = mixbart_pack(mbs_);
        bp.fields["nmix"] = nmix;
    } else if (drawdelta_) {
        bp.fields = List::create(Named("Deltadraw") = Deltadraw_,
                                  Named("nmix") = nmix);
    } else {
        bp.fields = List::create(Named("nmix") = nmix);
    }
    bp.caps = capabilities();
    return bp;
}

BackendCapabilities GlobalCovarianceBackend::capabilities() const {
    return {true, false, drawdelta_};
}
