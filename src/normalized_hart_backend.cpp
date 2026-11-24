#include "normalized_hart_backend.h"
#include "sigma_gibbs.h"
#include "bart_pack.h"

using namespace arma;
using namespace Rcpp;

NormalizedHartBackend::NormalizedHartBackend(
    mat const& theta, mat const& Z, vec const& mubar, mat const& Amu,
    double nu, mat const& Psi, int R, int keep, List const& bart_params)
    : HartBackend(theta, Z, R, keep, list_flag(bart_params, "store_trees", true)),
      mubar_(mubar), Amu_(Amu), Psi_(Psi), nu_(nu), bart_params_(bart_params),
      sparse_(list_flag(bart_params, "sparse", false)),
      burn_(sparse_ ? as<int>(bart_params["burn"]) : 0)
{
    if ((int)Amu_.n_rows != D_ || (int)Amu_.n_cols != D_)
        stop("NormalizedHartBackend: Amu must be nvar x nvar.");
    if ((int)mubar_.n_elem != D_)
        stop("NormalizedHartBackend: mubar must have length nvar.");

    // Start: mu at the coefficient means, Sigma = Psi / nu, stump trees.
    mu_ = n_ > 0 ? vec(mean(theta_, 0).t()) : mubar_;
    rooti_ = solve(trimatu(chol(Psi_ / nu_)), eye<mat>(D_, D_));
    L_ = rooti_.t();

    theta_std_.zeros(n_, D_);
    theta_std_cols_.resize(D_);
    for (int d = 0; d < D_; ++d) theta_std_cols_[d] = theta_std_.colptr(d);
    f_ = initializeBART(pZt(), n_, nz_, theta_std_cols_, bart_params_, gen_);
    fstar_.zeros(n_, D_);
    delta_Z_.zeros(n_, D_);

    probdraw_.ones(nkeep_, 1);
    compdraw_ = List(nkeep_);
    size_t m = as<size_t>(bart_params_["num_trees"]);
    treess_.resize(D_);
    for (int d = 0; d < D_; ++d) {
        treess_[d].precision(10);
        treess_[d] << nkeep_ << " " << m << " " << nz_ << std::endl;
    }
    varcount_.zeros(nz_, D_, nkeep_);
    varprob_.zeros(nz_, D_, nkeep_);
    zcache_.init(bart_params_, n_, D_, nkeep_, false, "NormalizedHartBackend");

    refresh_kernels();
}

void NormalizedHartBackend::draw_iteration(int rep) {
    draw_mu();
    draw_cholesky();
    draw_trees(rep);
    refresh_kernels();
}

void NormalizedHartBackend::draw_mu() {
    draw_mu_normalized(mu_, L_, theta_, fstar_, mubar_, Amu_, gen_);
}

void NormalizedHartBackend::draw_cholesky() {
    draw_cholesky_sweep(L_, theta_, fstar_, mu_, Psi_, nu_, gen_);
}

void NormalizedHartBackend::draw_trees(int rep) {
    // Whitened responses theta~_i = L (theta_i - mu), copied in place so the
    // ensembles' y pointers stay valid.
    mat ts = (theta_.each_row() - mu_.t()) * L_.t();
    std::copy(ts.begin(), ts.end(), theta_std_.begin());
    if (sparse_ && rep == burn_ + 1)
        for (int d = 0; d < D_; ++d) f_[d].startdart();
    for (int d = 0; d < D_; ++d) {
        f_[d].draw(1.0, gen_);
        for (int i = 0; i < n_; ++i) fstar_(i, d) = f_[d].f(i);
    }
}

void NormalizedHartBackend::refresh_kernels() {
    rooti_ = L_.t();
    delta_Z_ = solve(trimatl(L_), fstar_.t()).t();   // rows (L^{-1} f*_i)^T
    for (int i = 0; i < n_; ++i) prior_mean_[i] = mu_ + delta_Z_.row(i).t();
}

void NormalizedHartBackend::store(int mkeep) {
    compdraw_[mkeep - 1] = List::create(List::create(
        Named("mu")    = NumericVector(mu_.begin(), mu_.end()),
        Named("rooti") = rooti_));
    for (int d = 0; d < D_; ++d) {
        if (store_trees_)
            for (size_t j = 0; j < f_[d].getm(); ++j) treess_[d] << f_[d].gettree(j);
        varcount_.slice(mkeep - 1).col(d) = conv_to<vec>::from(f_[d].getnv());
        varprob_.slice(mkeep - 1).col(d)  = conv_to<vec>::from(f_[d].getpv());
    }
    zcache_.store_delta(mkeep, delta_Z_);
}

BackendPackPayload NormalizedHartBackend::pack() {
    List nmix = List::create(Named("probdraw") = probdraw_,
                             Named("zdraw")    = R_NilValue,
                             Named("compdraw") = compdraw_);
    List out = List::create(Named("varcount") = varcount_,
                            Named("varprob")  = varprob_,
                            Named("nmix")     = nmix);
    if (store_trees_) {
        std::vector<std::stringstream> tss(D_);
        for (int d = 0; d < D_; ++d) tss[d] << treess_[d].str();
        out["bart_models"] = pack_bart_list_(D_, f_, tss);
    }
    zcache_.pack(out);
    BackendPackPayload bp;
    bp.fields = out;
    bp.caps = capabilities();
    return bp;
}

List NormalizedHartBackend::checkpoint() {
    List mean(D_);
    for (int d = 0; d < D_; ++d) mean[d] = pack_bart_state_(f_[d]);
    return List::create(Named("mu") = NumericVector(mu_.begin(), mu_.end()),
                        Named("L") = L_, Named("mean") = mean);
}

void NormalizedHartBackend::restore(List const& state) {
    mu_ = as<vec>(state["mu"]);
    L_ = as<mat>(state["L"]);
    List mean = state["mean"];
    if ((int)mu_.n_elem != D_ || (int)L_.n_rows != D_ ||
        (int)L_.n_cols != D_ || mean.size() != D_)
        stop("NormalizedHartBackend: state dimension mismatch.");
    for (int d = 0; d < D_; ++d) {
        restore_bart_state_(f_[d], as<List>(mean[d]));
        for (int i = 0; i < n_; ++i) fstar_(i, d) = f_[d].f(i);
    }
    refresh_kernels();
}
