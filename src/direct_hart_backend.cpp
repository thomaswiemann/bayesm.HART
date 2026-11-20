#include "direct_hart_backend.h"
#include "bart_pack.h"

using namespace arma;
using namespace Rcpp;

DirectHartBackend::DirectHartBackend(
    mat const& theta, mat const& Z, vec const& mubar, mat const& Amu,
    int R, int keep, List const& bart_params, std::unique_ptr<DirectCov> cov)
    : HartBackend(theta, Z, R, keep, list_flag(bart_params, "store_trees", true)),
      mubar_(mubar), Amu_(Amu),
      a_tau_(as<double>(bart_params["a_tau"])),
      Hf_(as<size_t>(bart_params["num_trees"])),
      sparse_f_(as<bool>(bart_params["sparse"])),
      burn_f_(sparse_f_ ? as<int>(bart_params["burn"]) : 0),
      cov_(std::move(cov))
{
    if ((int)Amu_.n_rows != D_ || (int)Amu_.n_cols != D_)
        stop("DirectHartBackend: Amu must be nvar x nvar.");
    if ((int)mubar_.n_elem != D_)
        stop("DirectHartBackend: mubar must have length nvar.");
    vec tau0 = as<vec>(bart_params["tau0"]);
    if ((int)tau0.n_elem != D_)
        stop("DirectHartBackend: tau0 must have length nvar.");
    if (a_tau_ <= 1.0)
        stop("DirectHartBackend: a_tau must exceed 1.");
    b_tau_ = (a_tau_ - 1.0) * square(tau0);

    // Starting state (app:hs_sweep): stumps at zero,
    // tau_d^2 = b_tau_d / (a_tau + 1), mu at the coefficient means.
    mu_   = n_ > 0 ? vec(mean(theta_, 0).t()) : mubar_;
    tau2_ = b_tau_ / (a_tau_ + 1.0);
    E_.zeros(n_, D_);

    f_ = vecbart(Hf_, D_);
    f_.setprior(as<double>(bart_params["base"]), as<double>(bart_params["power"]));
    f_.setleafguard(as<size_t>(bart_params["nmin"]));
    f_.setdata(nz_, n_, pZt(), as<size_t>(bart_params["numcut"]));
    f_.setscale(tau2_, Hf_);
    setdart_from(f_, bart_params);

    mu_draw_.zeros(nkeep_, D_);
    scale_draw_.zeros(nkeep_, D_);
    varcount_.zeros(nz_, nkeep_);
    varprob_.zeros(nz_, nkeep_);
    treess_.precision(10);
    treess_ << nkeep_ << " " << Hf_ << " " << nz_ << " " << D_ << std::endl;

    cov_->init(n_, D_, nz_, pZt(), nkeep_);
    zcache_.init(bart_params, n_, D_, nkeep_, cov_->covariate_sigma(),
                 "DirectHartBackend");
    refresh_residuals();
    cov_->start(E_);
    refresh_kernels();
}

void DirectHartBackend::draw_iteration(int rep) {
    if (sparse_f_ && rep == burn_f_ + 1) f_.startdart();
    cov_->start_dart(rep);
    sweep(true, true, true, true);
}

void DirectHartBackend::sweep(bool mu, bool cov, bool f, bool tau) {
    if (mu) draw_mu();
    refresh_residuals();
    if (cov) cov_->draw(E_, gen_);
    if (f) draw_mean_trees();
    if (tau) draw_mean_scales();
    refresh_kernels();
}

void DirectHartBackend::draw_mu() {
    std::vector<mat> const& Q = cov_->precisions();
    mat const& f = f_.allfit();
    mat P = Amu_;
    vec b = Amu_ * mubar_;
    for (int i = 0; i < n_; ++i) {
        P += Q[i];
        b += Q[i] * (theta_.row(i).t() - f.row(i).t());
    }
    mat U = chol(P);                                        // U^T U = P
    vec m = solve(trimatu(U), solve(trimatl(U.t()), b));
    vec z(D_);
    for (int d = 0; d < D_; ++d) z(d) = gen_.normal();
    mu_ = m + solve(trimatu(U), z);                         // cov U^{-1} U^{-T}
}

void DirectHartBackend::refresh_residuals() {
    E_ = theta_ - f_.allfit();
    E_.each_row() -= mu_.t();
}

void DirectHartBackend::draw_mean_trees() {
    f_.setscale(tau2_, Hf_);
    mat y = theta_;
    y.each_row() -= mu_.t();
    f_.draw(y, cov_->precisions(), gen_);
}

void DirectHartBackend::draw_mean_scales() {
    double shape = a_tau_ + 0.5 * (double)f_.nleaves();
    vec ss = f_.leaf_sumsq();
    for (int d = 0; d < D_; ++d)
        tau2_(d) = 1.0 / gen_.gamma(shape, b_tau_(d) + 0.5 * (double)Hf_ * ss(d));
}

void DirectHartBackend::refresh_kernels() {
    mat const& f = f_.allfit();
    for (int i = 0; i < n_; ++i) prior_mean_[i] = mu_ + f.row(i).t();
}

void DirectHartBackend::store(int mkeep) {
    int s = mkeep - 1;
    mu_draw_.row(s)    = mu_.t();
    scale_draw_.row(s) = tau2_.t();
    varcount_.col(s)   = conv_to<vec>::from(f_.getnv());
    varprob_.col(s)    = conv_to<vec>::from(f_.getpv());
    if (store_trees_)
        for (size_t h = 0; h < Hf_; ++h) treess_ << f_.gettree(h);
    zcache_.store_delta(mkeep, f_.allfit());
    cov_->store(mkeep, mu_, store_trees_, zcache_);
}

BackendPackPayload DirectHartBackend::pack() {
    List out = List::create(Named("mu_draw")         = mu_draw_,
                            Named("mean_scale_draw") = scale_draw_,
                            Named("varcount")        = varcount_,
                            Named("varprob")         = varprob_);
    if (store_trees_) out["bart_models"] = pack_vecbart_entry_(f_, treess_);
    cov_->pack(out, store_trees_);
    zcache_.pack(out);
    BackendPackPayload bp;
    bp.fields = out;
    bp.caps = capabilities();
    return bp;
}

List DirectHartBackend::checkpoint() {
    return List::create(Named("mu") = NumericVector(mu_.begin(), mu_.end()),
                        Named("tau2") = NumericVector(tau2_.begin(),
                                                       tau2_.end()),
                        Named("mean") = pack_bart_state_(f_),
                        Named("covariance") = cov_->checkpoint());
}

void DirectHartBackend::restore(List const& state) {
    mu_ = as<vec>(state["mu"]);
    tau2_ = as<vec>(state["tau2"]);
    if ((int)mu_.n_elem != D_ || (int)tau2_.n_elem != D_)
        stop("DirectHartBackend: state dimension mismatch.");
    restore_bart_state_(f_, as<List>(state["mean"]));
    f_.setscale(tau2_, Hf_);
    cov_->restore(as<List>(state["covariance"]));
    refresh_residuals();
    refresh_kernels();
}
