#include "direct_hart_backend.h"
#include "bart_pack.h"

using namespace arma;
using namespace Rcpp;

DirectHartBackend::DirectHartBackend(
    mat const& theta, mat const& Z, vec const& mubar, mat const& Amu,
    int R, int keep, List const& bart_params, List const& var_params,
    List const& phi_params)
    : HartBackend(theta, Z, R, keep, list_flag(bart_params, "store_trees", true)),
      mubar_(mubar), Amu_(Amu),
      a_tau_(as<double>(bart_params["a_tau"])),
      Hf_(as<size_t>(bart_params["num_trees"])),
      Hv_(as<size_t>(var_params["num_trees"])),
      Hphi_(phi_params.size() > 0 ? as<size_t>(phi_params["num_trees"]) : 0),
      full_cov_(phi_params.size() > 0),
      sparse_f_(as<bool>(bart_params["sparse"])),
      sparse_v_(as<bool>(var_params["sparse"])),
      sparse_phi_(full_cov_ && as<bool>(phi_params["sparse"])),
      burn_f_(sparse_f_ ? as<int>(bart_params["burn"]) : 0),
      burn_v_(sparse_v_ ? as<int>(var_params["burn"]) : 0),
      burn_phi_(sparse_phi_ ? as<int>(phi_params["burn"]) : 0)
{
    if ((int)Amu_.n_rows != D_ || (int)Amu_.n_cols != D_)
        stop("DirectHartBackend: Amu must be nvar x nvar.");
    if ((int)mubar_.n_elem != D_)
        stop("DirectHartBackend: mubar must have length nvar.");
    if (full_cov_ != (D_ > 1))
        stop("DirectHartBackend: phi_params are required iff nvar > 1.");
    vec tau0 = as<vec>(bart_params["tau0"]);
    vec v0   = as<vec>(var_params["v0"]);
    double a_v = as<double>(var_params["a_v"]);
    if ((int)tau0.n_elem != D_ || (int)v0.n_elem != D_)
        stop("DirectHartBackend: tau0 and v0 must have length nvar.");
    if (a_tau_ <= 1.0 || a_v <= 1.0)
        stop("DirectHartBackend: a_tau and a_v must exceed 1.");
    b_tau_ = (a_tau_ - 1.0) * square(tau0);

    // Starting state (app:hs_sweep): stumps at zero, s_jh = v0_j^{1/H_v},
    // tau_d^2 = b_tau_d / (a_tau + 1), mu at the coefficient means.
    mu_   = n_ > 0 ? vec(mean(theta_, 0).t()) : mubar_;
    tau2_ = b_tau_ / (a_tau_ + 1.0);
    e_.zeros(n_, D_);
    r_sq_.zeros(n_, D_);
    vmat_.zeros(n_, D_);
    Y_phi_.zeros(n_, D_ * (D_ - 1) / 2);
    Q_.assign(n_, mat(D_, D_, fill::zeros));
    rootpi_.assign(n_, mat(D_, D_, fill::zeros));

    f_ = vecbart(Hf_, D_);
    f_.setprior(as<double>(bart_params["base"]), as<double>(bart_params["power"]));
    f_.setleafguard(as<size_t>(bart_params["nmin"]));
    f_.setdata(nz_, n_, pZt(), as<size_t>(bart_params["numcut"]));
    f_.setscale(tau2_, Hf_);
    setdart_from(f_, bart_params);

    // IG(a_v, b_vj) leaves as Inv-chi^2(2 a_v, b_vj / a_v).
    v_.resize(D_);
    for (int j = 0; j < D_; ++j) {
        double leaf0 = std::pow(v0(j), 1.0 / (double)Hv_);
        double b_vj  = (a_v - 1.0) * leaf0;
        v_[j] = varbart(Hv_);
        v_[j].setprior(as<double>(var_params["base"]),
                       as<double>(var_params["power"]), 1.0);
        v_[j].setleafguard(as<size_t>(var_params["nmin"]));
        v_[j].setdata(nz_, n_, pZt(), r_sq_.colptr(j),
                      as<size_t>(var_params["numcut"]));
        v_[j].setvarprior(2.0 * a_v, b_vj / a_v, leaf0);
        setdart_from(v_[j], var_params);
    }

    if (full_cov_) {
        double sd_phi = as<double>(phi_params["tau"]) / std::sqrt((double)Hphi_);
        phi_.resize(D_);
        for (int j = 1; j < D_; ++j) {
            phi_[j] = std::vector<vcbart>(j);
            for (int k = 0; k < j; ++k) {
                vcbart& m = phi_[j][k];
                m = vcbart(Hphi_);
                m.setprior(as<double>(phi_params["base"]),
                           as<double>(phi_params["power"]), sd_phi);
                m.setleafguard(as<size_t>(phi_params["nmin"]));
                m.setdata(nz_, n_, pZt(), Y_phi_.colptr(pair_index(j, k)),
                          as<size_t>(phi_params["numcut"]));
                setdart_from(m, phi_params);
            }
        }
    }

    mu_draw_.zeros(nkeep_, D_);
    scale_draw_.zeros(nkeep_, D_);
    varcount_.zeros(nz_, nkeep_);
    varprob_.zeros(nz_, nkeep_);
    var_varcount_.zeros(nz_, D_, nkeep_);
    var_varprob_.zeros(nz_, D_, nkeep_);
    treess_.precision(10);
    treess_ << nkeep_ << " " << Hf_ << " " << nz_ << " " << D_ << std::endl;
    var_treess_.resize(D_);
    for (int j = 0; j < D_; ++j) {
        var_treess_[j].precision(10);
        var_treess_[j] << nkeep_ << " " << Hv_ << " " << nz_ << std::endl;
    }
    if (full_cov_) {
        phi_treess_.resize(D_);
        for (int j = 1; j < D_; ++j) {
            phi_treess_[j] = std::vector<std::stringstream>(j);
            for (int k = 0; k < j; ++k) {
                phi_treess_[j][k].precision(10);
                phi_treess_[j][k] << nkeep_ << " " << Hphi_ << " " << nz_
                                  << std::endl;
            }
        }
    }
    zcache_.init(bart_params, n_, D_, nkeep_, true, "DirectHartBackend");

    for (int j = 0; j < D_; ++j)
        for (int i = 0; i < n_; ++i) vmat_(i, j) = v_[j].f(i);
    refresh_residuals();
    refresh_precision();
    refresh_kernels();
}

void DirectHartBackend::draw_iteration(int rep) {
    if (sparse_f_ && rep == burn_f_ + 1) f_.startdart();
    if (sparse_v_ && rep == burn_v_ + 1)
        for (int j = 0; j < D_; ++j) v_[j].startdart();
    if (sparse_phi_ && rep == burn_phi_ + 1)
        for (int j = 1; j < D_; ++j)
            for (int k = 0; k < j; ++k) phi_[j][k].startdart();
    sweep(true, true, true, true, true);
}

void DirectHartBackend::sweep(bool mu, bool v, bool phi, bool f, bool tau) {
    if (mu) draw_mu();
    refresh_residuals();
    if (v) draw_variance_trees();
    if (phi) draw_phi_trees();
    refresh_precision();
    if (f) draw_mean_trees();
    if (tau) draw_mean_scales();
    refresh_kernels();
}

void DirectHartBackend::draw_mu() {
    mat const& f = f_.allfit();
    mat P = Amu_;
    vec b = Amu_ * mubar_;
    for (int i = 0; i < n_; ++i) {
        P += Q_[i];
        b += Q_[i] * (theta_.row(i).t() - f.row(i).t());
    }
    mat U = chol(P);                                        // U^T U = P
    vec m = solve(trimatu(U), solve(trimatl(U.t()), b));
    vec z(D_);
    for (int d = 0; d < D_; ++d) z(d) = gen_.normal();
    mu_ = m + solve(trimatu(U), z);                         // cov U^{-1} U^{-T}
}

void DirectHartBackend::refresh_residuals() {
    e_ = theta_ - f_.allfit();
    e_.each_row() -= mu_.t();
    for (int i = 0; i < n_; ++i) {
        vec r = e_.row(i).t();
        cov::apply_L(r, eval(i));                           // r = T_i e_i
        r_sq_.row(i) = square(r).t();
    }
}

void DirectHartBackend::draw_variance_trees() {
    for (int j = 0; j < D_; ++j) {
        v_[j].draw(gen_);
        for (int i = 0; i < n_; ++i) vmat_(i, j) = v_[j].f(i);
    }
}

// y^{(jk)}_i = e_ij - sum_{l<j, l!=k} phi_jl(Z_i) e_il, regressor e_ik,
// variance v_j(Z_i); phi_jk is refreshed by its ensemble after each pair.
void DirectHartBackend::draw_phi_trees() {
    for (int j = 1; j < D_; ++j) {
        for (int k = 0; k < j; ++k) {
            double* y = Y_phi_.colptr(pair_index(j, k));
            for (int i = 0; i < n_; ++i) {
                double yi = e_(i, j);
                for (int l = 0; l < j; ++l)
                    if (l != k) yi -= phi_[j][l].f(i) * e_(i, l);
                y[i] = yi;
            }
            phi_[j][k].draw(e_.colptr(k), vmat_.colptr(j), gen_);
        }
    }
}

void DirectHartBackend::refresh_precision() {
    for (int i = 0; i < n_; ++i) {
        rootpi_[i] = cov::rootpi(eval(i));                  // R R^T = Q_i
        Q_[i] = rootpi_[i] * rootpi_[i].t();
    }
}

void DirectHartBackend::draw_mean_trees() {
    f_.setscale(tau2_, Hf_);
    mat y = theta_;
    y.each_row() -= mu_.t();
    f_.draw(y, Q_, gen_);
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
    for (int j = 0; j < D_; ++j) {
        var_varcount_.slice(s).col(j) = conv_to<vec>::from(v_[j].getnv());
        var_varprob_.slice(s).col(j)  = conv_to<vec>::from(v_[j].getpv());
    }
    if (store_trees_) {
        for (size_t h = 0; h < Hf_; ++h) treess_ << f_.gettree(h);
        for (int j = 0; j < D_; ++j)
            for (size_t h = 0; h < Hv_; ++h) var_treess_[j] << v_[j].gettree(h);
        for (int j = 1; j < D_; ++j)
            for (int k = 0; k < j; ++k)
                for (size_t h = 0; h < Hphi_; ++h)
                    phi_treess_[j][k] << phi_[j][k].gettree(h);
    }
    zcache_.store_delta(mkeep, f_.allfit());
    for (uword u = 0; u < zcache_.n_unique(); ++u) {
        mat A = solve(trimatu(rootpi_[zcache_.first_idx[u]]), eye<mat>(D_, D_));
        zcache_.store_sigma(mkeep, u, A.t() * A);          // Sigma = R^{-T} R^{-1}
    }
}

BackendPackPayload DirectHartBackend::pack() {
    List out = List::create(Named("mu_draw")         = mu_draw_,
                            Named("mean_scale_draw") = scale_draw_,
                            Named("varcount")        = varcount_,
                            Named("varprob")         = varprob_,
                            Named("var_varcount")    = var_varcount_,
                            Named("var_varprob")     = var_varprob_);
    if (store_trees_) {
        out["bart_models"] = pack_vecbart_entry_(f_, treess_);
        out["var_models"]  = pack_bart_list_(D_, v_, var_treess_);
        out["phi_models"]  = full_cov_
            ? (SEXP) pack_bart_jagged_(D_, phi_, phi_treess_) : R_NilValue;
    }
    zcache_.pack(out);
    BackendPackPayload bp;
    bp.fields = out;
    bp.caps = capabilities();
    return bp;
}
