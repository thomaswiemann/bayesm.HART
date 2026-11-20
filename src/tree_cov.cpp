#include "tree_cov.h"
#include "bart_pack.h"
#include "hart_backend.h"

using namespace arma;
using namespace Rcpp;

TreeCov::TreeCov(List const& var_params, List const& phi_params)
    : var_params_(var_params), phi_params_(phi_params),
      n_(0), D_(0), nz_(0),
      full_cov_(phi_params.size() > 0),
      Hv_(as<size_t>(var_params["num_trees"])),
      Hphi_(full_cov_ ? as<size_t>(phi_params["num_trees"]) : 0),
      sparse_v_(as<bool>(var_params["sparse"])),
      sparse_phi_(full_cov_ && as<bool>(phi_params["sparse"])),
      burn_v_(sparse_v_ ? as<int>(var_params["burn"]) : 0),
      burn_phi_(sparse_phi_ ? as<int>(phi_params["burn"]) : 0)
{}

void TreeCov::init(int n, int D, int nz, double* pZt, int nkeep) {
    n_ = n; D_ = D; nz_ = nz;
    if (full_cov_ != (D_ > 1))
        stop("TreeCov: phi_params are required iff nvar > 1.");
    vec v0     = as<vec>(var_params_["v0"]);
    double a_v = as<double>(var_params_["a_v"]);
    if ((int)v0.n_elem != D_) stop("TreeCov: v0 must have length nvar.");
    if (a_v <= 1.0) stop("TreeCov: a_v must exceed 1.");

    e_.zeros(n_, D_);
    r_sq_.zeros(n_, D_);
    vmat_.zeros(n_, D_);
    Y_phi_.zeros(n_, D_ * (D_ - 1) / 2);
    Q_.assign(n_, mat(D_, D_, fill::zeros));
    rootpi_.assign(n_, mat(D_, D_, fill::zeros));

    // IG(a_v, b_vj) leaves as Inv-chi^2(2 a_v, b_vj / a_v).
    v_.resize(D_);
    for (int j = 0; j < D_; ++j) {
        double leaf0 = std::pow(v0(j), 1.0 / (double)Hv_);
        double b_vj  = (a_v - 1.0) * leaf0;
        v_[j] = varbart(Hv_);
        v_[j].setprior(as<double>(var_params_["base"]),
                       as<double>(var_params_["power"]), 1.0);
        v_[j].setleafguard(as<size_t>(var_params_["nmin"]));
        v_[j].setdata(nz_, n_, pZt, r_sq_.colptr(j),
                      as<size_t>(var_params_["numcut"]));
        v_[j].setvarprior(2.0 * a_v, b_vj / a_v, leaf0);
        setdart_from(v_[j], var_params_);
    }

    if (full_cov_) {
        double sd_phi = as<double>(phi_params_["tau"]) / std::sqrt((double)Hphi_);
        phi_.resize(D_);
        for (int j = 1; j < D_; ++j) {
            phi_[j] = std::vector<vcbart>(j);
            for (int k = 0; k < j; ++k) {
                vcbart& m = phi_[j][k];
                m = vcbart(Hphi_);
                m.setprior(as<double>(phi_params_["base"]),
                           as<double>(phi_params_["power"]), sd_phi);
                m.setleafguard(as<size_t>(phi_params_["nmin"]));
                m.setdata(nz_, n_, pZt, Y_phi_.colptr(pair_index(j, k)),
                          as<size_t>(phi_params_["numcut"]));
                setdart_from(m, phi_params_);
            }
        }
    }

    var_varcount_.zeros(nz_, D_, nkeep);
    var_varprob_.zeros(nz_, D_, nkeep);
    var_treess_.resize(D_);
    for (int j = 0; j < D_; ++j) {
        var_treess_[j].precision(10);
        var_treess_[j] << nkeep << " " << Hv_ << " " << nz_ << std::endl;
    }
    if (full_cov_) {
        phi_treess_.resize(D_);
        for (int j = 1; j < D_; ++j) {
            phi_treess_[j] = std::vector<std::stringstream>(j);
            for (int k = 0; k < j; ++k) {
                phi_treess_[j][k].precision(10);
                phi_treess_[j][k] << nkeep << " " << Hphi_ << " " << nz_
                                  << std::endl;
            }
        }
    }
}

// Stumps at s_jh = v0_j^{1/H_v}, phi = 0 (app:hs_sweep initialization).
void TreeCov::start(mat const& E) {
    for (int j = 0; j < D_; ++j)
        for (int i = 0; i < n_; ++i) vmat_(i, j) = v_[j].f(i);
    refresh_residuals(E);
    refresh_precision();
}

void TreeCov::start_dart(int rep) {
    if (sparse_v_ && rep == burn_v_ + 1)
        for (int j = 0; j < D_; ++j) v_[j].startdart();
    if (sparse_phi_ && rep == burn_phi_ + 1)
        for (int j = 1; j < D_; ++j)
            for (int k = 0; k < j; ++k) phi_[j][k].startdart();
}

void TreeCov::draw(mat const& E, rn& gen) {
    refresh_residuals(E);
    draw_variance_trees(gen);
    draw_phi_trees(gen);
    refresh_precision();
}

void TreeCov::refresh_residuals(mat const& E) {
    e_ = E;
    for (int i = 0; i < n_; ++i) {
        vec r = e_.row(i).t();
        cov::apply_L(r, eval(i));                           // r = T_i e_i
        r_sq_.row(i) = square(r).t();
    }
}

void TreeCov::draw_variance_trees(rn& gen) {
    for (int j = 0; j < D_; ++j) {
        v_[j].draw(gen);
        for (int i = 0; i < n_; ++i) vmat_(i, j) = v_[j].f(i);
    }
}

// y^{(jk)}_i = e_ij - sum_{l<j, l!=k} phi_jl(Z_i) e_il, regressor e_ik,
// variance v_j(Z_i); phi_jk is refreshed by its ensemble after each pair.
void TreeCov::draw_phi_trees(rn& gen) {
    for (int j = 1; j < D_; ++j) {
        for (int k = 0; k < j; ++k) {
            double* y = Y_phi_.colptr(pair_index(j, k));
            for (int i = 0; i < n_; ++i) {
                double yi = e_(i, j);
                for (int l = 0; l < j; ++l)
                    if (l != k) yi -= phi_[j][l].f(i) * e_(i, l);
                y[i] = yi;
            }
            phi_[j][k].draw(e_.colptr(k), vmat_.colptr(j), gen);
        }
    }
}

void TreeCov::refresh_precision() {
    for (int i = 0; i < n_; ++i) {
        rootpi_[i] = cov::rootpi(eval(i));                  // R R^T = Q_i
        Q_[i] = rootpi_[i] * rootpi_[i].t();
    }
}

void TreeCov::store(int mkeep, vec const& /*mu*/, bool store_trees,
                    UniqueZCache& zc) {
    int s = mkeep - 1;
    for (int j = 0; j < D_; ++j) {
        var_varcount_.slice(s).col(j) = conv_to<vec>::from(v_[j].getnv());
        var_varprob_.slice(s).col(j)  = conv_to<vec>::from(v_[j].getpv());
    }
    if (store_trees) {
        for (int j = 0; j < D_; ++j)
            for (size_t h = 0; h < Hv_; ++h) var_treess_[j] << v_[j].gettree(h);
        for (int j = 1; j < D_; ++j)
            for (int k = 0; k < j; ++k)
                for (size_t h = 0; h < Hphi_; ++h)
                    phi_treess_[j][k] << phi_[j][k].gettree(h);
    }
    for (uword u = 0; u < zc.n_unique(); ++u) {
        mat A = solve(trimatu(rootpi_[zc.first_idx[u]]), eye<mat>(D_, D_));
        zc.store_sigma(mkeep, u, A.t() * A);                // Sigma = R^{-T} R^{-1}
    }
}

void TreeCov::pack(List& out, bool store_trees) {
    out["var_varcount"] = var_varcount_;
    out["var_varprob"]  = var_varprob_;
    if (store_trees) {
        out["var_models"] = pack_bart_list_(D_, v_, var_treess_);
        out["phi_models"] = full_cov_
            ? (SEXP) pack_bart_jagged_(D_, phi_, phi_treess_) : R_NilValue;
    }
}
