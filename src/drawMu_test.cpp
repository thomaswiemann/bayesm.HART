// TEMPORARY SHIM: drives drawMuHeterCov on synthetic D and Phi values so we
// can verify the empirical posterior matches the analytic
//   mu | . ~ N( V_post (Amu mubar0 + sum_i Sigma(Z_i)^{-1} theta_i),
//               V_post )
// with V_post^{-1} = Amu + sum_i Sigma(Z_i)^{-1}.
//
// Same construction trick as cov_helpers_test: encode each d_j(Z_i) and
// phi_jk(Z_i) as a 1-tree stump whose leaf is the prescribed value, so all
// observations see the same Sigma (independent of i).

#include "bayesm.HART.h"

namespace {

void configure_const_varbart(varbart& vb, double value,
                             arma::mat const& Zt,
                             std::vector<double>& Zt_buf,
                             std::vector<double>& y_buf) {
    size_t nz = Zt.n_rows;
    size_t n  = Zt.n_cols;
    std::copy(Zt.memptr(), Zt.memptr() + Zt.n_elem, Zt_buf.begin());
    std::fill(y_buf.begin(), y_buf.end(), 0.0);

    vb = varbart(/*m=*/1);
    vb.setprior(0.95, 2.0, 1.0);
    vb.setdata(nz, n, Zt_buf.data(), y_buf.data(), 100);
    vb.setvarprior(/*nu*/ 3.0, /*lambda*/ value);
}

void configure_const_heterbart(heterbart& hb, double value,
                               arma::mat const& Zt,
                               std::vector<double>& Zt_buf,
                               std::vector<double>& y_buf) {
    size_t nz = Zt.n_rows;
    size_t n  = Zt.n_cols;
    std::copy(Zt.memptr(), Zt.memptr() + Zt.n_elem, Zt_buf.begin());
    std::fill(y_buf.begin(), y_buf.end(), 0.0);

    hb = heterbart(/*m=*/1);
    hb.setprior(0.95, 2.0, 1.0);
    hb.setdata(nz, n, Zt_buf.data(), y_buf.data(), 100);
    hb.gettree(0).tonull();
    hb.gettree(0).settheta(value);
    hb.setm(1);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List drawMuHeterCov_test(arma::mat const& oldbetas,
                               arma::vec const& mubar0,
                               arma::mat const& Amu,
                               arma::vec const& D_vals,
                               arma::mat const& Phi_vals,
                               bool diagonal,
                               int  n_draws)
{
    const size_t D = D_vals.n_elem;
    const size_t n = oldbetas.n_rows;
    if (oldbetas.n_cols != D)         Rcpp::stop("oldbetas n_cols mismatch");
    if (Phi_vals.n_rows != D || Phi_vals.n_cols != D) Rcpp::stop("Phi dims mismatch");

    arma::mat Zt(/*nz=*/1, /*n=*/(arma::uword)n);
    for (size_t i = 0; i < n; i++) Zt(0, (arma::uword)i) = static_cast<double>(i);

    const size_t n_phi   = diagonal ? 0 : (D * (D - 1)) / 2;
    const size_t n_trees = D + n_phi;
    std::vector<std::vector<double>> Zt_buffers(n_trees, std::vector<double>(n));
    std::vector<std::vector<double>> y_buffers (n_trees, std::vector<double>(n));

    std::vector<varbart> d_models(D);
    for (size_t j = 0; j < D; j++) {
        configure_const_varbart(d_models[j], D_vals[j], Zt,
                                Zt_buffers[j], y_buffers[j]);
    }

    std::vector<std::vector<heterbart>> phi_models;
    if (!diagonal) {
        phi_models.resize(D);
        size_t buf_idx = D;
        for (size_t j = 1; j < D; j++) {
            phi_models[j].resize(j);
            for (size_t k = 0; k < j; k++) {
                configure_const_heterbart(
                    phi_models[j][k], Phi_vals(j, k), Zt,
                    Zt_buffers[buf_idx], y_buffers[buf_idx]);
                ++buf_idx;
            }
        }
    }

    arn gen;
    arma::mat draws(n_draws, D);
    for (int s = 0; s < n_draws; s++) {
        arma::vec mu = drawMuHeterCov(oldbetas, mubar0, Amu, d_models, phi_models, gen);
        draws.row(s) = mu.t();
    }

    return Rcpp::List::create(
        Rcpp::Named("draws") = draws
    );
}
