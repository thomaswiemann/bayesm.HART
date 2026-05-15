// TEMPORARY SHIM (to be removed once Phase 1 integration tests subsume it).
// Drives cov_helpers primitives on synthetic D and phi values so we can
// compare against a brute-force armadillo reference in R.
//
// Follows the bart ownership invariants documented in inst/include/bart.h:
//   * std::vector<bart-subclass> is sized once via default construction.
//   * Each entry is configured via op= (assigning from rvalue), then setdata.
//   * The vector is never resized after configuration.
//   * Backing buffers (Z transpose, y) are owned by std::vector<double> whose
//     storage is allocated up front and never reallocated.

#include "bayesm.HART.h"

namespace {

// Configure a freshly default-constructed varbart to be a 1-tree ensemble
// whose single stump leaf equals `value` (so cov_eval::d_j returns `value`
// for every i).  Buffers `Zt_buf` and `y_buf` MUST be pre-sized; the tree
// stores raw pointers into them.
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
    // m = 1, so lambda^(1/m) = lambda; setvarprior writes the single stump
    // leaf to `value` and refreshes allfit.
    vb.setvarprior(/*nu*/ 3.0, /*lambda*/ value);
}

// Same idea for heterbart.  bart::setm(1) re-runs predict so allfit reflects
// the manually set leaf.
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
    hb.setm(1);   // retriggers predict, refreshing allfit to the new leaf
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cov_helpers_test(arma::vec const& D_vals,
                            arma::mat const& Phi_vals,   // strict lower triangle used
                            arma::vec const& v_in,
                            arma::vec const& mu_in,
                            bool diagonal)
{
    const size_t D = D_vals.n_elem;
    if (Phi_vals.n_rows != D || Phi_vals.n_cols != D)
        Rcpp::stop("Phi_vals must be D x D");

    // Synthetic Z with two distinct rows so makexinfo has a non-degenerate
    // range; cov_eval only ever reads i = 0 in this shim.
    arma::mat Zt(/*nz=*/1, /*n=*/2);
    Zt(0, 0) = 0.0;
    Zt(0, 1) = 1.0;
    const size_t n = 2;

    // Backing buffers (one Zt + one y per tree).  Pre-sized; never reallocated.
    const size_t n_phi = diagonal ? 0 : (D * (D - 1)) / 2;
    const size_t n_trees = D + n_phi;
    std::vector<std::vector<double>> Zt_buffers(n_trees, std::vector<double>(n));
    std::vector<std::vector<double>> y_buffers (n_trees, std::vector<double>(n));

    // d_models: configured in-place via assign-then-setdata.
    std::vector<varbart> d_models(D);
    for (size_t j = 0; j < D; ++j) {
        configure_const_varbart(d_models[j], D_vals[j], Zt,
                                Zt_buffers[j], y_buffers[j]);
    }

    // phi_models: jagged, sized once, never resized.
    std::vector<std::vector<heterbart>> phi_models;
    if (!diagonal) {
        phi_models.resize(D);
        size_t buf_idx = D;
        for (size_t j = 1; j < D; ++j) {
            phi_models[j].resize(j);   // default-construct j heterbarts
            for (size_t k = 0; k < j; ++k) {
                configure_const_heterbart(
                    phi_models[j][k], Phi_vals(j, k), Zt,
                    Zt_buffers[buf_idx], y_buffers[buf_idx]);
                ++buf_idx;
            }
        }
    }

    cov::cov_eval ev{d_models, phi_models, /*i=*/0};

    arma::vec v_L     = v_in;  cov::apply_L         (v_L,    ev);
    arma::vec v_Linv  = v_in;  cov::apply_L_inv     (v_Linv, ev);
    arma::vec v_Dinv  = v_in;  cov::scale_D_inv_sqrt(v_Dinv, ev);
    arma::vec v_Dpos  = v_in;  cov::scale_D_sqrt    (v_Dpos, ev);
    arma::vec mhl     = cov::mahalanobis     (v_in, mu_in, ev);
    arma::vec sgs     = cov::sigma_sqrt_apply(v_in,         ev);
    arma::mat R       = cov::rootpi          (              ev);

    return Rcpp::List::create(
        Rcpp::Named("apply_L")     = v_L,
        Rcpp::Named("apply_L_inv") = v_Linv,
        Rcpp::Named("scale_Dinv")  = v_Dinv,
        Rcpp::Named("scale_Dpos")  = v_Dpos,
        Rcpp::Named("mahalanobis") = mhl,
        Rcpp::Named("sigma_sqrt")  = sgs,
        Rcpp::Named("rootpi")      = R
    );
}
