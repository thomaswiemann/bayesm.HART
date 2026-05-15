// TEMPORARY SMOKE TEST SHIM (to be removed after Phase 0 verification).
// Exposes varbart::draw on a fixed (Z, eta_sq) dataset so we can confirm
// in R that the product-of-trees ensemble recovers a known d(z) function.

#include "bayesm.HART.h"

// [[Rcpp::export]]
Rcpp::List varbart_smoke(arma::mat const& Z,
                         arma::vec const& eta_sq,
                         size_t m,
                         double nu,        // per-tree calibrated dof
                         double lambda,    // per-tree calibrated scale
                         double power,
                         double base,
                         int R_iter,
                         int burn,
                         int keep_every = 1)
{
    // Z arrives n x p; varbart expects column-major p x n (var-major-by-obs).
    arma::mat Zt = Z.t();
    size_t n  = Z.n_rows;
    size_t nz = Z.n_cols;

    // Mutable copies that own their storage so the C++ pointers stay valid.
    std::vector<double> Zt_buf(Zt.memptr(), Zt.memptr() + Zt.n_elem);
    std::vector<double> y_buf (eta_sq.memptr(), eta_sq.memptr() + eta_sq.n_elem);

    varbart vb(m);
    vb.setprior(base, power, /*tau (unused)*/ 1.0);
    vb.setdata(nz, n, Zt_buf.data(), y_buf.data(), /*numcut*/ 100);
    vb.setvarprior(nu, lambda);

    arn gen;

    int n_keep = (R_iter - burn + keep_every - 1) / keep_every;
    arma::mat dhat(n_keep, n, arma::fill::zeros);
    int row = 0;

    for (int it = 0; it < R_iter; it++) {
        vb.draw(gen);
        if (it >= burn && ((it - burn) % keep_every == 0) && row < n_keep) {
            for (size_t i = 0; i < n; i++) dhat(row, i) = vb.f(i);
            row++;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("d_draws") = dhat
    );
}
