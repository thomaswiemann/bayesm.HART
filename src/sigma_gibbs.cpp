// Conjugate (mu, L) blocks of the normalized HART model.  See sigma_gibbs.h.

#include "sigma_gibbs.h"

using namespace arma;

namespace {

// Statistics of the covariance block: V = Psi + sum x x^T, M = sum f* x^T.
void cholesky_stats(mat const& theta, mat const& fstar, vec const& mu,
                    mat const& Psi, mat& V, mat& M) {
    const uword D = Psi.n_rows;
    V = Psi;
    M.zeros(D, D);
    for (uword i = 0; i < theta.n_rows; ++i) {
        vec x = theta.row(i).t() - mu;
        V += x * x.t();
        M += fstar.row(i).t() * x.t();
    }
}

// Statistics of the mean block: P = Amu + n L^T L,
// b = Amu mubar + L^T sum_i (L theta_i - f*_i).
void mu_stats(mat const& L, mat const& theta, mat const& fstar,
              vec const& mubar, mat const& Amu, mat& P, vec& b) {
    const uword n = theta.n_rows, D = L.n_rows;
    P = Amu + (double)n * (L.t() * L);
    vec s(D, fill::zeros);
    for (uword i = 0; i < n; ++i) s += L * theta.row(i).t() - fstar.row(i).t();
    b = Amu * mubar + L.t() * s;
}

// Gaussian kernel of L_jk (k < j) given the other entries of row j.
void offdiag_kernel(mat const& V, mat const& M, mat const& L,
                    uword j, uword k, double& mean, double& prec) {
    double s = M(j, k);
    for (uword c = 0; c <= j; ++c) if (c != k) s -= V(k, c) * L(j, c);
    prec = V(k, k);
    mean = s / prec;
}

// Tilted kernel of L_jj: L^a exp(-b L^2 / 2 - c L) on L > 0 (j 0-based).
void diag_kernel(mat const& V, mat const& M, mat const& L, uword j,
                 double nu, uword n, uword D, double& a, double& b, double& c) {
    a = nu + (double)n - (double)D + (double)j;
    b = V(j, j);
    c = -M(j, j);
    for (uword cc = 0; cc < j; ++cc) c += V(j, cc) * L(j, cc);
}

// Slice sampler on u = log L_jj for
//   pi(u) propto exp{(a + 1) u - b e^{2u} / 2 - c e^u}.
// The density is unimodal for a > -1 (b t^2 + c t - (a + 1) = 0 has exactly
// one positive root t = e^u), so stepping out with a randomized budget plus
// shrinkage is valid.  Shrinkage terminates with probability one because u0
// lies inside the slice; the iteration cap turns a numerical breakdown into
// an explicit error rather than a silent rejection.
double slice_log_diag(double u0, double a, double b, double c,
                      double w, int max_steps, rn& gen) {
    auto logp = [&](double u) -> double {
        double eu = std::exp(u), e2u = eu * eu;
        if (!std::isfinite(e2u)) return -INFINITY;   // beyond any support
        return (a + 1.0) * u - 0.5 * b * e2u - c * eu;
    };
    double log_y = logp(u0) - gen.exp();
    double lo = u0 - w * gen.uniform(), hi = lo + w;
    int J = (int)std::floor(max_steps * gen.uniform());
    int K = (max_steps - 1) - J;
    while (J > 0 && logp(lo) > log_y) { lo -= w; --J; }
    while (K > 0 && logp(hi) > log_y) { hi += w; --K; }
    for (int it = 0; it < 100000; ++it) {
        double u = lo + (hi - lo) * gen.uniform();
        if (logp(u) > log_y) return u;
        if (u < u0) lo = u; else hi = u;
    }
    Rcpp::stop("slice_log_diag: shrinkage did not terminate (a=%g, b=%g, c=%g).",
               a, b, c);
}

}  // namespace

void draw_mu_normalized(vec& mu, mat const& L, mat const& theta,
                        mat const& fstar, vec const& mubar, mat const& Amu,
                        rn& gen) {
    mat P; vec b;
    mu_stats(L, theta, fstar, mubar, Amu, P, b);
    mat U = chol(P);                                        // U^T U = P
    vec m = solve(trimatu(U), solve(trimatl(U.t()), b));
    vec z(L.n_rows);
    for (uword d = 0; d < z.n_elem; ++d) z(d) = gen.normal();
    mu = m + solve(trimatu(U), z);                          // cov U^{-1} U^{-T}
}

void draw_cholesky_sweep(mat& L, mat const& theta, mat const& fstar,
                         vec const& mu, mat const& Psi, double nu, rn& gen) {
    const uword D = L.n_rows, n = theta.n_rows;
    mat V, M;
    cholesky_stats(theta, fstar, mu, Psi, V, M);
    for (uword j = 0; j < D; ++j) {
        for (uword k = 0; k < j; ++k) {
            double mean, prec;
            offdiag_kernel(V, M, L, j, k, mean, prec);
            L(j, k) = mean + gen.normal() / std::sqrt(prec);
        }
        double a, b, c;
        diag_kernel(V, M, L, j, nu, n, D, a, b, c);
        L(j, j) = std::exp(slice_log_diag(std::log(L(j, j)), a, b, c,
                                          1.0, 50, gen));
    }
}

// ---------------------------------------------------------------------------
// R-facing test harness (tests/testthat/test-sigma-gibbs-cholesky.R).
// ---------------------------------------------------------------------------

// Kernel parameters at the supplied state: for every (j, k < j) the Gaussian
// mean and precision of L_jk holding the rest of row j fixed; for every j the
// (a, b, c) of L_jj; and (P, b) of the mu block.
// [[Rcpp::export]]
Rcpp::List normalized_kernel_R(arma::mat const& L, arma::mat const& theta,
                               arma::mat const& fstar, arma::vec const& mu,
                               arma::mat const& Psi, double nu,
                               arma::vec const& mubar, arma::mat const& Amu) {
    const uword D = L.n_rows, n = theta.n_rows;
    mat V, M, P; vec b;
    cholesky_stats(theta, fstar, mu, Psi, V, M);
    mu_stats(L, theta, fstar, mubar, Amu, P, b);
    mat off_mean(D, D, fill::zeros), off_prec(D, D, fill::zeros);
    vec da(D), db(D), dc(D);
    for (uword j = 0; j < D; ++j) {
        for (uword k = 0; k < j; ++k)
            offdiag_kernel(V, M, L, j, k, off_mean(j, k), off_prec(j, k));
        diag_kernel(V, M, L, j, nu, n, D, da(j), db(j), dc(j));
    }
    return Rcpp::List::create(
        Rcpp::Named("offdiag_mean") = off_mean,
        Rcpp::Named("offdiag_prec") = off_prec,
        Rcpp::Named("diag_a") = da, Rcpp::Named("diag_b") = db,
        Rcpp::Named("diag_c") = dc,
        Rcpp::Named("mu_prec") = P, Rcpp::Named("mu_b") = b);
}

// `nsweep` iterations of (mu, L); either block can be pinned.
// [[Rcpp::export]]
Rcpp::List normalized_gibbs_R(arma::mat L, arma::vec mu,
                              arma::mat const& theta, arma::mat const& fstar,
                              arma::mat const& Psi, double nu,
                              arma::vec const& mubar, arma::mat const& Amu,
                              int nsweep = 1, bool update_mu = true,
                              bool update_L = true) {
    arn gen;
    for (int s = 0; s < nsweep; ++s) {
        if (update_mu) draw_mu_normalized(mu, L, theta, fstar, mubar, Amu, gen);
        if (update_L)  draw_cholesky_sweep(L, theta, fstar, mu, Psi, nu, gen);
    }
    return Rcpp::List::create(Rcpp::Named("L") = L, Rcpp::Named("mu") = mu);
}
