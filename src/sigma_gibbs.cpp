// Cholesky-Gibbs block update for (mu, Sigma) in the HART model.
//
// Math: discussions/2026-05-14-cholesky-gibbs-sigma.md
// Plan: discussions/2026-05-14-cholesky-gibbs-implementation-plan.md

#include "sigma_gibbs.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

namespace {

// Slice sample on u = log(L_jj).  Density (up to constant):
//   pi(u) propto exp((a+1) u - 0.5 b exp(2u) - c exp(u))
// Strictly log-concave for a > 0, so stepping-out + shrinkage is robust.
double slice_log_diag(double u0, double a, double b, double c,
                      double w, int max_steps, rn& gen) {
    auto logp = [&](double u) {
        double eu  = std::exp(u);
        double e2u = eu * eu;
        return (a + 1.0) * u - 0.5 * b * e2u - c * eu;
    };

    double f0 = logp(u0);
    double log_y = f0 - gen.exp();

    double L_pt = u0 - w * gen.uniform();
    double R_pt = L_pt + w;
    int J = (int)std::floor(max_steps * gen.uniform());
    int K = (max_steps - 1) - J;

    while (J > 0 && logp(L_pt) > log_y) {
        L_pt -= w;
        --J;
    }
    while (K > 0 && logp(R_pt) > log_y) {
        R_pt += w;
        --K;
    }

    for (int it = 0; it < 200; ++it) {
        double u_new = L_pt + (R_pt - L_pt) * gen.uniform();
        if (logp(u_new) > log_y) return u_new;
        if (u_new < u0) L_pt = u_new; else R_pt = u_new;
    }
    return u0;
}

}

void sigma_mu_block_gibbs(
    mat& L,
    vec& mu,
    mat const& theta,
    mat const& delta,
    mat const& Psi,
    double nu,
    vec const& mubar0,
    double Amu_scalar,
    rn& gen) {

    const uword D = L.n_rows;
    const uword n = theta.n_rows;

    // Sufficient statistics.  Handle n = 0 explicitly: with no data the
    // posterior reduces to the prior (IW(nu, Psi) on Sigma, N(mubar0,
    // Sigma/Amu) on mu); S, M, S_delta, theta_bar all collapse to zero.
    mat S(D, D, fill::zeros);
    mat M(D, D, fill::zeros);
    vec S_delta(D, fill::zeros);
    vec theta_bar(D, fill::zeros);
    if (n > 0) {
        mat X = theta;
        X.each_row() -= mu.t();
        S = X.t() * X;
        M = delta.t() * X;
        S_delta  = sum(delta, 0).t();
        theta_bar = mean(theta, 0).t();
    }
    mat V = Psi + S;

    // Cholesky sweep
    for (uword j = 0; j < D; ++j) {
        // Off-diagonal entries
        for (uword k = 0; k < j; ++k) {
            double b = V(k, k);
            double c_lin = -M(j, k);
            for (uword c = 0; c <= j; ++c) {
                if (c == k) continue;
                c_lin += V(k, c) * L(j, c);
            }
            double sigma2 = 1.0 / b;
            double mu_jk  = -c_lin * sigma2;
            L(j, k) = mu_jk + std::sqrt(sigma2) * gen.normal();
        }

        // Diagonal: slice on log scale
        double a = nu + (double)n - (double)D + (double)(j + 1) - 1.0;
        double b = V(j, j);
        double c_lin = -M(j, j);
        for (uword c = 0; c < j; ++c) {
            c_lin += V(j, c) * L(j, c);
        }
        double u0 = std::log(L(j, j));
        double u_new = slice_log_diag(u0, a, b, c_lin, 1.0, 50, gen);
        L(j, j) = std::exp(u_new);
    }

    // mu draw using the new L
    vec shift   = solve(trimatl(L), S_delta);
    double w    = 1.0 / ((double)n + Amu_scalar);
    vec mu_hat  = w * ((double)n * theta_bar + Amu_scalar * mubar0 - shift);

    vec z(D);
    for (uword i = 0; i < D; ++i) z(i) = gen.normal();
    vec mu_noise = solve(trimatl(L), z);
    mu = mu_hat + std::sqrt(w) * mu_noise;
}

// R-facing test wrapper (used by tests/testthat/test-sigma-gibbs-cholesky.R)
// [[Rcpp::export]]
Rcpp::List sigma_mu_block_gibbs_R(
    arma::mat L,
    arma::vec mu,
    arma::mat const& theta,
    arma::mat const& delta,
    arma::mat const& Psi,
    double nu,
    arma::vec const& mubar0,
    double Amu_scalar) {
    arn gen;
    sigma_mu_block_gibbs(L, mu, theta, delta, Psi, nu, mubar0, Amu_scalar, gen);
    return Rcpp::List::create(
        Rcpp::Named("L") = L,
        Rcpp::Named("mu") = mu);
}
