#ifndef SIGMA_GIBBS_H
#define SIGMA_GIBBS_H

#include <RcppArmadillo.h>
#include "rn.h"

// Conjugate blocks of the normalized HART model (appendix_models,
// app:models_computation):
//
//   theta_i = mu + L^{-1} { f*_i + eps_i },   eps_i ~ N(0, I),
//   L lower triangular with L^T L = Sigma^{-1},
//   mu ~ N(mubar, Amu^{-1}),   Sigma ~ IW(nu, Psi).
//
// Convention: L is bayesm's rooti^T; callers store rooti = L^T.
//
// draw_mu_normalized   mu | L, f*, theta        (eq:models_mu_update)
//   P = Amu + n L^T L,   b = Amu mubar + L^T sum_i (L theta_i - f*_i).
//
// draw_cholesky_sweep  L | mu, f*, theta        (eq:models_cholesky_offdiag,
//   eq:models_cholesky_diag) with x_i = theta_i - mu,
//   V = Psi + sum_i x_i x_i^T,  M = sum_i f*_i x_i^T.
//   Row j: L_jk ~ N((M_jk - sum_{c<=j,c!=k} V_kc L_jc) / V_kk, 1 / V_kk) for
//   k < j, then L_jj from L^a exp(-b L^2/2 - c L) with
//   a = nu + n - D + j - 1 (j 1-based), b = V_jj, c = sum_{c<j} V_jc L_jc - M_jj,
//   slice sampled on log L_jj with the Jacobian included.
//
// Cost: O(n D^2) for the statistics, O(D^3) for the sweep.

void draw_mu_normalized(arma::vec& mu, arma::mat const& L,
                        arma::mat const& theta, arma::mat const& fstar,
                        arma::vec const& mubar, arma::mat const& Amu,
                        rn& gen);

void draw_cholesky_sweep(arma::mat& L, arma::mat const& theta,
                         arma::mat const& fstar, arma::vec const& mu,
                         arma::mat const& Psi, double nu, rn& gen);

#endif
