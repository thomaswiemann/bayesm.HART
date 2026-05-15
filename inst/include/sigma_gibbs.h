#ifndef SIGMA_GIBBS_H
#define SIGMA_GIBBS_H

#include <RcppArmadillo.h>
#include "rn.h"

// Componentwise Cholesky-Gibbs block update for (mu, Sigma) under the HART
// model.  Replaces the joint IW + MN draw in `rmultireg` for callers where
// theta_i | mu, Sigma, delta_i ~ N(mu + L^{-1} delta_i, Sigma) with
// L^T L = Sigma^{-1}, L lower-triangular.
//
// Convention: L is bayesm's `rooti^T` -- lower-triangular with
// L^T L = Sigma^{-1}.  Caller stores `rooti = L^T` (upper-triangular) in
// component lists.
//
// Inputs:
//   L_inout:    D x D lower-tri, L^T L = Sigma^{-1} on entry & exit.
//   mu_inout:   length-D mean vector, updated in place.
//   theta:      n x D rows of unit-level parameters for this component.
//   delta:      n x D rows of BART mean shifts (whitened-scale).
//   Psi:        D x D prior scale matrix (IW(nu, Psi)).
//   nu:         IW degrees of freedom.
//   mubar0:     length-D prior mean for mu.
//   Amu_scalar: scalar prior precision (mu | Sigma ~ N(mubar0, Sigma/Amu)).
//   gen:        random number generator.
//
// Cost: O(n D^2) for sufficient statistics, O(D^3) for the Cholesky sweep.
void sigma_mu_block_gibbs(
    arma::mat& L_inout,
    arma::vec& mu_inout,
    arma::mat const& theta,
    arma::mat const& delta,
    arma::mat const& Psi,
    double nu,
    arma::vec const& mubar0,
    double Amu_scalar,
    rn& gen);

#endif
