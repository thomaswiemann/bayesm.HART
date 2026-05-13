/*
 *  cov_helpers — modified-Cholesky linear-algebra primitives for HART.
 *
 *  Encapsulates every Sigma(Z_i) operation needed by the MCMC loops so the
 *  loops stay thin orchestrators and each piece is unit-testable in isolation.
 *
 *  Convention (uniform across HART):
 *
 *      Sigma(Z_i)^{-1} = L(Z_i)^T D(Z_i)^{-1} L(Z_i)
 *
 *  with L(Z_i) unit lower-triangular, L_{jk} = -phi_{jk}(Z_i) for k < j,
 *  and D(Z_i) = diag(d_1(Z_i), ..., d_D(Z_i)) with d_j > 0.
 *
 *  The diagonal-only case (Phase 1) is signalled by an empty `phi` vector;
 *  every primitive degenerates to the diagonal arithmetic, no branching needed
 *  by the caller.
 */
#ifndef GUARD_cov_helpers_h
#define GUARD_cov_helpers_h

#include <RcppArmadillo.h>

#include "varbart.h"
#include "heterbart.h"

namespace cov {

//--------------------------------------------------
// Aggregating view of (L, D) at a single observation index i.
// Holds non-owning references to the variance and regression-coefficient
// ensembles; the caller is responsible for keeping them alive.
struct cov_eval {
    std::vector<varbart> const&                       d;     // size D
    std::vector<std::vector<heterbart>> const&        phi;   // jagged: phi[j][k] for k<j
    size_t                                            i;     // observation index

    // Number of dimensions (D) is given by d.size().
    size_t dim() const { return d.size(); }

    // Innovation variance d_j(Z_i).
    double d_j(size_t j) const { return d[j].f(i); }

    // Regression coefficient phi_{jk}(Z_i) for k < j.
    // Returns 0 in the diagonal case (phi empty).
    double phi_jk(size_t j, size_t k) const {
        return phi.empty() ? 0.0 : phi[j][k].f(i);
    }

    bool diagonal() const { return phi.empty(); }
};

//--------------------------------------------------
// Primitives — all in-place where it is safe to be in-place.
//
// apply_L:        v <- L v        =>  v_j <- v_j - sum_{k<j} phi_jk(Z_i) v_k.
//                 NOTE: writes must not contaminate later reads, so the sweep
//                 runs from j = dim()-1 down to 0 (each v_j depends only on
//                 v_0,...,v_{j-1}, which the sweep has not yet touched).
//
// apply_L_inv:    solve L y = v in place.  L is unit lower-triangular, so
//                 forward substitution: y_j = v_j + sum_{k<j} phi_jk(Z_i) y_k.
//                 Sweep runs j = 0,...,dim()-1; reads use already-overwritten
//                 entries, which by then equal the solution components.
//
// scale_D_inv_sqrt / scale_D_sqrt:  v_j <- v_j / sqrt(d_j(Z_i))   /   v_j * sqrt(d_j(Z_i)).
void apply_L         (arma::vec& v, cov_eval const& cov);
void apply_L_inv     (arma::vec& v, cov_eval const& cov);
void scale_D_inv_sqrt(arma::vec& v, cov_eval const& cov);
void scale_D_sqrt    (arma::vec& v, cov_eval const& cov);

//--------------------------------------------------
// Composite operations that the MCMC loops call directly.
//
//   mahalanobis(theta_i, mu, cov)  =  D^{-1/2} L (theta_i - mu)
//                                     ~ N(delta(Z_i), I) when (theta_i - mu)
//                                     ~ N(Sigma^{1/2} delta, Sigma).
//
//   sigma_sqrt_apply(pred, cov)    =  Sigma^{1/2} pred  =  L^{-1} D^{1/2} pred
//                                     (used to recover delta_Z from per-dim
//                                     mean-tree predictions).
//
//   rootpi(cov)                    =  L^T D^{-1/2}, an upper-triangular matrix
//                                     R satisfying R R^T = Sigma^{-1}.  Matches
//                                     the existing bayesm convention used by
//                                     mnlMetropOnce_con / lpostbeta.
arma::vec mahalanobis     (arma::vec const& theta_i, arma::vec const& mu,
                           cov_eval const& cov);
arma::vec sigma_sqrt_apply(arma::vec const& pred,
                           cov_eval const& cov);
arma::mat rootpi          (cov_eval const& cov);

}  // namespace cov

#endif
