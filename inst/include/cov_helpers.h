/*
 *  cov_helpers: modified-Cholesky primitives for
 *
 *      Sigma(Z_i)^{-1} = L(Z_i)^T D(Z_i)^{-1} L(Z_i),
 *
 *  L unit lower triangular with L_jk = -phi_jk(Z_i) for k < j, and
 *  D = diag(d_1(Z_i), ..., d_D(Z_i)).  An empty `phi` vector is the diagonal
 *  case.
 */
#ifndef GUARD_cov_helpers_h
#define GUARD_cov_helpers_h

#include <RcppArmadillo.h>
#include "varbart.h"
#include "vcbart.h"

namespace cov {

// View of (L, D) at unit i; non-owning references to the ensembles.
struct cov_eval {
    std::vector<varbart> const&             d;    // D variance ensembles
    std::vector<std::vector<vcbart>> const& phi;  // phi[j][k], k < j
    size_t                                  i;

    size_t dim() const { return d.size(); }
    double d_j(size_t j) const { return d[j].f(i); }
    double phi_jk(size_t j, size_t k) const {
        return phi.empty() ? 0.0 : phi[j][k].f(i);
    }
    bool diagonal() const { return phi.empty(); }
};

// v <- L v:  v_j <- v_j - sum_{k<j} phi_jk(Z_i) v_k.
void apply_L(arma::vec& v, cov_eval const& cov);

// R = L^T D^{-1/2}, upper triangular with R R^T = Sigma(Z_i)^{-1} (bayesm's
// rooti convention).
arma::mat rootpi(cov_eval const& cov);

}  // namespace cov

#endif
