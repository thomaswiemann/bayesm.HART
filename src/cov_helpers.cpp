/*
 *  cov_helpers — see inst/include/cov_helpers.h for conventions.
 */
#include "cov_helpers.h"

#include <cmath>

namespace cov {

//--------------------------------------------------
// v <- L v.  See header for sweep direction.
void apply_L(arma::vec& v, cov_eval const& cov)
{
    if (cov.diagonal()) return;          // L = I in diagonal case
    const size_t D = cov.dim();
    // Sweep from j = D-1 downward so each write to v[j] does not contaminate
    // subsequent reads of v[k] for k < j (which are still original entries).
    for (size_t j = D; j-- > 0; ) {
        double acc = v[j];
        for (size_t k = 0; k < j; ++k) {
            acc -= cov.phi_jk(j, k) * v[k];
        }
        v[j] = acc;
    }
}

//--------------------------------------------------
// Solve L y = v in place by forward substitution.  L unit lower-triangular,
// L_{jk} = -phi_{jk}, so y_j = v_j + sum_{k<j} phi_{jk} y_k.
void apply_L_inv(arma::vec& v, cov_eval const& cov)
{
    if (cov.diagonal()) return;          // L^{-1} = I in diagonal case
    const size_t D = cov.dim();
    for (size_t j = 0; j < D; ++j) {
        double acc = v[j];
        for (size_t k = 0; k < j; ++k) {
            acc += cov.phi_jk(j, k) * v[k];
        }
        v[j] = acc;
    }
}

//--------------------------------------------------
// v_j <- v_j / sqrt(d_j(Z_i)).
void scale_D_inv_sqrt(arma::vec& v, cov_eval const& cov)
{
    const size_t D = cov.dim();
    for (size_t j = 0; j < D; ++j) {
        v[j] /= std::sqrt(cov.d_j(j));
    }
}

//--------------------------------------------------
// v_j <- v_j * sqrt(d_j(Z_i)).
void scale_D_sqrt(arma::vec& v, cov_eval const& cov)
{
    const size_t D = cov.dim();
    for (size_t j = 0; j < D; ++j) {
        v[j] *= std::sqrt(cov.d_j(j));
    }
}

//--------------------------------------------------
// D^{-1/2} L (theta_i - mu).
arma::vec mahalanobis(arma::vec const& theta_i, arma::vec const& mu,
                      cov_eval const& cov)
{
    arma::vec out = theta_i - mu;
    apply_L         (out, cov);
    scale_D_inv_sqrt(out, cov);
    return out;
}

//--------------------------------------------------
// Sigma(Z_i)^{1/2} pred = L^{-1} D^{1/2} pred.
arma::vec sigma_sqrt_apply(arma::vec const& pred, cov_eval const& cov)
{
    arma::vec out = pred;
    scale_D_sqrt(out, cov);
    apply_L_inv (out, cov);
    return out;
}

//--------------------------------------------------
// R = L^T D^{-1/2}, upper-triangular, with R R^T = Sigma(Z_i)^{-1}.
//
//   R_{kk} =  d_k(Z_i)^{-1/2}
//   R_{jk} = -phi_{kj}(Z_i) * d_k(Z_i)^{-1/2}    for j < k
//   R_{jk} =  0                                   for j > k
arma::mat rootpi(cov_eval const& cov)
{
    const size_t D = cov.dim();
    arma::mat R(D, D, arma::fill::zeros);
    for (size_t k = 0; k < D; ++k) {
        double inv_sqrt_dk = 1.0 / std::sqrt(cov.d_j(k));
        R(k, k) = inv_sqrt_dk;
        for (size_t j = 0; j < k; ++j) {
            R(j, k) = -cov.phi_jk(k, j) * inv_sqrt_dk;
        }
    }
    return R;
}

}  // namespace cov
