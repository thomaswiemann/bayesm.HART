/*
 *  cov_helpers: see cov_helpers.h for the conventions.
 */
#include "cov_helpers.h"

#include <cmath>

namespace cov {

//--------------------------------------------------
// v <- L v, swept from j = D-1 down so each v_j reads untouched v_k, k < j.
void apply_L(arma::vec& v, cov_eval const& cov)
{
    if (cov.diagonal()) return;
    const size_t D = cov.dim();
    for (size_t j = D; j-- > 0; ) {
        double acc = v[j];
        for (size_t k = 0; k < j; ++k) acc -= cov.phi_jk(j, k) * v[k];
        v[j] = acc;
    }
}

//--------------------------------------------------
// R_kk = d_k^{-1/2},  R_jk = -phi_kj d_k^{-1/2} (j < k),  R_jk = 0 (j > k).
arma::mat rootpi(cov_eval const& cov)
{
    const size_t D = cov.dim();
    arma::mat R(D, D, arma::fill::zeros);
    for (size_t k = 0; k < D; ++k) {
        double inv_sqrt_dk = 1.0 / std::sqrt(cov.d_j(k));
        R(k, k) = inv_sqrt_dk;
        for (size_t j = 0; j < k; ++j) R(j, k) = -cov.phi_jk(k, j) * inv_sqrt_dk;
    }
    return R;
}

}  // namespace cov
