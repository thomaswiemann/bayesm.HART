/*
 *  varbart: product-of-trees variance ensemble.  See varbart.h.
 */
#include "varbart.h"

#include <cmath>

//--------------------------------------------------
void varbart::pr()
{
   cout << "+++++varbart object:" << std::endl;
   cout << "nu, lambda: " << nu << ", " << lambda << std::endl;
   bart::pr();
}

//--------------------------------------------------
// Set prior hyperparameters and reset all trees to stumps with leaf lambda^{1/m}
// so that the initial product equals lambda.  If data is already wired in
// (xi, allfit, p, n, x), refresh allfit accordingly.
void varbart::setvarprior(double nu_, double lambda_)
{
   nu = nu_;
   lambda = lambda_;

   double leaf_init = std::pow(lambda, 1.0 / static_cast<double>(m));
   for (size_t j = 0; j < m; j++) {
      t[j].tonull();
      t[j].settheta(leaf_init);
   }
   if (allfit && (xi.size() == p) && x) {
      predict_prod(p, n, x, allfit);
   }
}

//--------------------------------------------------
// Product-of-trees prediction: fp[i] = prod_l h_l(x_i).
void varbart::predict_prod(size_t p_, size_t n_, double *x_, double *fp)
{
   double *fptemp = new double[n_];
   for (size_t k = 0; k < n_; k++) fp[k] = 1.0;
   for (size_t j = 0; j < m; j++) {
      fit(t[j], xi, p_, n_, x_, fptemp);
      for (size_t k = 0; k < n_; k++) fp[k] *= fptemp[k];
   }
   delete[] fptemp;
}

//--------------------------------------------------
// Multiplicative backfitting sweep.
//
// Invariant maintained across the loop:
//   allfit[i] = prod_l h_l(x_i)  (current product over all trees)
//
// For each tree l:
//   1. Refit h_l to ftemp.
//   2. Divide ftemp out of allfit -> partial product s2_{-l}(x_i).
//   3. r[i] = y[i] / allfit[i]  is the variance residual e_i^2
//      (y holds the pre-tree squared innovation eta_i^2).
//   4. MH birth/death + chi^-2 leaf draw on tree l.
//   5. Refit new h_l to ftemp and multiply back into allfit.
void varbart::draw(rn& gen)
{
   for (size_t j = 0; j < m; j++) {
      fit(t[j], xi, p, n, x, ftemp);
      for (size_t k = 0; k < n; k++) {
         allfit[k] = allfit[k] / ftemp[k];   // s2_{-l}(x_k)
         r[k]      = y[k]      / allfit[k];  // e_k^2 = eta_k^2 / s2_{-l}
      }
      varbd(t[j], xi, di, pi, nu, lambda, nv, pv, aug, gen);
      vardrmu(t[j], xi, di, pi, nu, lambda, gen);
      fit(t[j], xi, p, n, x, ftemp);
      for (size_t k = 0; k < n; k++) allfit[k] *= ftemp[k];
   }

   if (dartOn) {
      draw_s(nv, lpv, theta, gen);
      draw_theta0(const_theta, theta, lpv, a, b, rho, gen);
      for (size_t j = 0; j < p; j++) pv[j] = ::exp(lpv[j]);
   }
}
