/*
 *  varbart: product-of-trees ensemble for positive scalar functions
 *  (e.g., innovation variances d_j(Z_i) in the modified Cholesky model).
 *
 *  Design (Pratola et al., 2020):
 *      f(x) = prod_{l=1}^{m'} h_l(x;  T'_l, M'_l),         h_l(x) > 0
 *      h_l(x) ~ Inv-chi^2(nu, lambda)                       (per-tree leaf prior;
 *                                                           caller supplies the
 *                                                           per-tree calibrated
 *                                                           hyperparameters)
 *
 *  Inherits all bookkeeping from `bart` (xinfo, dinfo, allfit, r, ftemp, dart).
 *  Differences from `bart` / `heterbart`:
 *   1. Backfitting is **multiplicative**:   allfit /= ftemp;  ...;  allfit *= ftemp.
 *   2. Tree update uses chi^-2 conjugacy:   varbd  + vardrmu.
 *   3. Leaf draws return s^2 > 0 directly (no scaling argument).
 *   4. Stumps are initialized to lambda^{1/m} so the initial product equals lambda.
 *
 *  Caller responsibilities:
 *   * Call setdata(...) first (sets up xinfo, dinfo, allfit, r, ftemp).
 *   * Then call setvarprior(nu, lambda) to fix prior hyperparameters and
 *     reinitialize stumps + allfit to a consistent state.
 *   * Each MCMC iteration: write the squared innovations eta_i^2 into y via
 *     updateY(...) (or hold the y pointer fixed and update its contents in place
 *     between draws), then call draw(gen).
 */
#ifndef GUARD_varbart_h
#define GUARD_varbart_h

#include "bart.h"
#include "varbartfuns.h"
#include "varbd.h"

class varbart : public bart {
public:
   varbart(): bart(),    nu(3.0), lambda(1.0) {}
   varbart(size_t im): bart(im), nu(3.0), lambda(1.0) {}

   double getnu()     const { return nu; }
   double getlambda() const { return lambda; }

   // Set chi^-2 prior, collapse all trees to stumps with leaf lambda^{1/m},
   // and recompute allfit (= product of stumps = lambda) if data is set.
   void setvarprior(double nu_, double lambda_);

   // Multiplicative backfitting sweep: one MH step per tree, then chi^-2 leaf draw.
   // No sigma argument — variance residuals e_i^2 = y[i] / allfit_{-l}(i)
   // are constructed internally from y (which holds eta_i^2) and the partial product.
   void draw(rn& gen);

   void pr();

protected:
   double nu;      // chi^-2 prior dof (per-tree, calibrated)
   double lambda;  // chi^-2 prior scale (per-tree, calibrated)

   // Multiplicative analogue of bart::predict; computes prod_l h_l(x_i).
   void predict_prod(size_t p_, size_t n_, double *x_, double *fp);
};

#endif
