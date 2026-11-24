/*
 *  vcbart: varying-coefficient sum of trees,
 *
 *      y_i = a_i f(x_i) + e_i,   e_i ~ N(0, var_i),   f(x) = sum_h g_h(x),
 *
 *  used for the modified-Cholesky off-diagonals phi_jk(Z) where the regressor
 *  a_i is another coordinate's residual and var_i its innovation variance
 *  (appendix_sampler, app:hs_phi).  Subclass of `bart` like varbart: same
 *  data and ownership rules; draw() takes the regressor and variance arrays
 *  per sweep.  f(i) returns f(x_i), not a_i f(x_i).
 */
#ifndef GUARD_vcbart_h
#define GUARD_vcbart_h
#include "bart.h"
#include "vcbartfuns.h"
#include "vcbd.h"

class vcbart : public bart
{
public:
   vcbart():bart() {}
   vcbart(size_t m):bart(m) {}
   void pr();
   void draw(double* reg, double* var, rn& gen);
};

#endif
