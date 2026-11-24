#ifndef GUARD_vcbartfuns_h
#define GUARD_vcbartfuns_h

#include <cmath>
#include "tree.h"
#include "treefuns.h"
#include "info.h"

// Varying-coefficient leaf model (appendix_sampler, eq:hs_phi_leaf):
//
//   y_i | lambda ~ N(a_i lambda, var_i),    lambda ~ N(0, tau^2),
//
// with sufficient statistics P = sum a_i^2 / var_i and b = sum a_i y_i / var_i.
// A leaf with P = 0 has score one and draws from its prior.

// log integral of the leaf prior against the likelihood, up to -log(tau),
// which cancels in the birth/death ratio.
inline double vclh(double P, double b, double tau) {
   double t2 = tau * tau;
   double k  = P * t2 + 1.0;
   return -0.5 * std::log(k) + 0.5 * b * b * t2 / k;
}
double vcdrawnodemu(double P, double b, double tau, rn& gen);

// Birth proposal: statistics of the two children of nx under split (v, c).
void vcgetsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di,
               double* a, double* var,
               size_t& nl, double& Pl, double& bl, size_t& nr, double& Pr, double& br);
// Death proposal: statistics of the existing children l and r.
void vcgetsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di,
               double* a, double* var,
               size_t& nl, double& Pl, double& bl, size_t& nr, double& Pr, double& br);
// All bottom nodes in one pass.
void vcallsuff(tree& x, xinfo& xi, dinfo& di, double* a, double* var,
               tree::npv& bnv, std::vector<double>& Pv, std::vector<double>& bv);
// Redraw every leaf from its conditional.
void vcdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double* a, double* var, rn& gen);

#endif
