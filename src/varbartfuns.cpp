/*
 *  Variance-tree sufficient statistics, integrated likelihood, and leaf draws.
 *  See varbartfuns.h for the conjugate model.
 */
#include "varbartfuns.h"

#include <map>

//--------------------------------------------------
// Log integrated marginal likelihood for chi^-2 leaf prior.
//
//   p(e_{leaf} | nu, lambda)
//      = int prod_i N(e_i | 0, s2_l) * Inv-chi^2(s2_l | nu, lambda) ds2_l
//      = pi^{-n/2} * (nu*lambda)^{nu/2} / Gamma(nu/2)
//                  * Gamma((n+nu)/2)    / (S + nu*lambda)^{(n+nu)/2}.
//
// All terms are kept.  The leaf-prior normalization
//      C1(nu, lambda) = (nu/2) log(nu*lambda) - lgamma(nu/2)
// scales with nu and contributes a +/- C1 offset to MH birth/death ratios
// (one extra leaf is integrated out under birth, one fewer under death), so
// it cannot be dropped without distorting the tree-structure posterior.  The
// -n/2 log(pi) term cancels exactly (n_l + n_r = n_t) but is included for
// auditability.
//--------------------------------------------------
double varlh(double S, double n, double nu, double lambda)
{
   double half    = 0.5 * (n + nu);
   double half_nu = 0.5 * nu;
   return  -0.5 * n * log(M_PI)
         + half_nu * log(nu * lambda)
         - lgamma(half_nu)
         + lgamma(half)
         - half * log(S + nu * lambda);
}

//--------------------------------------------------
// Birth proposal: count and sum-of-squares for the proposed
// left and right children of node nx splitting on (v, c).
//--------------------------------------------------
void vargetsuff(tree& x, tree::tree_p nx, size_t v, size_t c,
                xinfo& xi, dinfo& di,
                size_t& nl, double& Sl,
                size_t& nr, double& Sr)
{
   double *xx;
   nl = 0; nr = 0; Sl = 0.0; Sr = 0.0;

   for (size_t i = 0; i < di.n; i++) {
      xx = di.x + i * di.p;
      if (nx == x.bn(xx, xi)) {
         if (xx[v] < xi[v][c]) {
            nl += 1;
            Sl += di.y[i];
         } else {
            nr += 1;
            Sr += di.y[i];
         }
      }
   }
}

//--------------------------------------------------
// Death proposal: count and sum-of-squares for the existing
// left and right bots being collapsed.
//--------------------------------------------------
void vargetsuff(tree& x, tree::tree_p l, tree::tree_p r,
                xinfo& xi, dinfo& di,
                size_t& nl, double& Sl,
                size_t& nr, double& Sr)
{
   double *xx;
   nl = 0; nr = 0; Sl = 0.0; Sr = 0.0;

   for (size_t i = 0; i < di.n; i++) {
      xx = di.x + i * di.p;
      tree::tree_cp bn = x.bn(xx, xi);
      if (bn == l) {
         nl += 1;
         Sl += di.y[i];
      } else if (bn == r) {
         nr += 1;
         Sr += di.y[i];
      }
   }
}

//--------------------------------------------------
// Posterior draw of s^2_lk:
//   s^2_lk | data ~ Inv-chi^2(nu + n, (nu*lambda + S)/(nu + n))
//                = (nu*lambda + S) / chi^2_{nu + n}.
//--------------------------------------------------
double vardrawnodemu(double S, double n, double nu, double lambda, rn& gen)
{
   double dof = nu + n;
   double scale_sum = nu * lambda + S;
   return scale_sum / gen.chi_square(dof);
}

//--------------------------------------------------
// All-bottom-node sufficient statistics in one pass.
//--------------------------------------------------
void varallsuff(tree& x, xinfo& xi, dinfo& di,
                tree::npv& bnv,
                std::vector<size_t>& nv_leaf,
                std::vector<double>& Sv)
{
   tree::tree_cp tbn;
   size_t ni;
   double *xx;

   bnv.clear();
   x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   nv_leaf.assign(nb, 0);
   Sv.assign(nb, 0.0);

   std::map<tree::tree_cp, size_t> bnmap;
   for (bvsz i = 0; i < nb; i++) bnmap[bnv[i]] = i;

   for (size_t i = 0; i < di.n; i++) {
      xx = di.x + i * di.p;
      tbn = x.bn(xx, xi);
      ni = bnmap[tbn];
      nv_leaf[ni] += 1;
      Sv[ni] += di.y[i];
   }
}

//--------------------------------------------------
// Sweep all leaves and draw new s^2 values from their full conditionals.
//--------------------------------------------------
void vardrmu(tree& t, xinfo& xi, dinfo& di, pinfo& /*pi*/,
             double nu, double lambda, rn& gen)
{
   tree::npv bnv;
   std::vector<size_t> nv_leaf;
   std::vector<double> Sv;
   varallsuff(t, xi, di, bnv, nv_leaf, Sv);
   for (tree::npv::size_type i = 0; i < bnv.size(); i++) {
      bnv[i]->settheta(
         vardrawnodemu(Sv[i], static_cast<double>(nv_leaf[i]), nu, lambda, gen));
   }
}
