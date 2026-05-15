/*
 *  Birth/death Metropolis-Hastings step for variance trees.
 *  Structure mirrors heterbd; differences are confined to the leaf
 *  likelihood (varlh) and leaf-draw (vardrawnodemu) calls.
 */
#include "varbd.h"

bool varbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi,
           double nu, double lambda,
           std::vector<size_t>& nv, std::vector<double>& pv,
           bool aug, rn& gen)
{
   tree::npv goodbots;
   double PBx = getpb(x, xi, pi, goodbots);

   if (gen.uniform() < PBx) {
      //--------------------------------------------------
      // BIRTH proposal
      tree::tree_p nx;
      size_t v, c;
      double pr;
      bprop(x, xi, pi, goodbots, PBx, nx, v, c, pr, nv, pv, aug, gen);

      size_t nl, nr;
      double Sl, Sr;
      vargetsuff(x, nx, v, c, xi, di, nl, Sl, nr, Sr);

      double alpha = 0.0, lalpha = 0.0;
      if ((nl >= 5) && (nr >= 5)) {
         double n_l = static_cast<double>(nl);
         double n_r = static_cast<double>(nr);
         double lhl = varlh(Sl,        n_l,           nu, lambda);
         double lhr = varlh(Sr,        n_r,           nu, lambda);
         double lht = varlh(Sl + Sr,   n_l + n_r,     nu, lambda);

         alpha  = 1.0;
         lalpha = log(pr) + (lhl + lhr - lht);
         lalpha = std::min(0.0, lalpha);
      }

      if ((alpha > 0) && (log(gen.uniform()) < lalpha)) {
         double mul = vardrawnodemu(Sl, static_cast<double>(nl), nu, lambda, gen);
         double mur = vardrawnodemu(Sr, static_cast<double>(nr), nu, lambda, gen);
         x.birthp(nx, v, c, mul, mur);
         nv[v]++;
         return true;
      }
      return false;
   } else {
      //--------------------------------------------------
      // DEATH proposal
      double pr;
      tree::tree_p nx;
      dprop(x, xi, pi, goodbots, PBx, nx, pr, gen);

      size_t nl, nr;
      double Sl, Sr;
      vargetsuff(x, nx->getl(), nx->getr(), xi, di, nl, Sl, nr, Sr);

      double n_l = static_cast<double>(nl);
      double n_r = static_cast<double>(nr);
      double lhl = varlh(Sl,        n_l,         nu, lambda);
      double lhr = varlh(Sr,        n_r,         nu, lambda);
      double lht = varlh(Sl + Sr,   n_l + n_r,   nu, lambda);

      double lalpha = log(pr) + (lht - lhl - lhr);
      lalpha = std::min(0.0, lalpha);

      if (log(gen.uniform()) < lalpha) {
         double mu = vardrawnodemu(Sl + Sr, n_l + n_r, nu, lambda, gen);
         nv[nx->getv()]--;
         x.deathp(nx, mu);
         return true;
      }
      return false;
   }
}
