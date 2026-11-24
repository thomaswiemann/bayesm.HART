/*
 *  Vector-leaf sufficient statistics, integrated likelihood, and leaf draws.
 *  See vecbartfuns.h for the leaf model.
 */
#include "vecbartfuns.h"

#include <map>

using namespace arma;

//--------------------------------------------------
double veclh(vecsuff const& s, mat const& Cinv, double logdetC)
{
   mat Pt = Cinv + s.P;
   mat U  = chol(Pt);                              // U'U = Pt
   vec w  = solve(trimatl(U.t()), s.b);            // b'Pt^{-1}b = w'w
   return -0.5 * logdetC - sum(log(U.diag())) + 0.5 * dot(w, w);
}
//--------------------------------------------------
vec vecdrawnodemu(vecsuff const& s, mat const& Cinv, rn& gen)
{
   mat Pt = Cinv + s.P;
   mat U  = chol(Pt);
   vec m  = solve(trimatu(U), solve(trimatl(U.t()), s.b));
   vec z(s.b.n_elem);
   for (uword d = 0; d < z.n_elem; d++) z(d) = gen.normal();
   return m + solve(trimatu(U), z);                // cov U^{-1}U^{-T} = Pt^{-1}
}
//--------------------------------------------------
//statistics for left and right given bot and v,c
void vecgetsuff(vectree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi,
                dinfo& di, mat const& y, std::vector<mat> const& Q,
                vecsuff& sl, vecsuff& sr)
{
   double *xx;
   sl = vecsuff(y.n_cols); sr = vecsuff(y.n_cols);
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.topo().bn(xx,xi)) {
         if(xx[v] < xi[v][c]) sl.add(Q[i], y.row(i));
         else                 sr.add(Q[i], y.row(i));
      }
   }
}
//--------------------------------------------------
//statistics for left and right bots
void vecgetsuff(vectree& x, tree::tree_p l, tree::tree_p r, xinfo& xi,
                dinfo& di, mat const& y, std::vector<mat> const& Q,
                vecsuff& sl, vecsuff& sr)
{
   double *xx;
   sl = vecsuff(y.n_cols); sr = vecsuff(y.n_cols);
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.topo().bn(xx,xi);
      if(bn==l) sl.add(Q[i], y.row(i));
      if(bn==r) sr.add(Q[i], y.row(i));
   }
}
//--------------------------------------------------
//statistics for all bottom nodes in one pass
void vecallsuff(vectree& x, xinfo& xi, dinfo& di, mat const& y,
                std::vector<mat> const& Q, tree::npv& bnv, std::vector<vecsuff>& sv)
{
   bnv.clear();
   x.topo().getbots(bnv);
   sv.assign(bnv.size(), vecsuff(y.n_cols));
   std::map<tree::tree_cp,size_t> bnmap;
   for(size_t i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;
   double *xx;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      sv[bnmap[x.topo().bn(xx,xi)]].add(Q[i], y.row(i));
   }
}
//--------------------------------------------------
//draw all the bottom node leaves
void vecdrmu(vectree& t, xinfo& xi, dinfo& di, mat const& y,
             std::vector<mat> const& Q, mat const& Cinv, rn& gen)
{
   tree::npv bnv;
   std::vector<vecsuff> sv;
   vecallsuff(t,xi,di,y,Q,bnv,sv);
   for(size_t i=0;i!=bnv.size();i++)
      t.setleaf(bnv[i], vecdrawnodemu(sv[i], Cinv, gen));
}
