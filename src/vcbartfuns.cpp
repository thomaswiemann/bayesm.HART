/*
 *  Varying-coefficient sufficient statistics and leaf draws.
 *  See vcbartfuns.h for the leaf model; the score is heterlh.
 */
#include "vcbartfuns.h"

#include <map>

//--------------------------------------------------
//draw one leaf from its conditional (P = 0 gives a prior draw)
double vcdrawnodemu(double P, double b, double tau, rn& gen)
{
   double prec = 1.0/(tau*tau) + P;
   return b/prec + gen.normal()/sqrt(prec);
}
//--------------------------------------------------
//compute n, P, b for left and right given bot and v,c
void vcgetsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di,
               double* a, double* var,
               size_t& nl, double& Pl, double& bl, size_t& nr, double& Pr, double& br)
{
   double *xx;
   nl=0; Pl=0.0; bl=0.0; nr=0; Pr=0.0; br=0.0;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) {
         double w = a[i]/var[i];
         if(xx[v] < xi[v][c]) { nl++; Pl += w*a[i]; bl += w*di.y[i]; }
         else                 { nr++; Pr += w*a[i]; br += w*di.y[i]; }
      }
   }
}
//--------------------------------------------------
//compute n, P, b for left and right bots
void vcgetsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di,
               double* a, double* var,
               size_t& nl, double& Pl, double& bl, size_t& nr, double& Pr, double& br)
{
   double *xx;
   nl=0; Pl=0.0; bl=0.0; nr=0; Pr=0.0; br=0.0;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
      double w = a[i]/var[i];
      if(bn==l) { nl++; Pl += w*a[i]; bl += w*di.y[i]; }
      if(bn==r) { nr++; Pr += w*a[i]; br += w*di.y[i]; }
   }
}
//--------------------------------------------------
//statistics for all bottom nodes in one pass
void vcallsuff(tree& x, xinfo& xi, dinfo& di, double* a, double* var,
               tree::npv& bnv, std::vector<double>& Pv, std::vector<double>& bv)
{
   bnv.clear();
   x.getbots(bnv);
   Pv.assign(bnv.size(),0.0);
   bv.assign(bnv.size(),0.0);
   std::map<tree::tree_cp,size_t> bnmap;
   for(size_t i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;
   double *xx;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      size_t ni = bnmap[x.bn(xx,xi)];
      double w = a[i]/var[i];
      Pv[ni] += w*a[i];
      bv[ni] += w*di.y[i];
   }
}
//--------------------------------------------------
//draw all the bottom node leaves
void vcdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double* a, double* var, rn& gen)
{
   tree::npv bnv;
   std::vector<double> Pv, bv;
   vcallsuff(t,xi,di,a,var,bnv,Pv,bv);
   for(size_t i=0;i!=bnv.size();i++)
      bnv[i]->settheta(vcdrawnodemu(Pv[i],bv[i],pi.tau,gen));
}
