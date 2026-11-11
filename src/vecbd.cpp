/*
 *  vecbd: birth/death step for vector-leaf trees.  Mirrors bd; see vecbd.h.
 */
#include "vecbd.h"

bool vecbd(vectree& x, xinfo& xi, dinfo& di, pinfo& pi,
           arma::mat const& y, std::vector<arma::mat> const& Q,
           arma::mat const& Cinv, double logdetC,
           std::vector<size_t>& nv, std::vector<double>& pv,
           bool aug, rn& gen)
{
   tree& t = x.topo();
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(t,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(t,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      vecsuff sl, sr;
      vecgetsuff(x,nx,v,c,xi,di,y,Q,sl,sr);

      //--------------------------------------------------
      //compute alpha (admissibility floor as in bd)
      double alpha=0.0, lalpha=0.0;
      if((sl.n>=pi.nmin) && (sr.n>=pi.nmin)) {
         double lhl = veclh(sl,Cinv,logdetC);
         double lhr = veclh(sr,Cinv,logdetC);
         double lht = veclh(sl+sr,Cinv,logdetC);
         alpha=1.0;
         lalpha = std::min(0.0, log(pr) + (lhl+lhr-lht));
      }

      //--------------------------------------------------
      //try metrop
      if((alpha > 0) && (log(gen.uniform()) < lalpha)) {
         arma::vec ml = vecdrawnodemu(sl,Cinv,gen);
         arma::vec mr = vecdrawnodemu(sr,Cinv,gen);
         x.birth(nx,v,c,ml,mr);
         nv[v]++;
         return true;
      }
      return false;
   } else {
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree::tree_p nx; //nog node to death at
      dprop(t,xi,pi,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      vecsuff sl, sr;
      vecgetsuff(x,nx->getl(),nx->getr(),xi,di,y,Q,sl,sr);

      //--------------------------------------------------
      //compute alpha
      double lhl = veclh(sl,Cinv,logdetC);
      double lhr = veclh(sr,Cinv,logdetC);
      double lht = veclh(sl+sr,Cinv,logdetC);
      double lalpha = std::min(0.0, log(pr) + (lht-lhl-lhr));

      //--------------------------------------------------
      //try metrop
      if(log(gen.uniform()) < lalpha) {
         arma::vec m = vecdrawnodemu(sl+sr,Cinv,gen);
         nv[nx->getv()]--;
         x.death(nx,m);
         return true;
      }
      return false;
   }
}
