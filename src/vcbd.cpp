/*
 *  vcbd: birth/death step for varying-coefficient trees.  Mirrors bd; see
 *  vcbd.h.
 */
#include "vcbd.h"

bool vcbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double* a, double* var,
          std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen)
{
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      size_t nl,nr;
      double Pl,bl,Pr,br;
      vcgetsuff(x,nx,v,c,xi,di,a,var,nl,Pl,bl,nr,Pr,br);

      //--------------------------------------------------
      //compute alpha (admissibility floor as in bd)
      double alpha=0.0, lalpha=0.0;
      if((nl>=pi.nmin) && (nr>=pi.nmin)) {
         double lhl = vclh(Pl,bl,pi.tau);
         double lhr = vclh(Pr,br,pi.tau);
         double lht = vclh(Pl+Pr,bl+br,pi.tau);
         alpha=1.0;
         lalpha = std::min(0.0, log(pr) + (lhl+lhr-lht));
      }

      //--------------------------------------------------
      //try metrop
      if((alpha > 0) && (log(gen.uniform()) < lalpha)) {
         double mul = vcdrawnodemu(Pl,bl,pi.tau,gen);
         double mur = vcdrawnodemu(Pr,br,pi.tau,gen);
         x.birthp(nx,v,c,mul,mur);
         nv[v]++;
         return true;
      }
      return false;
   } else {
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree::tree_p nx; //nog node to death at
      dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      size_t nl,nr;
      double Pl,bl,Pr,br;
      vcgetsuff(x,nx->getl(),nx->getr(),xi,di,a,var,nl,Pl,bl,nr,Pr,br);

      //--------------------------------------------------
      //compute alpha
      double lhl = vclh(Pl,bl,pi.tau);
      double lhr = vclh(Pr,br,pi.tau);
      double lht = vclh(Pl+Pr,bl+br,pi.tau);
      double lalpha = std::min(0.0, log(pr) + (lht-lhl-lhr));

      //--------------------------------------------------
      //try metrop
      if(log(gen.uniform()) < lalpha) {
         double mu = vcdrawnodemu(Pl+Pr,bl+br,pi.tau,gen);
         nv[nx->getv()]--;
         x.deathp(nx,mu);
         return true;
      }
      return false;
   }
}
