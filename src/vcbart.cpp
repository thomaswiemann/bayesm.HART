/*
 *  vcbart: varying-coefficient sum of trees.  See vcbart.h.
 */
#include "vcbart.h"

//--------------------------------------------------
void vcbart::pr()
{
   cout << "+++++vcbart object:\n";
   bart::pr();
}
//--------------------------------------------------
// Backfitting sweep (as bart::draw).  The partial residual removes the
// other trees' contribution scaled by the regressor:
//   r_i = y_i - a_i * sum_{b != j} g_b(x_i).
void vcbart::draw(double* reg, double* var, rn& gen)
{
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-reg[k]*allfit[k];
      }
      vcbd(t[j],xi,di,pi,reg,var,nv,pv,aug,gen);
      vcdrmu(t[j],xi,di,pi,reg,var,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
   }
   if(dartOn) {
     draw_s(nv,lpv,theta,gen);
     draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
     for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
   }
}
