/*
 *  vecbart: sum of vector-leaf trees.  See vecbart.h.
 */
#include "vecbart.h"

//--------------------------------------------------
//constructor
vecbart::vecbart():m(200),D(0),t(m),pi(),p(0),n(0),x(0),xi(),di(),logdetC(0.0),
   dart(false),dartOn(false),aug(false),const_theta(false),
   a(0),b(0),rho(0),theta(1.),omega(1.) {}
vecbart::vecbart(size_t im, size_t iD):m(im),D(iD),t(m,vectree(iD)),pi(),
   p(0),n(0),x(0),xi(),di(),Cinv(arma::eye<arma::mat>(iD,iD)),logdetC(0.0),
   dart(false),dartOn(false),aug(false),const_theta(false),
   a(0),b(0),rho(0),theta(1.),omega(1.) {}
//--------------------------------------------------
void vecbart::setdata(size_t p, size_t n, double *x, size_t nc)
{
   this->p=p; this->n=n; this->x=x;
   if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);
   di.n=n; di.p=p; di.x=&x[0]; di.y=0;
   predict(p,n,x,fit_);
   r.zeros(n,D);
   ftemp.zeros(n,D);
   nv.assign(p,0);
   pv.assign(p,1.0/(double)p);
   lpv.assign(p,std::log(1.0/(double)p));
}
//--------------------------------------------------
void vecbart::setleafcov(arma::mat const& C)
{
   Cinv    = arma::inv_sympd(C);
   logdetC = arma::log_det_sympd(C);
}
//--------------------------------------------------
void vecbart::predict(size_t p, size_t n, double *x, arma::mat& fp)
//uses: m,t,xi
{
   arma::mat fptemp;
   fp.zeros(n,D);
   for(size_t j=0;j<m;j++) {
      t[j].fit(xi,p,n,x,fptemp);
      fp += fptemp;
   }
}
//--------------------------------------------------
// Backfitting sweep (as bart::draw): remove tree j from the fit, form the
// partial residual, update structure and leaves, add the tree back.
void vecbart::draw(arma::mat const& y, std::vector<arma::mat> const& Q, rn& gen)
{
   for(size_t j=0;j<m;j++) {
      t[j].fit(xi,p,n,x,ftemp);
      fit_ -= ftemp;
      r = y - fit_;
      vecbd(t[j],xi,di,pi,r,Q,Cinv,logdetC,nv,pv,aug,gen);
      vecdrmu(t[j],xi,di,r,Q,Cinv,gen);
      t[j].fit(xi,p,n,x,ftemp);
      fit_ += ftemp;
   }
   if(dartOn) {
     draw_s(nv,lpv,theta,gen);
     draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
     for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
   }
}
//--------------------------------------------------
size_t vecbart::nleaves() const
{
   size_t L=0;
   for(size_t j=0;j<m;j++) L += t[j].nleaves();
   return L;
}
//--------------------------------------------------
arma::vec vecbart::leaf_sumsq() const
{
   arma::vec s(D,arma::fill::zeros);
   for(size_t j=0;j<m;j++)
      for(auto const& kv : t[j].leaves()) s += arma::square(kv.second);
   return s;
}

//--------------------------------------------------
void vecbart::setxinfo(xinfo& _xi)
{
   size_t p=_xi.size();
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     size_t nc=_xi[i].size();
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = _xi[i][j];
   }
}
