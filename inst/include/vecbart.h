/*
 *  vecbart: sum of vector-leaf trees, one partition per tree,
 *
 *      f(x) = sum_h g(x; T_h, M_h),   g(x) = m_{h,leaf(x)} in R^D,
 *      m_{hg} ~ N_D(0, C),            y_i | f ~ N_D(f(x_i), Q_i^{-1}).
 *
 *  Mirrors `bart` with matrix responses and vector leaves.  Trees are value
 *  types (vectree) and the working buffers are Armadillo matrices, so the
 *  ensemble is copyable; only the predictor stack x (p x n) is caller-owned.
 */
#ifndef GUARD_vecbart_h
#define GUARD_vecbart_h

#include "vectree.h"
#include "info.h"
#include "bartfuns.h"
#include "vecbartfuns.h"
#include "vecbd.h"

class vecbart {
public:
   vecbart();
   vecbart(size_t m, size_t D);
   //------------------------------
   //get,set
   size_t getm() const { return m; }
   size_t getD() const { return D; }
   size_t getp() const { return p; }
   void setdata(size_t p, size_t n, double *x, size_t nc=100);
   void setprior(double alpha, double beta) { pi.alpha=alpha; pi.mybeta=beta; }
   void setleafguard(size_t nmin) { pi.nmin=nmin; }
   // Leaf prior covariance C; stores C^{-1} and log|C| for the leaf kernels.
   void setleafcov(arma::mat const& C);
   // The appendix's leaf prior N_D(0, B_tau / H) with B_tau = diag(tau2).
   void setscale(arma::vec const& tau2, size_t H)
      { setleafcov(arma::diagmat(tau2) / (double)H); }
   void setdart(double _a, double _b, double _rho, bool _aug, bool _dart,
                double _theta=0., double _omega=1.) {
     this->a=_a; this->b=_b; this->rho=_rho; this->aug=_aug;
     this->dart=_dart; this->omega=_omega;
     if(_theta==0.){ this->const_theta=false; this->theta=1.; }
     else          { this->const_theta=true;  this->theta=_theta; }
   }
   void startdart() { this->dartOn=!(this->dartOn); }
   vectree& gettree(size_t i) { return t[i]; }
   xinfo& getxinfo() { return xi; }
   void setxinfo(xinfo& _xi);
   std::vector<size_t>& getnv() { return nv; }
   std::vector<double>& getpv() { return pv; }
   std::vector<double>& getlpv() {return lpv;}
   bool getdartOn() const {return dartOn;}
   void setdartstate(double value, bool on) {theta=value; dartOn=on;}
   double gettheta() const {return theta;}
   void rebuildfit() {predict(p,n,x,fit_);}
   //------------------------------
   //public methods
   void predict(size_t p, size_t n, double *x, arma::mat& fp);
   // One backfitting sweep given the n x D response y and per-unit precisions.
   void draw(arma::mat const& y, std::vector<arma::mat> const& Q, rn& gen);
   arma::mat const& allfit() const { return fit_; }
   arma::rowvec f(size_t i) const { return fit_.row(i); }
   size_t nleaves() const;          // vector leaves across all trees
   arma::vec leaf_sumsq() const;    // sum_{h,g} m_{hg,d}^2, by coordinate d
protected:
   size_t m;  //number of trees
   size_t D;  //leaf dimension
   std::vector<vectree> t; //the trees
   pinfo pi; //prior and mcmc info
   //data
   size_t p,n; //x has dim p, n observations
   double *x;  //x is column stack, pxn
   xinfo xi; //cutpoint info
   dinfo di;
   //working
   arma::mat fit_, r, ftemp;   // n x D: current f(x_i), partial residual, one tree
   arma::mat Cinv;
   double logdetC;
   bool dart,dartOn,aug,const_theta;
   double a,b,rho,theta,omega;
   std::vector<size_t> nv;
   std::vector<double> pv, lpv;
};

#endif
