/*
 *  cpwvarbart — product-of-trees prediction for varbart ensembles.
 *
 *  Mirrors cpwbart (sum-of-trees) exactly EXCEPT:
 *    - Per-draw accumulator yhat(i, k) is initialized to 1.0 (multiplicative
 *      identity) instead of 0.0.
 *    - Within a draw, individual tree contributions are multiplied into
 *      yhat instead of summed.
 *
 *  Each leaf of a varbart tree carries s^2_{lk} (a positive variance scale),
 *  and the ensemble represents d_j(Z) = prod_l s^2_{l, leaf_l(Z)}.  This file
 *  is the prediction-time analog of varbart::draw's multiplicative
 *  backfitting (see src/varbart.cpp).
 */
#include "tree.h"
#include "treefuns.h"

typedef std::vector<tree> vtree;

#ifdef _OPENMP
#include <omp.h>
static void local_getpred_var(size_t nd, size_t p, size_t m, size_t np,
                              xinfo& xi, std::vector<vtree>& tmat,
                              double* px, Rcpp::NumericMatrix& yhat);
#endif

static void getpred_var(int beg, int end, size_t p, size_t m, size_t np,
                        xinfo& xi, std::vector<vtree>& tmat,
                        double* px, Rcpp::NumericMatrix& yhat);

RcppExport SEXP cpwvarbart(
   SEXP _itrees,    //treedraws list (cutpoints + serialized trees)
   SEXP _ix,        //x matrix to predict at (transposed: ncov x npred)
   SEXP _itc,       //thread count
   SEXP _verbose    //verbosity flag
)
{
   bool verbose = Rcpp::as<bool>(_verbose);
   if (verbose) Rprintf("*****In main of C++ for varbart prediction (product-of-trees)\n");

   int tc = Rcpp::as<int>(_itc);
   if (verbose) Rcpp::Rcout << "tc (threadcount): " << tc << std::endl;

   //--------------------------------------------------
   //process trees
   Rcpp::List trees(_itrees);
   Rcpp::CharacterVector itrees(Rcpp::wrap(trees["trees"]));
   std::string itv(itrees[0]);
   std::stringstream ttss(itv);

   size_t nd, m, p;
   ttss >> nd >> m >> p;
   if (verbose) {
      Rcpp::Rcout << "number of varbart draws:    " << nd << std::endl;
      Rcpp::Rcout << "number of trees per ensemble: " << m << std::endl;
      Rcpp::Rcout << "number of x columns:          " << p  << std::endl;
   }

   //--------------------------------------------------
   //process cutpoints
   Rcpp::List ixi(Rcpp::wrap(trees["cutpoints"]));
   size_t pp = ixi.size();
   if (p != pp && verbose)
      Rcpp::Rcout << "WARNING: p from trees and p from x don't agree\n";
   xinfo xi;
   xi.resize(p);
   for (size_t i = 0; i < p; i++) {
      Rcpp::NumericVector cutv(ixi[i]);
      xi[i].resize(cutv.size());
      std::copy(cutv.begin(), cutv.end(), xi[i].begin());
   }

   //--------------------------------------------------
   //process x (R-side transposes, so ncol(xpred) == npred)
   Rcpp::NumericMatrix xpred(_ix);
   size_t np = xpred.ncol();
   if (verbose)
      Rcpp::Rcout << "from x, npred, ncov: "
                  << xpred.nrow() << ", " << xpred.ncol() << std::endl;

   //--------------------------------------------------
   //read in trees
   std::vector<vtree> tmat(nd);
   for (size_t i = 0; i < nd; i++) tmat[i].resize(m);
   for (size_t i = 0; i < nd; i++) {
      for (size_t j = 0; j < m; j++) ttss >> tmat[i][j];
   }

   //--------------------------------------------------
   //get predictions: yhat starts at 1 (multiplicative identity)
   Rcpp::NumericMatrix yhat(nd, np);
   std::fill(yhat.begin(), yhat.end(), 1.0);
   double* px = &xpred(0, 0);

#ifndef _OPENMP
   if (verbose) Rcpp::Rcout << "***using serial code\n";
   getpred_var(0, nd - 1, p, m, np, xi, tmat, px, yhat);
#else
   if (tc == 1) {
      if (verbose) Rcpp::Rcout << "***using serial code\n";
      getpred_var(0, nd - 1, p, m, np, xi, tmat, px, yhat);
   } else {
      if (verbose) Rcpp::Rcout << "***using parallel code\n";
#pragma omp parallel num_threads(tc)
      local_getpred_var(nd, p, m, np, xi, tmat, px, yhat);
   }
#endif

   Rcpp::List ret;
   ret["yhat.test"] = yhat;
   return ret;
}

static void getpred_var(int beg, int end, size_t p, size_t m, size_t np,
                        xinfo& xi, std::vector<vtree>& tmat,
                        double* px, Rcpp::NumericMatrix& yhat)
{
   double* fptemp = new double[np];
   for (int i = beg; i <= end; i++) {
      for (size_t j = 0; j < m; j++) {
         fit(tmat[i][j], xi, p, np, px, fptemp);
         for (size_t k = 0; k < np; k++) yhat(i, k) *= fptemp[k];
      }
   }
   delete[] fptemp;
}

#ifdef _OPENMP
static void local_getpred_var(size_t nd, size_t p, size_t m, size_t np,
                              xinfo& xi, std::vector<vtree>& tmat,
                              double* px, Rcpp::NumericMatrix& yhat)
{
   int my_rank      = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int h   = nd / thread_count;
   int beg = my_rank * h;
   int end = beg + h - 1;
   if (my_rank == thread_count - 1) end = nd - 1;
   getpred_var(beg, end, p, m, np, xi, tmat, px, yhat);
}
#endif
