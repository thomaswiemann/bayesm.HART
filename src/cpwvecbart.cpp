/*
 *  Prediction for vector-leaf ensembles (vecbart).  Mirrors cpwbart: read the
 *  serialized draws, evaluate every tree at the prediction units, sum.
 *  Stream header: "nd m p D".  Returns an npred x D x nd cube.
 */
#include <RcppArmadillo.h>
#include <sstream>
#include "vectree.h"

// [[Rcpp::export]]
arma::cube cpwvecbart(Rcpp::List const& treedraws, arma::mat const& xt)
{
   // xt: p x npred column stack, one column per prediction unit.
   Rcpp::CharacterVector itrees = treedraws["trees"];
   std::string itv = Rcpp::as<std::string>(itrees[0]);
   std::stringstream ttss(itv);

   size_t nd, m, p, D;
   ttss >> nd >> m >> p >> D;
   if (!ttss) Rcpp::stop("cpwvecbart: malformed tree stream header.");

   Rcpp::List ixi = treedraws["cutpoints"];
   if ((size_t)ixi.size() != p)
      Rcpp::stop("cpwvecbart: cutpoints do not match the tree stream.");
   xinfo xi(p);
   for (size_t i = 0; i < p; i++) {
      Rcpp::NumericVector cutv(ixi[i]);
      xi[i].assign(cutv.begin(), cutv.end());
   }
   if (xt.n_rows != p)
      Rcpp::stop("cpwvecbart: x must have %d rows (one column per unit).", (int)p);

   size_t np = xt.n_cols;
   arma::cube yhat(np, D, nd, arma::fill::zeros);
   double* px = const_cast<double*>(xt.memptr());
   arma::mat ftmp;
   vectree t(D);
   for (size_t i = 0; i < nd; i++) {
      for (size_t j = 0; j < m; j++) {
         ttss >> t;
         if (!ttss)
            Rcpp::stop("cpwvecbart: malformed tree stream at draw %d, tree %d.",
                       (int)i + 1, (int)j + 1);
         t.fit(xi, p, np, px, ftmp);
         yhat.slice(i) += ftmp;
      }
   }
   return yhat;
}
