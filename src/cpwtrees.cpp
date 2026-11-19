/*
 *  Prediction for serialized tree ensembles, scalar or vector leaves.
 *  Reads the draws, evaluates every tree at the prediction units and
 *  accumulates: by sum for mean and off-diagonal ensembles, by product for
 *  the varbart variance ensembles (whose leaves carry positive scales
 *  s^2_{lk}), and by sum over vector leaves for vecbart.
 *  Stream header: "nd m p" with scalar leaves, "nd m p D" with vector leaves.
 *
 *  Derived from BART's cpwbart (McCulloch and Sparapani, GPL-2).
 */

#include <RcppArmadillo.h>
#include <sstream>
#include "tree.h"
#include "treefuns.h"
#include "vectree.h"

namespace {

// Leaf arity and the per-tree evaluation, which is all that differs.
struct scalar_leaves {
    typedef tree tree_t;
    static tree_t make(size_t) { return tree_t(); }
    static void fit_one(tree_t& t, xinfo& xi, size_t p, size_t np, double* px,
                        arma::mat& fv) {
        fit(t, xi, p, np, px, fv.memptr());
    }
};

struct vector_leaves {
    typedef vectree tree_t;
    static tree_t make(size_t D) { return tree_t(D); }
    static void fit_one(tree_t& t, xinfo& xi, size_t p, size_t np, double* px,
                        arma::mat& fv) {
        t.fit(xi, p, np, px, fv);
    }
};

template <class Leaf>
void getpred(int beg, int end, size_t p, size_t m, size_t np, size_t D,
             bool product, xinfo& xi,
             std::vector<std::vector<typename Leaf::tree_t> >& tmat,
             double* px, arma::cube& yhat)
{
    arma::mat fv(np, D);
    for (int i = beg; i <= end; i++) {
        for (size_t j = 0; j < m; j++) {
            Leaf::fit_one(tmat[i][j], xi, p, np, px, fv);
            if (product) yhat.slice(i) %= fv;
            else         yhat.slice(i) += fv;
        }
    }
}

#ifdef _OPENMP
// Shard the draws across the team; the last thread takes the remainder.
template <class Leaf>
void local_getpred(size_t nd, size_t p, size_t m, size_t np, size_t D,
                   bool product, xinfo& xi,
                   std::vector<std::vector<typename Leaf::tree_t> >& tmat,
                   double* px, arma::cube& yhat)
{
    int my_rank = omp_get_thread_num();
    int threads = omp_get_num_threads();
    int h   = nd / threads;
    int beg = my_rank * h;
    int end = (my_rank == threads - 1) ? (int)nd - 1 : beg + h - 1;
    getpred<Leaf>(beg, end, p, m, np, D, product, xi, tmat, px, yhat);
}
#endif

template <class Leaf>
arma::cube predict_trees(Rcpp::List const& treedraws,
                         Rcpp::NumericMatrix const& xt, bool vector_leaves,
                         bool product, int mc_cores, bool verbose)
{
    Rcpp::CharacterVector itrees(Rcpp::wrap(treedraws["trees"]));
    std::string itv(itrees[0]);
    std::stringstream ttss(itv);

    size_t nd, m, p, D = 1;
    ttss >> nd >> m >> p;
    if (vector_leaves) ttss >> D;
    if (!ttss) Rcpp::stop("cpwtrees: malformed tree stream header.");
    if (verbose)
        Rcpp::Rcout << "draws: " << nd << ", trees: " << m << ", x columns: "
                    << p << ", leaf dim: " << D << std::endl;

    Rcpp::List ixi(Rcpp::wrap(treedraws["cutpoints"]));
    if (p != (size_t)ixi.size())
        Rcpp::stop("cpwtrees: cutpoints do not match the tree stream.");
    xinfo xi;
    xi.resize(p);
    for (size_t i = 0; i < p; i++) {
        Rcpp::NumericVector cutv(ixi[i]);
        xi[i].assign(cutv.begin(), cutv.end());
    }

    // xt is p x npred: the R side transposes, so one column per unit.
    if ((size_t)xt.nrow() != p)
        Rcpp::stop("cpwtrees: x must have %d rows (one column per unit).",
                   (int)p);
    size_t np = xt.ncol();

    typedef typename Leaf::tree_t tree_t;
    std::vector<std::vector<tree_t> > tmat(nd);
    for (size_t i = 0; i < nd; i++) {
        tmat[i].assign(m, Leaf::make(D));
        for (size_t j = 0; j < m; j++) {
            ttss >> tmat[i][j];
            if (!ttss)
                Rcpp::stop("cpwtrees: malformed tree stream at draw %d, tree %d.",
                           (int)i + 1, (int)j + 1);
        }
    }

    arma::cube yhat(np, D, nd);
    yhat.fill(product ? 1.0 : 0.0);
    double* px = const_cast<double*>(&xt(0, 0));

#ifndef _OPENMP
    getpred<Leaf>(0, nd - 1, p, m, np, D, product, xi, tmat, px, yhat);
#else
    if (mc_cores == 1) {
        getpred<Leaf>(0, nd - 1, p, m, np, D, product, xi, tmat, px, yhat);
    } else {
#pragma omp parallel num_threads(mc_cores)
        local_getpred<Leaf>(nd, p, m, np, D, product, xi, tmat, px, yhat);
    }
#endif
    return yhat;
}

}  // namespace

// npred x leaf_dim x ndraws; `leaf_dim > 1` reads vector leaves.
// [[Rcpp::export]]
arma::cube cpwtrees(Rcpp::List const& treedraws,
                    Rcpp::NumericMatrix const& xt,
                    int leaf_dim, int mc_cores, bool product, bool verbose)
{
    if (leaf_dim > 1)
        return predict_trees<vector_leaves>(treedraws, xt, true, product,
                                            mc_cores, verbose);
    return predict_trees<scalar_leaves>(treedraws, xt, false, product,
                                        mc_cores, verbose);
}
