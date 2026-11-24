#ifndef GUARD_vecbartfuns_h
#define GUARD_vecbartfuns_h

#include <vector>
#include "vectree.h"
#include "info.h"

// Leaf model of a vector-leaf tree (appendix_sampler, eq:hs_mean_leaf):
//
//   y_i | m ~ N_D(m, Q_i^{-1}),   m ~ N_D(0, C),
//
// so a leaf's sufficient statistics are n, P = sum Q_i, b = sum Q_i y_i, and
//
//   m | . ~ N_D((C^{-1} + P)^{-1} b, (C^{-1} + P)^{-1}),
//   log M  = -1/2 log|C| - 1/2 log|C^{-1} + P| + 1/2 b'(C^{-1} + P)^{-1} b
//
// up to constants common to every partition (eq:hs_mean_score).  These are
// the vector analogues of getsuff / lh / drawnodemu / allsuff / drmu.

struct vecsuff {
    size_t    n;
    arma::mat P;
    arma::vec b;
    explicit vecsuff(size_t D = 0)
        : n(0), P(D, D, arma::fill::zeros), b(D, arma::fill::zeros) {}
    void add(arma::mat const& Qi, arma::rowvec const& yi) {
        n++; P += Qi; b += Qi * yi.t();
    }
    vecsuff operator+(vecsuff const& o) const {
        vecsuff s(*this); s.n += o.n; s.P += o.P; s.b += o.b; return s;
    }
};

double veclh(vecsuff const& s, arma::mat const& Cinv, double logdetC);
arma::vec vecdrawnodemu(vecsuff const& s, arma::mat const& Cinv, rn& gen);

// Birth proposal: statistics of the two children of nx under split (v, c).
void vecgetsuff(vectree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi,
                dinfo& di, arma::mat const& y, std::vector<arma::mat> const& Q,
                vecsuff& sl, vecsuff& sr);
// Death proposal: statistics of the existing children l and r.
void vecgetsuff(vectree& x, tree::tree_p l, tree::tree_p r, xinfo& xi,
                dinfo& di, arma::mat const& y, std::vector<arma::mat> const& Q,
                vecsuff& sl, vecsuff& sr);
// All bottom nodes in one pass.
void vecallsuff(vectree& x, xinfo& xi, dinfo& di, arma::mat const& y,
                std::vector<arma::mat> const& Q,
                tree::npv& bnv, std::vector<vecsuff>& sv);
// Redraw every leaf from its conditional.
void vecdrmu(vectree& t, xinfo& xi, dinfo& di, arma::mat const& y,
             std::vector<arma::mat> const& Q, arma::mat const& Cinv, rn& gen);

#endif
