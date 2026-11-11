/*
 *  Birth/death Metropolis-Hastings step for vector-leaf trees.
 *
 *  Mirrors bd: tree-structure prior, proposal mechanics (bprop / dprop /
 *  getpb), and DART bookkeeping are inherited unchanged from bartfuns; the
 *  leaf marginal (veclh) and leaf draws (vecdrawnodemu) come from
 *  vecbartfuns.  veclh keeps the per-leaf prior normalisation, so no extra
 *  log(sigma) term is needed in the acceptance ratio.
 */
#ifndef GUARD_vecbd_h
#define GUARD_vecbd_h

#include <RcppArmadillo.h>   // before common.h (Rcpp.h must not come first)
#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "vecbartfuns.h"

bool vecbd(vectree& x, xinfo& xi, dinfo& di, pinfo& pi,
           arma::mat const& y, std::vector<arma::mat> const& Q,
           arma::mat const& Cinv, double logdetC,
           std::vector<size_t>& nv, std::vector<double>& pv,
           bool aug, rn& gen);

#endif
