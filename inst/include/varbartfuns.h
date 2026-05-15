/*
 *  Variance-tree sufficient statistics, integrated likelihood, and leaf draws.
 *
 *  Counterpart to heterbartfuns.* for the chi^-2 conjugate model used by
 *  product-of-trees variance ensembles (varbart).
 *
 *  Conjugate model (per leaf k of tree l):
 *    e_i      | s2_lk ~ N(0, s2_lk)            (variance residual)
 *    s2_lk    ~ Inv-chi^2(nu, lambda)          (per-tree prior; calibrated
 *                                              from baseline (nu, lambda)
 *                                              via Pratola et al.)
 *  Sufficient statistics per leaf:
 *    n_k = |k|, S_k = sum_{i in k} e_i^2.
 */
#ifndef GUARD_varbartfuns_h
#define GUARD_varbartfuns_h

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "rn.h"

// Log integrated marginal likelihood for one leaf with sufficient stats (S, n).
// Includes the full leaf-prior normalization (nu/2)*log(nu*lambda) - lgamma(nu/2),
// because in MH birth/death this term contributes a +/- offset that scales with
// nu and therefore *does not* cancel.  (Contrast with heterlh, where the
// analogous dropped term is -log(tau) and is harmless.)
double varlh(double S, double n, double nu, double lambda);

// (n, S) for proposed left/right children of a birth at node nx.
void vargetsuff(tree& x, tree::tree_p nx, size_t v, size_t c,
                xinfo& xi, dinfo& di,
                size_t& nl, double& Sl,
                size_t& nr, double& Sr);

// (n, S) for left/right bots of a death candidate.
void vargetsuff(tree& x, tree::tree_p l, tree::tree_p r,
                xinfo& xi, dinfo& di,
                size_t& nl, double& Sl,
                size_t& nr, double& Sr);

// Draw s^2_lk from posterior Inv-chi^2(nu + n, (nu*lambda + S)/(nu + n)).
double vardrawnodemu(double S, double n, double nu, double lambda, rn& gen);

// (n, S) for all bottom nodes of t (single pass over the data).
void varallsuff(tree& x, xinfo& xi, dinfo& di,
                tree::npv& bnv,
                std::vector<size_t>& nv_leaf,
                std::vector<double>& Sv);

// Draw all leaf parameters of t from their chi^-2 full conditionals.
void vardrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi,
             double nu, double lambda, rn& gen);

#endif
