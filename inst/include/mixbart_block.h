#ifndef MIXBART_BLOCK_H
#define MIXBART_BLOCK_H

#include <RcppArmadillo.h>
#include <vector>
#include <sstream>
#include "bayesm.HART.h"

using namespace Rcpp;
using namespace arma;

struct MixBartState {
    // ---------- Inputs / shape parameters ----------
    int  nlgt, nz, nvar, R, keep;
    bool sparse;
    int  burn;

    // ---------- Tree ensembles (owned) ----------
    std::vector<bart>              bart_models;
    std::vector<std::stringstream> treess;
    arma::cube                     varcount, varprob;

    // ---------- Per-iter outputs consumed by Step E ----------
    arma::mat   delta_Z;        // nlgt x nvar
    Rcpp::List  oldcomp;        // mixture components (current draw)
    arma::vec   oldprob, ind;   // current mixture probs / labels

};

// One-time setup BEFORE the MCMC loop.  Allocates buffers and tree streams.
void mixbart_init(MixBartState& s,
                  int nlgt, int nz, int nvar, int R, int keep,
                  Rcpp::List const& bart_params,
                  arma::vec const& oldprob, arma::vec const& ind);

// Per iteration: rmixGibbs (with first-iter init dispatch) + DART startup +
// update_stdoldbetas + bart::draw + delta_Z recompute via L^{-1} from
// component-specific rooti.  Updates s.{oldcomp, oldprob, ind, delta_Z}.
//
// `oldbetas` is read-only here; the model-specific Step E updates it after
// this function returns.
//
// OWNERSHIP CONTRACT — `pstd_oldbetas_cols`:
//   Caller MUST allocate a writable working buffer (e.g.
//   `mat std_oldbetas = oldbetas;`) and populate `pstd_oldbetas_cols[j]` to
//   point at column `j` of that buffer BEFORE the first call.  This function
//   writes standardized betas into that storage via `update_stdoldbetas` and
//   reads them back through `s.bart_models[j].setdata` (set up in
//   `initializeBART` on rep == 0).  The MixBartState does NOT own this
//   buffer; freeing it before the loop ends is undefined behavior.  See
//   `discussions/2026-05-13-step0-recovery.md` §Lessons for context.
void mixbart_draw_iter(MixBartState& s,
                       arma::mat const& oldbetas,
                       arma::mat const& mubar, arma::mat const& Amu,
                       double nu, arma::mat const& V, arma::vec const& a_mix,
                       double* pZt,
                       std::vector<double*>& pstd_oldbetas_cols,
                       Rcpp::List const& bart_params,
                       int rep, arn& gen);

// Per-`keep` storage: serialize trees, copy varcount/varprob slices.
void mixbart_store(MixBartState& s, int mkeep);

// After the loop: returns a List with `bart_models`, `varcount`, `varprob`.
Rcpp::List mixbart_pack(MixBartState& s);

#endif
