/*
 *  Birth/death Metropolis-Hastings step for varying-coefficient trees.
 *
 *  Mirrors bd: tree-structure prior, proposal mechanics (bprop / dprop /
 *  getpb), and DART bookkeeping come unchanged from bartfuns; the leaf
 *  statistics are vcgetsuff and the score is vclh.  Admissibility is the raw
 *  count floor pi.nmin only; it depends on Z, not on the regressors.
 */
#ifndef GUARD_vcbd_h
#define GUARD_vcbd_h

#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "vcbartfuns.h"

bool vcbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double* a, double* var,
          std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen);

#endif
