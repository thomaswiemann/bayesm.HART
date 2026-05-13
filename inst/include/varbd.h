/*
 *  Birth/death Metropolis-Hastings step for variance trees.
 *
 *  Mirrors heterbd, but uses the chi^-2 conjugate marginal likelihood (varlh)
 *  and the chi^-2 leaf draws (vardrawnodemu) from varbartfuns.
 *
 *  Tree-structure prior, proposal mechanics (bprop / dprop / getpb), and
 *  optional DART support are inherited unchanged from bartfuns.
 */
#ifndef GUARD_varbd_h
#define GUARD_varbd_h

#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "varbartfuns.h"

bool varbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi,
           double nu, double lambda,
           std::vector<size_t>& nv, std::vector<double>& pv,
           bool aug, rn& gen);

#endif
