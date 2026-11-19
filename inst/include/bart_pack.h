#ifndef BART_PACK_H
#define BART_PACK_H

#include <RcppArmadillo.h>
#include <sstream>
#include <vector>
#include <limits>
#include "bart.h"
#include "vecbart.h"

using namespace Rcpp;

// Package a single bart-like model (any class exposing getxinfo() and
// gettree()) into the {treedraws: {cutpoints, trees}} List that the existing
// R-side post-processing expects.
template <class B>
static List pack_bart_entry_(B& model, std::stringstream& tss) {
    xinfo& xi = model.getxinfo();
    Rcpp::List xiret(xi.size());
    for (size_t j = 0; j < xi.size(); j++) {
        Rcpp::NumericVector vtemp(xi[j].size());
        std::copy(xi[j].begin(), xi[j].end(), vtemp.begin());
        xiret[j] = vtemp;
    }
    return List::create(
        Named("treedraws") = List::create(
            Named("cutpoints") = xiret,
            Named("trees")     = Rcpp::CharacterVector(tss.str())
        )
    );
}

// Vector-leaf ensemble (vecbart): pack_bart_entry_ plus the leaf dimension
// (stream header "nkeep m p D").
template <class B>
static List pack_vecbart_entry_(B& model, std::stringstream& tss) {
    List out = pack_bart_entry_(model, tss);
    List td  = out["treedraws"];
    td["leaf_dim"] = (int)model.getD();
    out["treedraws"] = td;
    return out;
}

// Flat List(nvar) of per-dimension model entries (bart, varbart).
template <class B>
static List pack_bart_list_(int nvar,
                            std::vector<B>& models,
                            std::vector<std::stringstream>& tss) {
    List out(nvar);
    for (int i = 0; i < nvar; i++) out[i] = pack_bart_entry_(models[i], tss[i]);
    return out;
}

// Jagged List(D) of List(j) entries for phi-trees: out[0] is NULL, out[j>0]
// has length j with entries phi[j][k] for k = 0, ..., j-1.
template <class B>
static List pack_bart_jagged_(int nvar,
                              std::vector<std::vector<B>>& models,
                              std::vector<std::vector<std::stringstream>>& tss) {
    List out(nvar);
    out[0] = R_NilValue;
    for (int j = 1; j < nvar; j++) {
        List inner(j);
        for (int k = 0; k < j; k++) {
            inner[k] = pack_bart_entry_(models[j][k], tss[j][k]);
        }
        out[j] = inner;
    }
    return out;
}

// Scalar streams omit the leaf dimension; vector streams carry D.
static size_t state_leaf_dim_(bart&) { return 0; }
static size_t state_leaf_dim_(vecbart& model) { return model.getD(); }

// One terminal ensemble, with the existing single-draw stream header.
template <class B>
static List pack_bart_state_(B& model) {
    std::stringstream ss;
    ss.precision(std::numeric_limits<double>::max_digits10);
    ss << 1 << " " << model.getm() << " " << model.getp();
    size_t D = state_leaf_dim_(model);
    if (D) ss << " " << D;
    ss << std::endl;
    for (size_t j = 0; j < model.getm(); j++) ss << model.gettree(j);
    List xi(model.getxinfo().size());
    for (size_t j = 0; j < model.getxinfo().size(); j++)
        xi[j] = wrap(model.getxinfo()[j]);
    return List::create(
        Named("cutpoints") = xi, Named("trees") = ss.str(),
        Named("nv") = NumericVector(model.getnv().begin(), model.getnv().end()),
        Named("pv") = wrap(model.getpv()), Named("lpv") = wrap(model.getlpv()),
        Named("theta") = model.gettheta(), Named("dartOn") = model.getdartOn());
}

// Target prior and data are allocated before terminal state is loaded.
template <class B>
static void restore_bart_state_(B& model, List const& state) {
    List cuts = state["cutpoints"];
    if ((size_t)cuts.size() != model.getp()) stop("cutpoint dimension mismatch");
    std::stringstream ss(as<std::string>(state["trees"]));
    std::string line;
    std::getline(ss, line);
    std::istringstream header(line);
    size_t nkeep, m, p, D = state_leaf_dim_(model), saved_D = 0;
    if (!(header >> nkeep >> m >> p)) stop("failed tree stream read");
    if (nkeep != 1 || m != model.getm()) stop("tree count mismatch");
    if (p != model.getp()) stop("cutpoint dimension mismatch");
    if (D && (!(header >> saved_D) || saved_D != D))
        stop("leaf dimension mismatch");
    header >> std::ws;
    if (!header.eof()) stop("leaf dimension mismatch");
    xinfo xi(p);
    for (size_t j = 0; j < p; j++) xi[j] = as<std::vector<double>>(cuts[j]);
    model.setxinfo(xi);
    for (size_t j = 0; j < m; j++)
        if (!(ss >> model.gettree(j))) stop("failed tree stream read");
    ss >> std::ws;
    if (!ss.eof()) stop("tree count mismatch");
    std::vector<size_t> nv = as<std::vector<size_t>>(state["nv"]);
    std::vector<double> pv = as<std::vector<double>>(state["pv"]);
    std::vector<double> lpv = as<std::vector<double>>(state["lpv"]);
    if (nv.size() != p || pv.size() != p || (!lpv.empty() && lpv.size() != p))
        stop("DART dimension mismatch");
    model.getnv() = nv;
    model.getpv() = pv;
    model.getlpv() = lpv;
    model.setdartstate(as<double>(state["theta"]), as<bool>(state["dartOn"]));
    model.rebuildfit();
}

#endif
