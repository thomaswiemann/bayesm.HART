#ifndef BART_PACK_H
#define BART_PACK_H

#include <RcppArmadillo.h>
#include <sstream>
#include <vector>

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

#endif
