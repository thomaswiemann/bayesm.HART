#include "mixbart_block.h"
#include "bart_pack.h"

using namespace Rcpp;
using namespace arma;

void mixbart_init(MixBartState& s,
                  int nlgt, int nz, int nvar, int R, int keep,
                  Rcpp::List const& bart_params,
                  arma::vec const& oldprob, arma::vec const& ind) {
    s.nlgt = nlgt;
    s.nz = nz;
    s.nvar = nvar;
    s.R = R;
    s.keep = keep;
    s.oldprob = oldprob;
    s.ind = ind;
    
    s.delta_Z.zeros(nlgt, nvar);
    
    s.treess.resize(nvar);
    size_t num_trees = bart_params["num_trees"];
    for (int k = 0; k < nvar; k++) {
        s.treess[k].precision(10);
        s.treess[k] << R / keep << " " << num_trees << " " << nz << std::endl;
    }
    
    s.varcount.set_size(nz, nvar, R / keep);
    s.varprob.set_size(nz, nvar, R / keep);
    s.varcount.zeros(nz, nvar, R / keep);
    s.varprob.zeros(nz, nvar, R / keep); 
    s.sparse = bart_params["sparse"];
    s.burn = 0;
    if (s.sparse) s.burn = bart_params["burn"];
}

void mixbart_draw_iter(MixBartState& s,
                       arma::mat const& oldbetas,
                       arma::mat const& mubar, arma::mat const& Amu,
                       double nu, arma::mat const& V, arma::vec const& a_mix,
                       double* pZt,
                       std::vector<double*>& pstd_oldbetas_cols,
                       Rcpp::List const& bart_params,
                       int rep, arn& gen) {
    
    List mgout;
    if (rep == 0) {
        mgout = rmixGibbs(oldbetas, mubar, Amu, nu, V, a_mix, s.oldprob, s.ind);
        s.oldcomp = mgout["comps"];
        s.oldprob = as<vec>(mgout["p"]);
        s.ind = as<vec>(mgout["z"]);
        update_stdoldbetas(oldbetas, pstd_oldbetas_cols, s.ind, s.oldcomp);
        s.bart_models = initializeBART(pZt, s.nlgt, s.nz, pstd_oldbetas_cols, bart_params, gen);
    } else {
        mgout = rmixGibbs(oldbetas - s.delta_Z, mubar, Amu, nu, V, a_mix, s.oldprob, s.ind);
        s.oldcomp = mgout["comps"];
        s.oldprob = as<vec>(mgout["p"]);
        s.ind = as<vec>(mgout["z"]);
    }
    
    if (s.sparse && rep == s.burn + 1) {
        for (int i = 0; i < s.nvar; i++) {
            s.bart_models[i].startdart();
        }
    }
    
    // Normalize oldbetas
    update_stdoldbetas(oldbetas, pstd_oldbetas_cols, s.ind, s.oldcomp);
    
    // Draw bart
    for (int i = 0; i < s.nvar; i++) {
        s.bart_models[i].draw(1.0, gen);
    }
    
    // Compute delta_Z
    vec pred_i(s.nvar);
    vec scaled_pred(s.nvar);
    for (int i = 0; i < s.nlgt; i++) {
        for (int j = 0; j < s.nvar; j++) {
            pred_i(j) = s.bart_models[j].f(i);
        }
        List comp_i = s.oldcomp[s.ind[i] - 1];
        mat rootii = trans(as<mat>(comp_i[1])); 
        scaled_pred = solve(trimatl(rootii), pred_i);
        s.delta_Z.row(i) = trans(scaled_pred);
    }
}

void mixbart_store(MixBartState& s, int mkeep) {
    for (int i = 0; i < s.nvar; i++) {
        for (size_t j = 0; j < s.bart_models[i].getm(); j++) {
            s.treess[i] << s.bart_models[i].gettree(j);
        }
        s.varcount.slice(mkeep - 1).col(i) = conv_to<vec>::from(s.bart_models[i].getnv());
        s.varprob.slice(mkeep - 1).col(i) = conv_to<vec>::from(s.bart_models[i].getpv());
    }
}

Rcpp::List mixbart_pack(MixBartState& s) {
    std::vector<std::stringstream> tss(s.nvar);
    for(int i = 0; i < s.nvar; ++i) tss[i] << s.treess[i].str();
    
    return List::create(
        Named("bart_models") = pack_bart_list_(s.nvar, s.bart_models, tss),
        Named("varcount") = s.varcount,
        Named("varprob") = s.varprob
    );
}
