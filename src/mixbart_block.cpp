#include "mixbart_block.h"
#include "bart_pack.h"
#include "sigma_gibbs.h"

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
    s.store_trees = bart_params.containsElementNamed("store_trees")
        ? as<bool>(bart_params["store_trees"])
        : true;
    s.has_unique_map = false;
    
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

    if (bart_params.containsElementNamed("z_index")) {
        IntegerVector z_index = bart_params["z_index"];
        if ((int)z_index.size() != nlgt)
            stop("mixbart_init: z_index length must match nreg.");
        int n_unique = bart_params.containsElementNamed("n_unique")
            ? as<int>(bart_params["n_unique"])
            : 0;
        if (n_unique <= 0)
            n_unique = Rcpp::max(z_index);
        if (n_unique <= 0)
            stop("mixbart_init: n_unique must be positive when z_index is supplied.");

        std::vector<int> first_idx(n_unique, -1);
        for (int i = 0; i < nlgt; i++) {
            int u = z_index[i] - 1; // 1-based from R
            if (u < 0 || u >= n_unique)
                stop("mixbart_init: z_index values must be in 1..n_unique.");
            if (first_idx[u] < 0) first_idx[u] = i;
        }
        s.unique_first_idx.set_size(n_unique);
        for (int u = 0; u < n_unique; u++) {
            if (first_idx[u] < 0)
                stop("mixbart_init: every unique index must appear at least once.");
            s.unique_first_idx[u] = static_cast<uword>(first_idx[u]);
        }
        s.delta_z_unique_draws.set_size(n_unique, nvar, R / keep);
        s.delta_z_unique_draws.zeros();
        s.has_unique_map = true;
    }
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
        // Cholesky-Gibbs (mu, Sigma) update for each mixture component.
        // Order matches rmixGibbs:
        //   1) update components given OLD labels
        //   2) draw labels given NEW components
        //   3) draw probs given NEW labels

        mat delta(s.nlgt, s.nvar);
        for (int i = 0; i < s.nlgt; i++) {
            for (int d = 0; d < s.nvar; d++) {
                delta(i, d) = s.bart_models[d].f(i);
            }
        }

        vec mubar0 = mubar.row(0).t();
        double Amu_scalar = Amu(0, 0);
        int ncomp = s.oldcomp.size();
        List new_oldcomp(ncomp);

        for (int k = 0; k < ncomp; k++) {
            uvec idx_k = find(s.ind == (k + 1)); // ind is 1-indexed
            int n_k = idx_k.n_elem;

            if (n_k > 0) {
                List comp_k = s.oldcomp[k];
                vec mu_k = as<vec>(comp_k[0]);
                mat rooti_k = as<mat>(comp_k[1]);
                mat L_k = rooti_k.t();

                mat theta_k = oldbetas.rows(idx_k);
                mat delta_k = delta.rows(idx_k);

                sigma_mu_block_gibbs(L_k, mu_k, theta_k, delta_k,
                                     V, nu, mubar0, Amu_scalar, gen);

                new_oldcomp[k] = List::create(
                    Named("mu") = NumericVector(mu_k.begin(), mu_k.end()),
                    Named("rooti") = L_k.t());
            } else {
                // Empty component fallback: prior draw (matches rmixGibbs path).
                mat S = solve(trimatu(chol(V)), eye(s.nvar, s.nvar));
                S = S * trans(S);
                List rw = rwishart(nu, S);
                mat IW = as<mat>(rw["IW"]);
                mat CI = as<mat>(rw["CI"]);
                mat rooti_k = solve(trimatu(chol(IW)), eye(s.nvar, s.nvar));
                vec r(s.nvar);
                for (int d = 0; d < s.nvar; d++) r(d) = gen.normal();
                vec mu_k = mubar0 + (CI * r) / std::sqrt(Amu_scalar);

                new_oldcomp[k] = List::create(
                    Named("mu") = NumericVector(mu_k.begin(), mu_k.end()),
                    Named("rooti") = rooti_k);
            }
        }
        s.oldcomp = new_oldcomp;

        s.ind = drawLabelsFromComps(oldbetas - s.delta_Z, s.oldprob, s.oldcomp);
        s.oldprob = drawPFromLabels(a_mix, s.ind);
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
        if (s.store_trees) {
            for (size_t j = 0; j < s.bart_models[i].getm(); j++) {
                s.treess[i] << s.bart_models[i].gettree(j);
            }
        }
        s.varcount.slice(mkeep - 1).col(i) = conv_to<vec>::from(s.bart_models[i].getnv());
        s.varprob.slice(mkeep - 1).col(i) = conv_to<vec>::from(s.bart_models[i].getpv());
    }
    if (s.has_unique_map) {
        s.delta_z_unique_draws.slice(mkeep - 1) = s.delta_Z.rows(s.unique_first_idx);
    }
}

Rcpp::List mixbart_pack(MixBartState& s) {
    List out = List::create(
        Named("varcount") = s.varcount,
        Named("varprob") = s.varprob
    );
    if (s.store_trees) {
        std::vector<std::stringstream> tss(s.nvar);
        for(int i = 0; i < s.nvar; ++i) tss[i] << s.treess[i].str();
        out["bart_models"] = pack_bart_list_(s.nvar, s.bart_models, tss);
    }
    if (s.has_unique_map) {
        out["DeltaZ_unique_draws"] = s.delta_z_unique_draws;
    }
    return out;
}
