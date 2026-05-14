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

        // ----------------------------------------------------------------
        // MH correction for the Sigma draw (ncomp == 1 only).
        //
        // rmixGibbs draws (mu_new, Sigma_new) from the IW conjugate update
        // using residuals theta_i - Delta*(Z_i).  In the (Sigma, delta)
        // parameterization, the exact full conditional for Sigma includes a
        // cross term from Sigma^{1/2} appearing in the mean Delta(Z_i).
        // The IW draw ignores this.  We correct with an MH accept/reject.
        //
        // The log MH ratio is:
        //   log alpha = -sum_i a_i^T b_i - 0.5 * sum_i ||b_i||^2
        // where:
        //   L_new  = rootpi_new^T   (lower-tri Cholesky of Sigma_new^{-1})
        //   L_old  = rootpi_old^T
        //   M      = L_new * L_old^{-1}
        //   b_i    = (M - I) * delta_i
        //   a_i    = L_new * epsilon_i
        //   epsilon_i = theta_i - mu_new - Delta*_old(Z_i)
        //
        // On reject: keep mu_new, revert rootpi to rootpi_prev.
        // ----------------------------------------------------------------
        if (s.oldcomp.size() == 1 && s.rootpi_prev.n_rows > 0) {
            List comp0 = s.oldcomp[0];
            vec  mu_new     = as<vec>(comp0[0]);
            mat  rootpi_new = trimatu(as<mat>(comp0[1]));

            // M = L_new * L_old^{-1} where L = rootpi^T (lower-tri).
            // M^T = L_old^{-T} * L_new^T = rootpi_old^{-1} * rootpi_new
            mat Mt = solve(trimatu(s.rootpi_prev), rootpi_new);
            mat MmI = Mt.t() - eye<mat>(s.nvar, s.nvar);  // M - I

            double log_alpha = 0.0;
            for (int i = 0; i < s.nlgt; i++) {
                vec delta_i(s.nvar), eps_i(s.nvar);
                for (int d = 0; d < s.nvar; d++) {
                    delta_i(d) = s.bart_models[d].f(i);
                    eps_i(d)   = oldbetas(i, d) - mu_new(d) - s.delta_Z(i, d);
                }
                vec b_i = MmI * delta_i;
                vec a_i = rootpi_new.t() * eps_i;  // L_new * eps_i
                log_alpha -= dot(a_i, b_i) + 0.5 * dot(b_i, b_i);
            }

            s.sigma_mh_total++;
            if (log(R::runif(0.0, 1.0)) > std::min(0.0, log_alpha)) {
                // Reject Sigma_new: keep mu_new, revert rootpi.
                List comp_reverted = List::create(
                    Named("mu") = mu_new,
                    Named("rooti") = wrap(s.rootpi_prev));
                s.oldcomp[0] = comp_reverted;
            } else {
                s.sigma_mh_accept++;
            }
        }

        // Store current rootpi for next iteration's MH comparison.
        if (s.oldcomp.size() == 1) {
            List comp0 = s.oldcomp[0];
            s.rootpi_prev = trimatu(as<mat>(comp0[1]));
        }
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
        Named("varprob") = s.varprob,
        Named("sigma_mh_accept_rate") = s.sigma_mh_total > 0
            ? (double)s.sigma_mh_accept / s.sigma_mh_total : NA_REAL
    );
}
