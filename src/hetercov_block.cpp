#include "hetercov_block.h"
#include "bart_pack.h"

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// hetercov_init
// ---------------------------------------------------------------------------
void hetercov_init(HeterCovState& s,
                   int nlgt, int nz, int nvar, int R, int keep,
                   List const& bart_params,
                   List const& var_params,
                   List const& phi_params) {
    s.nlgt = nlgt; s.nz = nz; s.nvar = nvar; s.R = R; s.keep = keep;
    s.useFullCov = (phi_params.size() > 0);
    s.n_pairs    = s.useFullCov ? (nvar * (nvar - 1)) / 2 : 0;

    // ---- DART configuration --------------------------------------------------
    // Each ensemble (mean / variance / phi) is independently sparse-eligible.
    // Defaults: dense (sparse = false, burn = 0).  startdart() is called once
    // per ensemble at iteration `burn + 1` inside hetercov_draw_iter.
    s.sparse     = false;  s.burn     = 0;
    s.var_sparse = false;  s.var_burn = 0;
    s.phi_sparse = false;  s.phi_burn = 0;
    if (bart_params.containsElementNamed("sparse")) {
        s.sparse = as<bool>(bart_params["sparse"]);
        if (s.sparse && bart_params.containsElementNamed("burn"))
            s.burn = as<int>(bart_params["burn"]);
    }
    if (var_params.containsElementNamed("sparse")) {
        s.var_sparse = as<bool>(var_params["sparse"]);
        if (s.var_sparse && var_params.containsElementNamed("burn"))
            s.var_burn = as<int>(var_params["burn"]);
    }
    if (s.useFullCov && phi_params.containsElementNamed("sparse")) {
        s.phi_sparse = as<bool>(phi_params["sparse"]);
        if (s.phi_sparse && phi_params.containsElementNamed("burn"))
            s.phi_burn = as<int>(phi_params["burn"]);
    }

    // ---- Output cubes / matrices --------------------------------------------
    s.varcount    .set_size(nz, nvar, R / keep); s.varcount    .zeros();
    s.varprob     .set_size(nz, nvar, R / keep); s.varprob     .zeros();
    s.var_varcount.set_size(nz, nvar, R / keep); s.var_varcount.zeros();
    s.var_varprob .set_size(nz, nvar, R / keep); s.var_varprob .zeros();
    s.mu_draw     .zeros(R / keep, nvar);

    // ---- Per-iter state -----------------------------------------------------
    s.mu_post.zeros(nvar);
    s.delta_Z.zeros(nlgt, nvar);

    // ---- Mean-tree streams (headers) ----------------------------------------
    s.treess.resize(nvar);
    size_t mean_num_trees = as<size_t>(bart_params["num_trees"]);
    for (int k = 0; k < nvar; k++) {
        s.treess[k].precision(10);
        s.treess[k] << R / keep << " " << mean_num_trees << " " << nz << std::endl;
    }

    // ---- Variance-tree streams + buffers -----------------------------------
    s.var_treess.resize(nvar);
    size_t var_num_trees = as<size_t>(var_params["num_trees"]);
    for (int j = 0; j < nvar; j++) {
        s.var_treess[j].precision(10);
        s.var_treess[j] << R / keep << " " << var_num_trees << " " << nz << std::endl;
    }
    s.eta_sq_buf.zeros(nlgt, nvar);
    s.peta_sq_cols.resize(nvar);
    for (int j = 0; j < nvar; j++) s.peta_sq_cols[j] = s.eta_sq_buf.colptr(j);

    // ---- Phi-tree streams + buffers (Phase 2 only) -------------------------
    if (s.useFullCov) {
        s.Y_phi.zeros(nlgt, s.n_pairs);
        s.sigma_phi_buf.zeros(nlgt);
        s.py_phi_cols.resize(nvar);
        size_t pair_idx = 0;
        for (int j = 1; j < nvar; j++) {
            s.py_phi_cols[j].resize(j);
            for (int k = 0; k < j; k++) {
                s.py_phi_cols[j][k] = s.Y_phi.colptr(pair_idx++);
            }
        }
        s.phi_treess.resize(nvar);
        size_t phi_num_trees = as<size_t>(phi_params["num_trees"]);
        for (int j = 1; j < nvar; j++) {
            s.phi_treess[j].resize(j);
            for (int k = 0; k < j; k++) {
                s.phi_treess[j][k].precision(10);
                s.phi_treess[j][k] << R / keep << " " << phi_num_trees << " " << nz << std::endl;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// hetercov_draw_iter
// ---------------------------------------------------------------------------
void hetercov_draw_iter(HeterCovState& s,
                        mat const& oldbetas,
                        vec const& mubar0,
                        mat const& Amu,
                        double* pZt,
                        std::vector<double*>& pstd_oldbetas_cols,
                        List const& bart_params,
                        List const& var_params,
                        List const& phi_params,
                        int rep, arn& gen) {

    const int nlgt = s.nlgt;
    const int nvar = s.nvar;
    const int nz   = s.nz;

    // ====================================================================
    // Step A: draw mu | {theta_i}, Sigma(.).
    //   rep == 0: baseline initialization (mu = colMeans(theta), eta_sq from
    //             (theta - mu)^2, var/phi/mean trees brought up).
    //   rep >= 1: drawMuHeterCov from conjugate normal posterior.
    // ====================================================================
    if (rep == 0) {
        s.mu_post = (sum(oldbetas, 0) / static_cast<double>(nlgt)).t();
        for (int j = 0; j < nvar; j++) {
            for (int i = 0; i < nlgt; i++) {
                s.eta_sq_buf(i, j) = std::pow(oldbetas(i, j) - s.mu_post(j), 2.0);
            }
        }
        s.var_models = initializeVarBART(pZt, nlgt, nz, s.peta_sq_cols, var_params, gen);
        if (s.useFullCov) {
            // phi_models initialized BEFORE update_stdoldbetas_het so the
            // Mahalanobis transform sees the full Cholesky from the start.
            // Y buffers are zero at init -> phi_models predict 0 (their prior
            // mean), the correct starting point.
            s.phi_models = initializePhiBART(pZt, nlgt, nz, s.py_phi_cols, phi_params, gen);
        }
        update_stdoldbetas_het(oldbetas, pstd_oldbetas_cols,
                               s.mu_post, s.var_models, s.phi_models);
        s.bart_models = initializeBART(pZt, nlgt, nz, pstd_oldbetas_cols, bart_params, gen);
    } else {
        s.mu_post = drawMuHeterCov(oldbetas, mubar0, Amu,
                                   s.var_models, s.phi_models, gen);
    }

    // ====================================================================
    // DART burn-in (mean and variance trees independently configured).
    // ====================================================================
    if (s.sparse && rep == s.burn + 1) {
        for (int j = 0; j < nvar; j++) s.bart_models[j].startdart();
    }
    if (s.var_sparse && rep == s.var_burn + 1) {
        for (int j = 0; j < nvar; j++) s.var_models[j].startdart();
    }
    if (s.useFullCov && s.phi_sparse && rep == s.phi_burn + 1) {
        for (int j = 1; j < nvar; j++)
            for (int k = 0; k < j; k++)
                s.phi_models[j][k].startdart();
    }

    // ====================================================================
    // Step C: variance trees d_j | mu, theta, current delta, phi.
    //   Compute eta_i = L (theta_i - mu) - D^{1/2} delta_i once per unit
    //   using the CURRENT bart_models (delta) and phi_models, then draw.
    //   Working in the (delta, d, phi) tuple of fundamental tree outputs
    //   (not delta_Z) avoids dropping cross terms when both delta_Z and phi
    //   are non-zero -- see discussions/2026-05-12-cholesky-trees-implementation.md
    //   "subtle pitfalls" section.
    // ====================================================================
    for (int i = 0; i < nlgt; i++) {
        cov::cov_eval ev{s.var_models, s.phi_models, (size_t)i};
        vec eta_i = oldbetas.row(i).t() - s.mu_post;
        cov::apply_L(eta_i, ev);                            // eta_i := L (theta_i - mu)
        for (int j = 0; j < nvar; j++) {
            double sqrt_dj = std::sqrt(s.var_models[j].f(i));
            double delta_j = s.bart_models[j].f(i);
            double e_ij    = eta_i[j] - sqrt_dj * delta_j;
            s.eta_sq_buf(i, j) = e_ij * e_ij;
        }
    }
    for (int j = 0; j < nvar; j++) s.var_models[j].draw(gen);

    // ====================================================================
    // Step D' (Phase 2 only): phi_jk | mu, theta, current delta, NEW d, OTHER phi.
    //   Per discussions/2026-05-12-phibart-vs-heterbart.md, fit each phi_jk
    //   via the varying-coefficient -> heterbart transform
    //
    //     r_i = phi_jk(Z_i) * c_i^{(k)} + eta_i^{(j)}     (data model)
    //     Y_i = r_i / c_i^{(k)},  sigma_i = sqrt(d_j) / |c_i^{(k)}|.
    //
    //     r_i = (theta_j - mu_j) - sqrt(d_j) delta_j
    //           - sum_{m<j, m != k} phi_jm (theta_m - mu_m)
    //
    //   c_i^{(k)} clamp: |c| >= 1e-8.  ESS guard inside heterbd ensures
    //   clamped observations contribute negligibly to splits.
    // ====================================================================
    if (s.useFullCov) {
        const double eps = 1e-8;
        for (int j = 1; j < nvar; j++) {
            for (int k = 0; k < j; k++) {
                double* y_buf = s.py_phi_cols[j][k];
                for (int i = 0; i < nlgt; i++) {
                    double c_ik = oldbetas(i, k) - s.mu_post(k);
                    if (std::abs(c_ik) < eps)
                        c_ik = (c_ik >= 0 ? eps : -eps);

                    double sqrt_dj = std::sqrt(s.var_models[j].f(i));
                    double delta_j = s.bart_models[j].f(i);
                    double r_i = oldbetas(i, j) - s.mu_post(j) - sqrt_dj * delta_j;
                    for (int m = 0; m < j; m++) {
                        if (m == k) continue;
                        r_i -= s.phi_models[j][m].f(i)
                             * (oldbetas(i, m) - s.mu_post(m));
                    }

                    y_buf[i]            = r_i / c_ik;
                    s.sigma_phi_buf[i]  = sqrt_dj / std::abs(c_ik);
                }
                s.phi_models[j][k].draw(s.sigma_phi_buf.memptr(), gen);
            }
        }
    }

    // ====================================================================
    // Step B: mean trees delta | mu, NEW d, NEW phi, theta.
    //   Standardize tilde_theta_i = D(Z_i)^{-1/2} L(Z_i) (theta_i - mu)
    //   using the FRESHLY updated d and phi, then backfit.  Mirrors Wiemann
    //   (2025) Step 4 (BART | mu, Sigma, theta).
    // ====================================================================
    update_stdoldbetas_het(oldbetas, pstd_oldbetas_cols,
                           s.mu_post, s.var_models, s.phi_models);
    for (int j = 0; j < nvar; j++) s.bart_models[j].draw(1.0, gen);

    // ====================================================================
    // delta_Z(i, .) = Sigma(Z_i)^{1/2} delta(Z_i) = L^{-1} D^{1/2} delta,
    // evaluated with CURRENT mu, d, phi, delta.  Step E consumes this
    // directly (no recompute).
    // ====================================================================
    vec pred_i(nvar);
    for (int i = 0; i < nlgt; i++) {
        for (int j = 0; j < nvar; j++) pred_i(j) = s.bart_models[j].f(i);
        cov::cov_eval ev{s.var_models, s.phi_models, (size_t)i};
        vec d_i = cov::sigma_sqrt_apply(pred_i, ev);
        s.delta_Z.row(i) = d_i.t();
    }
}

// ---------------------------------------------------------------------------
// hetercov_store
// ---------------------------------------------------------------------------
void hetercov_store(HeterCovState& s, int mkeep) {
    s.mu_draw.row(mkeep - 1) = s.mu_post.t();
    for (int j = 0; j < s.nvar; j++) {
        for (size_t k = 0; k < s.bart_models[j].getm(); k++)
            s.treess[j]     << s.bart_models[j].gettree(k);
        for (size_t k = 0; k < s.var_models[j].getm();  k++)
            s.var_treess[j] << s.var_models[j].gettree(k);
        s.varcount    .slice(mkeep - 1).col(j) = conv_to<vec>::from(s.bart_models[j].getnv());
        s.varprob     .slice(mkeep - 1).col(j) = conv_to<vec>::from(s.bart_models[j].getpv());
        s.var_varcount.slice(mkeep - 1).col(j) = conv_to<vec>::from(s.var_models[j].getnv());
        s.var_varprob .slice(mkeep - 1).col(j) = conv_to<vec>::from(s.var_models[j].getpv());
    }
    if (s.useFullCov) {
        for (int j = 1; j < s.nvar; j++) {
            for (int k = 0; k < j; k++) {
                for (size_t t = 0; t < s.phi_models[j][k].getm(); t++) {
                    s.phi_treess[j][k] << s.phi_models[j][k].gettree(t);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// hetercov_pack
// ---------------------------------------------------------------------------
List hetercov_pack(HeterCovState& s) {
    SEXP phi_models_pkg = s.useFullCov
        ? (SEXP) pack_bart_jagged_(s.nvar, s.phi_models, s.phi_treess)
        : R_NilValue;

    return List::create(
        Named("bart_models")  = pack_bart_list_(s.nvar, s.bart_models, s.treess),
        Named("var_models")   = pack_bart_list_(s.nvar, s.var_models,  s.var_treess),
        Named("phi_models")   = phi_models_pkg,
        Named("mu_draw")      = s.mu_draw,
        Named("varcount")     = s.varcount,
        Named("varprob")      = s.varprob,
        Named("var_varcount") = s.var_varcount,
        Named("var_varprob")  = s.var_varprob);
}
