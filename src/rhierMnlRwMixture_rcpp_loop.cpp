#include "bayesm.HART.h"

// Package a single bart-like model (any class exposing getxinfo() and
// gettree()) into the {treedraws: {cutpoints, trees}} List that the existing
// R-side post-processing expects.  Used directly for jagged phi ensembles and
// indirectly (via pack_bart_list_) for flat bart / varbart ensembles.
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

//FUNCTION SPECIFIC TO MAIN FUNCTION------------------------------------------------------
//[[Rcpp::export]]
double llmnl_con(vec const& betastar, vec const& y, mat const& X, vec const& SignRes = NumericVector::create(0)) {

    // Wayne Taylor 7/8/2016

    // Evaluates log-likelihood for the multinomial logit model WITH SIGN CONSTRAINTS
    // NOTE: this is exported only because it is used in the shell .R function, it will not be available to users

    //Reparameterize betastar to beta to allow for sign restrictions
    vec beta = betastar;

    //The default SignRes vector is a single element vector containing a zero
    //any() returns true if any elements of SignRes are non-zero
    if (any(SignRes)) {
        uvec signInd = find(SignRes != 0);
        beta.elem(signInd) = SignRes.elem(signInd) % exp(beta.elem(signInd));  //% performs element-wise multiplication
    }

    int n = y.size();
    int j = X.n_rows / n;
    mat Xbeta = X * beta;

    vec xby = zeros<vec>(n);
    vec denom = zeros<vec>(n);

    for (int i = 0; i < n;i++) {
        for (int p = 0;p < j;p++) denom[i] = denom[i] + exp(Xbeta[i * j + p]);
        xby[i] = Xbeta[i * j + y[i] - 1];
    }

    return(sum(xby - log(denom)));
}

//mnlRwMetropOnce=
//function(y,X,oldbeta,oldll,s,inc.root,betabar,rootpi){ 
//#
//# function to execute rw metropolis for the MNL
//# y is n vector with element = 1,...,j indicating which alt chosen
//# X is nj x k matrix of xvalues for each of j alt on each of n occasions
//# RW increments are N(0,s^2*t(inc.root)%*%inc.root)
//# prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
//#  inc.root, rootpi are upper triangular
//#  this means that we are using the UL decomp of Sigma^-1 for prior 
//# oldbeta is the current
//     stay=0
//     betac=oldbeta + s*t(inc.root)%*%(matrix(rnorm(ncol(X)),ncol=1))
//     cll=llmnl(betac,y,X)
//     clpost=cll+lndMvn(betac,betabar,rootpi)
//     ldiff=clpost-oldll-lndMvn(oldbeta,betabar,rootpi)
//     alpha=min(1,exp(ldiff))
//     if(alpha < 1) {unif=runif(1)} else {unif=0}
//     if (unif <= alpha)
//             {betadraw=betac; oldll=cll}
//           else
//             {betadraw=oldbeta; stay=1}
//return(list(betadraw=betadraw,stay=stay,oldll=oldll))
//}

mnlMetropOnceOut mnlMetropOnce_con(vec const& y, mat const& X, vec const& oldbeta,
    double oldll, double s, mat const& incroot,
    vec const& betabar, mat const& rootpi, vec const& SignRes = NumericVector::create(2)) {
    // Wayne Taylor 10/01/2014

    // function to execute rw metropolis for the MNL
    // y is n vector with element = 1,...,j indicating which alt chosen
    // X is nj x k matrix of xvalues for each of j alt on each of n occasions
    // RW increments are N(0,s^2*t(inc.root)%*%inc.root)
    // prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
    //  inc.root, rootpi are upper triangular
    //  this means that we are using the UL decomp of Sigma^-1 for prior 
    // oldbeta is the current


    mnlMetropOnceOut out_struct;

    double unif;
    vec betadraw, alphaminv(2);

    int stay = 0;
    vec betac = oldbeta + s * trans(incroot) * as<vec>(rnorm(X.n_cols));
    double cll = llmnl_con(betac, y, X, SignRes);
    double clpost = cll + lndMvn(betac, betabar, rootpi);
    double ldiff = clpost - oldll - lndMvn(oldbeta, betabar, rootpi);
    //alphaminv << 1 << exp(ldiff);
    alphaminv = { 1, exp(ldiff) };
    double alpha = min(alphaminv);

    if (alpha < 1) {
        unif = as_scalar(vec(runif(1)));
    }
    else {
        unif = 0;
    }
    if (unif <= alpha) {
        betadraw = betac;
        oldll = cll;
    }
    else {
        betadraw = oldbeta;
        stay = 1;
    }

    out_struct.betadraw = betadraw;
    out_struct.stay = stay;
    out_struct.oldll = oldll;

    return (out_struct);
}

//MAIN FUNCTION-------------------------------------------------------------------------------------

//[[Rcpp::export]]
List rhierMnlRwMixture_rcpp_loop(List const& lgtdata, mat const& Z,
    vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,
    double nu, mat const& V, double s,
    int R, int keep, int nprint, bool drawdelta,
    mat olddelta, vec const& a, vec oldprob, mat oldbetas, vec ind, vec const& SignRes,
    bool useBART, List const& bart_params,
    bool useHeterCov = false,
    List const& var_params = List::create(),
    List const& phi_params = List::create()) {

    // Wayne Taylor 10/01/2014
    // modified by Thomas Wiemann 2025
    // useHeterCov branch: heteroscedastic Sigma(Z_i) via modified-Cholesky tree
    //   ensembles.  Mutually exclusive with the mixture path.
    //   Phase 1 (phi_params empty): diagonal Sigma(Z_i) = D(Z_i).
    //   Phase 2 (phi_params non-empty): full modified-Cholesky
    //     Sigma(Z_i)^{-1} = L(Z_i)^T D(Z_i)^{-1} L(Z_i) with phi_jk(Z_i) trees.

    int nlgt = lgtdata.size();
    int nvar = V.n_cols;
    int nz = Z.n_cols;

    // Validation for the heter-cov path.
    if (useHeterCov && !useBART)    stop("useHeterCov requires useBART = TRUE");
    if (useHeterCov && !drawdelta)  stop("useHeterCov requires drawdelta = TRUE");

    mat rootpi, betabar, ucholinv, incroot;
    int mkeep;
    mnlMetropOnceOut metropout_struct;
    List lgtdatai, nmix;

    mat delta_Z;
    vec pred_i(nvar);
    vec scaled_pred(nvar);

    // convert List to std::vector of struct
    std::vector<moments> lgtdata_vector;
    moments lgtdatai_struct;
    for (int lgt = 0; lgt < nlgt; lgt++) {
        lgtdatai = lgtdata[lgt];

        lgtdatai_struct.y = as<vec>(lgtdatai["y"]);
        lgtdatai_struct.X = as<mat>(lgtdatai["X"]);
        lgtdatai_struct.hess = as<mat>(lgtdatai["hess"]);
        lgtdata_vector.push_back(lgtdatai_struct);
    }

    // allocate space for draws
    vec oldll = zeros<vec>(nlgt);
    cube betadraw(nlgt, nvar, R / keep);
    //cube DeltaZdraw(nlgt, nvar, R / keep);
    mat probdraw(R / keep, oldprob.size());
    vec loglike(R / keep);
    mat Deltadraw = drawdelta && !useBART ? zeros<mat>(R / keep, nz * nvar) : zeros<mat>(1, 1); //enlarge Deltadraw only if the space is required
    List compdraw(R / keep);
	if (drawdelta) delta_Z.zeros(nlgt, nvar);

    // Initialize BART-specific storage
    std::vector<bart> bart_models;
    std::vector<std::stringstream> treess;
    cube varcount, varprob;
    if (useBART && drawdelta) {
        treess.resize(nvar);
        size_t num_trees = bart_params["num_trees"];

        for (int k = 0; k < nvar; k++) {
            treess[k].precision(10);
            treess[k] << R / keep << " " << num_trees << " " << nz << endl;
        }

        // Initialize variable usage tracking
        varcount.set_size(nz, nvar, R / keep);
        varprob.set_size(nz, nvar, R / keep);
        varcount.zeros();
        varprob.zeros();
    }

    // ===== Heter-cov storage =====
    //   Phase 1: phi_models stays empty, cov_helpers takes the diagonal path.
    //   Phase 2: phi_models is populated by initializePhiBART; same code reads
    //   the full Cholesky via cov::cov_eval / cov::rootpi / sigma_sqrt_apply.
    std::vector<varbart> var_models;
    std::vector<std::vector<heterbart>> phi_models;
    mat eta_sq_buf;                                         // nlgt x nvar
    std::vector<double*> peta_sq_cols(nvar);
    vec mu_post(nvar);                                      // single posterior mean
    mat mu_draw;                                            // R/keep x nvar
    std::vector<std::stringstream> var_treess;
    cube var_varcount, var_varprob;
    int var_burn = 0;
    bool var_sparse = false;

    // Phase-2-only storage (full Cholesky via phi-trees).
    const bool useFullCov = useHeterCov && (phi_params.size() > 0);
    const size_t n_pairs  = useFullCov ? (nvar * (nvar - 1)) / 2 : 0;
    mat Y_phi;                                              // nlgt x n_pairs (per-(j,k) Y buffers)
    vec sigma_phi_buf;                                      // nlgt
    std::vector<std::vector<double*>> py_phi_cols;          // [j][k] -> Y_phi colptr
    int n_phi_treess = 0;
    std::vector<std::vector<std::stringstream>> phi_treess; // [j][k] tree streams
    if (useFullCov) {
        Y_phi.zeros(nlgt, n_pairs);
        sigma_phi_buf.zeros(nlgt);
        py_phi_cols.resize(nvar);
        size_t pair_idx = 0;
        for (int j = 1; j < nvar; j++) {
            py_phi_cols[j].resize(j);
            for (int k = 0; k < j; k++) {
                py_phi_cols[j][k] = Y_phi.colptr(pair_idx++);
            }
        }
        // Per-(j,k) tree-string streams (jagged, lower triangular).
        phi_treess.resize(nvar);
        size_t phi_num_trees = as<size_t>(phi_params["num_trees"]);
        for (int j = 1; j < nvar; j++) {
            phi_treess[j].resize(j);
            for (int k = 0; k < j; k++) {
                phi_treess[j][k].precision(10);
                phi_treess[j][k] << R / keep << " " << phi_num_trees << " " << nz << endl;
            }
            n_phi_treess += j;
        }
    }
    if (useHeterCov) {
        eta_sq_buf.zeros(nlgt, nvar);
        for (int j = 0; j < nvar; j++) peta_sq_cols[j] = eta_sq_buf.colptr(j);
        mu_post.zeros();
        mu_draw.zeros(R / keep, nvar);

        var_treess.resize(nvar);
        size_t var_num_trees = as<size_t>(var_params["num_trees"]);
        for (int j = 0; j < nvar; j++) {
            var_treess[j].precision(10);
            var_treess[j] << R / keep << " " << var_num_trees << " " << nz << endl;
        }
        var_varcount.set_size(nz, nvar, R / keep); var_varcount.zeros();
        var_varprob .set_size(nz, nvar, R / keep); var_varprob .zeros();

        if (var_params.containsElementNamed("sparse")) {
            var_sparse = as<bool>(var_params["sparse"]);
            if (var_sparse && var_params.containsElementNamed("burn"))
                var_burn = as<int>(var_params["burn"]);
        }
    }

    // Create column-major Z for BART
    mat Zt = Z.t();
    double* pZt = Zt.memptr(); // pointer to (transpose of) Z
    mat std_oldbetas = oldbetas; // standardized oldbetas, initialized at initial values
    std::vector<double*> pstd_oldbetas_cols(nvar); // pointers to columns
    for (int i = 0; i < nvar; i++) {
        pstd_oldbetas_cols[i] = std_oldbetas.colptr(i);
    }

    // Random number generator for BART
    arn gen;

    // Add burn in for DART
    int burn = 0;
    bool sparse = false;
    if (useBART) {
        sparse = bart_params["sparse"];
        if (sparse) burn = bart_params["burn"];
    }

    // Start MCMC loop
    if (nprint > 0) startMcmcTimer();

    for (int rep = 0; rep < R; rep++) {

        // ====================================================================
        // Step A: draw mu | {theta_i}, Sigma(.)
        //
        // Heter-cov path: conjugate normal (single component, no mixture).
        // Mixture path: rmixGibbs as before.
        // ====================================================================
        List mgout;
        List oldcomp;
        if (useHeterCov) {
            if (rep == 0) {
                // Initial mu = column means of theta_i.
                mu_post = (sum(oldbetas, 0) / static_cast<double>(nlgt)).t();
                // Naive eta_sq baseline (theta - mu)^2 for var_models init.
                for (int j = 0; j < nvar; j++) {
                    for (int lgt = 0; lgt < nlgt; lgt++) {
                        eta_sq_buf(lgt, j) = std::pow(oldbetas(lgt, j) - mu_post(j), 2.0);
                    }
                }
                var_models = initializeVarBART(pZt, nlgt, nz, peta_sq_cols, var_params, gen);
                if (useFullCov) {
                    // phi_models initialized BEFORE update_stdoldbetas_het so the
                    // Mahalanobis transform sees the full Cholesky from the start.
                    // Y buffers are zero at init -> phi_models predict 0 (their
                    // prior mean), which is the correct starting point.
                    phi_models = initializePhiBART(pZt, nlgt, nz, py_phi_cols,
                                                   phi_params, gen);
                }
                // Standardize betas using the initial (mu, var_models, phi_models) and
                // initialize the mean trees on the standardized residuals.
                update_stdoldbetas_het(oldbetas, pstd_oldbetas_cols,
                                       mu_post, var_models, phi_models);
                bart_models = initializeBART(pZt, nlgt, nz, pstd_oldbetas_cols,
                                             bart_params, gen);
            } else {
                mu_post = drawMuHeterCov(oldbetas, mubar.row(0).t(), Amu,
                                         var_models, phi_models, gen);
            }
        }
        else if (useBART && drawdelta) {
            if (rep == 0) {
                mgout = rmixGibbs(oldbetas, mubar, Amu, nu, V, a, oldprob, ind);
                update_stdoldbetas(oldbetas, pstd_oldbetas_cols, ind, mgout["comps"]);
                bart_models = initializeBART(pZt, nlgt, nz, pstd_oldbetas_cols, bart_params, gen);
            }
            else {
                mgout = rmixGibbs(oldbetas - delta_Z, mubar, Amu, nu, V, a, oldprob, ind);
            }
        }
        else if (drawdelta) {
            olddelta.reshape(nvar, nz);
            mgout = rmixGibbs(oldbetas - delta_Z, mubar, Amu, nu, V, a, oldprob, ind); //mgout = rmixGibbs(oldbetas - Z * trans(olddelta), mubar, Amu, nu, V, a, oldprob, ind);
        }
        else {
            mgout = rmixGibbs(oldbetas, mubar, Amu, nu, V, a, oldprob, ind);
        }

        if (!useHeterCov) {
            oldcomp = mgout["comps"];
            oldprob = as<vec>(mgout["p"]);
            ind = as<vec>(mgout["z"]);
        }

        // ====================================================================
        // heter-cov Gibbs cycle:  C -> D' -> B  (Wiemann 2025 §3.1 ordering:
        // covariance components Sigma updated BEFORE BART, so BART always sees
        // the freshest Sigma when standardizing tilde_theta).  delta_Z is
        // computed at the natural end of Step B and consumed unchanged by
        // Step E (no staleness, no recompute).
        //
        // Throughout we work in the (delta, d, phi) tuple of "fundamental"
        // tree outputs rather than the derived delta_Z = L^{-1} D^{1/2} delta:
        //
        //   eta_i = L(theta_i - mu) - D^{1/2} delta_i
        //   eta_i^{(j)} = ((theta_i - mu) - sum_{k<j} phi_jk (theta_k-mu_k))
        //               - sqrt(d_j(Z_i)) * delta_j(Z_i)
        //
        // This avoids a subtle bug: expanding eta_i = L(theta_i - mu - delta_Z)
        // componentwise as (theta_j - mu_j - delta_Z_j) - sum phi_jk
        // (theta_k - mu_k) silently drops the cross terms sum phi_jk
        // delta_Z_k whenever both delta_Z and phi are non-zero.
        // ====================================================================
        if (useHeterCov) {
            // DART burn-in (mean and variance trees independently configured).
            if (sparse && rep == burn + 1) {
                for (int j = 0; j < nvar; j++) bart_models[j].startdart();
            }
            if (var_sparse && rep == var_burn + 1) {
                for (int j = 0; j < nvar; j++) var_models[j].startdart();
            }

            // ===== Step C: variance trees d_j | mu, theta, current delta, phi.
            // Compute eta_i = L(theta_i - mu) - D^{1/2} delta_i once per unit
            // using the CURRENT bart_models (delta) and phi_models, then draw.
            for (int i = 0; i < nlgt; i++) {
                cov::cov_eval ev{var_models, phi_models, (size_t)i};
                vec eta_i = oldbetas.row(i).t() - mu_post;
                cov::apply_L(eta_i, ev);                            // eta_i := L (theta_i - mu)
                for (int j = 0; j < nvar; j++) {
                    double sqrt_dj = std::sqrt(var_models[j].f(i));
                    double delta_j = bart_models[j].f(i);
                    double e_ij    = eta_i[j] - sqrt_dj * delta_j;
                    eta_sq_buf(i, j) = e_ij * e_ij;
                }
            }
            for (int j = 0; j < nvar; j++) var_models[j].draw(gen);

            // ===== Step D' (Phase 2 only): phi_jk | mu, theta, current
            //       delta, NEW d, OTHER phi.
            // Per discussions/2026-05-12-phibart-vs-heterbart.md, fit each
            // phi_jk via the varying-coefficient -> heterbart transform
            //
            //   r_i = phi_jk(Z_i) * c_i^{(k)} + eta_i^{(j)}     (data model)
            //   Y_i = r_i / c_i^{(k)},  sigma_i = sqrt(d_j) / |c_i^{(k)}|.
            //
            //   r_i = (theta_j - mu_j) - sqrt(d_j) delta_j
            //         - sum_{m<j, m != k} phi_jm (theta_m - mu_m)
            //
            // c_i^{(k)} clamp: |c| >= 1e-8.  ESS guard inside heterbd ensures
            // clamped observations contribute negligibly to splits.
            if (useFullCov) {
                const double eps = 1e-8;
                for (int j = 1; j < nvar; j++) {
                    for (int k = 0; k < j; k++) {
                        double* y_buf = py_phi_cols[j][k];
                        for (int i = 0; i < nlgt; i++) {
                            double c_ik = oldbetas(i, k) - mu_post(k);
                            if (std::abs(c_ik) < eps)
                                c_ik = (c_ik >= 0 ? eps : -eps);

                            double sqrt_dj = std::sqrt(var_models[j].f(i));
                            double delta_j = bart_models[j].f(i);
                            // r_i: L-row (j, .) applied to (theta - mu), minus
                            // sqrt(d_j) delta_j, minus the other (j, m) phi
                            // contributions for m != k.
                            double r_i = oldbetas(i, j) - mu_post(j) - sqrt_dj * delta_j;
                            for (int m = 0; m < j; m++) {
                                if (m == k) continue;
                                r_i -= phi_models[j][m].f(i)
                                     * (oldbetas(i, m) - mu_post(m));
                            }

                            y_buf[i]         = r_i / c_ik;
                            sigma_phi_buf[i] = sqrt_dj / std::abs(c_ik);
                        }
                        phi_models[j][k].draw(sigma_phi_buf.memptr(), gen);
                    }
                }
            }

            // ===== Step B: mean trees delta | mu, NEW d, NEW phi, theta.
            // Standardize tilde_theta_i = D(Z_i)^{-1/2} L(Z_i) (theta_i - mu)
            // using the FRESHLY updated d and phi, then backfit.  Mirrors
            // Wiemann 2025 Step 4 (BART | mu, Sigma, theta) directly.
            update_stdoldbetas_het(oldbetas, pstd_oldbetas_cols,
                                   mu_post, var_models, phi_models);
            for (int j = 0; j < nvar; j++) bart_models[j].draw(1.0, gen);

            // delta_Z(i, .) = Sigma(Z_i)^{1/2} delta(Z_i) = L^{-1} D^{1/2} delta,
            // evaluated with CURRENT mu, d, phi, delta.  Step E and the
            // diagnostic accumulator both consume this directly.
            for (int i = 0; i < nlgt; i++) {
                for (int j = 0; j < nvar; j++) pred_i(j) = bart_models[j].f(i);
                cov::cov_eval ev{var_models, phi_models, (size_t)i};
                vec d_i = cov::sigma_sqrt_apply(pred_i, ev);
                delta_Z.row(i) = d_i.t();
            }
        }
        else if (useBART && drawdelta) {
            // Start DART after burn-in
            if (sparse && rep == burn + 1) {
                for (int i = 0; i < nvar; i++) {
                    bart_models[i].startdart();
                }
            }
            // Normalize oldbetas
            update_stdoldbetas(oldbetas, pstd_oldbetas_cols, ind, oldcomp);
            // Draw bart
            for (int i = 0; i < nvar; i++) {
                bart_models[i].draw(1.0, gen);
            }
            // Compute delta_Z
            for (int i = 0; i < nlgt; i++) {
                // Get predictions
                for (int j = 0; j < nvar; j++) {
                    pred_i(j) = bart_models[j].f(i);
                }
                // Re-scale
                List comp_i = oldcomp[ind[i] - 1];
                mat rootii = trans(as<mat>(comp_i[1])); 
                scaled_pred = solve(trimatl(rootii), pred_i);
                delta_Z.row(i) = trans(scaled_pred);
            }
        }
        else if (drawdelta) {
            // Draw regression
            olddelta = drawDelta(Z, oldbetas, ind, oldcomp, deltabar, Ad);
            // Compute delta_Z
            olddelta.reshape(nvar, nz);  // Reshape for multiplication
            delta_Z = Z * trans(olddelta);
        }

        // ====================================================================
        // Step E: theta_i Metropolis | Delta(.), Sigma(.), data.
        // ====================================================================
        for (int lgt = 0; lgt < nlgt; lgt++) {
            if (useHeterCov) {
                cov::cov_eval ev{var_models, phi_models, (size_t)lgt};
                rootpi  = cov::rootpi(ev);                         // upper-tri
                betabar = mu_post + delta_Z.row(lgt).t();
            } else {
                List oldcomplgt = oldcomp[ind[lgt] - 1];
                rootpi = as<mat>(oldcomplgt[1]);
                if (drawdelta) {
                    betabar = as<vec>(oldcomplgt[0]) + delta_Z.row(lgt).t();
                } else {
                    betabar = as<vec>(oldcomplgt[0]);
                }
            }

            if (rep == 0) oldll[lgt] = llmnl_con(vectorise(oldbetas(lgt, span::all)), lgtdata_vector[lgt].y, lgtdata_vector[lgt].X, SignRes);

            //compute inc.root
            ucholinv = solve(trimatu(chol(lgtdata_vector[lgt].hess + rootpi * trans(rootpi))), eye(nvar, nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
            incroot = chol(ucholinv * trans(ucholinv));

            metropout_struct = mnlMetropOnce_con(lgtdata_vector[lgt].y, lgtdata_vector[lgt].X, vectorise(oldbetas(lgt, span::all)),
                oldll[lgt], s, incroot, betabar, rootpi, SignRes);

            oldbetas(lgt, span::all) = trans(metropout_struct.betadraw);
            oldll[lgt] = metropout_struct.oldll;
        }

        //print time to completion and draw # every nprint'th draw
        if (nprint > 0) if ((rep + 1) % nprint == 0) infoMcmcTimer(rep, R);

        if ((rep + 1) % keep == 0) {
            mkeep = (rep + 1) / keep;
            betadraw.slice(mkeep - 1) = oldbetas;
            loglike[mkeep - 1] = sum(oldll);
            if (useHeterCov) {
                // Store posterior mu and BOTH bart + var trees + var counts.
                mu_draw.row(mkeep - 1) = mu_post.t();
                for (int j = 0; j < nvar; j++) {
                    for (size_t k = 0; k < bart_models[j].getm(); k++) treess[j]    << bart_models[j].gettree(k);
                    for (size_t k = 0; k < var_models[j].getm();  k++) var_treess[j]<< var_models[j].gettree(k);
                    varcount    .slice(mkeep - 1).col(j) = conv_to<vec>::from(bart_models[j].getnv());
                    varprob     .slice(mkeep - 1).col(j) = conv_to<vec>::from(bart_models[j].getpv());
                    var_varcount.slice(mkeep - 1).col(j) = conv_to<vec>::from(var_models[j].getnv());
                    var_varprob .slice(mkeep - 1).col(j) = conv_to<vec>::from(var_models[j].getpv());
                }
                // Phase-2 phi-tree serialization (one tree-string per (j, k)).
                if (useFullCov) {
                    for (int j = 1; j < nvar; j++) {
                        for (int k = 0; k < j; k++) {
                            for (size_t t = 0; t < phi_models[j][k].getm(); t++) {
                                phi_treess[j][k] << phi_models[j][k].gettree(t);
                            }
                        }
                    }
                }
            }
            else {
                probdraw(mkeep - 1, span::all) = trans(oldprob);
                compdraw[mkeep - 1] = oldcomp;
                if (useBART && drawdelta) {
                    // Store BART and variable usage information
                    for (int i = 0; i < nvar; i++) {
                        for (size_t j = 0; j < bart_models[i].getm(); j++) {
                            treess[i] << bart_models[i].gettree(j);
                        }
                        varcount.slice(mkeep - 1).col(i) = conv_to<vec>::from(bart_models[i].getnv());
                        varprob.slice(mkeep - 1).col(i) = conv_to<vec>::from(bart_models[i].getpv());
                    }
                }
                else if (drawdelta) {
                    Deltadraw(mkeep - 1, span::all) = trans(vectorise(olddelta));
                }
            }
        }
    }

    if (nprint > 0) endMcmcTimer();

    nmix = List::create(Named("probdraw") = probdraw,
        Named("zdraw") = R_NilValue, //sets the value to NULL in R
        Named("compdraw") = compdraw);

    //ADDED FOR CONSTRAINTS
    //If there are sign constraints, return f(betadraws) as "betadraws"
    //conStatus will be set to true if SignRes has any non-zero elements
    bool conStatus = any(SignRes);

    if (conStatus) {
        int SignResSize = SignRes.size();

        //loop through each sign constraint
        for (int i = 0;i < SignResSize; i++) {

            //if there is a constraint loop through each slice of betadraw
            if (SignRes[i] != 0) {
                for (int s = 0;s < R / keep; s++) {
                    betadraw(span(), span(i), span(s)) = SignRes[i] * exp(betadraw(span(), span(i), span(s)));
                }
            }

        }//end loop through SignRes
    }

    // Heter-cov return path (Phase 1 and Phase 2).
    //
    // The R-side predict.rhierMnlRwMixture method dispatches on class
    // "rhierMnlRwMixtureHeterCov" and consumes (bart_models, var_models,
    // phi_models, mu_draw) to evaluate Sigma(Z*)^{1/2} delta(Z*) at any
    // user-supplied Z*.  No diagnostic running-sum accumulators are returned
    // from C++; they would duplicate work that predict.* now does cleanly.
    if (useHeterCov) {
        SEXP phi_models_pkg = useFullCov
            ? (SEXP) pack_bart_jagged_(nvar, phi_models, phi_treess)
            : R_NilValue;

        return(List::create(
            Named("bart_models")  = pack_bart_list_(nvar, bart_models, treess),
            Named("var_models")   = pack_bart_list_(nvar, var_models,  var_treess),
            Named("phi_models")   = phi_models_pkg,
            Named("mu_draw")      = mu_draw,
            Named("varcount")     = varcount,
            Named("varprob")      = varprob,
            Named("var_varcount") = var_varcount,
            Named("var_varprob")  = var_varprob,
            Named("betadraw")     = betadraw,
            Named("loglike")      = loglike,
            Named("SignRes")      = SignRes));
    }

    if (useBART && drawdelta) {
        return(List::create(
            Named("bart_models") = pack_bart_list_(nvar, bart_models, treess),
            Named("varcount") = varcount,
            Named("varprob") = varprob,
			//Named("DeltaZdraw") = DeltaZdraw,
            Named("betadraw") = betadraw,
            Named("nmix") = nmix,
            Named("loglike") = loglike,
            Named("SignRes") = SignRes));
    }
    else if (drawdelta) {
        return(List::create(
			Named("Deltadraw") = Deltadraw,
			//Named("DeltaZdraw") = DeltaZdraw,
            Named("betadraw") = betadraw,
            Named("nmix") = nmix,
            Named("loglike") = loglike,
            Named("SignRes") = SignRes));
    }
    else {
        return(List::create(
            Named("betadraw") = betadraw,
            Named("nmix") = nmix,
            Named("loglike") = loglike,
            Named("SignRes") = SignRes));
    }

}