#include "bayesm.HART.h"
#include "mixbart_block.h"    // MixBartState + mixbart_* helpers (BART path)
#include "hetercov_block.h"   // HeterCovState + hetercov_* helpers

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
    bool useHeterCov,
    List const& var_params,
    List const& phi_params) {

    // Wayne Taylor 10/01/2014
    // modified by Thomas Wiemann 2025
    //
    // Three execution modes (mutually exclusive):
    //   * Heter-cov: Sigma(Z_i) via modified-Cholesky tree ensembles.
    //                Implementation lives in src/hetercov_block.cpp.
    //   * Mix-BART : sum-of-trees Delta(Z) on top of mixture-of-normals u_i.
    //                Implementation lives in src/mixbart_block.cpp; shared
    //                with rhierLinearMixture and rhierNegbinRw.
    //   * Linear   : original bayesm linear-Delta + rmixGibbs path.

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
    int nacceptbeta = 0;
    cube betadraw(nlgt, nvar, R / keep);
    mat probdraw(R / keep, oldprob.size());
    vec loglike(R / keep);
    mat Deltadraw = drawdelta && !useBART ? zeros<mat>(R / keep, nz * nvar) : zeros<mat>(1, 1); //enlarge Deltadraw only if the space is required
    List compdraw(R / keep);
	if (drawdelta) delta_Z.zeros(nlgt, nvar);

    // ===== Mix-BART state (only when useBART && drawdelta && !useHeterCov) =
    // Owns bart_models, treess, varcount, varprob, oldcomp, oldprob, ind,
    // delta_Z, plus DART burn-in scheduling.  See inst/include/mixbart_block.h.
    MixBartState mbs;
    if (useBART && drawdelta && !useHeterCov)
        mixbart_init(mbs, nlgt, nz, nvar, R, keep, bart_params, oldprob, ind);

    // ===== Heter-cov state ===============================================
    HeterCovState hcs;
    if (useHeterCov)
        hetercov_init(hcs, nlgt, nz, nvar, R, keep,
                      bart_params, var_params, phi_params);

    // BART operates on standardized betas; std_oldbetas is the caller-owned
    // working buffer that update_stdoldbetas* writes into and that BART reads
    // residuals from.  The mixbart_block / hetercov_block helpers do NOT own
    // this storage.
    mat Zt = Z.t();
    double* pZt = Zt.memptr();
    mat std_oldbetas = oldbetas;
    std::vector<double*> pstd_oldbetas_cols(nvar);
    for (int i = 0; i < nvar; i++) pstd_oldbetas_cols[i] = std_oldbetas.colptr(i);

    // Random number generator for BART
    arn gen;

    // Start MCMC loop
    if (nprint > 0) startMcmcTimer();

    for (int rep = 0; rep < R; rep++) {

        // ====================================================================
        // Steps A-D: heter-cov branch is fully encapsulated; mixture-BART
        // backfitting is fully delegated to mixbart_block.cpp.
        // ====================================================================
        List mgout;
        List oldcomp;
        if (useHeterCov) {
            // Steps A, B, C, D': Draw mu, variance trees (d_j), phi trees, and mean trees (Delta)
            hetercov_draw_iter(hcs, oldbetas, mubar.row(0).t(), Amu,
                               pZt, pstd_oldbetas_cols,
                               bart_params, var_params, phi_params,
                               rep, gen);
            delta_Z = hcs.delta_Z;
        }
        else if (useBART && drawdelta) {
            // Steps 2 & 3: Standardize unit-level betas and draw mean trees (Delta)
            mixbart_draw_iter(mbs, oldbetas, mubar, Amu, nu, V, a,
                              pZt, pstd_oldbetas_cols, bart_params, rep, gen);
            oldcomp = mbs.oldcomp;
            oldprob = mbs.oldprob;
            ind     = mbs.ind;
            delta_Z = mbs.delta_Z;
        }
        else if (drawdelta) {
            // Step 4: Draw mixture component assignments and update global mixture parameters
            olddelta.reshape(nvar, nz);
            mgout = rmixGibbs(oldbetas - delta_Z, mubar, Amu, nu, V, a, oldprob, ind);
            oldcomp = mgout["comps"];
            oldprob = as<vec>(mgout["p"]);
            ind = as<vec>(mgout["z"]);

            olddelta = drawDelta(Z, oldbetas, ind, oldcomp, deltabar, Ad);
            olddelta.reshape(nvar, nz);
            delta_Z = Z * trans(olddelta);
        }
        else {
            // Step 4: Draw mixture component assignments and update global mixture parameters
            mgout = rmixGibbs(oldbetas, mubar, Amu, nu, V, a, oldprob, ind);
            oldcomp = mgout["comps"];
            oldprob = as<vec>(mgout["p"]);
            ind = as<vec>(mgout["z"]);
        }

        // ====================================================================
        // Step E: theta_i Metropolis | Delta(.), Sigma(.), data.
        // ====================================================================
        for (int lgt = 0; lgt < nlgt; lgt++) {
            if (useHeterCov) {
                cov::cov_eval ev{hcs.var_models, hcs.phi_models, (size_t)lgt};
                rootpi  = cov::rootpi(ev);
                betabar = hcs.mu_post + delta_Z.row(lgt).t();
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
            if (metropout_struct.stay == 0) nacceptbeta = nacceptbeta + 1;

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
                hetercov_store(hcs, mkeep);
            }
            else {
                probdraw(mkeep - 1, span::all) = trans(oldprob);
                compdraw[mkeep - 1] = oldcomp;
                if (useBART && drawdelta) {
                    mixbart_store(mbs, mkeep);
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
        for (int i = 0;i < SignResSize; i++) {
            if (SignRes[i] != 0) {
                for (int s = 0;s < R / keep; s++) {
                    betadraw(span(), span(i), span(s)) = SignRes[i] * exp(betadraw(span(), span(i), span(s)));
                }
            }
        }
    }

    double acceptrbeta = nacceptbeta / (R * nlgt * 1.0) * 100;

    // Heter-cov return path: hetercov_pack supplies the heter-cov outputs;
    // caller merges in betadraw / loglike / SignRes.  The R-side
    // predict.rhierMnlRwMixture method dispatches on class
    // "bayesm.HART.HeterCov" and consumes (bart_models, var_models,
    // phi_models, mu_draw) to evaluate Sigma(Z*)^{1/2} delta(Z*) at any
    // user-supplied Z*.
    if (useHeterCov) {
        List out = hetercov_pack(hcs);
        out["betadraw"] = betadraw;
        out["loglike"]  = loglike;
        out["SignRes"]  = SignRes;
        out["acceptrbeta"] = acceptrbeta;
        return out;
    }

    if (useBART && drawdelta) {
        List out = mixbart_pack(mbs);
        out["betadraw"] = betadraw;
        out["nmix"]     = nmix;
        out["loglike"]  = loglike;
        out["SignRes"]  = SignRes;
        out["acceptrbeta"] = acceptrbeta;
        return out;
    }
    else if (drawdelta) {
        return(List::create(
			Named("Deltadraw") = Deltadraw,
            Named("betadraw") = betadraw,
            Named("nmix") = nmix,
            Named("loglike") = loglike,
            Named("SignRes") = SignRes,
            Named("acceptrbeta") = acceptrbeta));
    }
    else {
        return(List::create(
            Named("betadraw") = betadraw,
            Named("nmix") = nmix,
            Named("loglike") = loglike,
            Named("SignRes") = SignRes,
            Named("acceptrbeta") = acceptrbeta));
    }

}
