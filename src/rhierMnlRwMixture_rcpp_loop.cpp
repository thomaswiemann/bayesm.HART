#include "bayesm.HART.h"
#include "hart_backend.h"
#include "mnl_theta_updater.h"
#include "unified_sampler_driver.h"
#include <memory>

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
    std::string const& hart_model, std::string const& hart_cov,
    List const& bart_params,
    List const& var_params, List const& phi_params,
    Rcpp::Nullable<Rcpp::List> state = R_NilValue) {

    // Wayne Taylor 10/01/2014
    // modified by Thomas Wiemann 2025
    //
    // Covariance backend by hart_model (make_backend in hart_backend.h).

    int nlgt = lgtdata.size();
    int nvar = V.n_cols;
    if (hart_model != "none" && !drawdelta)
        stop("hart_model = '%s' requires Z.", hart_model.c_str());

    List lgtdatai;

    std::vector<moments> lgtdata_vector;
    moments lgtdatai_struct;
    for (int lgt = 0; lgt < nlgt; lgt++) {
        lgtdatai = lgtdata[lgt];

        lgtdatai_struct.y = as<vec>(lgtdatai["y"]);
        lgtdatai_struct.X = as<mat>(lgtdatai["X"]);
        lgtdatai_struct.hess = as<mat>(lgtdatai["hess"]);
        lgtdata_vector.push_back(lgtdatai_struct);
    }

    cube betadraw(nlgt, nvar, R / keep);
    vec loglike(R / keep);

    int offset = 0;
    List saved;
    if (state.isNotNull()) {
        saved = List(state);
        oldbetas = as<mat>(saved["betas"]);
        offset = as<int>(saved["completed_sweeps"]);
    }
    if (oldbetas.n_rows != (uword)nlgt || oldbetas.n_cols != (uword)nvar
        || offset < 0)
        stop("Invalid sampler state dimensions or completed sweeps.");

    std::unique_ptr<CovarianceBackend> backend = make_backend(
        hart_model, hart_cov, oldbetas, Z, mubar, Amu, nu, V, a,
        deltabar, Ad, drawdelta, olddelta, oldprob, ind, R, keep,
        bart_params, var_params, phi_params);
    if (state.isNotNull()) backend->restore(as<List>(saved["backend"]));

    MnlThetaUpdaterAdapter updater(oldbetas, lgtdata_vector, SignRes, s, R);

    run_unified_sampler(R, keep, nprint, nlgt,
                        *backend, updater, oldbetas, betadraw, loglike, offset);

    List terminal = List::create(Named("completed_sweeps") = offset + R,
                                 Named("betas") = oldbetas,
                                 Named("backend") = backend->checkpoint());

    // Public draws: SignRes * exp(betastar) on restricted coordinates.
    bool conStatus = any(SignRes);
    if (conStatus) {
        int SignResSize = SignRes.size();
        for (int i = 0;i < SignResSize; i++) {
            if (SignRes[i] != 0) {
                for (int ss = 0;ss < R / keep; ss++) {
                    betadraw(span(), span(i), span(ss)) = SignRes[i] * exp(betadraw(span(), span(i), span(ss)));
                }
            }
        }
    }

    AcceptanceMetrics am = updater.acceptance_metrics();
    BackendPackPayload bp = backend->pack();

    List out = bp.fields;
    out["state"] = terminal;
    out["betadraw"]    = betadraw;
    out["loglike"]     = loglike;
    out["SignRes"]     = SignRes;
    out["acceptrbeta"] = am.acceptrbeta;
    return out;
}
