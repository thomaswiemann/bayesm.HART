#include "bayesm.HART.h"
#include "mixbart_block.h"
#include "hetercov_block.h"

//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
double llnegbinpooled(std::vector<moments> regdata_vector, mat oldbetas, double alpha){

// Wayne Taylor 12/01/2014

// "Unlists" the regdata and calculates the negative binomial loglikelihood using individual-level betas

  int nreg = regdata_vector.size();
  double ll = 0.0;

  for(int reg = 0; reg<nreg; reg++){
  vec lambda = exp(regdata_vector[reg].X*trans(oldbetas(reg,span::all)));
  ll = ll + llnegbin(regdata_vector[reg].y,lambda,alpha,TRUE);
  }

  return(ll);
}

// [[Rcpp::export]]
List rhierNegbinRw_rcpp_loop(List const& regdata, List const& hessdata, mat const& Z,
                             mat oldbetas, vec const& deltabar, mat const& Ad,
                             mat const& mubar, mat const& Amu,
                             double nu, mat const& V, double a, double b,
                             int R, int keep, double sbeta, double alphacroot, int nprint,
                             bool drawdelta, mat olddelta, vec const& a_mix,
                             vec oldprob, vec ind,
                             double alpha, bool fixalpha,
                             bool useBART = false, List const& bart_params = List::create(),
                             bool useHeterCov = false,
                             List const& var_params = List::create(),
                             List const& phi_params = List::create()){

// Wayne Taylor 12/01/2014
// Modified by Thomas Wiemann 2025 -- HART/BART support, mixture-of-normals,
//   plus heter-cov branch (modified-Cholesky tree ensembles for Sigma(Z_i)).
//
// Three execution modes (mutually exclusive):
//   * Heter-cov: hetercov_block ensembles for d_j(Z_i) [+ phi_jk(Z_i)].
//   * Mix-BART : sum-of-trees Delta(Z) on top of mixture-of-normals u_i.
//   * Linear   : original bayesm linear-Delta + rmixGibbs path.

  // Heter-cov path validation.
  if (useHeterCov && !useBART)   stop("useHeterCov requires useBART = TRUE");
  if (useHeterCov && !drawdelta) stop("useHeterCov requires drawdelta = TRUE");

  double ldiff, acc, unif, logalphac, oldlpostalpha, clpostalpha;
  int mkeep, rep;
  int nreg = regdata.size();
  int nz = Z.n_cols;
  int nvar = V.n_cols;
  int nacceptbeta = 0;
  int nacceptalpha = 0;

  // convert regdata and hessdata Lists to std::vector of struct
  std::vector<moments> regdata_vector;
  moments regdatai_struct;
  List regdatai,hessi;

  // store vector with struct
  for (int reg = 0; reg<nreg; reg++){
    regdatai = regdata[reg];
    hessi = hessdata[reg];

    regdatai_struct.y = as<vec>(regdatai["y"]);
    regdatai_struct.X = as<mat>(regdatai["X"]);
    regdatai_struct.hess = as<mat>(hessi["hess"]);
    regdata_vector.push_back(regdatai_struct);
  }

  // Allocate common storage
  vec oldlpostbeta = zeros<vec>(nreg);
  vec clpostbeta = zeros<vec>(nreg);
  cube betadraw = zeros<cube>(nreg, nvar, R/keep);
  vec alphadraw = zeros<vec>(R/keep);
  vec loglike = zeros<vec>(R/keep);

  // =========================================================================
  // HETER-COV PATH  (Sigma(Z_i) via modified-Cholesky tree ensembles)
  // =========================================================================
  if (useHeterCov) {

    HeterCovState hcs;
    hetercov_init(hcs, nreg, nz, nvar, R, keep,
                  bart_params, var_params, phi_params);

    // BART operates on standardized betas; std_oldbetas is the caller-owned
    // working buffer that update_stdoldbetas_het writes into and that BART
    // reads residuals from.  hetercov_block does NOT own this storage.
    mat std_oldbetas = oldbetas;
    std::vector<double*> pstd_oldbetas_cols(nvar);
    for (int i = 0; i < nvar; i++) pstd_oldbetas_cols[i] = std_oldbetas.colptr(i);
    mat Zt = Z.t();
    double* pZt = Zt.memptr();
    arn gen;

    if (nprint>0) startMcmcTimer();

    for (rep = 0; rep < R; rep++) {

      // Steps A, B, C, D': Draw mu, variance trees (d_j), phi trees, and mean trees (Delta)
      // fully encapsulated in hetercov_draw_iter; on return
      // hcs.{mu_post, var_models, phi_models, bart_models, delta_Z} reflect
      // the freshly drawn Sigma(Z_i)^{1/2} delta_i.
      hetercov_draw_iter(hcs, oldbetas, mubar.row(0).t(), Amu,
                         pZt, pstd_oldbetas_cols,
                         bart_params, var_params, phi_params,
                         rep, gen);

      // Step E: per-unit Negbin RW Metropolis with prior
      //   beta_i | . ~ N(mu + Sigma(Z_i)^{1/2} delta_i,  Sigma(Z_i)),
      // i.e. rootpi = cov::rootpi(ev) (upper-tri square root of Sigma^{-1}).
      for (int reg = 0; reg < nreg; reg++) {
        cov::cov_eval ev{hcs.var_models, hcs.phi_models, (size_t)reg};
        mat rootpi  = cov::rootpi(ev);
        vec betabari = hcs.mu_post + hcs.delta_Z.row(reg).t();

        mat Vbetainv_reg = rootpi * trans(rootpi);
        mat betacvar = sbeta * solve(regdata_vector[reg].hess + Vbetainv_reg, eye(nvar, nvar));
        mat betaroot = trans(chol(betacvar));
        vec betac = vectorise(oldbetas(reg, span::all)) + betaroot * vec(rnorm(nvar));

        oldlpostbeta[reg] = lpostbeta(alpha, trans(oldbetas(reg, span::all)),
                                       regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
        clpostbeta[reg]   = lpostbeta(alpha, betac,
                                       regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
        ldiff = clpostbeta[reg] - oldlpostbeta[reg];
        acc = exp(ldiff);
        if (acc > 1) acc = 1;
        if (acc < 1) { unif = runif(1)[0]; } else { unif = 0; }
        if (unif <= acc) {
          oldbetas(reg, span::all) = trans(betac);
          nacceptbeta = nacceptbeta + 1;
        }
      }

      // Alpha MH (identical across paths).
      // Symmetric proposal on log scale: log(alpha_c) ~ N(log(alpha), alphacroot^2).
      // Working on the log scale, the target density is
      //   pi_tilde(log alpha) = pi(alpha) * |d alpha / d log alpha|
      //                       = (alpha^{a-1} e^{-b alpha}) * alpha
      //                       = alpha^{a} e^{-b alpha}
      // so the log-target uses  a*log(alpha) - b*alpha + ll(alpha)  (NOT (a-1)*log(alpha)).
      // The Jacobian factor `log(alpha)` is alpha-dependent and does NOT cancel between
      // old and candidate.  See discussions/2026-05-13-rhier-audit-results.md item M4.4.
      if (!fixalpha) {
        logalphac = log(alpha) + alphacroot * rnorm(1)[0];
        oldlpostalpha = llnegbinpooled(regdata_vector, oldbetas, alpha) + a * log(alpha)  - b * alpha;
        clpostalpha   = llnegbinpooled(regdata_vector, oldbetas, exp(logalphac)) + a * logalphac - b * exp(logalphac);
        ldiff = clpostalpha - oldlpostalpha;
        acc = exp(ldiff);
        if (acc > 1) acc = 1;
        if (acc < 1) { unif = runif(1)[0]; } else { unif = 0; }
        if (unif <= acc) {
          alpha = exp(logalphac);
          nacceptalpha = nacceptalpha + 1;
        }
      }

      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);

      if ((rep+1)%keep==0) {
        mkeep = (rep+1)/keep;
        betadraw.slice(mkeep-1) = oldbetas;
        alphadraw[mkeep-1]      = alpha;
        loglike[mkeep-1]        = llnegbinpooled(regdata_vector, oldbetas, alpha);
        hetercov_store(hcs, mkeep);
      }
    }

    if (nprint>0) endMcmcTimer();

    // Heter-cov return: hetercov_pack supplies tree ensembles + cubes + mu_draw.
    List out = hetercov_pack(hcs);
    out["betadraw"]     = betadraw;
    out["alphadraw"]    = alphadraw;
    out["loglike"]      = loglike;
    out["acceptrbeta"]  = nacceptbeta / (R * nreg * 1.0) * 100;
    out["acceptralpha"] = nacceptalpha / (R * 1.0) * 100;
    return out;
  }

  // =========================================================================
  // BART PATH (sum-of-trees Delta(Z) over mixture-of-normals u_i)
  // =========================================================================
  if (useBART) {

    // BART-specific storage
    mat probdraw(R/keep, oldprob.size());
    List compdraw(R/keep);

    MixBartState mbs;
    mixbart_init(mbs, nreg, nz, nvar, R, keep, bart_params, oldprob, ind);

    // BART operates on standardized betas; std_oldbetas is the caller-owned working
    // buffer that update_stdoldbetas writes into and that BART reads residuals
    // from.  pstd_oldbetas_cols[j] points to column j; the mixbart helpers do
    // NOT own this storage.
    mat std_oldbetas = oldbetas;
    std::vector<double*> pstd_oldbetas_cols(nvar);
    for (int i = 0; i < nvar; i++) pstd_oldbetas_cols[i] = std_oldbetas.colptr(i);
    mat Zt = Z.t();
    double* pZt = Zt.memptr();
    arn gen;

    if (nprint>0) startMcmcTimer();

    // MCMC loop (BART path)
    for (rep = 0; rep < R; rep++) {

      // Steps 2 & 3: Standardize unit-level betas and draw mean trees (Delta)
      mixbart_draw_iter(mbs, oldbetas, mubar, Amu, nu, V, a_mix, pZt, pstd_oldbetas_cols, bart_params, rep, gen);
      List oldcomp = mbs.oldcomp;
      oldprob = mbs.oldprob;
      ind = mbs.ind;
      mat delta_Z = mbs.delta_Z;

      // Step E: Draw unit-level parameters beta_i | Sigma_k^{-1}, etc.
      for (int reg = 0; reg < nreg; reg++) {
        List oldcompreg = oldcomp[ind[reg] - 1];
        mat rootpi = as<mat>(oldcompreg[1]);
        vec betabari = as<vec>(oldcompreg[0]) + trans(delta_Z.row(reg));

        mat Vbetainv_reg = rootpi * trans(rootpi);
        mat betacvar = sbeta * solve(regdata_vector[reg].hess + Vbetainv_reg, eye(nvar, nvar));
        mat betaroot = trans(chol(betacvar));
        vec betac = vectorise(oldbetas(reg, span::all)) + betaroot * vec(rnorm(nvar));

        oldlpostbeta[reg] = lpostbeta(alpha, trans(oldbetas(reg, span::all)),
                                       regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
        clpostbeta[reg] = lpostbeta(alpha, betac,
                                     regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
        ldiff = clpostbeta[reg] - oldlpostbeta[reg];
        acc = exp(ldiff);
        if (acc > 1) acc = 1;
        if (acc < 1) { unif = runif(1)[0]; } else { unif = 0; }
        if (unif <= acc) {
          oldbetas(reg, span::all) = trans(betac);
          nacceptbeta = nacceptbeta + 1;
        }
      }

      // Draw alpha (same in both paths).  See M4.4 in audit results: log-scale Jacobian
      // included via `a*log(alpha)` (not `(a-1)*log(alpha)`).
      if (!fixalpha) {
        logalphac = log(alpha) + alphacroot * rnorm(1)[0];
        oldlpostalpha = llnegbinpooled(regdata_vector, oldbetas, alpha) + a * log(alpha) - b * alpha;
        clpostalpha = llnegbinpooled(regdata_vector, oldbetas, exp(logalphac)) + a * logalphac - b * exp(logalphac);
        ldiff = clpostalpha - oldlpostalpha;
        acc = exp(ldiff);
        if (acc > 1) acc = 1;
        if (acc < 1) { unif = runif(1)[0]; } else { unif = 0; }
        if (unif <= acc) {
          alpha = exp(logalphac);
          nacceptalpha = nacceptalpha + 1;
        }
      }

      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);

      if ((rep+1)%keep==0) {
        mkeep = (rep+1)/keep;
        betadraw.slice(mkeep-1) = oldbetas;
        alphadraw[mkeep-1] = alpha;
        probdraw(mkeep-1, span::all) = trans(oldprob);
        compdraw[mkeep-1] = oldcomp;
        loglike[mkeep-1] = llnegbinpooled(regdata_vector, oldbetas, alpha);

        mixbart_store(mbs, mkeep);
      }
    }

    if (nprint>0) endMcmcTimer();

    // Package BART output
    List nmix = List::create(Named("probdraw") = probdraw,
        Named("zdraw") = R_NilValue,
        Named("compdraw") = compdraw);

    List out = mixbart_pack(mbs);
    out["alphadraw"] = alphadraw;
    out["betadraw"] = betadraw;
    out["nmix"] = nmix;
    out["loglike"] = loglike;
    out["acceptrbeta"] = nacceptbeta / (R * nreg * 1.0) * 100;
    out["acceptralpha"] = nacceptalpha / (R * 1.0) * 100;
    return out;

  } else {

    // =========================================================================
    // NON-BART PATH -- mirrors MNL's rmixGibbs + drawDelta pattern
    // =========================================================================
    mat delta_Z;
    mat probdraw(R/keep, oldprob.size());
    mat Deltadraw = drawdelta ? zeros<mat>(R/keep, nz*nvar) : zeros<mat>(1,1);
    List compdraw(R/keep);
    if (drawdelta) delta_Z.zeros(nreg, nvar);

    if (nprint>0) startMcmcTimer();

    //  start main iteration loop
    for (rep = 0; rep < R; rep++){

      // Step 4: Draw mixture component assignments and update global mixture parameters
      // ind,p need initialization comps is drawn first in sub-Gibbs
      List mgout;
      if (drawdelta) {
          olddelta.reshape(nvar, nz);
          mgout = rmixGibbs(oldbetas - delta_Z, mubar, Amu, nu, V, a_mix, oldprob, ind);
      } else {
          mgout = rmixGibbs(oldbetas, mubar, Amu, nu, V, a_mix, oldprob, ind);
      }

      List oldcomp = mgout["comps"];
      oldprob = as<vec>(mgout["p"]);
      ind = as<vec>(mgout["z"]);

      //now draw delta | {beta_i}, ind, comps
      if (drawdelta) {
          olddelta = drawDelta(Z, oldbetas, ind, oldcomp, deltabar, Ad);
          olddelta.reshape(nvar, nz);
          delta_Z = Z * trans(olddelta);
      }

      // Step 6: Loop over all regression equations drawing beta_i | ind[i],comp[ind[i]]
      for(int reg = 0; reg<nreg; reg++){
          List oldcompreg = oldcomp[ind[reg]-1];
          mat rootpi = as<mat>(oldcompreg[1]);

          //note: beta_i = Delta(z_i) + u_i
          vec betabari;
          if (drawdelta) {
              betabari = as<vec>(oldcompreg[0]) + delta_Z.row(reg).t();
          } else {
              betabari = as<vec>(oldcompreg[0]);
          }

          //compute inc.root using component-specific prior (mirrors MNL)
          mat Abeta = rootpi * trans(rootpi);
          mat betacvar = sbeta*solve(regdata_vector[reg].hess+Abeta,eye(nvar,nvar));
          mat betaroot = trans(chol(betacvar));
          vec betac = vectorise(oldbetas(reg,span::all)) + betaroot*vec(rnorm(nvar));

          oldlpostbeta[reg] = lpostbeta(alpha, trans(oldbetas(reg,span::all)), regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
          clpostbeta[reg] = lpostbeta(alpha, betac, regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
          ldiff = clpostbeta[reg] - oldlpostbeta[reg];
          acc = exp(ldiff);
          if (acc > 1) acc = 1;
          if(acc < 1) {unif=runif(1)[0];} else {unif=0;}
          if (unif <= acc){
            oldbetas(reg,span::all) = trans(betac);
            nacceptbeta = nacceptbeta + 1;
          }
      }

      // Draw alpha.  See M4.4 in audit results: log-scale Jacobian included via
      // `a*log(alpha)` (not `(a-1)*log(alpha)`).
      if (!fixalpha){
        logalphac = log(alpha) + alphacroot*rnorm(1)[0];
        oldlpostalpha = llnegbinpooled(regdata_vector,oldbetas,alpha) + a*log(alpha) - b*alpha;
        clpostalpha = llnegbinpooled(regdata_vector,oldbetas,exp(logalphac)) + a*logalphac - b*exp(logalphac);
        ldiff = clpostalpha - oldlpostalpha;
        acc = exp(ldiff);
        if (acc > 1) acc = 1;
        if(acc < 1) {unif=runif(1)[0];} else {unif=0;}
        if (unif <= acc){
          alpha = exp(logalphac);
          nacceptalpha = nacceptalpha + 1;
        }
      }

      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);

      if((rep+1)%keep==0){
        mkeep = (rep+1)/keep;
        betadraw.slice(mkeep-1) = oldbetas;
        alphadraw[mkeep-1] = alpha;
        probdraw(mkeep-1, span::all) = trans(oldprob);
        if(drawdelta) Deltadraw(mkeep-1,span::all) = trans(vectorise(olddelta));
        compdraw[mkeep-1] = oldcomp;
        loglike[mkeep-1] = llnegbinpooled(regdata_vector,oldbetas,alpha);
        }
    }

    if (nprint>0) endMcmcTimer();

    List nmix = List::create(Named("probdraw") = probdraw,
              Named("zdraw") = R_NilValue,
              Named("compdraw") = compdraw);

    if(drawdelta){
      return List::create(
        Named("loglike") = loglike,
        Named("betadraw") = betadraw,
        Named("alphadraw") = alphadraw,
        Named("Deltadraw") = Deltadraw,
        Named("nmix") = nmix,
        Named("acceptrbeta") = nacceptbeta/(R*nreg*1.0)*100,
        Named("acceptralpha") = nacceptalpha/(R*1.0)*100);
    } else {
      return List::create(
        Named("loglike") = loglike,
        Named("betadraw") = betadraw,
        Named("alphadraw") = alphadraw,
        Named("nmix") = nmix,
        Named("acceptrbeta") = nacceptbeta/(R*nreg*1.0)*100,
        Named("acceptralpha") = nacceptalpha/(R*1.0)*100);
    }
  }
}
