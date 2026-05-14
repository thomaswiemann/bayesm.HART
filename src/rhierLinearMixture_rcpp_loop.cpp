#include "bayesm.HART.h"
#include "mixbart_block.h"
#include "hetercov_block.h"

//[[Rcpp::export]]
List rhierLinearMixture_rcpp_loop(List const& regdata, mat const& Z,
                                  vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,
                                  double nu, mat const& V, double nu_e, vec const& ssq,
                                  int R, int keep, int nprint, bool drawdelta,
                                  mat olddelta, vec const& a, vec oldprob, vec ind, vec tau,
                                  bool useBART, List const& bart_params,
                                  bool useHeterCov,
                                  List const& var_params,
                                  List const& phi_params,
                                  mat const& Beta_init){

// Wayne Taylor 10/02/2014
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

  int nreg = regdata.size();
  int nvar = V.n_cols;
  int nz = Z.n_cols;

  mat rootpi, betabar, Abeta, Abetabar;
  int mkeep;
  unireg runiregout_struct;
  List regdatai, nmix;

  // convert List to std::vector of type "moments"
  std::vector<moments> regdata_vector;
  moments regdatai_struct;

  // store vector with struct
  for (int reg = 0; reg<nreg; reg++){
    regdatai = regdata[reg];

    regdatai_struct.y = as<vec>(regdatai["y"]);
    regdatai_struct.X = as<mat>(regdatai["X"]);
    regdatai_struct.XpX = as<mat>(regdatai["XpX"]);
    regdatai_struct.Xpy = as<vec>(regdatai["Xpy"]);
    regdata_vector.push_back(regdatai_struct);
  }

  // allocate space for draws
  mat oldbetas = zeros<mat>(nreg,nvar);
  mat taudraw(R/keep, nreg);
  cube betadraw(nreg, nvar, R/keep);
  mat probdraw(R/keep, oldprob.size());
  mat Deltadraw(1,1); if(drawdelta) Deltadraw.zeros(R/keep, nz*nvar);
  List compdraw(R/keep);
  vec loglike(R/keep);

  // =========================================================================
  // HETER-COV PATH  (Sigma(Z_i) via modified-Cholesky tree ensembles)
  // =========================================================================
  if (useHeterCov) {

    // Warm-start at per-unit OLS slopes (passed from R).  This matters more
    // for heter-cov than for the BART path because the variance trees fit to
    // squared residuals: starting from zeros would inflate initial d_j(.)
    // estimates and waste burn-in iterations.
    if (Beta_init.n_rows == (uword)nreg && Beta_init.n_cols == (uword)nvar) {
      oldbetas = Beta_init;
    }

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

    for (int rep = 0; rep < R; rep++) {

      // Steps A, B, C, D': Draw mu, variance trees (d_j), phi trees, and mean trees (Delta)
      // fully encapsulated in hetercov_draw_iter; on return
      // hcs.{mu_post, var_models, phi_models, bart_models, delta_Z} reflect
      // the freshly drawn Sigma(Z_i)^{1/2} delta_i.
      hetercov_draw_iter(hcs, oldbetas, mubar.row(0).t(), Amu,
                         pZt, pstd_oldbetas_cols,
                         bart_params, var_params, phi_params,
                         rep, gen);

      // Step E: per-unit Gaussian conjugate update via runiregG with
      //   prior:  beta_i | . ~ N(mu + Sigma(Z_i)^{1/2} delta_i,  Sigma(Z_i))
      // (rootpi = cov::rootpi(ev) is upper-tri with rootpi * rootpi^T = Sigma(Z_i)^{-1};
      //  see src/cov_helpers.cpp:85-89 for the contract).
      for (int reg = 0; reg < nreg; reg++) {
        cov::cov_eval ev{hcs.var_models, hcs.phi_models, (size_t)reg};
        rootpi  = cov::rootpi(ev);
        // Prior precision = Sigma^{-1} = rootpi * rootpi^T (NOT rootpi^T * rootpi;
        //  see audit M3.2 for the convention bug fixed here).
        Abeta   = rootpi * trans(rootpi);
        vec betabari = hcs.mu_post + hcs.delta_Z.row(reg).t();
        Abetabar = Abeta * betabari;

        runiregout_struct = runiregG(regdata_vector[reg].y, regdata_vector[reg].X,
                                      regdata_vector[reg].XpX, regdata_vector[reg].Xpy,
                                      tau[reg], Abeta, Abetabar, nu_e, ssq[reg]);
        oldbetas(reg, span::all) = trans(runiregout_struct.beta);
        tau[reg] = runiregout_struct.sigmasq;
      }

      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);

      if ((rep+1)%keep==0) {
        mkeep = (rep+1)/keep;
        betadraw.slice(mkeep-1) = oldbetas;
        taudraw(mkeep-1, span::all) = trans(tau);

        double ll = 0.0;
        for (int r = 0; r < nreg; r++) {
          vec resid = regdata_vector[r].y - regdata_vector[r].X * trans(oldbetas.row(r));
          int ni = resid.n_elem;
          ll += -0.5 * ni * log(2.0 * M_PI * tau[r]) - 0.5 * dot(resid, resid) / tau[r];
        }
        loglike[mkeep-1] = ll;

        hetercov_store(hcs, mkeep);
      }
    }

    if (nprint>0) endMcmcTimer();

    // Heter-cov return: hetercov_pack supplies tree ensembles + cubes + mu_draw.
    List out = hetercov_pack(hcs);
    out["betadraw"] = betadraw;
    out["taudraw"]  = taudraw;
    out["loglike"]  = loglike;
    return out;
  }

  // =========================================================================
  // BART PATH (sum-of-trees Delta(Z) over mixture-of-normals u_i)
  // =========================================================================
  if (useBART && drawdelta) {

    MixBartState mbs;
    mixbart_init(mbs, nreg, nz, nvar, R, keep, bart_params, oldprob, ind);

    // BART operates on standardized betas; std_oldbetas is the caller-owned
    // working buffer that update_stdoldbetas writes into and that BART reads
    // residuals from.  pstd_oldbetas_cols[j] points to column j of this
    // buffer; the mixbart helpers do NOT own this storage.
    mat std_oldbetas = oldbetas;
    std::vector<double*> pstd_oldbetas_cols(nvar);
    for (int i = 0; i < nvar; i++) pstd_oldbetas_cols[i] = std_oldbetas.colptr(i);
    mat Zt = Z.t();
    double* pZt = Zt.memptr();
    arn gen;

    if (nprint>0) startMcmcTimer();

    // MCMC loop (BART path)
    for (int rep = 0; rep < R; rep++) {

      // Steps 2 & 3: Standardize unit-level betas and draw mean trees (Delta)
      mixbart_draw_iter(mbs, oldbetas, mubar, Amu, nu, V, a, pZt, pstd_oldbetas_cols, bart_params, rep, gen);
      List oldcomp = mbs.oldcomp;
      oldprob = mbs.oldprob;
      ind = mbs.ind;
      mat delta_Z = mbs.delta_Z;

      // Step E: Draw unit-level parameters beta_i | Sigma(Z_i)^{-1}, etc.
      for (int reg = 0; reg < nreg; reg++) {
        List oldcompreg = oldcomp[ind[reg] - 1];
        rootpi = as<mat>(oldcompreg[1]);   // bayesm rooti: rootpi * rootpi^T = Sigma_k^{-1}
        // Prior precision = Sigma_k^{-1} = rootpi * rootpi^T (NOT rootpi^T * rootpi;
        //  see audit M3.2 for the convention bug fixed here).
        Abeta = rootpi * trans(rootpi);
        vec betabari = as<vec>(oldcompreg[0]) + trans(delta_Z.row(reg));
        Abetabar = Abeta * betabari;

        runiregout_struct = runiregG(regdata_vector[reg].y, regdata_vector[reg].X,
                                      regdata_vector[reg].XpX, regdata_vector[reg].Xpy,
                                      tau[reg], Abeta, Abetabar, nu_e, ssq[reg]);
        oldbetas(reg, span::all) = trans(runiregout_struct.beta);
        tau[reg] = runiregout_struct.sigmasq;
      }

      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);

      if ((rep+1)%keep==0) {
        mkeep = (rep+1)/keep;
        betadraw.slice(mkeep-1) = oldbetas;
        taudraw(mkeep-1, span::all) = trans(tau);
        probdraw(mkeep-1, span::all) = trans(oldprob);
        compdraw[mkeep-1] = oldcomp;

        // Compute log-likelihood: sum of Gaussian log-densities across units
        double ll = 0.0;
        for (int r = 0; r < nreg; r++) {
          vec resid = regdata_vector[r].y - regdata_vector[r].X * trans(oldbetas.row(r));
          int ni = resid.n_elem;
          ll += -0.5 * ni * log(2.0 * M_PI * tau[r]) - 0.5 * dot(resid, resid) / tau[r];
        }
        loglike[mkeep-1] = ll;

        mixbart_store(mbs, mkeep);
      }
    }

    if (nprint>0) endMcmcTimer();

    // Package BART output
    nmix = List::create(Named("probdraw") = probdraw,
        Named("zdraw") = R_NilValue,
        Named("compdraw") = compdraw);

    List out = mixbart_pack(mbs);
    out["taudraw"] = taudraw;
    out["betadraw"] = betadraw;
    out["nmix"] = nmix;
    out["loglike"] = loglike;
    return out;

  } else {

    // =========================================================================
    // ORIGINAL (NON-BART) PATH -- from bayesm
    // =========================================================================
    if (nprint>0) startMcmcTimer();

    for (int rep = 0; rep<R; rep++){

      // Step 4: Draw mixture component assignments and update global mixture parameters
      // ind,p need initialization comps is drawn first in sub-Gibbs
      List mgout;
      if(drawdelta) {
        olddelta.reshape(nvar,nz);
        mgout = rmixGibbs(oldbetas-Z*trans(olddelta),mubar,Amu,nu,V,a,oldprob,ind);
      } else {
        mgout = rmixGibbs(oldbetas,mubar,Amu,nu,V,a,oldprob,ind);
      }

      List oldcomp = mgout["comps"];
      oldprob = as<vec>(mgout["p"]);
      ind = as<vec>(mgout["z"]);

      //now draw delta | {beta_i}, ind, comps
      if(drawdelta) olddelta = drawDelta(Z,oldbetas,ind,oldcomp,deltabar,Ad);

      // Step E: Draw unit-level parameters beta_i | Sigma_k^{-1}, etc.
      for(int reg = 0; reg<nreg; reg++){
        List oldcompreg = oldcomp[ind[reg]-1];
        rootpi = as<mat>(oldcompreg[1]);

        //note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
        if(drawdelta){
          olddelta.reshape(nvar,nz);
          betabar = as<vec>(oldcompreg[0])+olddelta*vectorise(Z(reg,span::all));
        } else {
          betabar = as<vec>(oldcompreg[0]);
        }

        // Prior precision = Sigma_k^{-1} = rootpi * rootpi^T (bayesm rooti convention;
        //  see audit M3.2 for the convention bug fixed here).
        Abeta = rootpi*trans(rootpi);
        Abetabar = Abeta*betabar;

        runiregout_struct = runiregG(regdata_vector[reg].y, regdata_vector[reg].X,
                                regdata_vector[reg].XpX, regdata_vector[reg].Xpy,
                                tau[reg], Abeta, Abetabar, nu_e, ssq[reg]);

        oldbetas(reg,span::all) = trans(runiregout_struct.beta);
        tau[reg] = runiregout_struct.sigmasq;
      }

      //print time to completion and draw # every nprint'th draw
      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);

      if((rep+1)%keep==0){
        mkeep = (rep+1)/keep;
        taudraw(mkeep-1, span::all) = trans(tau);
        betadraw.slice(mkeep-1) = oldbetas;
        probdraw(mkeep-1, span::all) = trans(oldprob);
        if(drawdelta) Deltadraw(mkeep-1, span::all) = trans(vectorise(olddelta));
        compdraw[mkeep-1] = oldcomp;

        // Compute log-likelihood: sum of Gaussian log-densities across units
        double ll = 0.0;
        for (int r = 0; r < nreg; r++) {
          vec resid = regdata_vector[r].y - regdata_vector[r].X * trans(oldbetas.row(r));
          int ni = resid.n_elem;
          ll += -0.5 * ni * log(2.0 * M_PI * tau[r]) - 0.5 * dot(resid, resid) / tau[r];
        }
        loglike[mkeep-1] = ll;
      }
    }

    if (nprint>0) endMcmcTimer();

    nmix = List::create(Named("probdraw") = probdraw,
              Named("zdraw") = R_NilValue,
              Named("compdraw") = compdraw);

    if(drawdelta){
      return(List::create(
        Named("taudraw") = taudraw,
        Named("Deltadraw") = Deltadraw,
        Named("betadraw") = betadraw,
        Named("loglike") = loglike,
        Named("nmix") = nmix));
    } else {
      return(List::create(
        Named("taudraw") = taudraw,
        Named("betadraw") = betadraw,
        Named("loglike") = loglike,
        Named("nmix") = nmix));
    }
  }
}
