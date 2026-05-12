#include "bayesm.HART.h"

//[[Rcpp::export]]
List rhierLinearMixture_rcpp_loop(List const& regdata, mat const& Z,
                                  vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,
                                  double nu, mat const& V, double nu_e, vec const& ssq,
                                  int R, int keep, int nprint, bool drawdelta,
                                  mat olddelta, vec const& a, vec oldprob, vec ind, vec tau,
                                  bool useBART = false, List const& bart_params = List::create()){

// Wayne Taylor 10/02/2014
// Modified by Thomas Wiemann 2025 -- added HART/BART support

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
  // BART PATH
  // =========================================================================
  if (useBART && drawdelta) {

    mat delta_Z = zeros<mat>(nreg, nvar);
    vec pred_i(nvar);
    vec scaled_pred(nvar);

    // Initialize BART tree storage
    std::vector<bart> bart_models;
    std::vector<std::stringstream> treess(nvar);
    cube varcount, varprob;

    size_t num_trees = bart_params["num_trees"];
    for (int k = 0; k < nvar; k++) {
      treess[k].precision(10);
      treess[k] << R / keep << " " << num_trees << " " << nz << endl;
    }
    varcount.set_size(nz, nvar, R/keep);
    varprob.set_size(nz, nvar, R/keep);
    varcount.zeros();
    varprob.zeros();

    // Create column-major Z for BART
    mat Zt = Z.t();
    double* pZt = Zt.memptr();
    mat std_oldbetas = oldbetas;
    std::vector<double*> pstd_oldbetas_cols(nvar);
    for (int i = 0; i < nvar; i++) {
      pstd_oldbetas_cols[i] = std_oldbetas.colptr(i);
    }

    // Random number generator for BART
    arn gen;

    // DART burn-in
    int burn = 0;
    bool sparse = bart_params["sparse"];
    if (sparse) burn = bart_params["burn"];

    if (nprint>0) startMcmcTimer();

    // MCMC loop (BART path)
    for (int rep = 0; rep < R; rep++) {

      // Draw comps, ind, p | {beta_i}, delta_Z
      List mgout;
      if (rep == 0) {
        mgout = rmixGibbs(oldbetas, mubar, Amu, nu, V, a, oldprob, ind);
        update_stdoldbetas(oldbetas, pstd_oldbetas_cols, ind, mgout["comps"]);
        bart_models = initializeBART(pZt, nreg, nz, pstd_oldbetas_cols, bart_params, gen);
      } else {
        mgout = rmixGibbs(oldbetas - delta_Z, mubar, Amu, nu, V, a, oldprob, ind);
      }

      List oldcomp = mgout["comps"];
      oldprob = as<vec>(mgout["p"]);
      ind = as<vec>(mgout["z"]);

      // Start DART after burn-in
      if (sparse && rep == burn + 1) {
        for (int i = 0; i < nvar; i++) {
          bart_models[i].startdart();
        }
      }

      // Normalize oldbetas and draw BART
      update_stdoldbetas(oldbetas, pstd_oldbetas_cols, ind, oldcomp);
      for (int i = 0; i < nvar; i++) {
        bart_models[i].draw(1.0, gen);
      }

      // Compute delta_Z
      for (int i = 0; i < nreg; i++) {
        for (int j = 0; j < nvar; j++) {
          pred_i(j) = bart_models[j].f(i);
        }
        List comp_i = oldcomp[ind[i] - 1];
        mat rootii = trans(as<mat>(comp_i[1]));
        scaled_pred = solve(trimatl(rootii), pred_i);
        delta_Z.row(i) = trans(scaled_pred);
      }

      // Draw beta_i | ind[i], comps, delta_Z, data
      for (int reg = 0; reg < nreg; reg++) {
        List oldcompreg = oldcomp[ind[reg] - 1];
        rootpi = as<mat>(oldcompreg[1]);
        Abeta = trans(rootpi) * rootpi;
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

        // Store BART trees and variable usage
        for (int i = 0; i < nvar; i++) {
          for (size_t j = 0; j < bart_models[i].getm(); j++) {
            treess[i] << bart_models[i].gettree(j);
          }
          varcount.slice(mkeep-1).col(i) = conv_to<vec>::from(bart_models[i].getnv());
          varprob.slice(mkeep-1).col(i) = conv_to<vec>::from(bart_models[i].getpv());
        }
      }
    }

    if (nprint>0) endMcmcTimer();

    // Package BART output
    nmix = List::create(Named("probdraw") = probdraw,
        Named("zdraw") = R_NilValue,
        Named("compdraw") = compdraw);

    List bartModels(nvar);
    for (int i = 0; i < nvar; i++) {
      xinfo& xi = bart_models[i].getxinfo();
      Rcpp::List xiret(xi.size());
      for (size_t j = 0; j < xi.size(); j++) {
        Rcpp::NumericVector vtemp(xi[j].size());
        std::copy(xi[j].begin(), xi[j].end(), vtemp.begin());
        xiret[j] = vtemp;
      }
      bartModels[i] = List::create(
        Named("treedraws") = List::create(
          Named("cutpoints") = xiret,
          Named("trees") = Rcpp::CharacterVector(treess[i].str())
        )
      );
    }

    return List::create(
      Named("bart_models") = bartModels,
      Named("varcount") = varcount,
      Named("varprob") = varprob,
      Named("betadraw") = betadraw,
      Named("taudraw") = taudraw,
      Named("loglike") = loglike,
      Named("nmix") = nmix);

  } else {

    // =========================================================================
    // ORIGINAL (NON-BART) PATH -- from bayesm
    // =========================================================================
    if (nprint>0) startMcmcTimer();

    for (int rep = 0; rep<R; rep++){
   
      //first draw comps,ind,p | {beta_i}, delta
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
     
      //loop over all regression equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
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
      
        Abeta = trans(rootpi)*rootpi;
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
