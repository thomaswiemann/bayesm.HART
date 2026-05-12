#include "bayesm.HART.h"
  
//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
double llnegbinpooled(std::vector<moments> regdata_vector, mat Beta, double alpha){
  
// Wayne Taylor 12/01/2014

// "Unlists" the regdata and calculates the negative binomial loglikelihood using individual-level betas
  
  int nreg = regdata_vector.size();
  double ll = 0.0;
  
  for(int reg = 0; reg<nreg; reg++){
  vec lambda = exp(regdata_vector[reg].X*trans(Beta(reg,span::all)));
  ll = ll + llnegbin(regdata_vector[reg].y,lambda,alpha,TRUE);
  }
  
  return(ll);
}

// [[Rcpp::export]]
List rhierNegbinRw_rcpp_loop(List const& regdata, List const& hessdata, mat const& Z,
                             mat Beta, vec const& deltabar, mat const& Ad,
                             mat const& mubar, mat const& Amu,
                             double nu, mat const& V, double a, double b,
                             int R, int keep, double sbeta, double alphacroot, int nprint,
                             bool drawdelta, mat olddelta, vec const& a_mix,
                             vec oldprob, vec ind,
                             double alpha, bool fixalpha,
                             bool useBART = false, List const& bart_params = List::create()){
                            
// Wayne Taylor 12/01/2014
// Modified by Thomas Wiemann 2025 -- added HART/BART support

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
  cube Betadraw = zeros<cube>(nreg, nvar, R/keep);
  vec alphadraw = zeros<vec>(R/keep);
  vec loglike = zeros<vec>(R/keep);

  // =========================================================================
  // BART PATH
  // =========================================================================
  if (useBART) {

    mat delta_Z = zeros<mat>(nreg, nvar);
    vec pred_i(nvar);
    vec scaled_pred(nvar);

    // BART-specific storage
    mat probdraw(R/keep, oldprob.size());
    List compdraw(R/keep);

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
    mat std_oldbetas = Beta;
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
    for (rep = 0; rep < R; rep++) {

      // Draw comps, ind, p | {beta_i}, delta_Z
      List mgout;
      if (rep == 0) {
        mgout = rmixGibbs(Beta, mubar, Amu, nu, V, a_mix, oldprob, ind);
        update_stdoldbetas(Beta, pstd_oldbetas_cols, ind, mgout["comps"]);
        bart_models = initializeBART(pZt, nreg, nz, pstd_oldbetas_cols, bart_params, gen);
      } else {
        mgout = rmixGibbs(Beta - delta_Z, mubar, Amu, nu, V, a_mix, oldprob, ind);
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
      update_stdoldbetas(Beta, pstd_oldbetas_cols, ind, oldcomp);
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

      // Draw beta_i via MH with per-unit mixture prior
      for (int reg = 0; reg < nreg; reg++) {
        List oldcompreg = oldcomp[ind[reg] - 1];
        mat rootpi = as<mat>(oldcompreg[1]);
        vec betabari = as<vec>(oldcompreg[0]) + trans(delta_Z.row(reg));

        mat Vbetainv_reg = rootpi * trans(rootpi);
        mat betacvar = sbeta * solve(regdata_vector[reg].hess + Vbetainv_reg, eye(nvar, nvar));
        mat betaroot = trans(chol(betacvar));
        vec betac = vectorise(Beta(reg, span::all)) + betaroot * vec(rnorm(nvar));

        oldlpostbeta[reg] = lpostbeta(alpha, trans(Beta(reg, span::all)),
                                       regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
        clpostbeta[reg] = lpostbeta(alpha, betac,
                                     regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
        ldiff = clpostbeta[reg] - oldlpostbeta[reg];
        acc = exp(ldiff);
        if (acc > 1) acc = 1;
        if (acc < 1) { unif = runif(1)[0]; } else { unif = 0; }
        if (unif <= acc) {
          Beta(reg, span::all) = trans(betac);
          nacceptbeta = nacceptbeta + 1;
        }
      }

      // Draw alpha (same in both paths)
      if (!fixalpha) {
        logalphac = log(alpha) + alphacroot * rnorm(1)[0];
        oldlpostalpha = llnegbinpooled(regdata_vector, Beta, alpha) + (a - 1) * log(alpha) - b * alpha;
        clpostalpha = llnegbinpooled(regdata_vector, Beta, exp(logalphac)) + (a - 1) * logalphac - b * exp(logalphac);
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
        Betadraw.slice(mkeep-1) = Beta;
        alphadraw[mkeep-1] = alpha;
        probdraw(mkeep-1, span::all) = trans(oldprob);
        compdraw[mkeep-1] = oldcomp;
        loglike[mkeep-1] = llnegbinpooled(regdata_vector, Beta, alpha);

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
    List nmix = List::create(Named("probdraw") = probdraw,
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
      Named("loglike") = loglike,
      Named("Betadraw") = Betadraw,
      Named("alphadraw") = alphadraw,
      Named("nmix") = nmix,
      Named("acceptrbeta") = nacceptbeta / (R * nreg * 1.0) * 100,
      Named("acceptralpha") = nacceptalpha / (R * 1.0) * 100);

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
      
      //first draw comps,ind,p | {beta_i}, delta
      // ind,p need initialization comps is drawn first in sub-Gibbs
      List mgout;
      if (drawdelta) {
          olddelta.reshape(nvar, nz);
          mgout = rmixGibbs(Beta - delta_Z, mubar, Amu, nu, V, a_mix, oldprob, ind);
      } else {
          mgout = rmixGibbs(Beta, mubar, Amu, nu, V, a_mix, oldprob, ind);
      }

      List oldcomp = mgout["comps"];
      oldprob = as<vec>(mgout["p"]);
      ind = as<vec>(mgout["z"]);

      //now draw delta | {beta_i}, ind, comps
      if (drawdelta) {
          olddelta = drawDelta(Z, Beta, ind, oldcomp, deltabar, Ad);
          olddelta.reshape(nvar, nz);
          delta_Z = Z * trans(olddelta);
      }

      //loop over all regression equations drawing beta_i | ind[i],comp[ind[i]]
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
          vec betac = vectorise(Beta(reg,span::all)) + betaroot*vec(rnorm(nvar));
         
          oldlpostbeta[reg] = lpostbeta(alpha, trans(Beta(reg,span::all)), regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
          clpostbeta[reg] = lpostbeta(alpha, betac, regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootpi);
          ldiff = clpostbeta[reg] - oldlpostbeta[reg];
          acc = exp(ldiff);
          if (acc > 1) acc = 1;    
          if(acc < 1) {unif=runif(1)[0];} else {unif=0;}
          if (unif <= acc){
            Beta(reg,span::all) = trans(betac);
            nacceptbeta = nacceptbeta + 1;
          }
      }
      
      // Draw alpha
      if (!fixalpha){
        logalphac = log(alpha) + alphacroot*rnorm(1)[0];
        oldlpostalpha = llnegbinpooled(regdata_vector,Beta,alpha)+(a-1)*log(alpha) - b*alpha;
        clpostalpha = llnegbinpooled(regdata_vector,Beta,exp(logalphac))+(a-1)*logalphac - b*exp(logalphac);
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
        Betadraw.slice(mkeep-1) = Beta;
        alphadraw[mkeep-1] = alpha;
        probdraw(mkeep-1, span::all) = trans(oldprob);
        if(drawdelta) Deltadraw(mkeep-1,span::all) = trans(vectorise(olddelta));
        compdraw[mkeep-1] = oldcomp;
        loglike[mkeep-1] = llnegbinpooled(regdata_vector,Beta,alpha);
        } 
    }
    
    if (nprint>0) endMcmcTimer();
    
    List nmix = List::create(Named("probdraw") = probdraw,
              Named("zdraw") = R_NilValue,
              Named("compdraw") = compdraw);

    if(drawdelta){
      return List::create(
        Named("loglike") = loglike,
        Named("Betadraw") = Betadraw,
        Named("alphadraw") = alphadraw,
        Named("Deltadraw") = Deltadraw,
        Named("nmix") = nmix,
        Named("acceptrbeta") = nacceptbeta/(R*nreg*1.0)*100,
        Named("acceptralpha") = nacceptalpha/(R*1.0)*100);
    } else {
      return List::create(
        Named("loglike") = loglike,
        Named("Betadraw") = Betadraw,
        Named("alphadraw") = alphadraw,
        Named("nmix") = nmix,
        Named("acceptrbeta") = nacceptbeta/(R*nreg*1.0)*100,
        Named("acceptralpha") = nacceptalpha/(R*1.0)*100);
    }
  }
}
