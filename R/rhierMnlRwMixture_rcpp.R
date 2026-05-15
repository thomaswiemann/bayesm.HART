#' Bayesian Multinomial Logit Model with HART Prior
#' @description
#' `rhierMnlRwMixture` implements an MCMC algorithm for a Bayesian multinomial 
#' logit model with a Hierarchical Additive Regression Trees (HART) prior. HART 
#' is a hierarchical nonparametric prior that allows for flexible modeling of the 
#' representative consumer as a function of potentially many observed characteristics.
#' 
#' @note Currently, the mixture component is not implemented. Please use 
#' `ncomp = 1` in the Prior specification.
#'
#' @param Data A list containing:
#'   - `p`: Number of choice alternatives (integer).
#'   - `lgtdata`: A list of length `nlgt`. Each element `lgtdata[[i]]` must be a list with:
#'     - `y`: `n_i x 1` vector of multinomial outcomes (1 to `p`).
#'     - `X`: `(n_i * p) x nvar` matrix of alternative-specific attributes.
#'   - `Z` (optional): `nlgt x nz` matrix of observed characteristics for each unit.
#'     Should NOT contain an intercept and should be centered.
#' @param Prior A list containing prior parameters:
#'   - `ncomp`: Number of mixture components (required).
#'   - `a` (optional): `ncomp x 1` vector of Dirichlet prior parameters for mixture weights `pvec` (default: `rep(5, ncomp)`).
#'   - `deltabar` (optional): `nz * nvar x 1` prior mean for `vec(Delta)` (default: 0). Ignored if HART is used.
#'   - `Ad` (optional): Prior precision matrix for `vec(Delta)` (default: `0.01 * I`). Ignored if HART is used.
#'   - `mubar` (optional): `nvar x 1` prior mean vector for mixture component means (default: 0 if unrestricted, 2 if restricted).
#'   - `Amu` (optional): Prior precision for mixture component means (default: 0.01 if unrestricted, 0.1 if restricted).
#'   - `nu` (optional): Degrees of freedom for IW prior on component `Sigma` (default: `nvar+3` if unrestricted, `nvar+15` if restricted).
#'   - `V` (optional): Location matrix for IW prior on component `Sigma` (default: `nu * I` or scaled based on restriction).
#'   - `SignRes` (optional): `nvar x 1` vector of sign restrictions. Must contain values of 0, -1, or 1. The value 0 means no restriction, -1 ensures the coefficient is negative, and 1 ensures the coefficient is positive. For example, `SignRes = c(0,1,-1)` means the first coefficient is unconstrained, the second will be positive, and the third will be negative. Default: `rep(0, nvar)`.
#'   - `bart` (optional): List of parameters for the HART prior. If specified, this models the representative consumer \eqn{\Delta(Z)} as a scaled sum-of-trees factor model. See Details.
#'   - `vartree` (optional, *experimental extension beyond Wiemann 2025*): List of parameters enabling heteroscedastic covariance \eqn{\Sigma(Z_i)} via product-of-trees variance models on the Modified Cholesky diagonal \eqn{d_j(\cdot)}. Requires `bart` to also be specified and `ncomp = 1`. Compatible with `SignRes` (see "Heteroscedastic Covariance" in Details). When `nvar > 1`, the package automatically promotes to the full-Cholesky structure described under `phitree` (diagonal-only \eqn{\Sigma(Z_i)} would force every conditional cross-correlation to zero, which is essentially never a believable posterior).
#'   - `phitree` (optional, *experimental extension beyond Wiemann 2025*): List of parameters enabling sum-of-trees regression models on the Modified Cholesky off-diagonals \eqn{\phi_{jk}(\cdot)}. Requires `vartree`. Adds the full-Cholesky structure on top of the diagonal `vartree` model. Auto-enabled when `vartree` is supplied and `nvar > 1`.
#' @param Mcmc A list containing MCMC parameters:
#'   - `R`: Number of MCMC iterations (required).
#'   - `keep` (optional): Thinning parameter (default: 1).
#'   - `nprint` (optional): Print progress every `nprint` draws (default: 100, 0 for none).
#'   - `s` (optional): RW Metropolis scaling parameter (default: `2.93 / sqrt(nvar)`).
#'   - `w` (optional): Fractional likelihood weighting parameter (default: 0.1).
#' @param r_verbose Logical. Print startup messages? Default TRUE.
#'
#' @details
#' ## Model Specification
#' \eqn{y_i \sim MNL(X_i, \beta_i)} for unit \(i = 1, \ldots, nlgt\).
#' The unit-level coefficients (part-worths) \eqn{\beta_i} are modeled as:
#' \deqn{\beta_i = \Delta(Z_i) + u_i}
#' where \eqn{\Delta(Z_i)} is the *representative consumer* component, which depends on observed characteristics \eqn{Z_i}, and \eqn{u_i} is the unobserved heterogeneity component.
#'
#' The representative consumer component is specified as:
#'   - If `Z` is provided and `Prior$bart` is `NULL`: \eqn{\Delta(Z_i) = Z_i \Delta} where \eqn{\Delta} is an `nz x nvar` matrix (linear hierarchical model).
#'   - If `Z` is provided and `Prior$bart` is a list: \eqn{\Delta(Z_i)} is modeled with a HART prior (scaled sum-of-trees factor model).
#'   - If `Z` is `NULL`: \eqn{\Delta(Z_i) = 0}.
#'
#' With `ncomp = 1` (currently required), the unobserved heterogeneity component follows:
#' \deqn{u_i \sim N(\mu_1, \Sigma_1)}
#'
#' ## Prior Specifications
#' - **Mixture weights**: \eqn{pvec \sim Dirichlet(a)}
#' - **Linear model**: \eqn{\delta = vec(\Delta) \sim N(deltabar, A_d^{-1})}
#' - **Mixture component means**: \eqn{\mu_j \sim N(mubar, \Sigma_j \otimes Amu^{-1})} (covariance scaled by \eqn{\Sigma_j})
#' - **Mixture component covariance**: \eqn{\Sigma_j \sim IW(\nu, V)}
#' - **HART model**: A sum-of-trees prior is placed on each factor of the scaled sum-of-trees model (see HART details below).
#'
#' ## HART Prior Details
#' If `Prior$bart` is a list, it specifies a HART prior for the representative consumer \eqn{\Delta(Z)}. 
#' This replaces the conventional linear hierarchical specification. The HART prior models the representative 
#' consumer using a scaled vector of `nvar` sum-of-trees models.
#'
#' **HART Parameters** (defaults used if not specified in `Prior$bart`):
#'   - `num_trees`: Number of trees H in each sum-of-trees model (default: 200).
#'   - `power`, `base`: Parameters for the tree structure prior. The probability of a node at depth `q` splitting is \eqn{\alpha(1+q)^{-\beta}}, where `base`=\eqn{\alpha} and `power`=\eqn{\beta}. Defaults are `base=0.95`, `power=2`, which strongly favors shallow trees.
#'   - `tau`: Parameter controlling the prior variance of terminal leaf coefficients. The default is \eqn{\tau = 1/\sqrt{H}} where \eqn{\lambda_{dhg} \sim N(0, \tau^2)} for terminal leaf coefficients.
#'   - `numcut`: Number of grid points for proposing splitting rules for continuous variables (default: 100).
#'   - `sparse`: If `TRUE`, use the Dirichlet HART prior to induce sparsity in variable selection (default: `FALSE`).
#'
#' **Dirichlet HART** (`sparse = TRUE`): The Dirichlet HART model augments the HART prior to induce sparsity in variable selection, following Linero (2018). Instead of uniform probability for selecting splitting variables, the selection probabilities \eqn{\tau = (\tau^{(1)}, \ldots, \tau^{(K)})} are given a sparse Dirichlet prior: \eqn{(\tau^{(1)}, \ldots, \tau^{(K)}) \sim Dirichlet(\theta/K, \ldots, \theta/K)}, where K is the number of characteristics. The concentration parameter \eqn{\theta} is given a hierarchical prior: \eqn{\theta/(\theta+\rho) \sim Beta(a,b)}.
#'   - `a`, `b`: Shape parameters for the Beta hyperprior. The default (`a=0.5, b=1`) induces sparsity where few variables have high selection probabilities.
#'   - `rho`: Parameter influencing sparsity. Default is the number of characteristics K. Reducing rho below K encourages greater sparsity.
#'   - `theta`: When used, sets Dirichlet concentration parameter without additional hyper-prior (default: 0.0).
#'   - `burn`: Number of internal burn-in iterations for the Dirichlet HART sampler before variable selection is allowed (default: 100).
#'
#' ## Heteroscedastic Covariance \eqn{\Sigma(Z_i)} (experimental extension)
#'
#' Wiemann (2025) treats \eqn{\Sigma} as global with an inverse-Wishart prior, used as a fixed factor-model loading \eqn{\Sigma^{1/2}}. When `Prior$vartree` is supplied, the package replaces this with unit-specific
#' \deqn{\Sigma(Z_i)^{-1} = L(Z_i)^\top D(Z_i)^{-1} L(Z_i)}
#' where \eqn{L(Z_i)} is unit lower-triangular with \eqn{L_{jk} = -\phi_{jk}(Z_i)} for \eqn{k < j} and \eqn{D(Z_i) = \mathrm{diag}(d_1(Z_i), \ldots, d_D(Z_i))}, \eqn{d_j > 0}. Each \eqn{d_j(\cdot)} is modeled as a product of trees with \eqn{\chi^{-2}} leaves (Pratola et al.); each \eqn{\phi_{jk}(\cdot)} (when `Prior$phitree` is supplied) is modeled as a sum of trees with \eqn{N(0, \tau^2)} leaves. Only the single-component path (`ncomp = 1`) is supported in this mode.
#'
#' **`Prior$vartree` parameters** (defaults shown):
#'   - `num_trees` (40): Number of trees per dimension \eqn{m'}.
#'   - `nu` (10), `lambda` (auto-calibrated from \eqn{\mathrm{var}(\theta_i)}): Baseline parameters of the \eqn{\chi^{-2}} prior. Pratola per-tree calibration is applied internally.
#'   - `power` (2), `base` (0.95), `numcut` (100): Tree-structure prior parameters.
#'   - DART hyperparameters (`sparse`, `a`, `b`, `rho`, `theta`, `omega`, `aug`, `burn`): same meaning as `Prior$bart`.
#'
#' **`Prior$phitree` parameters** (defaults shown):
#'   - `num_trees` (40): Number of trees per \eqn{(j,k)} pair \eqn{m''}.
#'   - `tau` (\eqn{1/\sqrt{m''}}): Prior standard deviation of leaf \eqn{\lambda_{jkh}}.
#'   - `power` (2), `base` (0.95), `numcut` (100): Tree-structure prior parameters.
#'   - `nmin` (2), `ess_min` (5): Leaf-admissibility floors. The default `ess_min = 5` requires the effective sample size \eqn{\sum_{i \in \ell} (\theta_i^{(k)} - \mu^{(k)})^2 / d_j(Z_i)} of each candidate leaf to exceed 5; `nmin = 2` is the corresponding raw-count floor.
#'   - DART hyperparameters: same meaning as `Prior$bart`.
#'
#' When this extension is active, the returned object additionally inherits class `"rhierMnlRwMixtureHeterCov"` and contains slots `var_models`, `phi_models` (jagged, lower-triangular; `NULL` if `nvar == 1`), and `mu_draw` (single-component posterior mean draws). `predict()` dispatches on this class to evaluate \eqn{\Sigma(Z^*)} at new \eqn{Z^*}.
#'
#' Sign restrictions are fully compatible with this extension. Under `SignRes != 0`, the variance/phi trees model \eqn{\Sigma(Z_i)} on the unconstrained scale \eqn{\beta^*_i = \log|\beta_i|/SignRes} (exactly as in the standard MNL HART path), and the post-loop transform \eqn{\beta_i = SignRes \cdot \exp(\beta^*_i)} is applied unconditionally to every kept draw.
#'
#' ## Sign Restrictions
#' If `SignRes[k]` is non-zero, the k-th coefficient \eqn{\beta_{ik}} is modeled as
#' \deqn{\beta_{ik} = SignRes[k] \cdot \exp(\beta^*_{ik}).}
#' The `betadraw` output contains the draws for \eqn{\beta_{ik}} (with the restriction applied).
#' The `nmix` output contains draws for the *unrestricted* mixture components.
#' 
#' **Note:** Care should be taken when selecting priors on any sign restricted coefficients.
#'
#' @return A list containing:
#'   - `Deltadraw`: If `Z` provided and `bart=NULL`, `(R/keep) x (nz * nvar)` matrix of `vec(Delta)` draws.
#'   - `betadraw`: `nlgt x nvar x (R/keep)` array of unit-level `beta_i` draws.
#'   - `nmix`: List containing mixture draws with components (omitted when `vartree` is supplied):
#'     - `probdraw`: `(R/keep) x ncomp` matrix of mixture component probabilities.
#'     - `zdraw`: `(R/keep) x nlgt` matrix of component assignments for each unit.
#'     - `compdraw`: `(R/keep)` list of `ncomp` lists. `compdraw[[r]][[j]] = list(mu, rooti)` contains the draw of \eqn{\mu_j} and \eqn{\Sigma_j^{-1/2}} for component `j` at kept draw `r`.
#'   - `loglike`: `(R/keep) x 1` vector of log-likelihood values at kept draws.
#'   - `SignRes`: `nvar x 1` vector of sign restrictions used.
#'   - `acceptrbeta`: Metropolis acceptance rate (percent) for the unit-level
#'     MNL random-walk updates of `beta_i`.
#'   - `bart_models`: If HART used, list of length `nvar` containing tree structures and related parameters for the mean trees \eqn{\delta_j(\cdot)}.
#'   - `var_models` (only with `vartree`): list of length `nvar` of variance-tree ensembles for \eqn{d_j(\cdot)} (product of trees with \eqn{\chi^{-2}} leaves).
#'   - `phi_models` (only with `vartree` + `phitree`): jagged list of length `nvar`. `phi_models[[1]]` is `NULL`; for `j > 1`, `phi_models[[j]]` is a list of length `j-1` containing the sum-of-trees ensemble for \eqn{\phi_{jk}(\cdot)}.
#'   - `mu_draw` (only with `vartree`): `(R/keep) x nvar` matrix of \eqn{\mu} draws (single-component path).
#'
#' @references
#' Chipman, Hugh A., Edward I. George, and Robert E. McCulloch (2010). "BART: Bayesian Additive Regression Trees." Annals of Applied Statistics 4.1.
#' 
#' Linero, Antonio R. (2018). "Bayesian regression trees for high-dimensional prediction and variable selection." Journal of the American Statistical Association 113.522, pp. 626-636.
#' 
#' Rossi, Peter E., Greg M. Allenby, and Robert McCulloch (2009). Bayesian Statistics and Marketing. Reprint. Wiley Series in Probability and Statistics. Chichester: Wiley.
#' 
#' Rossi, Peter (2023). bayesm: Bayesian Inference for Marketing/Micro-Econometrics. Comprehensive R Archive Network.
#' 
#' Wiemann, Thomas (2025). "Personalization with HART." Working paper.
#'
#' @seealso [predict.rhierMnlRwMixture()]
#'
#' @author Peter Rossi (original bayesm code), Thomas Wiemann (HART modifications).
#'
#' @export
rhierMnlRwMixture=function(Data,Prior,Mcmc, r_verbose = TRUE){
  #
  # revision history:
  #   12/04 changed by rossi to fix bug in drawdelta when there is zero/one unit in a mixture component
  #   09/05 added loglike output, changed to reflect new argument order in llmnl, mnlHess 
  #   12/05 changed weighting scheme to (1-w)logl_i + w*Lbar (normalized) 
  #   03/07 added classes
  #   09/08 changed Dirichlet a check
  #   04/15 by Wayne Taylor: added nprint option to MCMC argument
  #   07/16 by Wayne Taylor: added sign restrictions
  #   10/10 by Dan Yavorsky: changed default priors when sign restrictions imposed
  #   12/18 by Peter Rossi: print out vector of sign-restrictions
  #   12/18 by Peter Rossi: changed Hessian for constrained parameters to reflect
  #                         reparameterization
  #   7/19 by Peter Rossi: further fixes for reparameterization of constrained parms
  #                        fixed Hessian as well as problems with initial values to find
  #                        constrained optima used to tune Metropolis. 
  #                        switched to Nelder-Mead to find constrained pooled optimum
  #                         BFGS sometimes has trouble with reparameterized model 
  #   6/20 by Peter Rossi: fixed check on size of betapooled to correct indexing
  #   6/25 by Thomas Wiemann: added HART prior + roxygen2 docs
  #
  if(missing(Data)) {pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")}
  if(is.null(Data$p)) {pandterm("Requires Data element p (# choice alternatives)") }
  p=Data$p
  if(is.null(Data$lgtdata)) {pandterm("Requires Data element lgtdata (list of data for each unit)")}
  lgtdata=Data$lgtdata
  nlgt=length(lgtdata)
  drawdelta=TRUE
  if(is.null(Data$Z)) { cat("Z not specified",fill=TRUE); fsh() ; drawdelta=FALSE}
  else {if (!is.matrix(Data$Z)) {pandterm("Z must be a matrix")}
    else {if (nrow(Data$Z) != nlgt) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number logits ",nlgt))}
      else {Z=Data$Z}}}
  if(drawdelta) {
    nz=ncol(Z)
    colmeans=apply(Z,2,mean)
    if(sum(colmeans) > .00001) 
    {pandterm(paste("Z does not appear to be de-meaned: colmeans= ",colmeans))}
  }
  #
  # check lgtdata for validity
  #
  ypooled=NULL
  Xpooled=NULL
  if(!is.null(lgtdata[[1]]$X & is.matrix(lgtdata[[1]]$X))) {oldncol=ncol(lgtdata[[1]]$X)}
  for (i in 1:nlgt) 
  {
    if(is.null(lgtdata[[i]]$y)) {pandterm(paste0("Requires element y of lgtdata[[",i,"]]"))}
    if(is.null(lgtdata[[i]]$X)) {pandterm(paste0("Requires element X of lgtdata[[",i,"]]"))}
    if(!is.matrix(lgtdata[[i]]$X)) {pandterm(paste0("lgtdata[[",i,"]]$X must be a matrix"))}
    if(!is.vector(lgtdata[[i]]$y, mode = "numeric") & !is.vector(lgtdata[[i]]$y, mode = "logical") & !is.matrix(lgtdata[[i]]$y)) 
    {pandterm(paste0("lgtdata[[",i,"]]$y must be a numeric or logical vector or matrix"))}
    if(is.matrix(lgtdata[[i]]$y)) { if(ncol(lgtdata[[i]]$y)>1) { pandterm(paste0("lgtdata[[",i,"]]$y must be a vector or one-column matrix")) } }
    ypooled=c(ypooled,lgtdata[[i]]$y)
    nrowX=nrow(lgtdata[[i]]$X)
    if((nrowX/p) !=length(lgtdata[[i]]$y)) {pandterm(paste("nrow(X) ne p*length(yi); exception at unit",i))}
    newncol=ncol(lgtdata[[i]]$X)
    if(newncol != oldncol) {pandterm(paste("All X elements must have same # of cols; exception at unit",i))}
    Xpooled=rbind(Xpooled,lgtdata[[i]]$X)
    oldncol=newncol
  }
  nvar=ncol(Xpooled)
  levely=as.numeric(levels(as.factor(ypooled)))
  if(length(levely) != p) {pandterm(paste("y takes on ",length(levely)," values -- must be = p"))}
  bady=FALSE
  for (i in 1:p )
  {
    if(levely[i] != i) bady=TRUE
  }
  if (r_verbose) {
    cat("Table of Y values pooled over all units",fill=TRUE)
    print(table(ypooled))
  }
  if (bady) 
  {pandterm("Invalid Y")}
  #
  # check on prior
  #
  if(missing(Prior)) {pandterm("Requires Prior list argument (at least ncomp)")} 
  if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp (num of mixture components)")} else {ncomp=Prior$ncomp}
  if(ncomp != 1) {pandterm("Only ncomp = 1 is currently supported")}
  if(is.null(Prior$SignRes)) {SignRes=rep(0,nvar)} else {SignRes=Prior$SignRes}                     
  if(length(SignRes) != nvar) {pandterm("The length SignRes must be equal to the dimension of X")}
  if(sum(!(SignRes %in% c(-1,0,1))>0)) {pandterm("All elements of SignRes must be equal to -1, 0, or 1")}
  if(is.null(Prior$mubar) & sum(abs(SignRes))==0) {
    mubar=matrix(rep(0,nvar),nrow=1)
  } else { 
    if(is.null(Prior$mubar) & sum(abs(SignRes)) >0) {
      mubar=matrix(rep(0,nvar)+2*abs(SignRes),nrow=1)
    } else { 
      mubar=matrix(Prior$mubar,nrow=1) } }
  if(ncol(mubar) != nvar) {pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ",ncol(mubar)))}
  if(is.null(Prior$Amu) & sum(abs(SignRes))==0) {
    Amu=matrix(BayesmConstant.A,ncol=1)
  } else {
    if(is.null(Prior$Amu) & sum(abs(SignRes)) >0) {
      Amu=matrix(BayesmConstant.A*10,ncol=1)
    } else {Amu=matrix(Prior$Amu,ncol=1) } }
  if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Am must be a 1 x 1 array")}
  if(is.null(Prior$nu) & sum(abs(SignRes))==0) {
    nu=nvar+BayesmConstant.nuInc
  } else {
    if(is.null(Prior$nu) & sum(abs(SignRes)) >0) {
      nu=nvar+BayesmConstant.nuInc+12
    } else {
      nu=Prior$nu } }
  if(nu < 1) {pandterm("invalid nu value")}
  if(is.null(Prior$V) & sum(abs(SignRes))==0) {
    V=nu*diag(nvar)
  } else {
    if(is.null(Prior$V) & sum(abs(SignRes)) >0) {
      V=nu*diag(abs(SignRes)*0.1+(!abs(SignRes))*4)
    } else {
      V=Prior$V } }
  if(sum(dim(V)==c(nvar,nvar)) !=2) pandterm("Invalid V in prior")
  useBART <- ifelse(!is.null(Prior$bart) & drawdelta, TRUE, FALSE)
  bart_params = Ad = deltabar = NULL # Set to NULL if unused
  if (useBART) {
    # Bart parameter defaults
    bart_params <- .parse_bart_params(Prior$bart, nz)
  } else {
    if(is.null(Prior$Ad) & drawdelta) {Ad=BayesmConstant.A*diag(nvar*nz)} else {Ad=Prior$Ad}
    if(drawdelta) {if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}}
    if(is.null(Prior$deltabar)& drawdelta) {deltabar=rep(0,nz*nvar)} else {deltabar=Prior$deltabar}
    if(drawdelta) {if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}}
  }
  # ---------------------------------------------------------------------------
  # Heteroscedastic covariance Sigma(Z_i) via modified-Cholesky tree ensembles.
  # Additive extension beyond Wiemann (2025); auto-enabled by Prior$vartree.
  # See R/heter_cov_priors.R for shared parsing/validation/printing.
  # ---------------------------------------------------------------------------
  hcv_flags   <- .parse_heter_cov_flags(Prior, nvar)
  useHeterCov <- hcv_flags$useHeterCov
  useFullCov  <- hcv_flags$useFullCov
  auto_full   <- hcv_flags$auto_full
  var_params  <- list()    # placeholder forwarded to C++ when useHeterCov=FALSE
  phi_params  <- list()    # placeholder forwarded to C++ when useFullCov=FALSE
  .validate_heter_cov(useHeterCov, useBART, ncomp, drawdelta, SignRes)
  if(is.null(Prior$a)) { a=rep(BayesmConstant.a,ncomp)} else {a=Prior$a}
  if(length(a) != ncomp) {pandterm("Requires dim(a)= ncomp (no of components)")}
  bada=FALSE
  for(i in 1:ncomp) { if(a[i] < 0) bada=TRUE}
  if(bada) pandterm("invalid values in a vector")
  
  if(is.null(Prior$nu)&sum(abs(SignRes))>0) {nu = nvar+15}                                        
  if(is.null(Prior$Amu)&sum(abs(SignRes))>0) {Amu = matrix(.1)}                                   
  if(is.null(Prior$V)&sum(abs(SignRes))>0) {V = nu*(diag(nvar)-diag(abs(SignRes)>0)*.8)} 
  
  #
  # check on Mcmc
  #
  if(missing(Mcmc)) 
  {pandterm("Requires Mcmc list argument")}
  else 
  { 
    if(is.null(Mcmc$s)) {s=BayesmConstant.RRScaling/sqrt(nvar)} else {s=Mcmc$s}
    if(is.null(Mcmc$w)) {w=BayesmConstant.w}  else {w=Mcmc$w}
    if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {R=Mcmc$R}
    if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
    if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
  }
  #
  # print out problem
  #
  if (r_verbose) {
    cat(" ",fill=TRUE)
    cat("Starting MCMC Inference for Hierarchical Logit:",fill=TRUE)
    cat("   Normal Mixture with",ncomp,"components for first stage prior",fill=TRUE)
    cat(paste("  ",p," alternatives; ",nvar," variables in X"),fill=TRUE)
    cat(paste("   for ",nlgt," cross-sectional units"),fill=TRUE)
    cat(" ",fill=TRUE)
    cat("Prior Parms: ",fill=TRUE)
    cat("nu =",nu,fill=TRUE)
    cat("V ",fill=TRUE)
    print(V)
    cat("mubar ",fill=TRUE)
    print(mubar)
    cat("Amu ", fill=TRUE)
    print(Amu)
    cat("a ",fill=TRUE)
    print(a)
    if(drawdelta) {
      if (useBART) {
        cat("\nHART prior parameters:\n")
        cat(sprintf("num_trees = %d\n", bart_params$num_trees))
        cat(sprintf("power = %g\n", bart_params$power))
        cat(sprintf("base = %g\n", bart_params$base))
        cat(sprintf("tau = %g\n", bart_params$tau))
        cat(sprintf("numcut = %d\n", bart_params$numcut))
        cat(sprintf("sparse = %d\n", bart_params$sparse))
        if (bart_params$sparse) {
          cat(sprintf("theta = %g\n", bart_params$theta))
          cat(sprintf("omega = %g\n", bart_params$omega))
          cat(sprintf("a = %g\n", bart_params$a))
          cat(sprintf("b = %g\n", bart_params$b))
          cat(sprintf("rho = %g\n", bart_params$rho))
          cat(sprintf("aug = %g\n", bart_params$aug))
          cat(sprintf("burn = %d\n", bart_params$burn))
        }
      } else {
        cat("\nLinear regression prior parameters:\n")
        cat("deltabar",fill=TRUE)
        print(deltabar)
        cat("Ad",fill=TRUE)
        print(Ad)
      }
    }
    if(sum(abs(SignRes)) != 0){
      cat("Sign Restrictions Vector (0: unconstrained, 1: positive, -1: negative)",fill=TRUE)
      print(matrix(SignRes,ncol=1))
    }
    
    cat(" ",fill=TRUE)
    cat("MCMC Parms: ",fill=TRUE)
    cat(paste("s=",round(s,3)," w= ",w," R= ",R," keep= ",keep," nprint= ",nprint),fill=TRUE)
    cat("",fill=TRUE)
  } # End r_verbose block for initial setup printing
  
  oldbetas = matrix(double(nlgt * nvar), ncol = nvar)
  
  #--------------------------------------------------------------------------------------------------
  #
  #  create functions needed
  #
  llmnlFract=                                                           
    function(beta,y,X,betapooled,rootH,w,wgt,SignRes = rep(0,ncol(X))){ 
      z=as.vector(rootH%*%(beta-betapooled))
      return((1-w)*llmnl_con(beta,y,X,SignRes)+w*wgt*(-.5*(z%*%z)))
    }
  
  mnlHess_con=function (betastar, y, X, SignRes = rep(0,ncol(X))) {     
    
    #Reparameterize betastar to beta to allow for sign restrictions
    beta = betastar
    beta[SignRes!=0] = SignRes[SignRes!=0]*exp(betastar[SignRes!=0])
    
    n = length(y)
    j = nrow(X)/n
    k = ncol(X)
    Xbeta = X %*% beta
    Xbeta = matrix(Xbeta, byrow = T, ncol = j)
    Xbeta = exp(Xbeta)
    iota = c(rep(1, j))
    denom = Xbeta %*% iota
    Prob = Xbeta/as.vector(denom)
    Hess = matrix(double(k * k), ncol = k)
    for (i in 1:n) {
      p = as.vector(Prob[i, ])
      A = diag(p) - outer(p, p)
      Xt = X[(j * (i - 1) + 1):(j * i), ]
      Hess = Hess + crossprod(Xt, A) %*% Xt
    }
    # modify Hessian for reparameterization
    #  Hess above is the hessian in the constrained parms (beta)
    #  we must express obtain hessian in betastar (unconstrained parms)
    #  Hess_beta = J^t Hess_betastar J
    #  Hess_betastar = (J^-1)t Hess_beta J^-1
    #  J: jacobian from beta to betastar
    #  J^-1: jacobian from  betastar to beta  -- see notes
    lambda = c(rep(1,length(SignRes)))
    lambda[SignRes == 1]= beta[SignRes == 1]
    lambda[SignRes == -1]= - beta[SignRes == -1]
    Hess=Hess * crossprod(t(lambda))
    # hess[i,j] = hess[i,j]*lambda[i]*lambda[j]
    return(Hess)
  }
  
  #-------------------------------------------------------------------------------------------------------
  #
  # intialize compute quantities for Metropolis
  #
  if (r_verbose) {
    cat("initializing Metropolis candidate densities for ",nlgt," units ...",fill=TRUE)
    fsh()
  }
  #
  #  now go thru and computed fraction likelihood estimates and hessians
  #
  #       Lbar=log(pooled likelihood^(n_i/N))
  #
  #       fraction loglike = (1-w)*loglike_i + w*Lbar
  #
  betainit=c(rep(0,nvar))
  noRes=c(rep(0,nvar))
  # run unconstrainted opt first
  out=optim(betainit,llmnl_con,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-6), 
            X=Xpooled,y=ypooled,SignRes=noRes)
  
  betainit=out$par
  betainit[SignRes!=0] = 0  # set constrained terms to zero -- implies setting "beta" to either 1, -1
  #
  #  compute pooled optimum
  #
  # changed to default method - Nelder-Mead for more robust optimization- sometimes BFGS
  #  fails to find optimum using exponential reparameterization
  out=optim(betainit,llmnl_con,control=list( fnscale=-1,trace=0,reltol=1e-6), 
            X=Xpooled,y=ypooled,SignRes=SignRes)
  
  betapooled=out$par
  #
  # warn user if the constrained pooled model has unreasonably small/large coefficients
  #
  if(sum(abs(betapooled[as.logical(SignRes)])>10))
  {
    cat("In tuning Metropolis algorithm, constrained pooled parameter estimates contain very small/large values",
        fill=TRUE)
    print(cbind(betapooled,SignRes))
    cat("check any constrained values with absolute value > 10 above",fill=TRUE)
    cat("      - implies abs(beta) > exp(10) or abs(beta) < exp(-10)",fill=TRUE)
  }
  
  
  H=mnlHess_con(betapooled,ypooled,Xpooled,SignRes)                                         
  rootH=chol(H)
  for (i in 1:nlgt) 
  {
    wgt=length(lgtdata[[i]]$y)/length(ypooled)
    out=optim(betapooled,llmnlFract,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-4), 
              X=lgtdata[[i]]$X,y=lgtdata[[i]]$y,betapooled=betapooled,rootH=rootH,w=w,wgt=wgt,SignRes=SignRes)
    if(out$convergence == 0) { 
      hess=mnlHess_con(out$par,lgtdata[[i]]$y,lgtdata[[i]]$X,SignRes)                          
      lgtdata[[i]]=c(lgtdata[[i]],list(converge=1,betafmle=out$par,hess=hess)) }
    else
    { lgtdata[[i]]=c(lgtdata[[i]],list(converge=0,betafmle=c(rep(0,nvar)),
                                       hess=diag(nvar))) }
    oldbetas[i,]=lgtdata[[i]]$betafmle
    if(i%%50 ==0 && r_verbose) cat("  completed unit #",i,fill=TRUE)
    fsh()
  }
  #
  #  initialize values
  #
  # set initial values for the indicators
  #     ind is of length(nlgt) and indicates which mixture component this obs
  #     belongs to.
  #
  ind=NULL
  ninc=floor(nlgt/ncomp)
  for (i in 1:(ncomp-1)) {ind=c(ind,rep(i,ninc))}
  if(ncomp != 1) {ind = c(ind,rep(ncomp,nlgt-length(ind)))} else {ind=rep(1,nlgt)}
  #
  # initialize probs
  #
  oldprob=rep(1/ncomp,ncomp)
  #
  #initialize delta
  #
  if (drawdelta & useBART){
    olddelta = 0
    deltabar = 0
    Ad = matrix(0) 
  } else if (drawdelta) {
    olddelta = rep(0,nz*nvar)
    bart_params = list(0)
  } else { #send placeholders to the _loop function if there is no Z matrix
    olddelta = 0
    Z = matrix(0)
    deltabar = 0
    Ad = matrix(0)
    bart_params = list(0)
  }
  # ---------------------------------------------------------------------------
  # Heter-cov hyperparameter setup.  oldbetas (per-unit MLE betas) are now
  # available, so lambda for the variance-tree chi^{-2} prior can be
  # data-calibrated if the user did not supply one.  All parsing lives in
  # R/heter_cov_priors.R for reuse by Negbin and Linear.
  # ---------------------------------------------------------------------------
  if (useHeterCov) {
    var_params <- .parse_var_params(Prior$vartree, oldbetas, nz)
    if (useFullCov) phi_params <- .parse_phi_params(Prior$phitree, nz)
    if (r_verbose)
      .print_heter_cov_summary(useFullCov, var_params, phi_params,
                               lambda_user_supplied = !is.null(Prior$vartree$lambda),
                               auto_full = auto_full)
  }

  ###################################################################
  # Wayne Taylor
  # 09/22/2014
  ###################################################################
  draws =  rhierMnlRwMixture_rcpp_loop(lgtdata, Z,
                                       deltabar, Ad, mubar, Amu,
                                       nu, V, s,
                                       R, keep, nprint, drawdelta,
                                       as.matrix(olddelta), a, oldprob, oldbetas, ind, SignRes,
                                       useBART, bart_params,
                                       useHeterCov, var_params, phi_params)
  ####################################################################

  if (useHeterCov) {
    attributes(draws$betadraw)$class <- c("bayesm.hcoef")
    attributes(draws$mu_draw)$class  <- c("bayesm.mat", "mcmc")
    attributes(draws$mu_draw)$mcpar  <- c(1, R, keep)
    # "bayesm.HART.HeterCov" is the marker superclass shared by all
    # heteroscedastic-Sigma(Z) hierarchical models (MNL, Negbin, Linear); it
    # lets predict_helpers.R dispatch generically without hardcoding any
    # model-specific class name.
    class(draws) <- c("rhierMnlRwMixtureHeterCov",
                      "bayesm.HART.HeterCov",
                      "rhierMnlRwMixture")
  } else {
    if(!useBART & drawdelta){
      attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
      attributes(draws$Deltadraw)$mcpar=c(1,R,keep)}
    attributes(draws$betadraw)$class=c("bayesm.hcoef")
    attributes(draws$nmix)$class="bayesm.nmix"
    class(draws) <- "rhierMnlRwMixture"
  }

  return(draws)
}

# Methods ======================================================================

# --- Shared Helper Functions (in R/predict_helpers.R) ---
# .validate_predict_inputs()  - Input validation for predict methods
# .calculate_delta_z()        - Compute Delta(Z) via Deltadraw or pwbart
# .add_mu_component()         - Add mixture component mean mu

# --- MNL-Specific Helper Functions ---

# Helper 4: Predict Posterior Probabilities
.predict_posterior_probs <- function(object, newdata, kept_draws_indices, 
                                     r_verbose) {
  # Minimal validation assumes object$betadraw is valid from main validation
  dims_beta <- dim(object$betadraw)
  nlgt <- dims_beta[1]
  nvar <- dims_beta[2]
  ndraws_out <- length(kept_draws_indices)
  
  # Minimal newdata validation
  if (is.null(newdata$p)) stop("newdata$p must be provided.")
  p <- newdata$p
  if (is.null(newdata$nlgtdata)) stop("newdata$nlgtdata must be provided.")
  
  # Initialize Probability List
  prob_pred_list <- vector("list", nlgt)
  names(prob_pred_list) <- paste0("Unit_", 1:nlgt)
  
  betadraw <- object$betadraw # Alias for readability
  
  # Calculate Probabilities
  for (i in 1:nlgt) { # Loop through original units
    if (r_verbose) cat("Calculating posterior probabilities for Unit", i, "of", 
                       nlgt, "...\n")
    
    # Minimal check on this unit's X (assuming it exists and is matrix)
    X_i <- newdata$nlgtdata[[i]]$X
    # Assume dimensions are correct based on prior checks/usage
    
    T_i <- nrow(X_i) / p
    # Initialize 3D array for this unit's probabilities
    prob_array_i <- array(0.0, dim = c(T_i, p, ndraws_out))
    
    draw_count <- 0
    for (s_orig in kept_draws_indices) { # Loop through kept draws indices
      draw_count <- draw_count + 1
      beta_is <- betadraw[i, , s_orig, drop = TRUE] # Stored posterior beta
      
      # Call helper function from utils.R (using :::)
      # Assume bayesm.HART namespace is available or ::: is appropriate
      probs_mat_is <- calculate_mnl_probs_from_beta(X_i, beta_is, p)
      
      # Assign the T_i x p matrix to the s-th slice of the 3D array
      prob_array_i[, , draw_count] <- probs_mat_is
    }#FOR s_orig
    prob_pred_list[[i]] <- prob_array_i # Assign 3D array to unit i
  }#FOR i
  if (r_verbose) cat("Posterior probability calculation complete.\n")
  
  return(prob_pred_list)
}#PREDICT_POSTERIOR_PROBS

# Helper 5: Predict Prior Probabilities
#
# For each prediction unit i and kept draw s, samples nsim values of
#   beta_i = DeltaZ_is + eta_k
# where eta_k is drawn from the mixture / heter-cov heterogeneity prior, and
# averages choice probabilities over the nsim draws.  Two paths:
#
#   Mixture (BART or linear): eta_k ~ mixture(pvec_s, compdraw[[s]]).
#   Heter-cov:                eta_k = mu_s + Sigma(Z*_i)^{1/2} z_k,
#                             z_k ~ N(0, I).
#
# The heter-cov path consumes precomputed `hetercov_comps` (built once by
# predict.rhierMnlRwMixture); when called outside that dispatcher, it
# re-evaluates the relevant trees on the supplied newdata$Z.
.predict_prior_probs <- function(object, newdata, delta_z_array,
                                 kept_draws_indices_nmix, nsim, r_verbose,
                                 hetercov_comps = NULL) {
  is_heter <- inherits(object, "bayesm.HART.HeterCov")
  if (!is_heter && !requireNamespace("bayesm", quietly = TRUE)) {
    stop("The 'bayesm' package is required for type='prior_probs'. Please install it.")
  }

  dims_beta    <- dim(object$betadraw)
  nvar         <- dims_beta[2]
  ndraws_total <- dims_beta[3]
  npred        <- dim(delta_z_array)[1]
  ndraws_out   <- dim(delta_z_array)[3]

  if (is.null(newdata$p)) stop("newdata$p must be provided.")
  p <- newdata$p
  if (is.null(newdata$X)) stop("newdata$X must be provided as a list.")

  if (is_heter) {
    if (is.null(newdata$Z))
      stop("Heter-cov prior_probs requires newdata$Z (same Z as DeltaZ).")
    if (npred != nrow(newdata$Z))
      stop("npred (from delta_z_array) does not match nrow(newdata$Z).")
    if (is.null(object$mu_draw))
      stop("Heter-cov object missing mu_draw.")
    if (is.null(hetercov_comps))
      hetercov_comps <- .hetercov_components(object, newdata$Z, npred, nvar,
                                             ndraws_total, r_verbose)
    mu_kept <- object$mu_draw[kept_draws_indices_nmix, , drop = FALSE]  # ndraws_out x nvar
  } else {
    .assert_single_component_nmix(object, "prior_probs prediction")
  }

  prob_pred_list        <- vector("list", npred)
  names(prob_pred_list) <- paste0("PredUnit_", seq_len(npred))

  for (i in seq_len(npred)) {
    if (r_verbose) cat("Calculating prior probabilities for Prediction Unit", i,
                       "of", npred, "...\n")
    X_i <- newdata$X[[i]]
    T_i <- nrow(X_i) / p
    prob_array_i <- array(0.0, dim = c(T_i, p, ndraws_out))

    for (draw_idx in seq_len(ndraws_out)) {
      DeltaZ_is <- delta_z_array[i, , draw_idx, drop = TRUE]
      s_orig    <- kept_draws_indices_nmix[draw_idx]

      if (is_heter) {
        # eta = mu_s + Sigma(Z*_i)^{1/2} z, with z ~ N(0, I), nsim copies.
        z_mat <- matrix(rnorm(nsim * nvar), nrow = nsim, ncol = nvar)
        d_is  <- hetercov_comps$d_arr[i, , s_orig]
        x_mat <- z_mat * matrix(sqrt(d_is), nrow = nsim, ncol = nvar,
                                byrow = TRUE)
        if (hetercov_comps$use_full) {
          phi_is <- hetercov_comps$phi_arr[i, , , s_orig]   # nvar x nvar lower
          for (k in seq_len(nsim)) {
            x <- x_mat[k, ]
            for (j in 2:nvar) {
              x[j] <- x[j] + sum(phi_is[j, seq_len(j - 1)] * x[seq_len(j - 1)])
            }
            x_mat[k, ] <- x
          }
        }
        mu_s    <- mu_kept[draw_idx, ]
        eta_mat <- x_mat + matrix(mu_s, nrow = nsim, ncol = nvar, byrow = TRUE)

        prob_is_sum <- matrix(0.0, nrow = T_i, ncol = p)
        for (k in seq_len(nsim)) {
          beta_is_k    <- DeltaZ_is + eta_mat[k, ]
          probs_mat_is <-
            calculate_mnl_probs_from_beta(X_i, beta_is_k, p)
          prob_is_sum  <- prob_is_sum + probs_mat_is
        }
      } else {
        comps_s <- object$nmix$compdraw[[s_orig]]
        pvec_s  <- object$nmix$probdraw[s_orig, ]
        prob_is_sum <- matrix(0.0, nrow = T_i, ncol = p)
        for (k in seq_len(nsim)) {
          eta_is_k     <- as.vector(bayesm::rmixture(n = 1, pvec = pvec_s,
                                                     comps = comps_s)$x)
          beta_is_k    <- DeltaZ_is + eta_is_k
          probs_mat_is <-
            calculate_mnl_probs_from_beta(X_i, beta_is_k, p)
          prob_is_sum  <- prob_is_sum + probs_mat_is
        }
      }

      prob_array_i[, , draw_idx] <- prob_is_sum / nsim
    }
    prob_pred_list[[i]] <- prob_array_i
  }
  if (r_verbose) cat("Prior probability calculation complete.\n")

  return(prob_pred_list)
}#PREDICT_PRIOR_PROBS


# --- Refactored predict.rhierMnlRwMixture Method ---

#' Predict Method for rhierMnlRwMixture Objects
#' @param object A fitted rhierMnlRwMixture object.
#' @param newdata Optional list containing data for prediction. Structure depends
#'   on `type`:
#'   - For `type %in% c("DeltaZ", "DeltaZ+mu", "SigmaZ")`: Requires `newdata$Z`, a matrix
#'     with `npred` rows for prediction units (if model was fit with Z).
#'   - For `type = "posterior_probs"`: Requires `newdata$nlgtdata`, a list of
#'     length `nlgt` (original number of units). Each element `\\[[i]]` must
#'     contain `$X`, the design matrix `(T_i*p) x nvar` for unit `i`. Also
#'     requires `newdata$p`, the number of alternatives.
#'   - For `type = "prior_probs"`: Requires `newdata$Z` (if model fit with Z,
#'     determining `npred`), `newdata$p`, and `newdata$X` (a list of length
#'     `npred`, each element `\\[[i]]` having the design matrix `(T_i*p) x nvar`).
#' @param type Type of prediction:
#'   - `"DeltaZ"`: Expected part-worths of the representative consumer, \eqn{\Delta(Z)}.
#'   - `"DeltaZ+mu"`: Expected part-worths plus the mean of the unobserved heterogeneity component, \eqn{\Delta(Z) + \mu_1}. The package supports only `ncomp = 1`.
#'   - `"posterior_probs"`: Posterior predictive choice probabilities for the original
#'     estimation units using stored `betadraw`.
#'   - `"prior_probs"`: Prior predictive choice probabilities for new prediction units
#'     (based on `newdata$Z` or the overall mixture if no Z was used). Probabilities
#'     are averaged over `nsim` draws from the heterogeneity mixture distribution
#'     per posterior draw.
#'   - `"SigmaZ"`: Draws of the heteroscedastic covariance matrix \eqn{\Sigma(Z)}.
#'     Available only for models fit with `Prior$vartree` (class marker
#'     `"bayesm.HART.HeterCov"`).
#' @param burn Integer, number of initial MCMC draws to discard.
#' @param nsim Integer, number of draws from the heterogeneity mixture distribution
#'   per posterior draw for `type = "prior_probs"`.
#' @param r_verbose Logical, print progress updates?
#' @param ... Additional arguments passed to underlying prediction functions
#'   (e.g., `mc.cores`, `verbose` for BART `DeltaZ` predictions via `pwbart`).
#' @return Depends on `type`:
#'   - For `type %in% c("DeltaZ", "DeltaZ+mu")`: 3D array `[npred, nvar, ndraws_out]`
#'     of predicted expected part-worths.
#'   - For `type = "SigmaZ"`: 4D array `[npred, nvar, nvar, ndraws_out]` of
#'     covariance draws at each prediction unit.
#'   - For `type = "posterior_probs"`: List of length `nlgt`. Each element `\\[[i]]`
#'     is a 3D array `[T_i, p, ndraws_out]` of posterior predictive choice probabilities
#'     for unit `i`.
#'   - For `type = "prior_probs"`: List of length `npred`. Each element `\\[[i]]`
#'     is a 3D array `[T_i, p, ndraws_out]` of prior predictive choice
#'     probabilities for prediction unit `i`.
#' @keywords internal
#' @export
predict.rhierMnlRwMixture <- function(object, newdata = NULL,
                                      type = "DeltaZ+mu", burn = 0, nsim = 10,
                                      r_verbose = TRUE, ...) {
  mnl_types    <- c("DeltaZ", "DeltaZ+mu", "posterior_probs", "prior_probs", "SigmaZ")
  ndraws_total <- .validate_predict_inputs(object, newdata, type, burn, nsim,
                                           valid_types = mnl_types)
  kept_draws_indices <- if (burn > 0) (burn + 1):ndraws_total else 1:ndraws_total

  if (type == "posterior_probs")
    return(.predict_posterior_probs(object, newdata, kept_draws_indices,
                                    r_verbose))

  # Heter-cov: evaluate every BART / varBART / phi-BART ensemble at newdata$Z
  # exactly ONCE, then share the result with both .calculate_delta_z and
  # (when type == "prior_probs") .predict_prior_probs / .calculate_sigma_z.
  hetercov_comps <- NULL
  if (inherits(object, "bayesm.HART.HeterCov") &&
      type %in% c("DeltaZ", "DeltaZ+mu", "prior_probs", "SigmaZ")) {
    if (is.null(newdata$Z))
      stop("Heter-cov predictions require newdata$Z.")
    nvar  <- dim(object$betadraw)[2]
    npred <- nrow(newdata$Z)
    hetercov_comps <- .hetercov_components(object, newdata$Z, npred, nvar,
                                           ndraws_total, r_verbose)
  }
  if (type == "SigmaZ") {
    return(.calculate_sigma_z(object, newdata, burn, r_verbose,
                              hetercov_comps = hetercov_comps))
  }

  delta_z <- .calculate_delta_z(object, newdata, burn, r_verbose,
                                hetercov_comps = hetercov_comps, ...)
  if (type == "DeltaZ")    return(delta_z)
  if (type == "DeltaZ+mu") return(.add_mu_component(delta_z, object,
                                                    kept_draws_indices))

  # type == "prior_probs"
  .predict_prior_probs(object, newdata, delta_z, kept_draws_indices, nsim,
                       r_verbose, hetercov_comps = hetercov_comps)
}#PREDICT.RHIERMNLRWMIXTURE