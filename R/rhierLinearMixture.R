#' Bayesian Hierarchical Linear Model with HART Prior
#' @description
#' `rhierLinearMixture` implements an MCMC algorithm for a Bayesian hierarchical
#' linear regression model with a Hierarchical Additive Regression Trees (HART) prior. HART 
#' is a hierarchical nonparametric prior that allows for flexible modeling of the 
#' representative unit as a function of potentially many observed characteristics.
#' \code{Prior$bart} allows for modeling the conditional mean of the normal prior.
#' 
#' @note Currently, only \code{ncomp = 1} is supported.
#'
#' @param Data A list containing:
#'   - `regdata`: A list of length `nreg`. Each element `regdata[[i]]` must be a list with:
#'     - `y`: `n_i x 1` vector of responses.
#'     - `X`: `n_i x nvar` design matrix.
#'   - `Z` (optional): `nreg x nz` matrix of unit-level covariates. If omitted,
#'     `drawdelta` is set to FALSE and no second-stage covariates are used.
#' @param Prior A list containing prior parameters:
#'   - `ncomp` (required): Number of components. Must be 1.
#'   - `nu.e` (optional): Degrees of freedom for error variance prior (default: 3).
#'   - `ssq` (optional): Scale for error variance prior, `nreg x 1` vector (default: var(y_i)).
#'   - `mubar` (optional): `1 x nvar` prior mean vector for the normal prior (default: 0).
#'   - `Amu` (optional): Prior precision for the normal prior on \eqn{\mu} (default: 0.01). Scalar on the bayesm paths. Under a HART model the prior is the independent \eqn{N(\bar\mu, Amu^{-1})}, and `Amu` may be a scalar (expanded to `Amu * I`) or an `nvar x nvar` SPD matrix.
#'   - `nu` (optional): Degrees of freedom for IW prior on the normal prior covariance (default: nvar + 3).
#'   - `V` (optional): `nvar x nvar` location matrix for IW prior (default: nu * I).
#'   - `Ad` (optional): Prior precision for vec(Delta) (default: 0.01 * I). Only used when `drawdelta = TRUE` and `useBART = FALSE`.
#'   - `deltabar` (optional): Prior mean for vec(Delta) (default: 0). Only used when `drawdelta = TRUE` and `useBART = FALSE`.
#'   - `hart_model` (optional): Preference-model selector, one of `"none"`
#'     (default; bayesm constant or linear \eqn{\Delta(Z)}), `"normalized"`
#'     (HART prior on the representative unit), or `"direct"` (HART mean and
#'     covariance trees). Both HART models require `Z`. The tree lists below
#'     are optional under a HART model and are an error under `"none"`. See
#'     `rhierMnlRwMixture` for the model definitions.
#'   - `bart` (optional): List of hyperparameters for the HART mean trees
#'     (`hart_model = "normalized"` or `"direct"`). See Details.
#'   - `vartree` (optional, `hart_model = "direct"` only): List of
#'     hyperparameters for the product-of-trees variance models on the Modified
#'     Cholesky diagonal \eqn{d_j(\cdot)}. See "Heteroscedastic Covariance"
#'     in Details.
#'   - `phitree` (optional, `hart_model = "direct"` and `nvar > 1` only): List
#'     of hyperparameters for the sum-of-trees models on the Modified Cholesky
#'     off-diagonals \eqn{\phi_{jk}(\cdot)}.
#' @param Mcmc A list containing MCMC parameters:
#'   - `R`: Number of MCMC iterations (required).
#'   - `keep` (optional): Thinning parameter (default: 1).
#'   - `nprint` (optional): Print progress every `nprint` draws (default: 100, 0 for none).
#' @param r_verbose Logical. Print startup messages? Default TRUE.
#'
#' @details
#' ## Model Specification
#' \eqn{y_i = X_i \beta_i + \epsilon_i}, \eqn{\epsilon_i \sim N(0, \tau_i I)}
#'
#' The unit-level coefficients are modeled as:
#'   \eqn{\beta_i = Z_i \Delta + u_i} \cr
#'   \eqn{u_i \sim N(0, \Sigma_i)} \cr
#'
#' ## HART Prior Details
#' If `Prior$bart` is a list, \eqn{\Delta(Z_i)} is modeled via a sum-of-trees (BART) prior
#' instead of a linear hierarchical specification. See `rhierMnlRwMixture` for complete 
#' details and hyperparameter definitions.
#'
#' ## Heteroscedastic Covariance \eqn{\Sigma(Z_i)} (experimental extension)
#' When `Prior$vartree` is supplied, the homoscedastic prior covariance \eqn{\Sigma} is
#' replaced with a unit-specific modified Cholesky decomposition where the components
#' are modeled as tree ensembles. See `rhierMnlRwMixture` for complete details.
#'
#' When this extension is active, the returned object additionally inherits
#' classes `"rhierLinearMixtureHeterCov"` and `"bayesm.HART.HeterCov"`. `predict()` 
#' dispatches on these classes to evaluate \eqn{\Sigma(Z^*)} at any new \eqn{Z^*}.
#'
#' @return A list of class `"rhierLinearMixture"` containing:
#'   - `betadraw`: `nreg x nvar x (R/keep)` array of unit-level beta draws.
#'   - `taudraw`: `(R/keep) x nreg` matrix of error variance draws.
#'   - `loglike`: `(R/keep) x 1` vector of log-likelihoods.
#'   - `nmix`: Legacy list containing prior draws (omitted under heter-cov).
#'   - If `drawdelta` and non-BART: `Deltadraw`.
#'   - If BART: `bart_models`, `varcount`, `varprob`.
#'   - If heter-cov (additional class `"rhierLinearMixtureHeterCov"`):
#'     `var_models`, `phi_models` (jagged or `NULL`), `mu_draw`,
#'     `mean_scale_draw`, `var_varcount`, `var_varprob`.
#'
#' @seealso [predict.rhierLinearMixture()], [marginal_effects.rhierLinearMixture()]
#'
#' @author Peter Rossi, Wayne Taylor (original bayesm code), Thomas Wiemann (HART modifications).
#'
#' @examples
#' \donttest{
#' set.seed(20260513)
#' sim <- bayesm.HART::sim_hier_linear(
#'   nreg = 30, nobs = 12, nvar = 1, nz = 2,
#'   const = TRUE, het_observed = "linear",
#'   target_var_betabar = 1.0, target_var_eps = 0.5,
#'   sigma_sq = 0.5)
#' Prior <- list(ncomp = 1L,
#'               bart    = list(num_trees = 20),
#'               vartree = list(num_trees = 20))
#' Mcmc  <- list(R = 200L, keep = 1L, nprint = 0L)
#' fit   <- bayesm.HART::rhierLinearMixture(
#'   Data = list(regdata = sim$regdata, Z = sim$Z),
#'   Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
#' str(fit, max.level = 1)
#' }
#'
#' @export
rhierLinearMixture=
function(Data,Prior,Mcmc, r_verbose = TRUE)
{
#
# Revision History
#     P. Rossi / W. Taylor
#     6/25 T. Wiemann - added HART prior + roxygen2 docs
#
# Purpose:
#   run hiearchical linear regression model
#
# Arguments:
#   Data list of regdata,Z
#     regdata is a list of lists each list with members y, X
#        e.g. regdata[[i]]=list(y=y,X=X)
#     X has nvar columns
#     Z is nreg=length(regdata) x nz
#   Prior list of prior hyperparameters
#     ncomp, nu.e, ssq, mubar, Amu, nu, V, a, Ad, deltabar
#   Mcmc
#     list of Mcmc parameters
#     R is number of draws
#     keep is thining parameter -- keep every keepth draw
#     nprint - print estimated time remaining on every nprint'th draw
#
# Output:
#   list of
#   betadraw -- nreg x nvar x R/keep array of individual regression betas
#   taudraw -- R/keep x nreg  array of error variances for each regression
#   Deltadraw -- R/keep x nz x nvar array of Delta draws (if drawdelta)
#   nmix -- list with probdraw, zdraw, compdraw
#   loglike -- R/keep x 1 vector of log-likelihoods
#
# Model:
# nreg regression equations
#        y_i = X_ibeta_i + epsilon_i
#        epsilon_i ~ N(0,tau_i)
#             nvar X vars in each equation
#
# Priors:
#        tau_i ~ nu.e*ssq_i/chisq(nu.e)
#        beta_i ~ N(mu_j + Delta*z_i, Sigma_j)  j = ind[i]
#        normal prior on beta heterogeneity
#
#  create needed functions
#
#------------------------------------------------------------------------------
append=function(l) { l=c(l,list(XpX=crossprod(l$X),Xpy=crossprod(l$X,l$y)))}
#
getvar=function(l) {
     v=var(l$y)
     if(is.na(v)) return(1)
     if(v>0) return (v) else return (1)}
#
#------------------------------------------------------------------------------

#
# check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata, and (possibly) Z")}
if(is.null(Data$regdata)) {pandterm("Requires Data element regdata (list of data for each unit)")}
regdata=Data$regdata
nreg=length(regdata)
drawdelta=TRUE
if(is.null(Data$Z)) {
  if(r_verbose) cat("Z not specified",fill=TRUE)
  fsh()
  drawdelta=FALSE
}
else {
  if (!is.matrix(Data$Z)) {pandterm("Z must be a matrix")}
  else {
    if (nrow(Data$Z) != nreg) {pandterm(paste("Nrow(Z) ",nrow(Data$Z),"ne number regressions ",nreg))}
    else {Z=Data$Z}
  }
}
if(drawdelta) {
  nz=ncol(Z)
  colmeans=apply(Z,2,mean)
  if(sum(colmeans) > 1e-05) {
    if (r_verbose) cat(paste("Warning: Z does not appear to be de-meaned: colmeans= ",
                              paste(round(colmeans, 4), collapse=", ")),fill=TRUE)
  }
}
#
# check data for validity
#
for (i in 1:nreg) {
  if(!is.matrix(regdata[[i]]$X)) {pandterm(paste0("regdata[[",i,"]]$X must be a matrix"))}
  if(!is.vector(regdata[[i]]$y, mode = "numeric") & !is.vector(regdata[[i]]$y, mode = "logical") & !is.matrix(regdata[[i]]$y))
    {pandterm(paste0("regdata[[",i,"]]$y must be a numeric or logical vector or matrix"))}
  if(is.matrix(regdata[[i]]$y)){ if(ncol(regdata[[i]]$y)>1) {pandterm(paste0("regdata[[",i,"]]$y must be a vector or one-column matrix"))}}
}

dimfun=function(l) {c(length(l$y),dim(l$X))}
dims=sapply(regdata,dimfun)
dims=t(dims)
nvar=quantile(dims[,3],prob=.5)

for (i in 1:nreg)
{
   if(dims[i,1] != dims[i,2]  || dims[i,3] !=nvar)
      {pandterm(paste("Bad Data dimensions for unit ",i," dims(y,X) =",dims[i,]))}
}
#
# check for Prior
#
if(missing(Prior))
   { pandterm("Requires Prior list argument (at least ncomp)")}
if(is.null(Prior$nu.e)) {nu.e=BayesmConstant.nu.e}
   else {nu.e=Prior$nu.e}
if(is.null(Prior$ssq)) {ssq=sapply(regdata,getvar)}
   else {ssq=Prior$ssq}
if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp")}
   else {ncomp=Prior$ncomp}
if(ncomp != 1) {pandterm("Only ncomp = 1 is currently supported")}
# HART preference model (R/hart_priors.R): "none" keeps the bayesm paths.
cfg <- .hart_config(Prior, nvar = nvar, nz = if (drawdelta) nz else 0L,
                    drawdelta = drawdelta,
                    coef_names = colnames(regdata[[1]]$X))
if(is.null(Prior$mubar)) { mubar=matrix(rep(0,nvar),nrow=1) }
   else { mubar=matrix(Prior$mubar, nrow=1) }
if(ncol(mubar) != nvar) {pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ",ncol(mubar)))}
if (cfg$model == "none") {
  if(is.null(Prior$Amu)) { Amu=matrix(BayesmConstant.A,ncol=1) }
     else { Amu=matrix(Prior$Amu, ncol=1) }
  if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Am must be a 1 x 1 array")}
} else {
  # HART: independent prior mu ~ N(mubar, Amu^{-1}); scalar or nvar x nvar.
  Amu <- .hart_Amu(Prior$Amu, nvar, default = BayesmConstant.A)
}
cfg$Amu <- Amu
if(is.null(Prior$nu)) {nu=nvar+BayesmConstant.nuInc}
   else {nu=Prior$nu}
if(nu < 1) {pandterm("invalid nu value")}
if(is.null(Prior$V)) {V=nu*diag(nvar)}
   else {V=Prior$V}
if(sum(dim(V)==c(nvar,nvar)) != 2) pandterm("Invalid V in prior")

# ---- HART preference model (see R/hart_priors.R) ----
# Prior$hart_model selects "none" (bayesm paths), "normalized", or "direct".
useBART     <- cfg$model != "none"
useHeterCov <- cfg$model == "direct"
bart_params <- cfg$bart
var_params  <- cfg$vartree %||% list()
phi_params  <- cfg$phitree %||% list()
# ---- end HART preference model ----

if(is.null(Prior$Ad) & drawdelta) {Ad=BayesmConstant.A*diag(nvar*nz)}
   else {Ad=Prior$Ad}
if(drawdelta) {
  if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}
}
if(is.null(Prior$deltabar) & drawdelta) {deltabar=rep(0,nz*nvar)}
   else {deltabar=Prior$deltabar}
if(drawdelta) {
  if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}
}
if(is.null(Prior$a)) { a=rep(BayesmConstant.a,ncomp)} else {a=Prior$a}
if(length(a) != ncomp) {pandterm("Requires dim(a)= ncomp (no of components)")}
bada=FALSE
for (i in 1:ncomp) { if(a[i] < 0) bada=TRUE }
if(bada) pandterm("invalid values in a vector")

#
# check MCMC argument
#
if(missing(Mcmc)) {pandterm("Requires Mcmc list argument")}
else
   {
    if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {R=Mcmc$R}
    if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
    if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
   }
#
# print out problem
#
if (r_verbose) {
cat(" ", fill=TRUE)
cat("Starting MCMC Inference for Hierarchical Linear Model:",fill=TRUE)
cat("   Normal Mixture with",ncomp,"components for first stage prior",fill=TRUE)
cat(paste("   for ",nreg," cross-sectional units"),fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("nu.e =",nu.e,fill=TRUE)
cat("nu =",nu,fill=TRUE)
cat("V ",fill=TRUE)
print(V)
cat("mubar ",fill=TRUE)
print(mubar)
cat("Amu ",fill=TRUE)
print(Amu)
cat("a ",fill=TRUE)
print(a)
if (useBART) {
  .print_hart_prior(cfg)
  if (useHeterCov) .print_hart_covariance(cfg)
} else if(drawdelta) {
  cat("deltabar",fill=TRUE)
  print(deltabar)
  cat("Ad",fill=TRUE)
  print(Ad)
}
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat("",fill=TRUE)
} # end r_verbose

#
# initialize and call C++ loop
#
regdata=lapply(regdata,append)
tau=sapply(regdata,getvar)
ind=NULL
ninc=floor(nreg/ncomp)
for (i in 1:(ncomp-1)) {ind=c(ind,rep(i,ninc))}
if(ncomp != 1) {ind = c(ind,rep(ncomp,nreg-length(ind)))} else {ind=rep(1,nreg)}

# ---- HART initialization ----
if (drawdelta & useBART){
  olddelta = rep(0, 1)
  deltabar = 0
  Ad = matrix(0)
} else if(drawdelta) {
  olddelta=rep(0,nz*nvar)
} else {
  olddelta=0
  Z=matrix(0)
  deltabar=0
  Ad=matrix(0)
  bart_params = list(0)
}
if (!useBART) {
  bart_params = list(0)
}
# ---- end HART initialization ----

if (drawdelta && useBART) {
  z_map <- .build_unique_z_map(Z)
  bart_params$z_index <- z_map$z_index
  bart_params$n_unique <- nrow(z_map$Z_unique)
} else {
  z_map <- NULL
}

oldprob=rep(1/ncomp,ncomp)

# Per-unit OLS slopes as the warm start under the direct model; zero for a
# unit with singular XpX.
Beta_init <- matrix(0, nrow = nreg, ncol = nvar)
if (useHeterCov) {
  Beta_init <- t(vapply(regdata, function(rd)
    tryCatch(as.numeric(solve(rd$XpX, rd$Xpy)),
             error = function(e) rep(0, nvar)),
    numeric(nvar)))
}

###################################################################
# Wayne Taylor 10/02/2014
# Modified by Thomas Wiemann 2025 -- added HART support,
#   and heter-cov (Sigma(Z_i)) extension.
###################################################################
draws = rhierLinearMixture_rcpp_loop(regdata, Z, deltabar,
    Ad, mubar, Amu, nu, V, nu.e, ssq, R, keep, nprint, drawdelta,
    as.matrix(olddelta), a, oldprob, ind, tau,
    cfg$model, bart_params,
    var_params, phi_params, Beta_init)
###################################################################

#
# format output
#
attributes(draws$taudraw)$class=c("bayesm.mat","mcmc")
attributes(draws$taudraw)$mcpar=c(1,R,keep)
draws <- .finish_hart_fit(draws, cfg, z_map, R, keep,
                          base_class = "rhierLinearMixture")

return(draws)
}

# ==============================================================================
# Predict Method
# ==============================================================================

#' Predict Method for rhierLinearMixture Objects
#'
#' Computes coefficient-level or predictive response draws.
#'
#' @param object A fitted `rhierLinearMixture` object.
#' @param newdata A list whose required fields depend on `mode`:
#'   - For `mode = "coefficients"`: requires `newdata$Z` (if the model was fit with Z).
#'   - For `mode = "posterior"`: requires `newdata$regdata`, a list with per-unit design
#'     matrices `X` matching the fitted units.
#'   - For `mode = "prior"`: requires `newdata$X`, a list of design matrices for prediction
#'     units, and `newdata$Z` when the model was fit with Z.
#' @param type Character; interpretation depends on `mode`:
#'   - For `mode = "coefficients"`: one of `"DeltaZ"`, `"DeltaZ+mu"`, `"SigmaZ"`.
#'   - For `mode %in% c("prior","posterior")`: must be `"response"`.
#' @param mode Prediction mode: `"coefficients"` (default), `"posterior"`, or `"prior"`.
#' @param burn Integer, number of initial MCMC draws to discard.
#' @param nsim Integer, number of Monte Carlo draws for prior predictive response mixing.
#' @param r_verbose Logical, print progress updates?
#' @param ... Additional arguments passed to `pwbart` for BART models.
#'
#' @return Depends on `mode` / `type`:
#'   - For `mode = "coefficients"`, `type %in% c("DeltaZ", "DeltaZ+mu")`: 3D array
#'     `[npred, ncoef, ndraws_out]` of predicted betabar values.
#'   - For `mode = "coefficients"`, `type = "SigmaZ"`: 4D array `[npred, ncoef, ncoef, ndraws_out]` of
#'     covariance draws at each prediction unit.
#'   - For `mode %in% c("prior","posterior")`, `type = "response"`: list of matrices
#'     with simulated response draws per unit (`nobs_i x ndraws_out`).
#' @export
predict.rhierLinearMixture <- function(object, newdata = NULL,
                                        type = "DeltaZ+mu", burn = 0,
                                        mode = "coefficients", nsim = 10,
                                        r_verbose = TRUE,
                                        force_tree_eval = FALSE, ...) {
  mode <- .resolve_predict_mode(mode)
  if (mode %in% c("prior", "posterior")) {
    if (type %in% c("DeltaZ", "DeltaZ+mu", "SigmaZ")) {
      type <- "response"
    }
    if (!identical(type, "response")) {
      stop("For rhierLinearMixture predictive modes, use type='response'.")
    }
    return(.predictive_dispatch(
      object, newdata, mode = mode, type = type,
      burn = burn, nsim = nsim, r_verbose = r_verbose,
      force_tree_eval = force_tree_eval, ...
    ))
  }

  valid_types <- c("DeltaZ", "DeltaZ+mu", "SigmaZ")
  .validate_predict_inputs(object, newdata, type, burn, nsim = 1,
                           valid_types = valid_types)
  .predict_structural_common(
    object, newdata, type, burn, r_verbose,
    force_tree_eval = force_tree_eval, ...
  )
}
