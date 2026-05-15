#' Bayesian Hierarchical Linear Model with Normal Mixture and HART Prior
#' @description
#' `rhierLinearMixture` implements a Gibbs sampler for a Bayesian hierarchical
#' linear regression model with a mixture of normals prior on the unit-level
#' coefficients. Supports both standard linear hierarchical prior and HART
#' (sum-of-trees) prior. Structurally symmetric to `rhierMnlRwMixture`.
#'
#' @param Data A list containing:
#'   - `regdata`: A list of length `nreg`. Each element `regdata[[i]]` must be a list with:
#'     - `y`: `n_i x 1` vector of responses.
#'     - `X`: `n_i x nvar` design matrix.
#'   - `Z` (optional): `nreg x nz` matrix of unit-level covariates. If omitted,
#'     `drawdelta` is set to FALSE and no second-stage covariates are used.
#' @param Prior A list containing prior parameters:
#'   - `ncomp` (required): Number of normal mixture components.
#'   - `nu.e` (optional): Degrees of freedom for error variance prior (default: 3).
#'   - `ssq` (optional): Scale for error variance prior, `nreg x 1` vector (default: var(y_i)).
#'   - `mubar` (optional): `1 x nvar` prior mean for mixture component means (default: 0).
#'   - `Amu` (optional): `1 x 1` prior precision for mixture component means (default: 0.01).
#'   - `nu` (optional): Degrees of freedom for IW prior on mixture component covariances (default: nvar + 3).
#'   - `V` (optional): `nvar x nvar` location matrix for IW prior (default: nu * I).
#'   - `a` (optional): `ncomp x 1` Dirichlet prior parameters (default: rep(5, ncomp)).
#'   - `Ad` (optional): Prior precision for vec(Delta) (default: 0.01 * I). Only used when `drawdelta = TRUE` and `useBART = FALSE`.
#'   - `deltabar` (optional): Prior mean for vec(Delta) (default: 0). Only used when `drawdelta = TRUE` and `useBART = FALSE`.
#'   - `bart` (optional): List of HART prior parameters. See Details.
#'   - `vartree` (optional, *experimental extension beyond Wiemann 2025*): List of
#'     parameters enabling heteroscedastic covariance \eqn{\Sigma(Z_i)} via
#'     product-of-trees variance models on the Modified Cholesky diagonal
#'     \eqn{d_j(\cdot)}. Requires `bart` to also be specified and `ncomp = 1`.
#'     When `nvar > 1`, the package automatically promotes to the full-Cholesky
#'     structure described under `phitree`. See "Heteroscedastic Covariance"
#'     in Details.
#'   - `phitree` (optional, *experimental extension beyond Wiemann 2025*): List of
#'     parameters enabling sum-of-trees regression models on the Modified
#'     Cholesky off-diagonals \eqn{\phi_{jk}(\cdot)}. Requires `vartree`.
#'     Auto-enabled when `vartree` is supplied and `nvar > 1`.
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
#' \deqn{\beta_i = \Delta(Z_i) + u_i}
#' where \eqn{u_i \sim N(\mu_{k_i}, \Sigma_{k_i})} with \eqn{k_i} drawn from
#' a mixture of normals with `ncomp` components.
#'
#' ## HART Prior Details
#' If `Prior$bart` is a list, \eqn{\Delta(Z_i)} is modeled via a sum-of-trees (BART) prior.
#' Parameters (with defaults):
#'   - `num_trees`: Number of trees (default: 200).
#'   - `power`, `base`: Tree structure prior parameters (default: 2, 0.95).
#'   - `tau`: Terminal leaf coefficient prior variance (default: 1/sqrt(num_trees)).
#'   - `numcut`: Number of cutpoints (default: 100).
#'   - `sparse`: Use Dirichlet HART for sparsity (default: FALSE).
#'   - `theta`, `omega`, `a`, `b`, `rho`, `aug`, `burn`: Dirichlet HART parameters.
#'
#' ## Heteroscedastic Covariance \eqn{\Sigma(Z_i)} (experimental extension)
#'
#' When `Prior$vartree` is supplied, the homoscedastic mixture-of-normals
#' covariance \eqn{\Sigma_{k_i}} is replaced with a unit-specific
#' \deqn{\Sigma(Z_i)^{-1} = L(Z_i)^\top D(Z_i)^{-1} L(Z_i)}
#' where \eqn{L(Z_i)} is unit lower-triangular with \eqn{L_{jk} = -\phi_{jk}(Z_i)}
#' for \eqn{k < j} and \eqn{D(Z_i) = \mathrm{diag}(d_1(Z_i), \ldots, d_{nvar}(Z_i))},
#' \eqn{d_j > 0}. Each \eqn{d_j(\cdot)} is modeled as a product of trees with
#' \eqn{\chi^{-2}} leaves; each \eqn{\phi_{jk}(\cdot)} (when `Prior$phitree` is
#' supplied or auto-promoted) is modeled as a sum of trees with
#' \eqn{N(0, \tau^2)} leaves. Only the single-component path (`ncomp = 1`) is
#' supported in this mode. Per-unit \eqn{\beta_i} are drawn from the conjugate
#' Gaussian posterior using \eqn{\Sigma(Z_i)^{-1}} as the prior precision.
#'
#' **`Prior$vartree` parameters** (defaults shown):
#'   - `num_trees` (40): Number of trees per dimension.
#'   - `nu` (10), `lambda` (auto-calibrated from per-unit OLS): \eqn{\chi^{-2}}
#'     prior parameters.
#'   - `power` (2), `base` (0.95), `numcut` (100): Tree-structure prior parameters.
#'   - DART hyperparameters (`sparse`, `a`, `b`, `rho`, `theta`, `omega`, `aug`,
#'     `burn`): same meaning as `Prior$bart`.
#'
#' **`Prior$phitree` parameters** (defaults shown):
#'   - `num_trees` (40): Number of trees per \eqn{(j,k)} pair.
#'   - `tau` (1/sqrt(num_trees)): Prior standard deviation of leaf parameters.
#'   - `power` (2), `base` (0.95), `numcut` (100): Tree-structure prior parameters.
#'   - `nmin` (2), `ess_min` (5): Leaf-admissibility floors.
#'   - DART hyperparameters: same meaning as `Prior$bart`.
#'
#' When this extension is active, the returned object additionally inherits
#' classes `"rhierLinearMixtureHeterCov"` and `"bayesm.HART.HeterCov"` and
#' contains slots `var_models`, `phi_models` (jagged, lower-triangular; `NULL`
#' if `nvar == 1`), and `mu_draw`. `predict()` dispatches on these classes to
#' evaluate \eqn{\Sigma(Z^*)} at any new \eqn{Z^*}.
#'
#' @return A list of class `"rhierLinearMixture"` containing:
#'   - `betadraw`: `nreg x nvar x (R/keep)` array of unit-level beta draws.
#'   - `taudraw`: `(R/keep) x nreg` matrix of error variance draws.
#'   - `loglike`: `(R/keep) x 1` vector of log-likelihoods.
#'   - `nmix`: List with `probdraw`, `zdraw`, `compdraw` (omitted under heter-cov).
#'   - If `drawdelta` and non-BART: `Deltadraw`.
#'   - If BART: `bart_models`, `varcount`, `varprob`.
#'   - If heter-cov (additional class `"rhierLinearMixtureHeterCov"`):
#'     `var_models`, `phi_models` (jagged or `NULL`), `mu_draw`,
#'     `var_varcount`, `var_varprob`.
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
#   run hiearchical linear regression model with mixture of normals
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
#        mixture of normals prior on beta heterogeneity
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
if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp (num of mixture components)")}
   else {ncomp=Prior$ncomp}
if(ncomp != 1) {pandterm("Only ncomp = 1 is currently supported")}
if(is.null(Prior$mubar)) { mubar=matrix(rep(0,nvar),nrow=1) }
   else { mubar=matrix(Prior$mubar, nrow=1) }
if(ncol(mubar) != nvar) {pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ",ncol(mubar)))}
if(is.null(Prior$Amu)) { Amu=matrix(BayesmConstant.A,ncol=1) }
   else { Amu=matrix(Prior$Amu, ncol=1) }
if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Am must be a 1 x 1 array")}
if(is.null(Prior$nu)) {nu=nvar+BayesmConstant.nuInc}
   else {nu=Prior$nu}
if(nu < 1) {pandterm("invalid nu value")}
if(is.null(Prior$V)) {V=nu*diag(nvar)}
   else {V=Prior$V}
if(sum(dim(V)==c(nvar,nvar)) != 2) pandterm("Invalid V in prior")

# ---- HART prior modification ----
useBART <- ifelse(!is.null(Prior$bart) & drawdelta, TRUE, FALSE)
bart_params <- NULL
if (useBART) {
  bart_params <- .parse_bart_params(Prior$bart, nz)
}
# ---- end HART prior modification ----

# ---------------------------------------------------------------------------
# Heteroscedastic covariance Sigma(Z_i) via modified-Cholesky tree ensembles.
# Auto-enabled by Prior$vartree.  Requires Prior$bart, drawdelta, ncomp == 1.
# Linear lambda baseline + MCMC warm-start use per-unit OLS slopes (computed
# below once XpX/Xpy are available); see R/heter_cov_priors.R for the parsing.
# ---------------------------------------------------------------------------
hcv_flags   <- .parse_heter_cov_flags(Prior, nvar)
useHeterCov <- hcv_flags$useHeterCov
useFullCov  <- hcv_flags$useFullCov
auto_full   <- hcv_flags$auto_full
var_params  <- list()
phi_params  <- list()
.validate_heter_cov(useHeterCov, useBART, ncomp, drawdelta, SignRes = NULL)

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
  cat("\nHART prior parameters:\n")
  cat(sprintf("num_trees = %d\n", bart_params$num_trees))
  cat(sprintf("power = %g\n", bart_params$power))
  cat(sprintf("base = %g\n", bart_params$base))
  cat(sprintf("tau = %g\n", bart_params$tau))
  cat(sprintf("numcut = %d\n", bart_params$numcut))
  cat(sprintf("sparse = %d\n", bart_params$sparse))
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

oldprob=rep(1/ncomp,ncomp)

# Per-unit OLS slopes: used as the lambda baseline for vartree auto-calibration
# under heter-cov, and also as the C++ warm-start for `oldbetas` (avoids the
# zero-init pathology where d_j(.) trees fit to inflated squared residuals).
# Falls back to zeros if XpX is singular for a unit.
Beta_init <- matrix(0, nrow = nreg, ncol = nvar)
if (useHeterCov) {
  Beta_init <- t(vapply(regdata, function(rd)
    tryCatch(as.numeric(solve(rd$XpX, rd$Xpy)),
             error = function(e) rep(0, nvar)),
    numeric(nvar)))
  var_params <- .parse_var_params(Prior$vartree, Beta_init, nz)
  if (useFullCov) phi_params <- .parse_phi_params(Prior$phitree, nz)
  if (r_verbose)
    .print_heter_cov_summary(useFullCov, var_params, phi_params,
                             lambda_user_supplied = !is.null(Prior$vartree$lambda),
                             auto_full = auto_full)
}

###################################################################
# Wayne Taylor 10/02/2014
# Modified by Thomas Wiemann 2025 -- added HART support, mixture-of-normals,
#   and heter-cov (Sigma(Z_i)) extension.
###################################################################
draws = rhierLinearMixture_rcpp_loop(regdata, Z, deltabar,
    Ad, mubar, Amu, nu, V, nu.e, ssq, R, keep, nprint, drawdelta,
    as.matrix(olddelta), a, oldprob, ind, tau,
    useBART, bart_params,
    useHeterCov, var_params, phi_params, Beta_init)
###################################################################

#
# format output
#
attributes(draws$taudraw)$class=c("bayesm.mat","mcmc")
attributes(draws$taudraw)$mcpar=c(1,R,keep)
attributes(draws$betadraw)$class=c("bayesm.hcoef")

if (useHeterCov) {
  # "bayesm.HART.HeterCov" is the marker superclass shared by all
  # heteroscedastic-Sigma(Z) hierarchical models; it lets predict_helpers.R
  # dispatch generically without hardcoding any model-specific class name.
  attributes(draws$mu_draw)$class <- c("bayesm.mat", "mcmc")
  attributes(draws$mu_draw)$mcpar <- c(1, R, keep)
  class(draws) <- c("rhierLinearMixtureHeterCov",
                    "bayesm.HART.HeterCov",
                    "rhierLinearMixture")
} else {
  if(!useBART & drawdelta) {
    attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
    attributes(draws$Deltadraw)$mcpar=c(1,R,keep)
  }
  attributes(draws$nmix)$class="bayesm.nmix"
  class(draws) <- "rhierLinearMixture"
}

return(draws)
}

# ==============================================================================
# Predict Method
# ==============================================================================

#' Predict Method for rhierLinearMixture Objects
#'
#' Computes posterior draws of the systematic component Delta(Z) for new or
#' existing covariate values.
#'
#' @param object A fitted `rhierLinearMixture` object.
#' @param newdata A list with element `Z`: an `npred x nz` matrix of covariate
#'   values at which to predict.
#' @param type Character. `"DeltaZ"` for the systematic component only, or
#'   `"DeltaZ+mu"` to add the mixture component mean (BART models only), or
#'   `"SigmaZ"` for heteroscedastic covariance draws \eqn{\Sigma(Z)} when the
#'   model was fit with `Prior$vartree`.
#' @param burn Integer, number of initial MCMC draws to discard.
#' @param r_verbose Logical, print progress updates?
#' @param ... Additional arguments passed to `pwbart` for BART models.
#'
#' @return Depends on `type`:
#'   - For `type %in% c("DeltaZ", "DeltaZ+mu")`: 3D array
#'     `[npred, ncoef, ndraws_out]` of predicted betabar values.
#'   - For `type = "SigmaZ"`: 4D array `[npred, ncoef, ncoef, ndraws_out]` of
#'     covariance draws at each prediction unit.
#' @export
predict.rhierLinearMixture <- function(object, newdata = NULL,
                                        type = "DeltaZ+mu", burn = 0,
                                        r_verbose = TRUE, ...) {
  valid_types <- c("DeltaZ", "DeltaZ+mu", "SigmaZ")
  ndraws_total <- .validate_predict_inputs(object, newdata, type, burn, nsim = 1,
                                           valid_types = valid_types)
  kept_draws_indices <- if (burn > 0) (burn + 1):ndraws_total else 1:ndraws_total

  if (type == "SigmaZ") {
    return(.calculate_sigma_z(object, newdata, burn, r_verbose))
  }

  result <- .calculate_delta_z(object, newdata, burn, r_verbose, ...)

  if (type == "DeltaZ+mu") {
    result <- .add_mu_component(result, object, kept_draws_indices)
  }

  return(result)
}
