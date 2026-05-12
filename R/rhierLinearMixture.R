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
#' @return A list of class `"rhierLinearMixture"` containing:
#'   - `betadraw`: `nreg x nvar x (R/keep)` array of unit-level beta draws.
#'   - `taudraw`: `(R/keep) x nreg` matrix of error variance draws.
#'   - `loglike`: `(R/keep) x 1` vector of log-likelihoods.
#'   - `nmix`: List with `probdraw`, `zdraw`, `compdraw` (mixture draws).
#'   - If `drawdelta` and non-BART: `Deltadraw`.
#'   - If BART: `bart_models`, `varcount`, `varprob`.
#'
#' @seealso [predict.rhierLinearMixture()], [marginal_effects.rhierLinearMixture()]
#'
#' @author Peter Rossi, Wayne Taylor (original bayesm code), Thomas Wiemann (HART modifications).
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
bart_params = NULL
if (useBART) {
  bart_params <- list(
    num_trees = ifelse(is.null(Prior$bart$num_trees), 200, Prior$bart$num_trees),
    power = ifelse(is.null(Prior$bart$power), 2.0, Prior$bart$power),
    base = ifelse(is.null(Prior$bart$base), 0.95, Prior$bart$base),
    tau = ifelse(is.null(Prior$bart$tau), 1.0 / sqrt(
      ifelse(is.null(Prior$bart$num_trees), 200, Prior$bart$num_trees)
    ), Prior$bart$tau),
    numcut = ifelse(is.null(Prior$bart$numcut), 100, Prior$bart$numcut),
    sparse = ifelse(is.null(Prior$bart$sparse), FALSE, Prior$bart$sparse),
    theta = ifelse(is.null(Prior$bart$theta), 0.0, Prior$bart$theta),
    omega = ifelse(is.null(Prior$bart$omega), 1.0, Prior$bart$omega),
    a = ifelse(is.null(Prior$bart$a), 0.5, Prior$bart$a),
    b = ifelse(is.null(Prior$bart$b), 1.0, Prior$bart$b),
    rho = ifelse(is.null(Prior$bart$rho), nz, Prior$bart$rho),
    aug = ifelse(is.null(Prior$bart$aug), FALSE, Prior$bart$aug),
    burn = ifelse(is.null(Prior$bart$burn), 100, Prior$bart$burn)
  )
}
# ---- end HART prior modification ----

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

###################################################################
# Wayne Taylor 10/02/2014
# Modified by Thomas Wiemann 2025 -- added HART support
###################################################################
draws = rhierLinearMixture_rcpp_loop(regdata, Z, deltabar,
    Ad, mubar, Amu, nu, V, nu.e, ssq, R, keep, nprint, drawdelta,
    as.matrix(olddelta), a, oldprob, ind, tau,
    useBART, bart_params)
###################################################################

#
# format output
#
attributes(draws$taudraw)$class=c("bayesm.mat","mcmc")
attributes(draws$taudraw)$mcpar=c(1,R,keep)
if(!useBART & drawdelta) {
  attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
  attributes(draws$Deltadraw)$mcpar=c(1,R,keep)
}
attributes(draws$betadraw)$class=c("bayesm.hcoef")
attributes(draws$nmix)$class="bayesm.nmix"
class(draws) <- "rhierLinearMixture"

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
#'   `"DeltaZ+mu"` to add the mixture component mean (BART models only).
#' @param burn Integer, number of initial MCMC draws to discard.
#' @param r_verbose Logical, print progress updates?
#' @param ... Additional arguments passed to `pwbart` for BART models.
#'
#' @return A 3D array `[npred, ncoef, ndraws_out]` of predicted betabar values.
#' @export
predict.rhierLinearMixture <- function(object, newdata = NULL,
                                        type = "DeltaZ+mu", burn = 0,
                                        r_verbose = TRUE, ...) {
  ndraws_total <- .validate_predict_inputs(object, newdata, type, burn, nsim = 1)
  kept_draws_indices <- if (burn > 0) (burn + 1):ndraws_total else 1:ndraws_total

  result <- .calculate_delta_z(object, newdata, burn, r_verbose, ...)

  if (type == "DeltaZ+mu") {
    result <- .add_mu_component(result, object, kept_draws_indices)
  }

  return(result)
}
