#' Bayesian Hierarchical Negative Binomial Model with Normal Mixture and HART Prior
#' @description
#' `rhierNegbinRw` implements an MCMC algorithm for a Bayesian hierarchical
#' negative binomial regression model with a mixture of normals prior on the
#' unit-level coefficients. Per-unit \eqn{\beta_i} are drawn via Random Walk
#' Metropolis-Hastings; the dispersion \eqn{\alpha} is drawn on the log scale.
#' Supports both standard linear hierarchical prior and HART (sum-of-trees) prior.
#' Structurally symmetric to `rhierLinearMixture` and `rhierMnlRwMixture`.
#'
#' @param Data A list containing:
#'   - `regdata`: A list of length `nreg`. Each element `regdata[[i]]` must be a list with:
#'     - `y`: `n_i x 1` vector of count responses.
#'     - `X`: `n_i x nvar` design matrix (including intercept column).
#'   - `Z` (optional): `nreg x nz` matrix of unit-level covariates. If omitted,
#'     `drawdelta` is set to FALSE and no second-stage covariates are used.
#' @param Prior A list containing prior parameters:
#'   - `ncomp` (required): Number of normal mixture components.
#'   - `nu` (optional): Degrees of freedom for IW prior on mixture component covariances (default: nvar + 3).
#'   - `V` (optional): `nvar x nvar` location matrix for IW prior (default: nu * I).
#'   - `a` (optional): Shape parameter for gamma prior on alpha (default: 0.5). Note: This differs from `rhierLinearMixture` where `a` is the Dirichlet prior.
#'   - `b` (optional): Rate parameter for gamma prior on alpha (default: 0.1).
#'   - `mubar` (optional): `1 x nvar` prior mean for mixture component means (default: 0).
#'   - `Amu` (optional): `1 x 1` prior precision for mixture component means (default: 0.01).
#'   - `a_mix` (optional): `ncomp x 1` Dirichlet prior parameters for mixture weights (default: rep(5, ncomp)). Named `a_mix` here to avoid collision with the gamma prior shape parameter `a`.
#'   - `Ad` (optional): `nvar*nz x nvar*nz` prior precision for vec(Delta) (default: 0.01 * I). Only used when `drawdelta = TRUE` and `useBART = FALSE`.
#'   - `deltabar` (optional): `nvar*nz` prior mean for vec(Delta) (default: 0). Only used when `drawdelta = TRUE` and `useBART = FALSE`.
#'   - `bart` (optional): List of HART prior parameters. See `rhierLinearMixture` for details.
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
#'   - `s_beta` (optional): RW scaling for beta (default: 2.93/sqrt(nvar)).
#'   - `s_alpha` (optional): RW scaling for alpha (default: 2.93).
#'   - `w` (optional): Fractional likelihood weight used in candidate Hessian construction (default: 0.1).
#'   - `alpha` (optional): Fixed alpha value. If provided, alpha is not sampled (sets `fixalpha = TRUE`).
#'   - `fixalpha` (optional): Logical, whether to fix alpha (default: FALSE).
#'
#' @param r_verbose Logical. Print startup messages? Default TRUE.
#'
#' @details
#' ## Model Specification
#' \eqn{(y_i|\lambda_i, \alpha) \sim NegBin(\lambda_i, \alpha)},
#' \eqn{\ln(\lambda_i) = X_i \beta_i}
#'
#' The unit-level coefficients are modeled as:
#' \deqn{\beta_i = \Delta(Z_i) + u_i}
#' where \eqn{u_i \sim N(\mu_{k_i}, \Sigma_{k_i})} with \eqn{k_i} drawn from
#' a mixture of normals with `ncomp` components.
#'
#' ## HART Prior Details
#' Same as `rhierLinearMixture`. If `Prior$bart` is a list, \eqn{\Delta(Z_i)} is
#' modeled via a sum-of-trees (BART) prior instead of a linear hierarchical
#' specification.
#'
#' ## Heteroscedastic Covariance \eqn{\Sigma(Z_i)} (experimental extension)
#'
#' Same modified-Cholesky construction as `rhierLinearMixture`. When
#' `Prior$vartree` is supplied, the homoscedastic mixture covariance
#' \eqn{\Sigma_{k_i}} is replaced with a unit-specific
#' \deqn{\Sigma(Z_i)^{-1} = L(Z_i)^\top D(Z_i)^{-1} L(Z_i)}
#' where \eqn{d_j(\cdot)} is a product-of-trees ensemble with \eqn{\chi^{-2}}
#' leaves and \eqn{\phi_{jk}(\cdot)} is a sum-of-trees ensemble with
#' \eqn{N(0, \tau^2)} leaves. Per-unit \eqn{\beta_i} are drawn via a
#' Random-Walk Metropolis-Hastings step that uses \eqn{\Sigma(Z_i)^{-1}} as
#' the prior precision. Only the single-component path (`ncomp = 1`) is
#' supported in this mode.
#'
#' Hyperparameters for `Prior$vartree` and `Prior$phitree` mirror those
#' documented in `rhierLinearMixture`. The `lambda` parameter of `vartree`
#' is auto-calibrated from per-unit MLE \eqn{\hat\beta_i} when not supplied.
#'
#' When this extension is active, the returned object additionally inherits
#' classes `"rhierNegbinRwHeterCov"` and `"bayesm.HART.HeterCov"` and
#' contains slots `var_models`, `phi_models` (jagged, lower-triangular; `NULL`
#' if `nvar == 1`), and `mu_draw`. `predict()` dispatches on these classes to
#' evaluate \eqn{\Sigma(Z^*)} at any new \eqn{Z^*}.
#'
#' @return A list of class `"rhierNegbinRw"` containing:
#'   - `betadraw`: `nreg x nvar x (R/keep)` array of unit-level beta draws.
#'   - `alphadraw`: `(R/keep) x 1` vector of overdispersion draws.
#'   - `loglike`: `(R/keep) x 1` vector of log-likelihoods.
#'   - `nmix`: List with `probdraw`, `zdraw`, `compdraw` (omitted under heter-cov).
#'   - `acceptrbeta`, `acceptralpha`: Metropolis acceptance rates (percent).
#'   - If `drawdelta` and non-BART: `Deltadraw`.
#'   - If BART: `bart_models`, `varcount`, `varprob`.
#'   - If heter-cov (additional class `"rhierNegbinRwHeterCov"`):
#'     `var_models`, `phi_models` (jagged or `NULL`), `mu_draw`,
#'     `var_varcount`, `var_varprob`.
#'
#' @seealso [predict.rhierNegbinRw()], [marginal_effects.rhierNegbinRw()]
#'
#' @author Sridhar Narayanan, Peter Rossi (original bayesm code), Thomas Wiemann (HART + mixture modifications).
#'
#' @examples
#' \donttest{
#' set.seed(20260513)
#' sim <- bayesm.HART::sim_hier_negbin(
#'   nreg = 30, nobs = 12, nvar = 1, nz = 2,
#'   const = TRUE, het_observed = "linear",
#'   target_var_betabar = 0.5, target_var_eps = 0.25,
#'   alpha = 5.0)
#' Prior <- list(ncomp = 1L,
#'               bart    = list(num_trees = 20),
#'               vartree = list(num_trees = 20))
#' Mcmc  <- list(R = 200L, keep = 1L, nprint = 0L)
#' fit   <- bayesm.HART::rhierNegbinRw(
#'   Data = list(regdata = sim$regdata, Z = sim$Z),
#'   Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
#' str(fit, max.level = 1)
#' }
#'
#' @export
rhierNegbinRw= function(Data, Prior, Mcmc, r_verbose = TRUE) {
#   Revision History
#	  Sridhar Narayanan - 05/2005
#         P. Rossi 6/05
#         fixed error with nobs not specified and changed llnegbinFract 9/05
#         3/07 added classes
#         3/08 fixed fractional likelihood
#   W. Taylor 4/15 - added nprint option to MCMC argument
#   T. Wiemann 6/25 - added HART prior + roxygen2 docs
#
#   Model
#       (y_i|lambda_i,alpha) ~ Negative Binomial(Mean = lambda_i, Overdispersion par = alpha)
#
#       ln(lambda_i) =  X_i * beta_i
#
#       beta_i = Delta'*z_i + nu_i
#               nu_i~N(0,Vbeta)
#
#   Priors
#       vec(Delta|Vbeta) ~ N(vec(Deltabar), Vbeta (x) (Adelta^-1))
#       Vbeta ~ Inv Wishart(nu, V)
#       alpha ~ Gamma(a,b) where mean = a/b and variance = a/(b^2)

#
# Definitions of functions used within rhierNegbinRw (but outside of Rcpp loop)
#
llnegbinR = function(par,X,y, nvar) {
# Computes the log-likelihood
    beta = par[1:nvar]
    alpha = exp(par[nvar+1])+1.0e-50
    mean=exp(X%*%beta)
    prob=alpha/(alpha+mean)
    prob=ifelse(prob<1.0e-100,1.0e-100,prob)
     out=stats::dnbinom(y,size=alpha,prob=prob,log=TRUE)
     return(sum(out))
}

llnegbinFract = 
function(par,X,y,Xpooled, ypooled, w,wgt, nvar,lnalpha)  {
# Computes the fractional log-likelihood at the unit level
    theta = c(par,lnalpha)
    (1-w)*llnegbinR(theta,X,y,nvar) + w*wgt*llnegbinR(theta,Xpooled,ypooled, nvar) 
}

#
# Error Checking
#
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata and (possibly) Z")}

if(is.null(Data$regdata)) {
    pandterm("Requires Data element regdata -- list of data for each unit : y and X")
}
regdata=Data$regdata
nreg = length(regdata)

drawdelta=TRUE
if (is.null(Data$Z)) {
    if (r_verbose) cat("Z not specified - using a column of ones instead", fill = TRUE)
    fsh()
    drawdelta=FALSE
    Z = matrix(rep(1,nreg),ncol=1)
}
else {
  if (!is.matrix(Data$Z)) {
      pandterm("Z must be a matrix")
    }
    else {
      if (nrow(Data$Z) != nreg) {
          pandterm(paste("Nrow(Z) ", nrow(Data$Z), "ne number units ",nreg))
      }
      else {
          Z = Data$Z
      }
    }
}
nz = ncol(Z)
for (i in 1:nreg) {
  if(!is.matrix(regdata[[i]]$X)) {pandterm(paste0("regdata[[",i,"]]$X must be a matrix"))}
  if(!is.vector(regdata[[i]]$y, mode = "numeric") & !is.vector(regdata[[i]]$y, mode = "logical") & !is.matrix(regdata[[i]]$y)) 
    {pandterm(paste0("regdata[[",i,"]]$y must be a numeric or logical vector or matrix"))}
  if(is.matrix(regdata[[i]]$y)) { if(ncol(regdata[[i]]$y)>1) {pandterm(paste0("regdata[[",i,"]]$y must be a vector or one-column matrix"))}}
}
dimfun = function(l) {
    c(length(l$y),dim(l$X))
}
dims=sapply(regdata,dimfun)
dims = t(dims)
nvar = quantile(dims[,3],prob=0.5)
for (i in 1:nreg) {
        if (dims[i, 1] != dims[i, 2] || dims[i, 3] != nvar) {
            pandterm(paste("Bad Data dimensions for unit", i, 
                "dims(y,X) =", dims[i, ]))
        }
}

ypooled = NULL
Xpooled = NULL
for (i in 1:nreg) {
    ypooled = c(ypooled,regdata[[i]]$y)
    Xpooled = rbind(Xpooled,regdata[[i]]$X)
}
nobs= length(ypooled)

nvar=ncol(Xpooled)
#
# check for prior elements
#
if(missing(Prior)) {
    pandterm("Requires Prior list argument (at least ncomp)")
}
if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp (num of mixture components)")}
   else {ncomp=Prior$ncomp}
if(is.null(Prior[["nu"]])) {nu=nvar+BayesmConstant.nuInc} else {nu=Prior[["nu"]]}
if(is.null(Prior[["V"]])) {V=nu*diag(nvar)} else {V=Prior[["V"]]}
if(is.null(Prior[["a"]])) {a=BayesmConstant.agammaprior} else {a=Prior[["a"]]}
if(is.null(Prior[["b"]])) {b=BayesmConstant.bgammaprior} else {b=Prior[["b"]]}
if(nu < nvar) pandterm("invalid nu value")
if(sum(dim(V)==c(nvar,nvar)) != 2) pandterm("V is of incorrect dimension")
if((length(a) != 1) | (a <=0)) pandterm("a should be a positive number")
if((length(b) != 1) | (b <=0)) pandterm("b should be a positive number")

# Mixture prior parameters (mirrors MNL/linear)
if(is.null(Prior$mubar)) { mubar=matrix(rep(0,nvar),nrow=1) } else { mubar=matrix(Prior$mubar, nrow=1) }
if(ncol(mubar) != nvar) {pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ",ncol(mubar)))}
if(is.null(Prior$Amu)) { Amu=matrix(BayesmConstant.A,ncol=1) } else { Amu=matrix(Prior$Amu, ncol=1) }
if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Amu must be a 1 x 1 array")}
if(is.null(Prior$a_mix)) { a_mix=rep(BayesmConstant.a, ncomp) } else { a_mix=Prior$a_mix }
if(length(a_mix) != ncomp) {pandterm("Requires dim(a_mix)= ncomp (no of components)")}

# ---- HART prior modification ----
useBART <- ifelse(!is.null(Prior$bart) & drawdelta, TRUE, FALSE)
bart_params = NULL
if (useBART) {
  bart_params <- .parse_bart_params(Prior$bart, nz)
}

# Delta prior (for non-BART drawdelta path)
if(is.null(Prior$Ad) & drawdelta) {Ad=BayesmConstant.A*diag(nvar*nz)}
   else {Ad=Prior$Ad}
if(drawdelta & !useBART) {
  if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}
}
if(is.null(Prior$deltabar) & drawdelta) {deltabar=rep(0,nz*nvar)}
   else {deltabar=Prior$deltabar}
if(drawdelta & !useBART) {
  if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}
}
# ---- end HART prior modification ----

# ---------------------------------------------------------------------------
# Heteroscedastic covariance Sigma(Z_i) via modified-Cholesky tree ensembles.
# Auto-enabled by Prior$vartree.  Requires Prior$bart, drawdelta, ncomp == 1.
# Negbin lambda baseline uses per-unit MLE betas (Beta_init below); see
# R/heter_cov_priors.R for the parsing helpers.
# ---------------------------------------------------------------------------
hcv_flags   <- .parse_heter_cov_flags(Prior, nvar)
useHeterCov <- hcv_flags$useHeterCov
useFullCov  <- hcv_flags$useFullCov
auto_full   <- hcv_flags$auto_full
var_params  <- list()
phi_params  <- list()
.validate_heter_cov(useHeterCov, useBART, ncomp, drawdelta, SignRes = NULL)

#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
  if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
if(is.null(Mcmc$s_alpha)) { s_alpha=BayesmConstant.RRScaling} 
    else {s_alpha= Mcmc$s_alpha }
if(is.null(Mcmc$s_beta)) { s_beta=BayesmConstant.RRScaling/sqrt(nvar)} 
    else {s_beta=Mcmc$s_beta }
if(is.null(Mcmc$w)) { w=BayesmConstant.w} 
    else {w = Mcmc$w}
if(is.null(Mcmc$alpha)) {fixalpha=FALSE} else {fixalpha=TRUE; alpha=Mcmc$alpha}
if(is.null(Mcmc$fixalpha)) {} else {fixalpha=Mcmc$fixalpha}
if(fixalpha && alpha<=0) pandterm("alpha is not positive")

# print out problem
#
if (r_verbose) {
cat(" ",fill=TRUE)
cat("Starting Random Walk Metropolis Sampler for Hierarchical Negative Binomial Regression",fill=TRUE)
cat("  ",nobs," obs; ",nvar," covariates (including the intercept); ",fill=TRUE)
if(drawdelta) {
  cat("  ",nz," individual characteristics (including the intercept) ",fill=TRUE)
} else {
  cat("  ","Z not specified (intercept only) ",fill=TRUE)
}
cat(" ",fill=TRUE)
cat("Prior Parameters:",fill=TRUE)
if (useBART) {
  cat("\nHART prior parameters:\n")
  cat(sprintf("num_trees = %d\n", bart_params$num_trees))
  cat(sprintf("power = %g\n", bart_params$power))
  cat(sprintf("base = %g\n", bart_params$base))
  cat(sprintf("tau = %g\n", bart_params$tau))
  cat(sprintf("numcut = %d\n", bart_params$numcut))
  cat(sprintf("sparse = %d\n", bart_params$sparse))
} else {
  cat("deltabar",fill=TRUE)
  print(deltabar)
  cat("Ad",fill=TRUE)
  print(Ad)
}
cat("nu",fill=TRUE)
print(nu)
cat("V",fill=TRUE)
print(V)
if (!fixalpha) {
    cat("a",fill=TRUE)
    print(a)
    cat("b",fill=TRUE)
    print(b)
}
cat(" ",fill=TRUE)
cat("MCMC Parameters:",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat("s_alpha = ",s_alpha,fill=TRUE)
cat("s_beta = ",s_beta,fill=TRUE)
if (fixalpha) {
    cat("alpha",fill=TRUE)
    print(alpha)
}
cat("Fractional Likelihood Weight Parameter = ",w,fill=TRUE)
cat(" ",fill=TRUE)
} # end r_verbose

par = rep(0,(nvar+1))
if (r_verbose) {
  cat("initializing Metropolis candidate densities for ",nreg," units ...",fill=TRUE)
  fsh()
}
mle = optim(par,llnegbinR, X=Xpooled, y=ypooled, nvar=nvar, 
      method="L-BFGS-B", upper=c(rep(Inf,nvar),log(100000000)), hessian=TRUE, control=list(fnscale=-1))
fsh()
beta_mle=mle$par[1:nvar]
alpha_mle = exp(mle$par[nvar+1])
varcovinv = -mle$hessian
Beta = t(matrix(rep(beta_mle,nreg),ncol=nreg))
if(!fixalpha) {alpha = alpha_mle}
alphacvar = s_alpha/varcovinv[nvar+1,nvar+1]
alphacroot = sqrt(alphacvar)
if (r_verbose) {
  cat("beta_mle = ",beta_mle,fill=TRUE)
  cat("alpha_mle = ",alpha_mle, fill = TRUE)
  fsh()
}

hess_i=NULL
# Per-unit MLE betas (collected for heter-cov lambda calibration; harmless
# overhead for non-heter-cov paths).  Default to pooled MLE on convergence
# failure so lambda_default = mean(apply(Beta_init, 2, var)) stays finite.
Beta_init = matrix(rep(beta_mle, nreg), nrow = nreg, ncol = nvar, byrow = TRUE)
if(nobs > 1000){
  sind=sample(c(1:nobs),size=1000)
  ypooleds=ypooled[sind]
  Xpooleds=Xpooled[sind,]
  }
else{
	ypooleds=ypooled
	Xpooleds=Xpooled
}
# Find the individual candidate hessian
for (i in 1:nreg) {
    wgt = length(regdata[[i]]$y)/length(ypooleds)
    mle2 = optim(mle$par[1:nvar],llnegbinFract, X=regdata[[i]]$X, y=regdata[[i]]$y, Xpooled=Xpooleds, 
           ypooled=ypooleds, w=w,wgt=wgt, nvar=nvar, lnalpha=mle$par[nvar+1], 
           method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=0,reltol=1e-6))
    if (mle2$convergence==0) {
        hess_i[[i]] = list(hess=-mle2$hessian)
        Beta_init[i, ] = mle2$par
    } else {
        hess_i[[i]] = list(hess=diag(rep(1,nvar)))
        # Beta_init[i, ] keeps its pooled-MLE default.
    }
   if(i%%50 ==0 & r_verbose) cat("  completed unit #",i,fill=TRUE)	
   fsh()
}

# Heter-cov hyperparameter setup (Beta_init now available for lambda calibration).
if (useHeterCov) {
  var_params <- .parse_var_params(Prior$vartree, Beta_init, nz)
  if (useFullCov) phi_params <- .parse_phi_params(Prior$phitree, nz)
  if (r_verbose)
    .print_heter_cov_summary(useFullCov, var_params, phi_params,
                             lambda_user_supplied = !is.null(Prior$vartree$lambda),
                             auto_full = auto_full)
  # Warm-start MCMC at per-unit MLEs under heter-cov; the BART path keeps the
  # original pooled-MLE start for backward compatibility.
  Beta = Beta_init
}

# ---- Mixture initialization (mirrors MNL/linear) ----
ninc=floor(nreg/ncomp)
ind=NULL
for (i in 1:(ncomp-1)) {ind=c(ind,rep(i,ninc))}
if(ncomp != 1) {ind = c(ind,rep(ncomp,nreg-length(ind)))} else {ind=rep(1,nreg)}
oldprob=rep(1/ncomp,ncomp)

if (drawdelta & useBART){
  deltabar = 0
  Ad = matrix(0)
} else if (!drawdelta) {
  Z = matrix(0)
  deltabar = 0
  Ad = matrix(0)
  bart_params = list(0)
}
if (!useBART) {
  bart_params = list(0)
}

if (drawdelta) {
  olddelta = rep(0, nz*nvar)
} else {
  olddelta = 0
}
# ---- end initialization ----

###################################################################
# Wayne Taylor 12/01/2014
# Modified by Thomas Wiemann 2025 -- mixture-of-normals support
###################################################################
if (fixalpha) {alpha=Mcmc$alpha}
draws=rhierNegbinRw_rcpp_loop(regdata, hess_i, Z, Beta, deltabar,
                              Ad, mubar, Amu, nu, V, a, b,
                              R, keep, s_beta, alphacroot, nprint,
                              drawdelta, as.matrix(olddelta), a_mix,
                              oldprob, ind,
                              alpha, fixalpha,
                              useBART, bart_params,
                              useHeterCov, var_params, phi_params)
###################################################################

# === Format output ===
# Normalize Betadraw -> betadraw for consistency with shared predict helpers
# (Now handled directly in C++)
attributes(draws$alphadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$alphadraw)$mcpar=c(1,R,keep)
attributes(draws$betadraw)$class=c("bayesm.hcoef")

if (useHeterCov) {
  # "bayesm.HART.HeterCov" is the marker superclass shared by all
  # heteroscedastic-Sigma(Z) hierarchical models; it lets predict_helpers.R
  # dispatch generically without hardcoding any model-specific class name.
  attributes(draws$mu_draw)$class <- c("bayesm.mat", "mcmc")
  attributes(draws$mu_draw)$mcpar <- c(1, R, keep)
  class(draws) <- c("rhierNegbinRwHeterCov",
                    "bayesm.HART.HeterCov",
                    "rhierNegbinRw")
} else {
  if(!useBART & drawdelta) {
    attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
    attributes(draws$Deltadraw)$mcpar=c(1,R,keep)
  }
  attributes(draws$nmix)$class <- "bayesm.nmix"
  class(draws) <- "rhierNegbinRw"
}

return(draws)
}

# ==============================================================================
# Predict Method
# ==============================================================================

#' Predict Method for rhierNegbinRw Objects
#'
#' Computes posterior draws of the systematic component Delta(Z) for new or
#' existing covariate values.
#'
#' @param object A fitted `rhierNegbinRw` object.
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
predict.rhierNegbinRw <- function(object, newdata = NULL,
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
