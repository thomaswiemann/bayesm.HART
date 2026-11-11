#' Bayesian Hierarchical Negative Binomial Model with HART Prior
#' @description
#' `rhierNegbinRw` implements an MCMC algorithm for a Bayesian hierarchical
#' negative binomial regression model with a Hierarchical Additive Regression Trees (HART) prior. HART 
#' is a hierarchical nonparametric prior that allows for flexible modeling of the 
#' representative unit as a function of potentially many observed characteristics.
#' Per-unit \eqn{\beta_i} are drawn via Random Walk Metropolis-Hastings; the dispersion \eqn{\alpha} is drawn on the log scale.
#' \code{Prior$bart} allows for modeling the conditional mean of the normal prior.
#' 
#' @note Currently, only \code{ncomp = 1} is supported.
#'
#' @param Data A list containing:
#'   - `regdata`: A list of length `nreg`. Each element `regdata[[i]]` must be a list with:
#'     - `y`: `n_i x 1` vector of count responses.
#'     - `X`: `n_i x nvar` design matrix (including intercept column).
#'   - `Z` (optional): `nreg x nz` matrix of unit-level covariates. If omitted,
#'     `drawdelta` is set to FALSE and no second-stage covariates are used.
#' @param Prior A list containing prior parameters:
#'   - `ncomp` (required): Number of components. Must be 1.
#'   - `nu` (optional): Degrees of freedom for IW prior on the normal prior covariance (default: nvar + 3).
#'   - `V` (optional): `nvar x nvar` location matrix for IW prior (default: nu * I).
#'   - `a` (optional): Shape parameter for gamma prior on alpha (default: 0.5). Note: This differs from `rhierLinearMixture` where `a` is the Dirichlet prior.
#'   - `b` (optional): Rate parameter for gamma prior on alpha (default: 0.1).
#'   - `mubar` (optional): `1 x nvar` prior mean vector for the normal prior (default: 0).
#'   - `Amu` (optional): Prior precision for the normal prior on \eqn{\mu} (default: 0.01). Scalar on the bayesm paths. Under a HART model the prior is the independent \eqn{N(\bar\mu, Amu^{-1})}, and `Amu` may be a scalar (expanded to `Amu * I`) or an `nvar x nvar` SPD matrix.
#'   - `Ad` (optional): `nvar*nz x nvar*nz` prior precision for vec(Delta) (default: 0.01 * I). Only used when `drawdelta = TRUE` and `useBART = FALSE`.
#'   - `deltabar` (optional): `nvar*nz` prior mean for vec(Delta) (default: 0). Only used when `drawdelta = TRUE` and `useBART = FALSE`.
#'   - `hart_model` (optional): Preference-model selector, one of `"none"`
#'     (default; bayesm constant or linear \eqn{\Delta(Z)}), `"normalized"`
#'     (HART prior on the representative unit), or `"direct"` (HART mean and
#'     covariance trees). Both HART models require `Z`. The tree lists below
#'     are optional under a HART model and are an error under `"none"`. See
#'     `rhierMnlRwMixture` for the model definitions.
#'   - `bart` (optional): List of hyperparameters for the HART mean trees
#'     (`hart_model = "normalized"` or `"direct"`). See `rhierMnlRwMixture`.
#'   - `vartree` (optional, `hart_model = "direct"` only): List of
#'     hyperparameters for the product-of-trees variance models on the Modified
#'     Cholesky diagonal \eqn{d_j(\cdot)}.
#'   - `phitree` (optional, `hart_model = "direct"` and `nvar > 1` only): List
#'     of hyperparameters for the sum-of-trees models on the Modified Cholesky
#'     off-diagonals \eqn{\phi_{jk}(\cdot)}.
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
#'   \eqn{\beta_i = Z_i \Delta + u_i} \cr
#'   \eqn{u_i \sim N(0, \Sigma_i)} \cr
#'   Note: This documentation describes `ncomp = 1` behavior as higher-order components are not supported. \cr
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
#' classes `"rhierNegbinRwHeterCov"` and `"bayesm.HART.HeterCov"`. `predict()` 
#' dispatches on these classes to evaluate \eqn{\Sigma(Z^*)} at any new \eqn{Z^*}.
#'
#' @return A list of class `"rhierNegbinRw"` containing:
#'   - `betadraw`: `nreg x nvar x (R/keep)` array of unit-level beta draws.
#'   - `alphadraw`: `(R/keep) x 1` vector of overdispersion draws.
#'   - `loglike`: `(R/keep) x 1` vector of log-likelihoods.
#'   - `nmix`: Legacy list containing prior draws (omitted under heter-cov).
#'   - `acceptrbeta`, `acceptralpha`: Metropolis acceptance rates (percent).
#'   - If `drawdelta` and non-BART: `Deltadraw`.
#'   - If BART: `bart_models`, `varcount`, `varprob`.
#'   - If heter-cov (additional class `"rhierNegbinRwHeterCov"`):
#'     `var_models`, `phi_models` (jagged or `NULL`), `mu_draw`,
#'     `mean_scale_draw`, `var_varcount`, `var_varprob`.
#'
#' @seealso [predict.rhierNegbinRw()], [marginal_effects.rhierNegbinRw()]
#'
#' @author Sridhar Narayanan, Peter Rossi (original bayesm code), Thomas Wiemann (HART modifications).
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
if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp")}
   else {ncomp=Prior$ncomp}
if(ncomp != 1) {pandterm("Only ncomp = 1 is currently supported")}
# HART preference model (R/hart_priors.R): "none" keeps the bayesm paths.
cfg <- .hart_config(Prior, nvar = nvar, nz = if (drawdelta) nz else 0L,
                    drawdelta = drawdelta,
                    coef_names = colnames(regdata[[1]]$X))
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
if (cfg$model == "none") {
  if(is.null(Prior$Amu)) { Amu=matrix(BayesmConstant.A,ncol=1) } else { Amu=matrix(Prior$Amu, ncol=1) }
  if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Amu must be a 1 x 1 array")}
} else {
  # HART: independent prior mu ~ N(mubar, Amu^{-1}); scalar or nvar x nvar.
  Amu <- .hart_Amu(Prior$Amu, nvar, default = BayesmConstant.A)
}
cfg$Amu <- Amu
if(is.null(Prior$a_mix)) { a_mix=rep(BayesmConstant.a, ncomp) } else { a_mix=Prior$a_mix }
if(length(a_mix) != ncomp) {pandterm("Requires dim(a_mix)= ncomp (no of components)")}

# ---- HART preference model (see R/hart_priors.R) ----
# Prior$hart_model selects "none" (bayesm paths), "normalized", or "direct".
useBART     <- cfg$model != "none"
useHeterCov <- cfg$model == "direct"
bart_params <- cfg$bart
var_params  <- cfg$vartree %||% list()
phi_params  <- cfg$phitree %||% list()
# ---- end HART preference model ----

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
  .print_hart_prior(cfg)
  if (useHeterCov) .print_hart_covariance(cfg)
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

# Warm start at the per-unit MLEs under the direct model.
if (useHeterCov) Beta = Beta_init

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

if (drawdelta && useBART) {
  z_map <- .build_unique_z_map(Z)
  bart_params$z_index <- z_map$z_index
  bart_params$n_unique <- nrow(z_map$Z_unique)
} else {
  z_map <- NULL
}

###################################################################
# Wayne Taylor 12/01/2014
# Modified by Thomas Wiemann 2025 -- HART support
###################################################################
if (fixalpha) {alpha=Mcmc$alpha}
draws=rhierNegbinRw_rcpp_loop(regdata, hess_i, Z, Beta, deltabar,
                              Ad, mubar, Amu, nu, V, a, b,
                              R, keep, s_beta, alphacroot, nprint,
                              drawdelta, as.matrix(olddelta), a_mix,
                              oldprob, ind,
                              alpha, fixalpha,
                              cfg$model, bart_params,
                              var_params, phi_params)
###################################################################

# === Format output ===
attributes(draws$alphadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$alphadraw)$mcpar=c(1,R,keep)
draws <- .finish_hart_fit(draws, cfg, z_map, R, keep,
                          base_class = "rhierNegbinRw")

return(draws)
}

# ==============================================================================
# Predict Method
# ==============================================================================

#' Predict Method for rhierNegbinRw Objects
#'
#' Computes coefficient-level or predictive count-response draws.
#'
#' @param object A fitted `rhierNegbinRw` object.
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
#'     with simulated count draws per unit (`nobs_i x ndraws_out`).
#' @export
predict.rhierNegbinRw <- function(object, newdata = NULL,
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
      stop("For rhierNegbinRw predictive modes, use type='response'.")
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
