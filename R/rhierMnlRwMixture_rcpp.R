#' Bayesian Multinomial Logit Model with HART Prior
#' @description
#' `rhierMnlRwMixture` implements an MCMC algorithm for a Bayesian hierarchical multinomial 
#' logit model with a Hierarchical Additive Regression Trees (HART) prior. HART 
#' is a hierarchical nonparametric prior that allows for flexible modeling of the 
#' representative consumer as a function of potentially many observed characteristics.
#' \code{Prior$bart} allows for modeling the conditional mean of the normal prior.
#' 
#' @note Currently, only \code{ncomp = 1} is supported.
#'
#' @param Data A list containing:
#'   - `p`: Number of choice alternatives (integer).
#'   - `lgtdata`: A list of length `nlgt`. Each element `lgtdata[[i]]` must be a list with:
#'     - `y`: `n_i x 1` vector of multinomial outcomes (1 to `p`).
#'     - `X`: `(n_i * p) x nvar` matrix of alternative-specific attributes.
#'   - `Z` (optional): `nlgt x nz` matrix of observed characteristics for each unit.
#'     Should NOT contain an intercept and should be centered.
#'   - `unit_id` (optional): Unit IDs, unique and non-missing, in unit
#'     order. Default: names of the unit list if complete and unique,
#'     otherwise "1", "2", ...
#'   With `previous`, omit `Data` or use NULL for unchanged data. Otherwise
#'   `Data` is a batch and `unit_id` is required: observations append to matching
#'   IDs, new IDs append after old units, and units absent from the batch remain.
#'   X columns must match the saved design. If the fit has Z, supply one batch
#'   Z row per ID with the saved columns; old-unit rows must equal stored rows.
#'   Z is not recentered or rescaled, and cannot be added to a fit without Z.
#'   `p` is inherited when omitted and must match if supplied. Alternative
#'   blocks remain consecutive within each choice observation.
#' @param Prior A list containing prior parameters:
#'   - `ncomp` (required): Number of components. Must be 1.
#'   - `deltabar` (optional): `nz * nvar x 1` prior mean for `vec(Delta)` (default: 0). Ignored if HART is used.
#'   - `Ad` (optional): Prior precision matrix for `vec(Delta)` (default: `0.01 * I`). Ignored if HART is used.
#'   - `mubar` (optional): `nvar x 1` prior mean vector for the normal prior (default: 0 if unrestricted, 2 if restricted).
#'   - `Amu` (optional): Prior precision for the normal prior on \eqn{\mu} (default: 0.01 if unrestricted, 0.1 if restricted). Scalar on the bayesm paths. Under a HART model the prior is the independent \eqn{N(\bar\mu, Amu^{-1})}, and `Amu` may be a scalar (expanded to `Amu * I`) or an `nvar x nvar` SPD matrix.
#'   - `nu` (optional): Degrees of freedom for IW prior on component `Sigma` (default: `nvar+3` if unrestricted, `nvar+15` if restricted). Not used under `hart_model = "direct"` with `hart_cov = "heteroskedastic"`.
#'   - `V` (optional): Location matrix for IW prior on component `Sigma` (default: `nu * I` or scaled based on restriction). Not used under `hart_model = "direct"` with `hart_cov = "heteroskedastic"`.
#'   - `SignRes` (optional): `nvar x 1` vector of sign restrictions. Must contain values of 0, -1, or 1. The value 0 means no restriction, -1 ensures the coefficient is negative, and 1 ensures the coefficient is positive. For example, `SignRes = c(0,1,-1)` means the first coefficient is unconstrained, the second will be positive, and the third will be negative. Default: `rep(0, nvar)`.
#'   - `hart_model` (optional): Preference-model selector, one of `"none"` (default; the bayesm constant or linear \eqn{\Delta(Z)} paths), `"normalized"` (HART prior on the representative consumer, \eqn{\theta_i = \mu + S\{f^*(Z_i) + \varepsilon_i\}}), or `"direct"` (HART mean and covariance trees, \eqn{\theta_i = \mu + f(Z_i) + S_i \varepsilon_i}). Both HART models require `Z`. The tree lists below are optional under a HART model and are an error under `"none"`.
#'   - `bart` (optional): List of hyperparameters for the HART mean trees (`hart_model = "normalized"` or `"direct"`). See Details.
#'   - `hart_cov` (optional, `hart_model = "direct"` only): Covariance model, `"heteroskedastic"` (default; covariate-dependent \eqn{\Sigma(Z_i)} from the variance and off-diagonal ensembles) or `"constant"` (one \eqn{\Sigma \sim IW(\nu, V)} shared by all units, drawn by its conjugate update). Under `"constant"`, `vartree` and `phitree` are an error and the fit carries bayesm's `nmix` instead of the heteroscedastic outputs.
#'   - `vartree` (optional, `hart_model = "direct"` only): List of hyperparameters for the product-of-trees variance models on the Modified Cholesky diagonal \eqn{d_j(\cdot)}. Compatible with `SignRes`. See "Direct HART" in Details.
#'   - `phitree` (optional, `hart_model = "direct"` and `nvar > 1` only): List of hyperparameters for the sum-of-trees models on the Modified Cholesky off-diagonals \eqn{\phi_{jk}(\cdot)}. The off-diagonal ensembles always exist when `nvar > 1`; this list only sets their hyperparameters.
#'   With `previous`, prior settings are inherited; see Continuing a Fit.
#' @param Mcmc A list containing MCMC parameters:
#'   - `R`: Number of MCMC iterations (required).
#'   - `keep` (optional): Thinning parameter (default: 1).
#'   - `nprint` (optional): Print progress every `nprint` draws (default: 100, 0 for none).
#'   - `s` (optional): RW Metropolis scaling parameter (default: `2.93 / sqrt(nvar)`).
#'   - `w` (optional): Fractional likelihood weighting parameter (default: 0.1).
#'   With `previous`, only `R`, `keep` and `nprint` are accepted; `R` counts
#'   additional sweeps and the other controls are inherited.
#' @param r_verbose Logical. Print startup messages? Default TRUE.
#'
#' @param previous A fitted object with `restart`, or NULL for a fresh fit.
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
#' - **Linear model**: \eqn{\delta = vec(\Delta) \sim N(deltabar, A_d^{-1})}
#' - **Component means**: \eqn{\mu_1 \sim N(mubar, \Sigma_1 \otimes Amu^{-1})} (covariance scaled by \eqn{\Sigma_1}) on the bayesm paths. Under a HART model, \eqn{\mu \sim N(mubar, Amu^{-1})} independently of \eqn{\Sigma}.
#' - **Component covariance**: \eqn{\Sigma_1 \sim IW(\nu, V)}
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
#'   - `tau` (`hart_model = "normalized"`): Prior scale of each sum-of-trees factor on the standardized coordinates. Terminal leaf coefficients are \eqn{\lambda_{dhg} \sim N(0, \tau^2 / H)}, so \eqn{\tau} is the prior standard deviation of \eqn{f^*_d(Z)}. Default \eqn{\tau = 1}. Under `hart_model = "direct"` the scales are learned; see `a_tau`, `tau0` below.
#'   - `numcut`: Number of grid points for proposing splitting rules for continuous variables (default: 100).
#'   - `nmin`: Minimum number of units in each child of a split (default: 5).
#'   - `sparse`: If `TRUE`, use the Dirichlet HART prior to induce sparsity in variable selection (default: `FALSE`).
#'
#' **Dirichlet HART** (`sparse = TRUE`): The Dirichlet HART model augments the HART prior to induce sparsity in variable selection, following Linero (2018). Instead of uniform probability for selecting splitting variables, the selection probabilities \eqn{\tau = (\tau^{(1)}, \ldots, \tau^{(K)})} are given a sparse Dirichlet prior: \eqn{(\tau^{(1)}, \ldots, \tau^{(K)}) \sim Dirichlet(\theta/K, \ldots, \theta/K)}, where K is the number of characteristics. The concentration parameter \eqn{\theta} is given a hierarchical prior: \eqn{\theta/(\theta+\rho) \sim Beta(a,b)}.
#'   - `a`, `b`: Shape parameters for the Beta hyperprior. The default (`a=0.5, b=1`) induces sparsity where few variables have high selection probabilities.
#'   - `rho`: Parameter influencing sparsity. Default is the number of characteristics K. Reducing rho below K encourages greater sparsity.
#'   - `theta`: When used, sets Dirichlet concentration parameter without additional hyper-prior (default: 0.0).
#'   - `burn`: Number of internal burn-in iterations for the Dirichlet HART sampler before variable selection is allowed (default: 100).
#'
#' ## Direct HART (`hart_model = "direct"`)
#'
#' The direct model places the trees in coefficient units and lets the covariance vary with \eqn{Z}:
#' \deqn{\theta_i = \mu + f(Z_i) + S_i \varepsilon_i, \qquad \Sigma(Z_i) = S_i S_i^\top = L(Z_i)^{-1} D(Z_i) L(Z_i)^{-\top},}
#' where \eqn{f(Z)} is a sum of `num_trees` trees with vector leaves \eqn{m_{hg} \sim N(0, \mathrm{diag}(\tau_1^2, \ldots, \tau_D^2) / H)} on one shared partition, \eqn{L(Z_i)} is unit lower-triangular with \eqn{L_{jk} = -\phi_{jk}(Z_i)} for \eqn{k < j}, and \eqn{D(Z_i) = \mathrm{diag}(d_1(Z_i), \ldots, d_D(Z_i))}. Each \eqn{d_j(\cdot)} is a product of trees with inverse-gamma leaves (Pratola et al., 2020); each \eqn{\phi_{jk}(\cdot)} is a sum of trees with normal leaves. The prior is not invariant to the ordering of the columns of `X`, which is recorded in `hart$coef_names`. `Prior$nu` and `Prior$V` are not used.
#'
#' **`Prior$bart` under the direct model**: `num_trees` (200), `power` (2), `base` (0.95), `numcut` (100), `nmin` (5), the DART fields, and the scale hyperprior \eqn{\tau_d^2 \sim IG(a_\tau, (a_\tau - 1)\tau_{0d}^2)} through `a_tau` (3) and `tau0` (1; length 1 or `nvar`), so that \eqn{E[\tau_d^2] = \tau_{0d}^2}. `tau` is not used.
#'
#' **`Prior$vartree` parameters** (defaults shown):
#'   - `num_trees` (40), `power` (2), `base` (0.95), `numcut` (100), `nmin` (5): Tree-structure prior parameters.
#'   - `tau` (1), `v0` (1 unrestricted, 0.1 sign-restricted; length 1 or `nvar`): Ensemble scale and level. Leaves are \eqn{s \sim IG(a_v, b_{vj})} with \eqn{H_v \psi'(a_v) = \tau^2} and \eqn{b_{vj} = (a_v - 1) v_{0j}^{1/H_v}}, so that \eqn{\log d_j(Z)} has prior standard deviation \eqn{\tau} and \eqn{E[d_j(Z)] = v_{0j}} for any `num_trees` (\eqn{\psi'} is the trigamma function).
#'   - DART hyperparameters (`sparse`, `a`, `b`, `rho`, `theta`, `omega`, `aug`, `burn`): same meaning as `Prior$bart`.
#'
#' **`Prior$phitree` parameters** (defaults shown):
#'   - `num_trees` (40), `power` (2), `base` (0.95), `numcut` (100), `nmin` (5): Tree-structure prior parameters.
#'   - `tau` (1): Ensemble scale; leaves are \eqn{N(0, \tau^2 / H_\phi)}.
#'   - DART hyperparameters: same meaning as `Prior$bart`.
#'
#' The returned object additionally inherits class `"rhierMnlRwMixtureHeterCov"` and contains `mu_draw`, `mean_scale_draw`, `var_models`, and `phi_models` (jagged; `NULL` if `nvar == 1`). `predict()` dispatches on this class to evaluate \eqn{\Sigma(Z^*)} at new \eqn{Z^*}.
#'
#' ## Sign Restrictions
#' If `SignRes[k]` is non-zero, the k-th coefficient \eqn{\beta_{ik}} is modeled as
#' \deqn{\beta_{ik} = SignRes[k] \cdot \exp(\beta^*_{ik}).}
#' The `betadraw` output contains the draws for \eqn{\beta_{ik}} (with the restriction applied).
#' The `nmix` output contains draws for the *unrestricted* prior covariance and mean.
#' 
#' **Note:** Care should be taken when selecting priors on any sign restricted coefficients.
#'
#' ## Continuing a Fit
#' The result contains only the new run's draws. Old draws are not pooled
#' with them. Priors, tree cutpoints and sign restrictions are inherited
#' from the first fit. New units outside the original Z range are resolved
#' on the original grid. Settling sweeps may be needed after adding data.
#' Passing the same batch twice counts its observations twice. Fits retain
#' their data and are therefore larger. Fits made before restart support
#' cannot be continued. `Prior` is not accepted with `previous`.
#'
#' @return A list containing:
#'   - `Deltadraw`: If `Z` provided and `bart=NULL`, `(R/keep) x (nz * nvar)` matrix of `vec(Delta)` draws.
#'   - `betadraw`: `nlgt x nvar x (R/keep)` array of unit-level `beta_i` draws.
#'   - `nmix`: bayesm list of `(mu, rooti)` draws (omitted under `hart_model = "direct"` with `hart_cov = "heteroskedastic"`).
#'   - `loglike`: `(R/keep) x 1` vector of log-likelihood values at kept draws.
#'   - `SignRes`: `nvar x 1` vector of sign restrictions used.
#'   - `acceptrbeta`: Metropolis acceptance rate (percent) for the unit-level
#'     MNL random-walk updates of `beta_i`.
#'   - `bart_models`: Under `hart_model = "normalized"`, a list of length `nvar` of scalar tree ensembles for \eqn{f^*_d(\cdot)}. Under `hart_model = "direct"`, one vector-leaf ensemble (`treedraws` with `leaf_dim = nvar`) for \eqn{f(\cdot)}.
#'   - `varcount`, `varprob`: Split counts and selection probabilities per characteristic; `nz x nvar x (R/keep)` under `"normalized"`, `nz x (R/keep)` under `"direct"`.
#'   - `var_models` (`hart_model = "direct"`): list of length `nvar` of variance-tree ensembles for \eqn{d_j(\cdot)}, with `var_varcount`, `var_varprob`.
#'   - `phi_models` (`hart_model = "direct"`, `nvar > 1`): jagged list of length `nvar`. `phi_models[[1]]` is `NULL`; for `j > 1`, `phi_models[[j]]` is a list of length `j-1` containing the sum-of-trees ensemble for \eqn{\phi_{jk}(\cdot)}.
#'   - `mu_draw` (`hart_model = "direct"`): `(R/keep) x nvar` matrix of \eqn{\mu} draws.
#'   - `mean_scale_draw` (`hart_model = "direct"`): `(R/keep) x nvar` matrix of \eqn{\tau_d^2} draws.
#'   - `restart`: Retained data, resolved configuration, proposal preparation
#'     and terminal sampler state for continuation through `previous`.
#'
#' @references
#' Chipman, Hugh A., Edward I. George, and Robert E. McCulloch (2010). "BART: Bayesian Additive Regression Trees." Annals of Applied Statistics 4.1.
#' 
#' Linero, Antonio R. (2018). "Bayesian regression trees for high-dimensional prediction and variable selection." Journal of the American Statistical Association 113.522, pp. 626-636.
#' 
#' Pratola, M. T., Chipman, H. A., George, E. I., and McCulloch, R. E. (2020). "Heteroscedastic BART via Multiplicative Regression Trees." Journal of Computational and Graphical Statistics 29.2, pp. 405-417.
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
#' @examples
#' \donttest{
#' sim <- sim_hier_mnl(nlgt = 12, nT = 5, p = 3, nz = 2, nXa = 2,
#'                     nXd = 0, const = FALSE, seed = 11)
#' fit <- rhierMnlRwMixture(Data = sim, Prior = list(ncomp = 1),
#'                          Mcmc = list(R = 200, keep = 1, nprint = 0))
#' continued <- rhierMnlRwMixture(previous = fit, Mcmc = list(R = 200))
#' batch_units <- lapply(sim$lgtdata[1:2], function(u)
#'   list(y = sample.int(sim$p, 1), X = u$X[seq_len(sim$p), , drop = FALSE]))
#' batch <- list(lgtdata = batch_units, unit_id = c("1", "new"),
#'               Z = sim$Z[1:2, , drop = FALSE])
#' extended <- rhierMnlRwMixture(Data = batch, previous = continued,
#'                                 Mcmc = list(R = 200))
#' }
#'
#' @export
rhierMnlRwMixture=function(Data,Prior,Mcmc, r_verbose = TRUE, previous = NULL){
  #
  # revision history:
  #   12/04 changed by rossi to fix bug in drawdelta when there is zero/one unit in a prior component
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
  saved=NULL
  state=NULL
  batch=FALSE
  if(!is.null(previous)) {
    saved=.hart_previous(previous, "mnl")
    batch=!missing(Data) && !is.null(Data)
    if(!missing(Prior)) {
      pandterm("Prior is inherited with previous")
    }
    if(missing(Mcmc)) pandterm("Requires R in Mcmc for continuation")
    Mcmc=.hart_continue_controls(Mcmc, saved)
    if(batch) Data=.hart_merge_data(saved, Data) else Data=saved$data
    state=saved$state
  }
  
  if(missing(Data)) {pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")}
  if(is.null(Data$p)) {pandterm("Requires Data element p (# choice alternatives)") }
  p=Data$p
  if(is.null(Data$lgtdata)) {pandterm("Requires Data element lgtdata (list of data for each unit)")}
  lgtdata=Data$lgtdata
  restart_data=.hart_restart_data(Data, "lgtdata")
  nlgt=length(lgtdata)
  drawdelta=TRUE
  if(is.null(Data$Z)) { cat("Z not specified",fill=TRUE); fsh() ; drawdelta=FALSE}
  else {if (!is.matrix(Data$Z)) {pandterm("Z must be a matrix")}
    else {if (nrow(Data$Z) != nlgt) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number logits ",nlgt))}
      else {Z=Data$Z}}}
  if(drawdelta) {
    nz=ncol(Z)
    if(is.null(saved)) {
    colmeans=apply(Z,2,mean)
    if(sum(colmeans) > .00001) 
    {pandterm(paste("Z does not appear to be de-meaned: colmeans= ",colmeans))}
    }
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
  if(is.null(saved)) {
  if(missing(Prior)) {pandterm("Requires Prior list argument (at least ncomp)")} 
  if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp")} else {ncomp=Prior$ncomp}
  if(ncomp != 1) {pandterm("Only ncomp = 1 is currently supported")}
  if(is.null(Prior$SignRes)) {SignRes=rep(0,nvar)} else {SignRes=Prior$SignRes}                     
  if(length(SignRes) != nvar) {pandterm("The length SignRes must be equal to the dimension of X")}
  if(sum(!(SignRes %in% c(-1,0,1))>0)) {pandterm("All elements of SignRes must be equal to -1, 0, or 1")}
  # HART preference model (R/hart_priors.R): "none" keeps the bayesm paths.
  # Sign-restricted coordinates get the smaller default heterogeneity scale.
  cfg <- .hart_config(Prior, nvar = nvar, nz = if (drawdelta) nz else 0L,
                      drawdelta = drawdelta,
                      coef_names = colnames(lgtdata[[1]]$X),
                      v0_default = ifelse(SignRes != 0, 0.1, 1))
  if(is.null(Prior$mubar) & sum(abs(SignRes))==0) {
    mubar=matrix(rep(0,nvar),nrow=1)
  } else { 
    if(is.null(Prior$mubar) & sum(abs(SignRes)) >0) {
      mubar=matrix(rep(0,nvar)+2*abs(SignRes),nrow=1)
    } else { 
      mubar=matrix(Prior$mubar,nrow=1) } }
  if(ncol(mubar) != nvar) {pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ",ncol(mubar)))}
  if (cfg$model == "none") {
    if(is.null(Prior$Amu) & sum(abs(SignRes))==0) {
      Amu=matrix(BayesmConstant.A,ncol=1)
    } else {
      if(is.null(Prior$Amu) & sum(abs(SignRes)) >0) {
        Amu=matrix(BayesmConstant.A*10,ncol=1)
      } else {Amu=matrix(Prior$Amu,ncol=1) } }
    if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Am must be a 1 x 1 array")}
  } else {
    # HART: independent prior mu ~ N(mubar, Amu^{-1}); scalar or nvar x nvar.
    Amu <- .hart_Amu(Prior$Amu, nvar,
                     default = if (sum(abs(SignRes)) == 0) BayesmConstant.A
                               else BayesmConstant.A * 10)
  }
  cfg$Amu <- Amu
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
  # ---------------------------------------------------------------------------
  # HART preference model.  Prior$hart_model selects "none" (bayesm paths),
  # "normalized", or "direct"; see R/hart_priors.R.
  # ---------------------------------------------------------------------------
  useBART     <- cfg$model != "none"
  useHeterCov <- cfg$model == "direct" && identical(cfg$cov, "heteroskedastic")
  bart_params <- cfg$bart
  var_params  <- cfg$vartree %||% list()
  phi_params  <- cfg$phitree %||% list()
  Ad = deltabar = NULL
  if (!useBART) {
    if(is.null(Prior$Ad) & drawdelta) {Ad=BayesmConstant.A*diag(nvar*nz)} else {Ad=Prior$Ad}
    if(drawdelta) {if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}}
    if(is.null(Prior$deltabar)& drawdelta) {deltabar=rep(0,nz*nvar)} else {deltabar=Prior$deltabar}
    if(drawdelta) {if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}}
  }
  if(is.null(Prior$a)) { a=rep(BayesmConstant.a,ncomp)} else {a=Prior$a}
  if(length(a) != ncomp) {pandterm("Requires dim(a)= ncomp (no of components)")}
  bada=FALSE
  for(i in 1:ncomp) { if(a[i] < 0) bada=TRUE}
  if(bada) pandterm("invalid values in a vector")
  
  if(is.null(Prior$nu)&sum(abs(SignRes))>0) {nu = nvar+15}                                        
  if(cfg$model == "none" & is.null(Prior$Amu) & sum(abs(SignRes))>0) {Amu = matrix(.1)}                                   
  if(is.null(Prior$V)&sum(abs(SignRes))>0) {V = nu*(diag(nvar)-diag(abs(SignRes)>0)*.8)} 
  
  } else {
    cfg=.hart_previous_config(saved)
    ncomp=cfg$ncomp
    mubar=cfg$prior$mubar
    Amu=cfg$prior$Amu
    nu=cfg$prior$nu
    V=cfg$prior$V
    Ad=cfg$prior$Ad
    deltabar=cfg$prior$deltabar
    a=cfg$prior$a
    SignRes=cfg$prior$SignRes
    useBART=cfg$model != "none"
    useHeterCov=cfg$model == "direct" &&
                identical(cfg$cov, "heteroskedastic")
    bart_params=cfg$bart
    var_params=cfg$vartree %||% list()
    phi_params=cfg$phitree %||% list()
  }

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
        .print_hart_prior(cfg)
        if (useHeterCov) .print_hart_covariance(cfg)
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
  } # end r_verbose
  
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
  if(is.null(saved) || batch) {
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
  # BFGS also for nvar = 1; the default Nelder-Mead warns on scalars.
  out=optim(betainit,llmnl_con,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-6), 
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
  } else {
    for(i in 1:nlgt)
      lgtdata[[i]]$hess=saved$preparation$hess[[i]]
    oldbetas=state$betas
  }
  #
  #  initialize values
  #
  # set initial values for the indicators
  #     ind is of length(nlgt) and indicates which prior component this obs
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
  # Unique-row map for exact cache storage in C++ prediction payload.
  if (drawdelta && useBART) {
    z_map <- .build_unique_z_map(Z)
    bart_params$z_index <- z_map$z_index
    bart_params$n_unique <- nrow(z_map$Z_unique)
  } else {
    z_map <- NULL
  }

  ###################################################################
  # Wayne Taylor
  # 09/22/2014
  ###################################################################
  if(batch) state=.hart_extend_state(saved, oldbetas)

  draws =  rhierMnlRwMixture_rcpp_loop(lgtdata, Z,
                                       deltabar, Ad, mubar, Amu,
                                       nu, V, s,
                                       R, keep, nprint, drawdelta,
                                       as.matrix(olddelta), a, oldprob, oldbetas, ind, SignRes,
                                       cfg$model, cfg$cov %||% "heteroskedastic",
                                       bart_params,
                                       var_params, phi_params, state)
  ####################################################################

  draws <- .finish_hart_fit(draws, cfg, z_map, R, keep,
                            base_class = "rhierMnlRwMixture")

  draws$restart=.hart_restart("mnl", restart_data, cfg, ncomp, drawdelta,
      list(mubar=mubar, Amu=Amu, nu=nu, V=V, Ad=Ad, deltabar=deltabar,
           a=a, SignRes=SignRes),
      list(keep=keep, nprint=nprint, s=s, w=w),
      list(hess=lapply(lgtdata, function(l) l$hess)))

  draws$restart$state=draws$state
  draws$state=NULL

  return(draws)
}

# Methods ======================================================================

#' Predict Method for rhierMnlRwMixture Objects
#' @param object A fitted rhierMnlRwMixture object.
#' @param newdata Optional list containing data for prediction. Structure depends
#'   on `mode` and `type`:
#'   - For `mode = "coefficients"` with `type %in% c("DeltaZ", "DeltaZ+mu", "SigmaZ")`:
#'     requires `newdata$Z`, a matrix with `npred` rows for prediction units
#'     (if model was fit with Z).
#'   - For `mode = "posterior"` with `type = "choice_probs"`: requires
#'     `newdata$nlgtdata`, a list of length `nlgt` (original number of units).
#'     Each element `\\[[i]]` must contain `$X`, the design matrix `(T_i*p) x nvar`
#'     for unit `i`. Also requires `newdata$p`, the number of alternatives.
#'   - For `mode = "prior"` with `type = "choice_probs"`: requires `newdata$Z`
#'     (if model fit with Z, determining `npred`), `newdata$p`, and `newdata$X`
#'     (a list of length `npred`, each element `\\[[i]]` having design matrix
#'     `(T_i*p) x nvar`).
#' @param type Type of prediction within the selected `mode`:
#'   - `"DeltaZ"`: Expected part-worths of the representative consumer, \eqn{\Delta(Z)}.
#'   - `"DeltaZ+mu"`: Expected part-worths plus the mean of the unobserved heterogeneity component, \eqn{\Delta(Z) + \mu_1}. The package supports only `ncomp = 1`.
#'   - `"choice_probs"`: Predictive choice probabilities (use with
#'     `mode = "posterior"` or `mode = "prior"`).
#'   - `"SigmaZ"`: Draws of the covariance matrix \eqn{\Sigma(Z)}. Under
#'     `hart_model = "direct"` with `hart_cov = "heteroskedastic"` these vary
#'     with `Z`; otherwise the single \eqn{\Sigma} draw is repeated on every
#'     row of `newdata$Z`.
#' @param mode Prediction mode:
#'   - `"coefficients"`: coefficient-level outputs (`DeltaZ`, `DeltaZ+mu`, `SigmaZ`).
#'   - `"posterior"`: posterior predictive output (`type = "choice_probs"`).
#'   - `"prior"`: prior predictive output (`type = "choice_probs"`).
#' @param burn Integer, number of initial MCMC draws to discard.
#' @param nsim Integer, number of draws from the heterogeneity distribution
#'   per posterior draw for `mode = "prior"` and `type = "choice_probs"`.
#' @param force_tree_eval Logical. Evaluate the stored trees at `newdata$Z`
#'   even when the cached draws cover every row (default `FALSE`).
#' @param r_verbose Logical, print progress updates?
#' @param ... Additional arguments passed to underlying prediction functions
#'   (e.g., `mc.cores`, `verbose` for tree-based `DeltaZ` predictions).
#' @return Depends on `type`:
#'   - For `type %in% c("DeltaZ", "DeltaZ+mu")`: 3D array `[npred, nvar, ndraws_out]`
#'     of predicted expected part-worths.
#'   - For `type = "SigmaZ"`: 4D array `[npred, nvar, nvar, ndraws_out]` of
#'     covariance draws at each prediction unit.
#'   - For `mode = "posterior", type = "choice_probs"`: List of length `nlgt`.
#'     Each element `\\[[i]]`
#'     is a 3D array `[T_i, p, ndraws_out]` of posterior predictive choice probabilities
#'     for unit `i`.
#'   - For `mode = "prior", type = "choice_probs"`: List of length `npred`.
#'     Each element `\\[[i]]`
#'     is a 3D array `[T_i, p, ndraws_out]` of prior predictive choice
#'     probabilities for prediction unit `i`.
#' @keywords internal
#' @export
predict.rhierMnlRwMixture <- function(object, newdata = NULL,
                                      type = "DeltaZ+mu", burn = 0, nsim = 10,
                                      mode = "coefficients",
                                      r_verbose = TRUE,
                                      force_tree_eval = FALSE, ...) {
  mode <- .resolve_predict_mode(mode)

  if (mode %in% c("prior", "posterior")) {
    if (!identical(type, "choice_probs")) {
      stop("For rhierMnlRwMixture predictive modes, use type='choice_probs'.")
    }
    return(.predictive_dispatch(
      object, newdata, mode = mode, type = "choice_probs",
      burn = burn, nsim = nsim, r_verbose = r_verbose,
      force_tree_eval = force_tree_eval, ...
    ))
  }

  structural_types <- c("DeltaZ", "DeltaZ+mu", "SigmaZ")
  ndraws_total <- .validate_predict_inputs(object, newdata, type, burn, nsim,
                                           valid_types = structural_types)
  kept_draws_indices <- .kept_draw_indices(ndraws_total, burn)

  # Under heteroskedastic covariance the ensembles are evaluated once at
  # newdata$Z and shared by Delta(Z) and Sigma(Z).
  hetercov_comps <- NULL
  if (type == "SigmaZ") {
    return(.predict_structural_common(
      object, newdata, type, burn, r_verbose,
      force_tree_eval = force_tree_eval,
      hetercov_comps = hetercov_comps
    ))
  }

  if (type %in% c("DeltaZ", "DeltaZ+mu")) {
    return(.predict_structural_common(
      object, newdata, type, burn, r_verbose,
      force_tree_eval = force_tree_eval,
      hetercov_comps = hetercov_comps, ...
    ))
  }
}#PREDICT.RHIERMNLRWMIXTURE