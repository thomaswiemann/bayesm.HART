#' Bayesian Multinomial Logit Model with HART Prior
#' @description
#' `rhierMnlRwMixture` implements a MCMC algorithm for a Bayesian multinomial 
#' logit model with a Hierarchical Additive Regression Trees (HART) prior. The model 
#' allows for flexible modeling of the representative consumer and captures 
#' unobserved preference heterogeneity.
#'
#' @param Data A list containing:
#'   - `p`: Number of choice alternatives (integer).
#'   - `lgtdata`: A list of length `nlgt`. Each element `lgtdata[[i]]` must be a list with:
#'     - `y`: `n_i x 1` vector of multinomial outcomes (1 to `p`).
#'     - `X`: `(n_i * p) x nvar` matrix of alternative-specific attributes.
#'   - `Z` (optional): `nlgt x nz` matrix of observed characteristics for each unit.
#'     **Should NOT contain an intercept** and is typically centered.
#' @param Prior A list containing prior parameters:
#'   - `ncomp`: Number of mixture components (required).
#'   - `a` (optional): `ncomp x 1` vector of Dirichlet prior parameters for mixture weights `pvec` (default: `rep(5, ncomp)`).
#'   - `deltabar` (optional): `nz * nvar x 1` prior mean for `vec(Delta)` (default: 0). Ignored if BART is used.
#'   - `Ad` (optional): Prior precision matrix for `vec(Delta)` (default: `0.01 * I`). Ignored if BART is used.
#'   - `mubar` (optional): `nvar x 1` prior mean for component means (default: 0 if unrestricted, 2 if restricted).
#'   - `Amu` (optional): Prior precision for component means (default: 0.01 if unrestricted, 0.1 if restricted).
#'   - `nu` (optional): Degrees of freedom for IW prior on component `Sigma` (default: `nvar+3` if unrestricted, `nvar+15` if restricted).
#'   - `V` (optional): Location matrix for IW prior on component `Sigma` (default: `nu * I` or scaled based on restriction).
#'   - `SignRes` (optional): `nvar x 1` vector of sign restrictions (0=none, 1=pos, -1=neg). Default: `rep(0, nvar)`.
#'   - `bart` (optional): List of parameters for the HART prior. If specified, this models the representative consumer \eqn{\Delta(Z)} as a sum-of-trees. See Details.
#' @param Mcmc A list containing MCMC parameters:
#'   - `R`: Number of MCMC iterations (required).
#'   - `keep` (optional): Thinning parameter (default: 1).
#'   - `nprint` (optional): Print progress every `nprint` draws (default: 100, 0 for none).
#'   - `s` (optional): RW Metropolis scaling parameter (default: `2.93 / sqrt(nvar)`).
#'   - `w` (optional): Fractional likelihood weighting parameter (default: 0.1).
#' @param r_verbose Logical. Print startup messages? Default TRUE.
#'
#' @return A list containing:
#'   - `Deltadraw`: If `Z` provided and `bart=NULL`, `(R/keep) x (nz * nvar)` matrix of `vec(Delta)` draws.
#'   - `betadraw`: `nlgt x nvar x (R/keep)` array of unit-level `beta_i` draws.
#'   - `nmix`: List containing mixture draws (`probdraw`, `zdraw`, `compdraw`). See Details.
#'   - `loglike`: `(R/keep) x 1` vector of log-likelihood values at kept draws.
#'   - `SignRes`: `nvar x 1` vector of sign restrictions used.
#'   - `bart_trees`: If BART used, list containing tree structures.
#'
#' @details
#' ## Model and Priors
#' \eqn{y_i \sim MNL(X_i, \beta_i)} for unit \(i = 1, \ldots, nlgt\).
#' The unit-level coefficients (part-worths) \eqn{\beta_i} are modeled as:
#' \deqn{\beta_i = \Delta(Z_i) + u_i}
#' where \eqn{\Delta(Z_i)} is the *representative consumer* component, which depends on observed characteristics \eqn{Z_i}, and \eqn{u_i} is the unobserved heterogeneity component.
#'   - If `Z` is provided and `Prior$bart` is `NULL`: \eqn{\Delta(Z_i) = Z_i \Delta} (linear model).
#'   - If `Z` is provided and `Prior$bart` is a list: \eqn{\Delta(Z_i)} is modeled with a HART prior (sum-of-trees).
#'   - If `Z` is `NULL`: \eqn{\Delta(Z_i) = 0}.
#'
#' The unobserved heterogeneity component \(u_i\) follows a mixture of normals:
#' \deqn{u_i \sim \sum_{j=1}^{ncomp} p_j N(\mu_j, \Sigma_j)}
#'
#' **Priors:**
#' - For mixture weights: \eqn{pvec \sim Dirichlet(a)}
#' - For the linear model: \eqn{\delta = vec(\Delta) \sim N(deltabar, A_d^{-1})}
#' - For the HART model: A sum-of-trees prior is placed on each dimension of a standardized \eqn{\Delta(Z_i)}. See BART details below.
#' - For mixture component means: \eqn{\mu_j \sim N(mubar, \Sigma_j \otimes Amu^{-1})} (Note: Scaled by Sigma_j)
#' - For mixture component covariance: \eqn{\Sigma_j \sim IW(\nu, V)}
#'
#' ## Argument Details
#' **Data List:**
#'   - `p`: Number of alternatives.
#'   - `lgtdata`: List of `nlgt` lists. `lgtdata[[i]] = list(y, X)`.
#'   - `Z`: `nlgt x nz` matrix (optional). Centered, no intercept.
#'
#' **Prior List:**
#'   - `ncomp`: Number of mixture components.
#'   - See `@param` descriptions for defaults.
#'
#' **Mcmc List:**
#'   - `R`: Number of draws.
#'   - See `@param` descriptions for defaults.
#'
#' **HART Prior (`Prior$bart`):** 
#' If `Prior$bart` is a list, it specifies a HART prior for the representative consumer, \eqn{\Delta(Z)}. 
#' This replaces the conventional linear specification \eqn{Z \Delta}. The HART prior models each dimension of 
#' the (standardized) representative consumer as a sum-of-trees.
#' Relevant parameters (defaults used if not specified in `Prior$bart`):
#'   - `num_trees`: Number of trees in each sum-of-trees model (default: 200).
#'   - `power`, `base`: Parameters for the tree structure prior. The probability of a node at depth `q` splitting is \eqn{\alpha(1+q)^{-\beta}}, where `base`=`\alpha` and `power`=`\beta`. Defaults are `base=0.95`, `power=2`, which favors shallow trees.
#'   - `tau`: The standard deviation for the normal prior on terminal leaf coefficients, \eqn{\lambda_{dhg} \sim N(0, \tau^2)}. The default value is `1/sqrt(num_trees)`, which regularizes the model by shrinking individual tree contributions to be small.
#'   - `numcut`: Number of grid points for proposing splitting rules for continuous variables (default: 100).
#'   - `sparse`: If `TRUE`, use the Dirichlet HART prior to induce sparsity. See next section. (default: `FALSE`).
#'   - `burn`: Number of internal burn-in iterations for the BART-sampler within each MCMC iteration (default: 100).
#'
#' ## Dirichlet HART (`sparse = TRUE`)
#' The Dirichlet HART model augments the HART prior to induce sparsity in variable selection, following Linero (2018). The selection probabilities for splitting variables are given a `Dirichlet(zeta/K, ..., zeta/K)` prior, where `K` is the number of characteristics. The concentration parameter `zeta` is given a `Beta(a,b)` hyperprior on `zeta/(zeta+rho)`.
#'   - `a`, `b`: Shape parameters for the Beta hyperprior. The default (`a=0.5, b=1`) induces sparsity.
#'   - `rho`: A parameter that influences the number of selected variables. Default is the number of characteristics.
#'   - `theta`, `omega`: Additional parameters for the sparsity-inducing prior (defaults: 0.0, 1.0).
#'   - `aug`: Logical. For internal use, not relevant for the logit model.
#'
#' ## Sign Restrictions
#' If `SignRes[k]` is non-zero, the k-th coefficient \eqn{\beta_{ik}} is modeled as
#' \deqn{\beta_{ik} = SignRes[k] \cdot exp(\beta^*_{ik}).}
#' The `betadraw` output contains the draws for \eqn{\beta_{ik}} (with the restriction applied).
#' The `nmix` output contains draws for the *unrestricted* mixture components.
#'
#' ## `nmix` Details
#' `nmix` is a list: `list(probdraw, zdraw, compdraw)`
#'   - `probdraw`: `(R/keep) x ncomp` matrix of mixture component probabilities.
#'   - `zdraw`: `(R/keep) x nlgt` matrix of component *assignments* for each unit (i.e., which `j` for \eqn{u_i}).
#'   - `compdraw`: `(R/keep)` list of `ncomp` lists. `compdraw[[r]][[j]] = list(mu, rooti)` contains the draw of \eqn{\mu_j} and \eqn{\Sigma_j^{-1/2}} for component `j` at kept draw `r`.
#'
#' @references
#' Rossi, Peter E., Greg M. Allenby, and Robert McCulloch (2009). Bayesian Statistics and Marketing. Reprint. Wiley Series in Probability and Statistics. Chichester: Wiley.
#' 
#' Rossi, Peter (2023). bayesm: Bayesian Inference for Marketing/Micro-Econometrics. Comprehensive R Archive Network.
#' 
#' Wiemann, Thomas (2025). "Personalization with HART." Working paper.
#'
#' @seealso `rmnlIndepMetrop`
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
  #
  # purpose: run hierarchical mnl logit model with mixture of normals 
  #   using RW and cov(RW inc) = (hess_i + Vbeta^-1)^-1
  #   uses normal approximation to pooled likelihood
  #
  # Arguments:
  #   Data contains a list of (p,lgtdata, and possibly Z)
  #      p is number of choice alternatives
  #      lgtdata is a list of lists (one list per unit)
  #          lgtdata[[i]]=list(y,X)
  #             y is a vector indicating alternative chosen
  #               integers 1:p indicate alternative
  #             X is a length(y)*p x nvar matrix of values of
  #               X vars including intercepts
  #             Z is an length(lgtdata) x nz matrix of values of variables
  #               note: Z should NOT contain an intercept
  #   Prior contains a list of (deltabar,Ad,mubar,Amu,nu,V,ncomp,SignRes) 
  #      ncomp is the number of components in normal mixture
  #           if elements of Prior (other than ncomp) do not exist, defaults are used
  #      SignRes is a vector of sign restrictions
  #   Mcmc contains a list of (s,c,R,keep,nprint)
  #
  # Output:  as list containing
  #   Deltadraw R/keep  x nz*nvar matrix of draws of Delta, first row is initial value
  #   betadraw is nlgt x nvar x R/keep array of draws of betas
  #   probdraw is R/keep x ncomp matrix of draws of probs of mixture components
  #   compdraw is a list of list of lists (length R/keep)
  #      compdraw[[rep]] is the repth draw of components for mixtures
  #   loglike  log-likelikelhood at each kept draw
  #
  # Priors:
  #    beta_i = D %*% z[i,] + u_i
  #       u_i ~ N(mu_ind[i],Sigma_ind[i])
  #       ind[i] ~multinomial(p)
  #       p ~ dirichlet (a)
  #       D is a k x nz array
  #          delta= vec(D) ~ N(deltabar,A_d^-1)
  #    mu_j ~ N(mubar,A_mu^-1(x)Sigma_j)
  #    Sigma_j ~ IW(nu,V^-1)
  #    ncomp is number of components
  #
  # MCMC parameters
  #   s is the scaling parameter for the RW inc covariance matrix; s^2 Var is inc cov
  #      matrix
  #   w is parameter for weighting function in fractional likelihood
  #      w is the weight on the normalized pooled likelihood 
  #   R is number of draws
  #   keep is thinning parameter, keep every keepth draw
  #   nprint - print estimated time remaining on every nprint'th draw
  #
  #  check arguments
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
    bart_params <- list(
      num_trees = ifelse(is.null(Prior$bart$num_trees), 200, 
                         Prior$bart$num_trees),
      power = ifelse(is.null(Prior$bart$power), 2.0, Prior$bart$power),
      base = ifelse(is.null(Prior$bart$base), 0.95, Prior$bart$base),
      tau = ifelse(is.null(Prior$bart$tau), 1.0 / sqrt( # Default calculation
        ifelse(is.null(Prior$bart$num_trees), 200, Prior$bart$num_trees)
      ), Prior$bart$tau), # Use provided tau if exists
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
  } else {
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
        cat("\nBART prior parameters:\n")
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
  
  ###################################################################
  # Wayne Taylor
  # 09/22/2014
  ###################################################################
  draws =  rhierMnlRwMixture_rcpp_loop(lgtdata, Z,
                                       deltabar, Ad, mubar, Amu,
                                       nu, V, s,
                                       R, keep, nprint, drawdelta,
                                       as.matrix(olddelta), a, oldprob, oldbetas, ind, SignRes,
                                       useBART, bart_params)
  ####################################################################
  
  if(!useBART & drawdelta){
    attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
    attributes(draws$Deltadraw)$mcpar=c(1,R,keep)}
  attributes(draws$betadraw)$class=c("bayesm.hcoef")
  attributes(draws$nmix)$class="bayesm.nmix"
  class(draws) <- "rhierMnlRwMixture"
  
  return(draws)
}

# Methods ======================================================================

# --- Internal Helper Functions for predict.rhierMnlRwMixture ---

# Helper 1: Minimal Input Validation
.validate_predict_inputs <- function(object, newdata, type, burn, nsim) {
  # Minimal validation: Type and basic object structure
  valid_types <- c("DeltaZ", "DeltaZ+mu", "posterior_probs", "prior_probs")
  if (!(type %in% valid_types)) {
    stop(paste("Invalid type specified. Choose one of:", 
               paste(valid_types, collapse = ", ")))
  }
  if (is.null(object$betadraw) || !is.array(object$betadraw) || 
      length(dim(object$betadraw)) != 3) {
    stop("Invalid object: object$betadraw missing or not a 3D array.")
  }
  
  ndraws_total <- dim(object$betadraw)[3]
  
  # Minimal validation: burn and nsim (assuming numeric)
  if (burn < 0 || burn >= ndraws_total) {
    stop(paste0("burn must be >= 0 and < ndraws_total (", ndraws_total, ")."))
  }
  if (type == "prior_probs" && nsim <= 0) {
    stop("nsim must be > 0 for type='prior_probs'.")
  }
  
  return(ndraws_total)
}#VALIDATE_PREDICT_INPUTS

# Helper 2: Calculate DeltaZ component
.calculate_delta_z <- function(object, newdata, burn, r_verbose, ...) {
  nvar <- dim(object$betadraw)[2]
  ndraws_total <- dim(object$betadraw)[3]
  
  # Check if Z covariates were used in the original model fit
  model_has_Z <- !is.null(object$Deltadraw) || !is.null(object$bart_models)
  
  pred <- NULL # Initialize pred
  npred <- 1 # Default if no Z
  
  if (model_has_Z) {
    # Model with Z covariates - Assume newdata$Z is provided and is a matrix
    if (is.null(newdata$Z)) {
      stop("Model was fit with Z. 'newdata$Z' must be provided.")
    }
    npred <- nrow(newdata$Z)
    nz <- ncol(newdata$Z)
    
    pred <- raw_pred <- array(0, dim = c(npred, nvar, ndraws_total))
    
    # Calculate DeltaZ component (systematic part)
    if (!is.null(object$Deltadraw)) {
      # Linear model prediction case
      if (ncol(object$Deltadraw) != nz * nvar) {
        stop(paste("Dimension mismatch: ncol(Deltadraw) is", 
                   ncol(object$Deltadraw), "but expected nz*nvar =", nz * nvar))
      }
      for (i in 1:ndraws_total) {
        Delta_draw <- matrix(object$Deltadraw[i, ], nrow = nz, byrow = TRUE)
        pred[, , i] <- newdata$Z %*% Delta_draw
      }#FOR
    } else if (!is.null(object$bart_models)) {
      # BART model prediction case
      if (length(object$bart_models) != nvar) {
        stop(paste("Number of BART models (", length(object$bart_models), 
                   ") does not match nvar (", nvar, ")"))
      }
      # --- Prepare arguments for pwbart, handling ... ---
      passed_args <- list(...)
      pwbart_arg_names <- names(formals(bayesm.HART:::pwbart))
      valid_passed_args <- passed_args[names(passed_args) %in% pwbart_arg_names]
      base_args <- list(
        x.test = newdata$Z,
        mu = 0,
        transposed = FALSE,
        dodraws = TRUE
      )
      final_args_template <- utils::modifyList(base_args, valid_passed_args)
      final_args_template <- 
        final_args_template[names(final_args_template) %in% pwbart_arg_names]
      # --- End argument preparation ---
      
      for (j in 1:nvar) {
        if (is.null(object$bart_models[[j]]$treedraws)) {
          stop(paste("Missing treedraws for BART model of coefficient", j))
        }
        if (r_verbose) {
          cat("Predicting coefficient", j, "with BART model\n")
        }
        current_iter_args <- final_args_template
        current_iter_args$treedraws <- object$bart_models[[j]]$treedraws
        bart_pred <- do.call(bayesm.HART:::pwbart, current_iter_args)
        raw_pred[, j, ] <- t(bart_pred)
      }#FOR
      
      # Apply scaling using rooti (specific to this package's BART impl.)
      if (is.null(object$nmix) || is.null(object$nmix$compdraw)) {
        stop("Missing nmix$compdraw needed for BART scaling.")
      }
      for (s in 1:ndraws_total) {
        if (is.null(object$nmix$compdraw[[s]][[1]]$rooti)) {
          stop(paste("Missing nmix$compdraw[[ ", s, " ]][[1]]$rooti", sep = ""))
        }
        root_i_s <- object$nmix$compdraw[[s]][[1]]$rooti
        if (inherits(try(solve(root_i_s), silent = TRUE), "try-error")) {
          stop(paste("rooti matrix is singular for draw", s))
        }
        root_s_inv <- solve(root_i_s)
        pred[, , s] <- raw_pred[, , s] %*% root_s_inv
      }#FOR
    }#ELSEIF
  } else {
    # Model without Z covariates
    npred <- 1 # Predict only the baseline heterogeneity component
    pred <- array(0, dim = c(npred, nvar, ndraws_total))
  }#ELSE (model_has_Z)
  
  # Apply burn-in
  if (burn >= ndraws_total) stop("Burn-in period too large.") # Should be caught earlier
  if (burn > 0) {
    kept_draws <- (burn + 1):ndraws_total
    pred <- pred[, , kept_draws, drop = FALSE]
  }#IF
  
  return(pred)
}#CALCULATE_DELTA_Z

# Helper 3: Add mu component
.add_mu_component <- function(delta_z_array, object, kept_draws_indices) {
  # Check ncomp and issue warning if > 1
  ncomp <- 1 # Default
  if (!is.null(object$nmix$probdraw)) {
    ncomp <- ncol(object$nmix$probdraw)
  } else if (!is.null(object$nmix$compdraw)) {
    # Fallback: check length of first draw's components if probdraw missing
    ncomp <- length(object$nmix$compdraw[[1]]) 
  }
  if (ncomp > 1) {
    warning("DeltaZ+mu prediction currently only uses mu from the first mixture component.")
  }
  
  # Extract mu from the first component for kept draws
  mudraw <- tryCatch({
    sapply(object$nmix$compdraw[kept_draws_indices], function(x) {
      if (is.null(x[[1]]$mu)) stop("mu missing in component draw")
      x[[1]]$mu
    })
  }, error = function(e) {
    stop(paste("Error extracting mu from nmix$compdraw:", e$message))
  })
  
  # Get dimensions from inputs
  npred <- dim(delta_z_array)[1]
  nvar <- dim(delta_z_array)[2]
  ndraws_out <- dim(delta_z_array)[3]
  
  if (!is.matrix(mudraw) || nrow(mudraw) != nvar || ncol(mudraw) != ndraws_out) {
    stop("Extracted mudraw has incorrect dimensions.")
  }
  
  # Reshape mudraw and add to delta_z_array
  mudraw_array <- array(mudraw, dim = c(nvar, ndraws_out, npred))
  mudraw_array <- aperm(mudraw_array, c(3, 1, 2))
  
  return(delta_z_array + mudraw_array)
}#ADD_MU_COMPONENT

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
      probs_mat_is <- bayesm.HART:::calculate_mnl_probs_from_beta(X_i, beta_is, p)
      
      # Assign the T_i x p matrix to the s-th slice of the 3D array
      prob_array_i[, , draw_count] <- probs_mat_is
    }#FOR s_orig
    prob_pred_list[[i]] <- prob_array_i # Assign 3D array to unit i
  }#FOR i
  if (r_verbose) cat("Posterior probability calculation complete.\n")
  
  return(prob_pred_list)
}#PREDICT_POSTERIOR_PROBS

# Helper 5: Predict Prior Probabilities
.predict_prior_probs <- function(object, newdata, delta_z_array, 
                                 kept_draws_indices_nmix, nsim, r_verbose) {
  if (!requireNamespace("bayesm", quietly = TRUE)) {
    stop("The 'bayesm' package is required for type='prior_probs'. Please install it.")
  }
  
  # Minimal validation assumes object structure is valid
  dims_beta <- dim(object$betadraw) # Used only for nvar
  nvar <- dims_beta[2]
  
  npred <- dim(delta_z_array)[1]
  ndraws_out <- dim(delta_z_array)[3]
  
  # Minimal newdata validation
  if (is.null(newdata$p)) stop("newdata$p must be provided.")
  p <- newdata$p
  if (is.null(newdata$X)) stop("newdata$X must be provided as a list.")
  
  # Initialize Probability List
  prob_pred_list <- vector("list", npred)
  names(prob_pred_list) <- paste0("PredUnit_", 1:npred)
  
  # --- Calculate Probabilities ---
  for (i in 1:npred) { # Loop through prediction units
    if (r_verbose) cat("Calculating prior probabilities for Prediction Unit", i, 
                       "of", npred, "...\n")
    
    # Minimal check on this unit's X (assuming it exists and is matrix)
    X_i <- newdata$X[[i]]
    # Assume dimensions are correct based on prior checks/usage
    
    T_i <- nrow(X_i) / p
    # Initialize 3D array for this unit's probabilities
    prob_array_i <- array(0.0, dim = c(T_i, p, ndraws_out))
    
    for (draw_idx in 1:ndraws_out) { # Loop through kept posterior draws
      # Get systematic component for this unit/draw
      DeltaZ_is <- delta_z_array[i, , draw_idx, drop = TRUE] # Already post-burn
      
      # Get mixture components for this draw (original index)
      s_orig <- kept_draws_indices_nmix[draw_idx]
      comps_s <- object$nmix$compdraw[[s_orig]]
      pvec_s <- object$nmix$probdraw[s_orig, ]
      
      # Simulate eta, calculate beta, get probs, average over nsim
      prob_is_sum <- matrix(0.0, nrow = T_i, ncol = p)
      for (k in 1:nsim) {
        # Simulate heterogeneity component eta_is ~ N(mu_s, Sigma_s) mixture
        # bayesm::rmixture returns a list with $z (comp indicator) and $x (draw)
        eta_is_k <- bayesm::rmixture(n = 1, pvec = pvec_s, comps = comps_s)$x
        eta_is_k <- as.vector(eta_is_k) # Ensure it's a vector
        
        # Combine systematic and heterogeneity parts
        beta_is_k <- DeltaZ_is + eta_is_k
        
        # Calculate probabilities for this simulated beta
        # Assume bayesm.HART namespace is available or ::: is appropriate
        probs_mat_is_k <- 
          bayesm.HART:::calculate_mnl_probs_from_beta(X_i, beta_is_k, p)
        
        prob_is_sum <- prob_is_sum + probs_mat_is_k
      } # End nsim loop
      
      # Average probabilities over nsim simulations
      prob_array_i[, , draw_idx] <- prob_is_sum / nsim
      
    }#FOR draw_idx
    prob_pred_list[[i]] <- prob_array_i # Assign 3D array to prediction unit i
  }#FOR i
  if (r_verbose) cat("Prior probability calculation complete.\n")
  
  return(prob_pred_list)
}#PREDICT_PRIOR_PROBS


# --- Refactored predict.rhierMnlRwMixture Method ---

#' Predict Method for rhierMnlRwMixture Objects
#' @param object A fitted rhierMnlRwMixture object.
#' @param newdata Optional list containing data for prediction. Structure depends
#'   on `type`:
#'   - For `type %in% c("DeltaZ", "DeltaZ+mu")`: Requires `newdata$Z`, a matrix
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
#'   - `"DeltaZ+mu"`: Expected part-worths plus the mean of the unobserved heterogeneity component, \eqn{\Delta(Z) + \mu_j}. Note: for mixtures (`ncomp > 1`), this uses the mean \eqn{\mu_1} from the first component.
#'   - `"posterior_probs"`: Posterior predictive choice probabilities for the original
#'     estimation units using stored `betadraw`.
#'   - `"prior_probs"`: Prior predictive choice probabilities for new prediction units
#'     (based on `newdata$Z` or the overall mixture if no Z was used). Probabilities
#'     are averaged over `nsim` draws from the heterogeneity mixture distribution
#'     per posterior draw.
#' @param burn Integer, number of initial MCMC draws to discard.
#' @param nsim Integer, number of draws from the heterogeneity mixture distribution
#'   per posterior draw for `type = "prior_probs"`.
#' @param r_verbose Logical, print progress updates?
#' @param ... Additional arguments passed to underlying prediction functions
#'   (e.g., `mc.cores`, `verbose` for BART `DeltaZ` predictions via `pwbart`).
#' @return Depends on `type`:
#'   - For `type %in% c("DeltaZ", "DeltaZ+mu")`: 3D array `[npred, nvar, ndraws_out]`
#'     of predicted expected part-worths.
#'   - For `type = "posterior_probs"`: List of length `nlgt`. Each element `\\[[i]]`
#'     is a 3D array `[T_i, p, ndraws_out]` of posterior predictive choice probabilities
#'     for unit `i`.
#'   - For `type = "prior_probs"`: List of length `npred`. Each element `\\[[i]]`
#'     is a 3D array `[T_i, p, ndraws_out]` of prior predictive choice
#'     probabilities for prediction unit `i`.
#' @keywords internal
#' @export
predict.rhierMnlRwMixture <- function(object, newdata = NULL, type = "DeltaZ+mu", burn = 0, nsim = 10, r_verbose = TRUE, ...) {
  
  # 1. Initial Validation & Setup
  ndraws_total <- .validate_predict_inputs(object, newdata, type, burn, nsim)
  kept_draws_indices <- if (burn > 0) (burn + 1):ndraws_total else 1:ndraws_total
  # Note: ndraws_out is derived within helpers where needed
  
  # 2. Dispatch based on type
  if (type == "posterior_probs") {
    
    result <- .predict_posterior_probs(object, newdata, kept_draws_indices, 
                                       r_verbose)
    
  } else if (type == "prior_probs") {
    
    # Calculate the base DeltaZ component (handles burn-in internally)
    # Pass ellipsis for potential pwbart args
    delta_z_array <- .calculate_delta_z(object, newdata, burn, r_verbose, ...) 
    
    # Predict prior probabilities (uses kept_draws_indices for nmix access)
    result <- .predict_prior_probs(object, newdata, delta_z_array, 
                                   kept_draws_indices, nsim, r_verbose)
    
  } else if (type %in% c("DeltaZ", "DeltaZ+mu")) {
    
    # Calculate the base DeltaZ component (handles burn-in internally)
    # Pass ellipsis for potential pwbart args
    result <- .calculate_delta_z(object, newdata, burn, r_verbose, ...) 
    
    if (type == "DeltaZ+mu") {
      # Add mu component (uses kept_draws_indices for nmix access)
      result <- .add_mu_component(result, object, kept_draws_indices)
    }
    
  } else {
    # Should be caught by initial validation, but defensive programming
    stop("Internal error: Unhandled prediction type.") 
  }
  
  # 3. Return result
  return(result)
}#PREDICT.RHIERMNLRWMIXTURE