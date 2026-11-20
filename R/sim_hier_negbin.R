#' Simulate Hierarchical Negative Binomial Data
#'
#' Generates simulated data from a hierarchical negative binomial model with
#' optional observed heterogeneity via Z covariates.
#'
#' @param nreg Number of cross-sectional units.
#' @param nobs Number of observations per unit.
#' @param nvar Number of X variables (excluding intercept if const=TRUE).
#' @param nz Number of Z variables. Set to 0 for no observed heterogeneity.
#' @param const Logical. Include an intercept in X? Default TRUE.
#' @param het_observed Character. Functional form of observed heterogeneity.
#'   Options: "none", "linear", "step".
#' @param target_var_betabar Numeric. Target variance for the first coefficient
#'   of the observed component. Default 0.5.
#' @param target_var_eps Numeric. Target variance for each coefficient in the
#'   unobserved component eps_i. Default 0.25.
#' @param alpha Numeric. Overdispersion parameter for negative binomial. Default 5.
#' @param Sigma_func Function of one row of `Z` returning the `ncoef x ncoef`
#'   covariance of `eps_i`. Overrides `target_var_eps` with covariate-dependent
#'   heterogeneity; requires `Z`. Default `NULL` (constant covariance).
#' @param seed Integer. Optional random seed.
#'
#' @return A list containing:
#'   \item{regdata}{List of length nreg. Each element has y, X.}
#'   \item{hessdata}{List of length nreg. Each element has hess (approx Hessian).}
#'   \item{Z}{nreg x nz matrix of covariates (NULL if nz=0).}
#'   \item{true_values}{List with beta_true, betabar_true, alpha, Sigma_eps,
#'     and Sigma_true (nreg x ncoef x ncoef, with Sigma_func).}
#'
#' @importFrom stats rnorm rnbinom var
#' @export
sim_hier_negbin <- function(nreg = 100, nobs = 50, nvar = 2, nz = 3,
                             const = TRUE,
                             het_observed = c("none", "linear", "step"),
                             target_var_betabar = 0.5,
                             target_var_eps = 0.25,
                             alpha = 5.0,
                             Sigma_func = NULL,
                             seed = NULL) {

  het_observed <- match.arg(het_observed)
  if (!is.null(seed)) set.seed(seed)

  ncoef <- nvar + as.integer(const)

  # Z
  Z <- NULL
  if (nz > 0 && het_observed != "none") {
    Z <- matrix(rnorm(nreg * nz), nrow = nreg, ncol = nz)
    Z <- scale(Z, center = TRUE, scale = FALSE)
  }

  # betabar_i
  betabar_true <- matrix(0, nrow = nreg, ncol = ncoef)

  if (!is.null(Z) && het_observed != "none") {
    raw_k1 <- numeric(nreg)
    if (het_observed == "linear") {
      Delta_true <- matrix(0, nrow = nz, ncol = ncoef)
      nz_used <- min(nz, 5)
      fixed_coefs <- c(0.5, -0.5, 0.3, -0.3, 0.1)
      if (nz_used > 0) Delta_true[1:nz_used, 1] <- fixed_coefs[1:nz_used]
      betabar_true <- Z %*% Delta_true
      raw_k1 <- betabar_true[, 1]
    } else if (het_observed == "step") {
      for (i in 1:nreg) raw_k1[i] <- ifelse(Z[i, 1] > 0, 0.5, -0.5)
    }

    raw_var <- var(raw_k1)
    if (raw_var > 1e-10 && target_var_betabar > 0) {
      scale_factor <- sqrt(target_var_betabar / raw_var)
      raw_k1 <- (raw_k1 - mean(raw_k1)) * scale_factor
    } else {
      raw_k1 <- rep(0, nreg)
    }
    betabar_true[, 1] <- raw_k1

    # Small fixed values keep lambda stable.
    if (ncoef > 1) {
      fixed_remaining <- c(0.5, -0.3, 0.2, -0.1)
      n_remaining <- min(ncoef - 1, length(fixed_remaining))
      for (k in 2:(1 + n_remaining)) {
        betabar_true[, k] <- fixed_remaining[k - 1]
      }
    }
  } else {
    # Small constant betabar.
    if (ncoef >= 1) betabar_true[, 1] <- 0.5
    if (ncoef >= 2) betabar_true[, 2] <- -0.3
  }

  # eps_i
  Sigma_eps <- diag(target_var_eps, ncoef)
  Sigma_true <- NULL
  if (is.null(Sigma_func)) {
    eps_true <- matrix(rnorm(nreg * ncoef), nrow = nreg, ncol = ncoef)
    eps_true <- eps_true %*% chol(Sigma_eps)
  } else {
    if (is.null(Z))
      stop("Sigma_func requires Z (nz > 0, het_observed != 'none').")
    Sigma_true <- .sigma_true_array(Sigma_func, Z, ncoef)
    eps_true <- .draw_heteroskedastic_eps(Sigma_true)
  }

  beta_true <- betabar_true + eps_true

  # regdata
  regdata <- vector("list", nreg)
  hessdata <- vector("list", nreg)

  for (i in 1:nreg) {
    Xi <- matrix(rnorm(nobs * nvar, sd = 0.5), nrow = nobs, ncol = nvar)
    if (const) Xi <- cbind(1, Xi)

    lambda_i <- exp(Xi %*% beta_true[i, ])
    # Clamp lambda.
    lambda_i <- pmin(lambda_i, 1e6)
    prob_i <- alpha / (alpha + lambda_i)

    yi <- rnbinom(nobs, size = alpha, prob = prob_i)

    # Information matrix (positive definite) -- bayesm convention
    w_i <- as.vector(lambda_i * alpha / (alpha + lambda_i))
    hess_i <- crossprod(Xi * sqrt(w_i))

    regdata[[i]] <- list(y = as.double(yi), X = Xi)
    hessdata[[i]] <- list(hess = hess_i)
  }

  true_values <- list(
    beta_true = beta_true,
    betabar_true = betabar_true,
    eps_true = eps_true,
    alpha = alpha,
    Sigma_eps = Sigma_eps,
    dimensions = list(nreg = nreg, nobs = nobs, ncoef = ncoef, nz = nz)
  )
  if (het_observed == "linear" && exists("Delta_true")) {
    true_values$Delta <- Delta_true
  }
  if (!is.null(Sigma_true)) true_values$Sigma_true <- Sigma_true

  return(list(regdata = regdata, hessdata = hessdata, Z = Z, true_values = true_values))
}
