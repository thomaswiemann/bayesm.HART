#' Simulate Hierarchical Linear Model Data
#'
#' Generates simulated data from a hierarchical linear model with optional
#' observed heterogeneity (linear or nonlinear via Z covariates).
#'
#' @param nreg Number of cross-sectional units.
#' @param nobs Number of observations per unit.
#' @param nvar Number of X variables (excluding intercept if const=TRUE).
#' @param nz Number of Z variables. Set to 0 for no observed heterogeneity.
#' @param const Logical. Include an intercept in X? Default TRUE.
#' @param het_observed Character. Functional form of observed heterogeneity.
#'   Options: "none", "linear", "step", "friedman".
#' @param target_var_betabar Numeric. Target variance for the first coefficient
#'   of the observed component betabar_i = f(Z_i). Default 1.0.
#' @param target_var_eps Numeric. Target variance for each coefficient in the
#'   unobserved component eps_i. Default 0.5.
#' @param sigma_sq Numeric. Error variance for y_i = X_i * beta_i + e_i. Default 1.0.
#' @param Sigma_func Function of one row of `Z` returning the `ncoef x ncoef`
#'   covariance of `eps_i`. Overrides `target_var_eps` with covariate-dependent
#'   heterogeneity; requires `Z`. Default `NULL` (constant covariance).
#' @param seed Integer. Optional random seed.
#'
#' @return A list containing:
#'   \item{regdata}{List of length nreg. Each element is a list with y, X, XpX, Xpy.}
#'   \item{Z}{nreg x nz matrix of covariates (NULL if nz=0).}
#'   \item{true_values}{List with beta_true, betabar_true, eps_true, sigma_sq,
#'     Delta (if linear), and Sigma_true (nreg x ncoef x ncoef, with
#'     Sigma_func).}
#'
#' @importFrom stats rnorm runif
#' @export
sim_hier_linear <- function(nreg = 100, nobs = 10, nvar = 2, nz = 3,
                             const = TRUE,
                             het_observed = c("none", "linear", "step", "friedman"),
                             target_var_betabar = 1.0,
                             target_var_eps = 0.5,
                             sigma_sq = 1.0,
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

  # betabar_i = f(Z_i)
  betabar_true <- matrix(0, nrow = nreg, ncol = ncoef)

  if (!is.null(Z) && het_observed != "none") {
    raw_k1 <- numeric(nreg)

    if (het_observed == "linear") {
      # Delta: nz x ncoef, only first column is nonzero
      Delta_true <- matrix(0, nrow = nz, ncol = ncoef)
      nz_used <- min(nz, 5)
      fixed_coefs <- c(3, -3, 2, -2, 1)
      if (nz_used > 0) Delta_true[1:nz_used, 1] <- fixed_coefs[1:nz_used]
      betabar_true <- Z %*% Delta_true
      raw_k1 <- betabar_true[, 1]
    } else if (het_observed == "step") {
      for (i in 1:nreg) {
        raw_k1[i] <- ifelse(Z[i, 1] > 0, 1, -1)
      }
    } else if (het_observed == "friedman") {
      if (nz < 5) stop("Friedman function requires nz >= 5.")
      for (i in 1:nreg) {
        Zi <- Z[i, ]
        term1 <- 10 * sin(pi * Zi[1] * Zi[2])
        term2 <- 20 * (Zi[3] - 0.5)^2
        term3 <- 10 * Zi[4]
        term4 <- 5 * Zi[5]
        raw_k1[i] <- (term1 + term2 + term3 + term4) / 8
      }
    }

    # First coefficient scaled to the target variance and de-meaned.
    raw_var <- var(raw_k1)
    if (raw_var > 1e-10 && target_var_betabar > 0) {
      scale_factor <- sqrt(target_var_betabar / raw_var)
      raw_k1 <- (raw_k1 - mean(raw_k1)) * scale_factor
    } else {
      raw_k1 <- rep(0, nreg)
    }

    betabar_true[, 1] <- raw_k1
    # Remaining coefficients fixed.
    if (ncoef > 1) {
      fixed_remaining <- c(-1, 1, -1, 1)
      n_remaining <- min(ncoef - 1, length(fixed_remaining))
      for (k in 2:(1 + n_remaining)) {
        betabar_true[, k] <- fixed_remaining[k - 1]
      }
    }
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

  # beta_i = betabar_i + eps_i
  beta_true <- betabar_true + eps_true

  # regdata
  regdata <- vector("list", nreg)
  for (i in 1:nreg) {
    Xi <- matrix(rnorm(nobs * nvar), nrow = nobs, ncol = nvar)
    if (const) Xi <- cbind(1, Xi)
    yi <- as.vector(Xi %*% beta_true[i, ] + rnorm(nobs, sd = sqrt(sigma_sq)))
    regdata[[i]] <- list(
      y = yi,
      X = Xi,
      XpX = crossprod(Xi),
      Xpy = crossprod(Xi, yi)
    )
  }

  # true_values
  true_values <- list(
    beta_true = beta_true,
    betabar_true = betabar_true,
    eps_true = eps_true,
    sigma_sq = sigma_sq,
    Sigma_eps = Sigma_eps,
    dimensions = list(nreg = nreg, nobs = nobs, ncoef = ncoef, nz = nz)
  )
  if (het_observed == "linear" && exists("Delta_true")) {
    true_values$Delta <- Delta_true
  }
  if (!is.null(Sigma_true)) true_values$Sigma_true <- Sigma_true

  result <- list(
    regdata = regdata,
    Z = Z,
    true_values = true_values
  )

  return(result)
}
