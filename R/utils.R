# ==============================================================================
# General Helper Functions
# ==============================================================================

# Helper function to calculate MNL probabilities for one unit/draw
# @param X_i [(T_i*p) x nvar] matrix
# @param beta_is [nvar x 1] vector
# @param p scalar, number of alternatives
# @return [T_i x p] matrix of probabilities
# @keywords internal
calculate_mnl_probs_from_beta <- function(X_i, beta_is, p) {
    if (!is.matrix(X_i)) stop("Internal error: X_i must be a matrix")
    if (!is.vector(beta_is)) stop("Internal error: beta_is must be a vector")
    
    T_i <- nrow(X_i) / p
    if(abs(T_i - round(T_i)) > 1e-6 || T_i <= 0) {
        stop("Internal error: nrow(X_i) must be a positive multiple of p")
    }
    T_i <- round(T_i)
    
    utils_vec_is <- X_i %*% beta_is
    utils_mat_is <- matrix(utils_vec_is, nrow = T_i, ncol = p, byrow = TRUE)
    
    max_utils_is <- apply(utils_mat_is, 1, max) # Vector length T_i
    
    # Explicitly replicate the row maximums into a T_i x p matrix
    max_utils_rep_mat <- matrix(max_utils_is, nrow = T_i, ncol = p, byrow = FALSE)
    
    exp_utils_is <- exp(utils_mat_is - max_utils_rep_mat)
    row_sums_is <- rowSums(exp_utils_is)
    
    probs_mat_is <- exp_utils_is / row_sums_is
            
    return(probs_mat_is)
}#CALCULATE_MNL_PROBS_FROM_BETA

#' Null Coalescing Operator
#'
#' Provides a default value if the left-hand side expression evaluates to NULL.
#' Similar to `rlang::%||%`.
#'
#' @param a The expression to check.
#' @param b The default value if `a` is NULL.
#' @return `a` if not NULL, otherwise `b`.
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}#%||%


#' Analyze Inverse-Wishart Prior Parameters
#'
#' Calculates the expectation and variance of diagonal elements for an
#' Inverse-Wishart distribution W^-1(V, nu).
#' Uses defaults from `rhierMnlRwMixture` if `nu` or `V` are NULL.
#'
#' @param nvar Number of variables (dimension of the matrix).
#' @param nu Degrees of freedom (scalar > nvar - 1). Optional.
#' @param V A positive definite scale matrix (nvar x nvar). Optional.
#' @return A list containing the expectation matrix 'E_X' and a vector
#'   'Var_X_diag' with variances of the diagonal elements. Returns NA for
#'   values if conditions on 'nu' are not met.
#' @references \url{https://en.wikipedia.org/wiki/Inverse-Wishart_distribution}
#' @export
analyze_inv_wishart_prior <- function(nvar, nu = NULL, V = NULL) {

  # If nu or V are not provided, use default values from rhierMnlRwMixture
  nu <- nu %||% (nvar + 3)
  V <- V %||% (nu * diag(nvar)) # Default V depends on default nu

  if (!is.numeric(nvar) || length(nvar) != 1 || nvar <= 0 || floor(nvar) != nvar) {
      stop("nvar must be a positive integer.")
  }
  if (!is.matrix(V) || nrow(V) != ncol(V) || nrow(V) != nvar) {
    stop(paste("V must be a square matrix of dimension", nvar))
  }#IF
  if (!all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 1e-8)) {
    stop("V must be positive definite.")
  }#IF
  p <- nvar # Dimension parameter in IW formulas is nvar here
  if (!is.numeric(nu) || length(nu) != 1 || nu <= p - 1) {
    stop(paste("nu must be a scalar greater than p - 1 =", p - 1))
  }#IF

  # Expectation
  E_X <- NA
  if (nu > p + 1) {
    E_X <- V / (nu - p - 1)
  } else {
    warning("Expectation requires nu > p + 1.")
  }#ELSE

  # Variance of diagonal elements
  Var_X_diag <- rep(NA_real_, p)
  if (nu > p + 3) {
    nu_p_1 <- nu - p - 1
    nu_p_3 <- nu - p - 3
    # Denominator per Wikipedia page (Moments section) for Var(X_ij)
    denom <- (nu - p) * (nu - p - 1)^2 * (nu - p - 3)
    # Formula Var(X_ii) = [ (nu-p+1)*V_ii^2 + (nu-p-1)*sum_k(V_ik^2) ] / denom

    for (i in 1:p) {
       Var_X_diag[i] <- ( (nu - p + 1) * V[i,i]^2 + (nu - p - 1) * sum(V[i,]^2) ) / denom
    }#FOR

  } else {
    warning("Variance calculation requires nu > p + 3.")
  }#ELSE

  return(list(E_X = E_X, Var_X_diag = Var_X_diag))
}#ANALYZE_INV_WISHART_PRIOR
