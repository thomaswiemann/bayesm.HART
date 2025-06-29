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

