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

# Null Coalescing Operator (internal)
#
# Provides a default value if the left-hand side expression evaluates to NULL.
# Behaves like `rlang::%||%`.  Not exported and not documented via roxygen --
# `\name` cannot legally contain `%`, `|`, or `@`, so an Rd file would fail
# `R CMD check` (see WARN "checking Rd files").
#
# Usage: a %||% b  ->  a if !is.null(a) else b
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}#%||%

# Parse BART prior parameters, supplying defaults
# @param bart_list Prior$bart list
# @param nz number of covariates (used for default rho)
# @return fully populated list of parameters
# @keywords internal
.parse_bart_params <- function(bart_list, nz) {
  list(
    num_trees = bart_list$num_trees %||% 200,
    power     = bart_list$power %||% 2.0,
    base      = bart_list$base %||% 0.95,
    tau       = bart_list$tau %||% (1.0 / sqrt(bart_list$num_trees %||% 200)),
    numcut    = bart_list$numcut %||% 100,
    sparse    = bart_list$sparse %||% FALSE,
    theta     = bart_list$theta %||% 0.0,
    omega     = bart_list$omega %||% 1.0,
    a         = bart_list$a %||% 0.5,
    b         = bart_list$b %||% 1.0,
    rho       = bart_list$rho %||% nz,
    aug       = bart_list$aug %||% FALSE,
    burn      = bart_list$burn %||% 100
  )
}#.parse_bart_params

# Build exact unique-row mapping for Z matrices.
# Returns:
#   - Z_unique: unique rows in first-appearance order
#   - z_index:  length nrow(Z), 1-based map into Z_unique
#   - z_key:    stable row keys aligned with Z_unique (for fast matching)
.build_unique_z_map <- function(Z) {
  if (!is.matrix(Z)) stop("Internal error: Z must be a matrix.")
  if (nrow(Z) < 1L) stop("Internal error: Z must have at least one row.")

  # 17 significant digits is sufficient to preserve exact double identity for
  # values generated/retained within the same R session object.
  row_key <- apply(Z, 1, function(r) paste(sprintf("%.17g", r), collapse = "\r"))
  keep_unique <- !duplicated(row_key)

  Z_unique <- Z[keep_unique, , drop = FALSE]
  z_key <- row_key[keep_unique]
  z_index <- match(row_key, z_key)
  if (anyNA(z_index)) {
    stop("Internal error: failed to build z_index for unique-row map.")
  }

  list(
    Z_unique = Z_unique,
    z_index = as.integer(z_index),
    z_key = z_key
  )
}
