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
    
    max_utils_rep_mat <- matrix(max_utils_is, nrow = T_i, ncol = p, byrow = FALSE)
    
    exp_utils_is <- exp(utils_mat_is - max_utils_rep_mat)
    row_sums_is <- rowSums(exp_utils_is)
    
    probs_mat_is <- exp_utils_is / row_sums_is
            
    return(probs_mat_is)
}#CALCULATE_MNL_PROBS_FROM_BETA

# Null-coalescing operator: a if !is.null(a) else b.  No roxygen block because
# `%` is not a legal Rd name.
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}#%||%

# Build exact unique-row mapping for Z matrices.
# Returns:
#   - Z_unique: unique rows in first-appearance order
#   - z_index:  length nrow(Z), 1-based map into Z_unique
#   - z_key:    stable row keys aligned with Z_unique (for fast matching)
.build_unique_z_map <- function(Z) {
  if (!is.matrix(Z)) stop("Internal error: Z must be a matrix.")
  if (nrow(Z) < 1L) stop("Internal error: Z must have at least one row.")

  # 17 significant digits reproduce a double exactly.
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
    z_index  = as.integer(z_index),
    z_key    = z_key,
    z_counts = as.integer(tabulate(z_index, nbins = nrow(Z_unique)))
  )
}

# Move the DeltaZ/SigmaZ draws from the C++ output into the hart_cache
# attribute (SigmaZ reshaped to 4D) together with the z_map metadata.
.assemble_hart_cache <- function(draws, z_map) {
  cache <- attr(draws, "hart_cache", exact = TRUE)
  if (is.null(cache)) cache <- list()

  if (!is.null(z_map)) {
    cache$Z_unique <- z_map$Z_unique
    cache$z_index  <- z_map$z_index
    cache$z_key    <- z_map$z_key
    cache$z_counts <- z_map$z_counts
  }

  if (!is.null(draws$DeltaZ_unique_draws)) {
    cache$DeltaZ_unique_draws <- draws$DeltaZ_unique_draws
    draws$DeltaZ_unique_draws <- NULL
  }

  if (!is.null(draws$SigmaZ_unique_draws_flat) && !is.null(cache$Z_unique)) {
    ndraws   <- dim(draws$betadraw)[3]
    nvar_fit <- dim(draws$betadraw)[2]
    n_unique <- nrow(cache$Z_unique)
    cache$SigmaZ_unique_draws <- array(draws$SigmaZ_unique_draws_flat,
      dim = c(n_unique, nvar_fit, nvar_fit, ndraws))
    draws$SigmaZ_unique_draws_flat <- NULL
  }

  if (length(cache) > 0L) attr(draws, "hart_cache") <- cache
  draws
}#.ASSEMBLE_HART_CACHE

# Shared output finishing for the three samplers: bayesm classes and `mcpar`
# attributes, the HART metadata field, and the unique-Z cache attribute.
# Likelihood-specific attributes (taudraw, alphadraw) are set by the wrapper
# before calling this.
#
# `cfg` is the .hart_config() list (with vartree/phitree parsed), `z_map` the
# .build_unique_z_map() result or NULL, `base_class` the bayesm class name.
.finish_hart_fit <- function(draws, cfg, z_map, R, keep, base_class) {
  attributes(draws$betadraw)$class <- "bayesm.hcoef"
  if (cfg$model == "direct") {
    # "bayesm.HART.HeterCov" is the marker superclass shared by all
    # heteroscedastic-Sigma(Z) fits; predict_helpers.R dispatches on it.
    attributes(draws$mu_draw)$class <- c("bayesm.mat", "mcmc")
    attributes(draws$mu_draw)$mcpar <- c(1, R, keep)
    class(draws) <- c(paste0(base_class, "HeterCov"),
                      "bayesm.HART.HeterCov", base_class)
  } else {
    if (!is.null(draws$Deltadraw)) {
      attributes(draws$Deltadraw)$class <- c("bayesm.mat", "mcmc")
      attributes(draws$Deltadraw)$mcpar <- c(1, R, keep)
    }
    attributes(draws$nmix)$class <- "bayesm.nmix"
    class(draws) <- base_class
  }
  # HART metadata; bayesm-path fits carry none.
  if (cfg$model != "none") {
    draws$hart <- list(model      = cfg$model,
                       schema     = 1L,
                       nvar       = cfg$nvar,
                       nz         = cfg$nz,
                       coef_names = cfg$coef_names,
                       priors     = cfg[c("bart", "vartree", "phitree",
                                          "Amu", "store_trees")])
  }
  .assemble_hart_cache(draws, z_map)
}#.FINISH_HART_FIT
