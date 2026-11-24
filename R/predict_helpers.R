# ==============================================================================
# Shared Predict Helpers for Hierarchical Models
# ==============================================================================
# These helpers are model-agnostic and used by predict methods for
# rhierMnlRwMixture, rhierLinearMixture, and rhierNegbinRw.

# Checks type, betadraw, burn and nsim; returns the total number of draws.
.validate_predict_inputs <- function(object, newdata, type, burn, nsim,
                                      valid_types = c("DeltaZ", "DeltaZ+mu")) {
  if (!(type %in% valid_types)) {
    stop(paste("Invalid type specified. Choose one of:",
               paste(valid_types, collapse = ", ")))
  }
  if (is.null(object$betadraw) || !is.array(object$betadraw) ||
      length(dim(object$betadraw)) != 3) {
    stop("Invalid object: object$betadraw missing or not a 3D array.")
  }

  ndraws_total <- dim(object$betadraw)[3]

  if (burn < 0 || burn >= ndraws_total) {
    stop(paste0("burn must be >= 0 and < ndraws_total (", ndraws_total, ")."))
  }
  if (startsWith(type, "prior_") && nsim <= 0) {
    stop("nsim must be > 0 for prior predictive mode.")
  }

  return(ndraws_total)
}#VALIDATE_PREDICT_INPUTS

# Global-mixture prediction requires a single component.
.assert_single_component_nmix <- function(object, context) {
  if (.hart_cov_kind(object) == "heteroskedastic") return(invisible(NULL))
  if (is.null(object$nmix) || is.null(object$nmix$probdraw)) return(invisible(NULL))
  ncomp <- ncol(as.matrix(object$nmix$probdraw))
  if (is.na(ncomp) || ncomp != 1) {
    stop(sprintf("%s supports only ncomp = 1. Found ncomp = %s.",
                 context, as.character(ncomp)))
  }
  invisible(NULL)
}

.row_keys_from_Z <- function(Z) {
  apply(Z, 1, function(r) paste(sprintf("%.17g", r), collapse = "\r"))
}

.cache_get <- function(object, name) {
  if (!is.null(object[[name]])) return(object[[name]])
  cache <- attr(object, "hart_cache", exact = TRUE)
  if (!is.null(cache) && !is.null(cache[[name]])) return(cache[[name]])
  NULL
}

.match_cached_unique_rows <- function(object, Z) {
  Z_unique <- .cache_get(object, "Z_unique")
  z_index <- .cache_get(object, "z_index")
  if (is.null(Z_unique) || is.null(z_index)) return(NULL)
  z_key_unique <- .cache_get(object, "z_key")
  if (is.null(z_key_unique)) z_key_unique <- .row_keys_from_Z(Z_unique)
  idx <- match(.row_keys_from_Z(Z), z_key_unique)
  list(idx = idx, all_seen = all(!is.na(idx)))
}

.has_tree_payload <- function(object) {
  if (.hart_cov_kind(object) == "heteroskedastic") {
    return(!is.null(object$bart_models$treedraws) &&
             !is.null(object$var_models))
  }
  !is.null(object$bart_models)
}

# Returns a (npred x nvar x ndraws_kept) array of representative-consumer
# predictions Delta(Z*) at user-supplied newdata$Z; .hart_delta() knows where
# each fit keeps them.
.calculate_delta_z <- function(object, newdata, burn, r_verbose,
                               hetercov_comps = NULL,
                               force_tree_eval = FALSE, ...) {
  ndraws_total <- dim(object$betadraw)[3]
  if (burn >= ndraws_total) stop("Burn-in period too large.")
  kept <- .kept_draw_indices(ndraws_total, burn)
  .hart_delta(object, newdata$Z, kept, force_tree_eval = force_tree_eval,
              r_verbose = r_verbose, comps = hetercov_comps, ...)
}#CALCULATE_DELTA_Z

# Returns a [npred, nvar, nvar, ndraws_out] array.  Under a constant
# covariance the same matrix is repeated on every row of Z.
.calculate_sigma_z <- function(object, newdata, burn, r_verbose,
                               hetercov_comps = NULL,
                               force_tree_eval = FALSE) {
  ndraws_total <- dim(object$betadraw)[3]
  if (burn >= ndraws_total) stop("Burn-in period too large.")
  kept <- .kept_draw_indices(ndraws_total, burn)
  .hart_sigma(object, newdata$Z, kept, force_tree_eval = force_tree_eval,
              r_verbose = r_verbose, comps = hetercov_comps)
}#CALCULATE_SIGMA_Z

# Adds mu to Delta(Z): from nmix$compdraw on the mixture paths, from mu_draw
# under heteroskedastic covariance.
.add_mu_component <- function(delta_z_array, object, kept_draws_indices) {
  npred      <- dim(delta_z_array)[1]
  nvar       <- dim(delta_z_array)[2]
  ndraws_out <- dim(delta_z_array)[3]

  mudraw <- t(.hart_mu(object, kept_draws_indices))        # nvar x ndraws_out

  if (!is.matrix(mudraw) || nrow(mudraw) != nvar || ncol(mudraw) != ndraws_out)
    stop("Extracted mudraw has incorrect dimensions.")

  mudraw_array <- aperm(array(mudraw, dim = c(nvar, ndraws_out, npred)),
                        c(3, 1, 2))
  return(delta_z_array + mudraw_array)
}#ADD_MU_COMPONENT
