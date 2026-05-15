# ==============================================================================
# Shared Predictive Core
# ==============================================================================
# Internal utilities shared by model-specific predict methods.

# Normalize and validate predictive mode for unified predictive API work.
.resolve_predict_mode <- function(mode = "coefficients") {
  # Backward-compatible alias: "structural" -> "coefficients".
  if (identical(mode, "structural")) mode <- "coefficients"
  valid_modes <- c("coefficients", "prior", "posterior")
  if (!(mode %in% valid_modes)) {
    stop(paste("Invalid mode specified. Choose one of:",
               paste(valid_modes, collapse = ", ")))
  }
  mode
}

.kept_draw_indices <- function(ndraws_total, burn) {
  if (burn > 0) (burn + 1):ndraws_total else seq_len(ndraws_total)
}

# Shared structural prediction path for DeltaZ / DeltaZ+mu / SigmaZ.
.predict_structural_common <- function(object, newdata, type, burn, r_verbose,
                                       force_tree_eval = FALSE,
                                       hetercov_comps = NULL, ...) {
  ndraws_total <- dim(object$betadraw)[3]
  kept_draws_indices <- .kept_draw_indices(ndraws_total, burn)

  if (type == "SigmaZ") {
    return(.calculate_sigma_z(object, newdata, burn, r_verbose,
                              hetercov_comps = hetercov_comps,
                              force_tree_eval = force_tree_eval))
  }

  result <- .calculate_delta_z(object, newdata, burn, r_verbose,
                               hetercov_comps = hetercov_comps,
                               force_tree_eval = force_tree_eval, ...)

  if (type == "DeltaZ+mu") {
    result <- .add_mu_component(result, object, kept_draws_indices)
  }
  result
}
