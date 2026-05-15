# ==============================================================================
# Predictive API (Phase B scaffold)
# ==============================================================================
# Internal dispatcher for clean prior/posterior predictive interfaces.
# Public rollout can be done after all three model families are wired.

.predictive_dispatch <- function(object, newdata, mode, type,
                                 burn = 0, nsim = 10,
                                 r_verbose = TRUE,
                                 force_tree_eval = FALSE, ...) {
  if (inherits(object, "rhierMnlRwMixture")) {
    return(.predictive_dispatch_mnl(
      object, newdata, mode, type, burn, nsim,
      r_verbose, force_tree_eval, ...
    ))
  }
  if (inherits(object, "rhierLinearMixture")) {
    return(.predictive_dispatch_linear(
      object, newdata, mode, type, burn, nsim,
      r_verbose, force_tree_eval, ...
    ))
  }
  if (inherits(object, "rhierNegbinRw")) {
    return(.predictive_dispatch_negbin(
      object, newdata, mode, type, burn, nsim,
      r_verbose, force_tree_eval, ...
    ))
  }
  stop("No predictive dispatcher implemented for this model class yet.")
}

.predictive_dispatch_mnl <- function(object, newdata, mode, type,
                                     burn = 0, nsim = 10,
                                     r_verbose = TRUE,
                                     force_tree_eval = FALSE, ...) {
  mode <- .resolve_predict_mode(mode)
  if (!(mode %in% c("prior", "posterior"))) {
    stop("MNL predictive dispatcher supports mode in {'prior','posterior'}.")
  }
  if (!identical(type, "choice_probs")) {
    stop("MNL predictive dispatcher currently supports type='choice_probs' only.")
  }

  ndraws_total <- .validate_predict_inputs(
    object, newdata,
    type = if (mode == "prior") "prior_choice_probs" else "posterior_choice_probs",
    burn = burn, nsim = nsim,
    valid_types = c("prior_choice_probs", "posterior_choice_probs")
  )
  kept_draws_indices <- .kept_draw_indices(ndraws_total, burn)

  if (mode == "posterior") {
    return(.predict_posterior_probs(object, newdata, kept_draws_indices, r_verbose))
  }

  hetercov_comps <- NULL
  if (inherits(object, "bayesm.HART.HeterCov")) {
    if (is.null(newdata$Z))
      stop("Heter-cov predictions require newdata$Z.")
    map_seen <- .match_cached_unique_rows(object, newdata$Z)
    seen_all <- !is.null(map_seen) && map_seen$all_seen
    has_sigma_cache <- !is.null(.cache_get(object, "SigmaZ_unique_draws"))

    need_tree_eval <- .has_tree_payload(object) && (
      force_tree_eval || !(has_sigma_cache && seen_all)
    )
    if (!.has_tree_payload(object) && !(has_sigma_cache && seen_all)) {
      stop("prior_probs prediction requires stored tree objects for unseen Z; refit with store_trees=TRUE.")
    }
    if (need_tree_eval) {
      nvar <- dim(object$betadraw)[2]
      npred <- nrow(newdata$Z)
      hetercov_comps <- .hetercov_components(object, newdata$Z, npred, nvar,
                                             ndraws_total, r_verbose)
    }
  }

  delta_z <- .predict_structural_common(
    object, newdata, "DeltaZ", burn, r_verbose,
    force_tree_eval = force_tree_eval,
    hetercov_comps = hetercov_comps, ...
  )
  .predict_prior_probs(object, newdata, delta_z, kept_draws_indices, nsim,
                       r_verbose, hetercov_comps = hetercov_comps)
}

.predictive_dispatch_linear <- function(object, newdata, mode, type,
                                        burn = 0, nsim = 10,
                                        r_verbose = TRUE,
                                        force_tree_eval = FALSE, ...) {
  mode <- .resolve_predict_mode(mode)
  if (!(mode %in% c("prior", "posterior"))) {
    stop("Linear predictive dispatcher supports mode in {'prior','posterior'}.")
  }
  if (!identical(type, "response")) {
    stop("Linear predictive dispatcher currently supports type='response' only.")
  }

  ndraws_total <- .validate_predict_inputs(
    object, newdata,
    type = if (mode == "prior") "prior_linear_response" else "posterior_linear_response",
    burn = burn, nsim = nsim,
    valid_types = c("prior_linear_response", "posterior_linear_response")
  )
  kept_draws_indices <- .kept_draw_indices(ndraws_total, burn)

  if (mode == "posterior") {
    return(.predict_linear_posterior_response(
      object, newdata, kept_draws_indices, r_verbose
    ))
  }

  if (is.null(newdata$X) || !is.list(newdata$X)) {
    stop("Linear prior prediction requires newdata$X as a list of design matrices.")
  }

  hetercov_comps <- NULL
  if (inherits(object, "bayesm.HART.HeterCov")) {
    if (is.null(newdata$Z))
      stop("Heter-cov linear prior prediction requires newdata$Z.")
    map_seen <- .match_cached_unique_rows(object, newdata$Z)
    seen_all <- !is.null(map_seen) && map_seen$all_seen
    has_sigma_cache <- !is.null(.cache_get(object, "SigmaZ_unique_draws"))

    need_tree_eval <- .has_tree_payload(object) && (
      force_tree_eval || !(has_sigma_cache && seen_all)
    )
    if (!.has_tree_payload(object) && !(has_sigma_cache && seen_all)) {
      stop("Linear prior response prediction requires stored tree objects for unseen Z; refit with store_trees=TRUE.")
    }
    if (need_tree_eval) {
      nvar <- dim(object$betadraw)[2]
      npred <- nrow(newdata$Z)
      hetercov_comps <- .hetercov_components(object, newdata$Z, npred, nvar,
                                             ndraws_total, r_verbose)
    }
  }

  has_model_z <- !is.null(object$Deltadraw) ||
    !is.null(object$bart_models) ||
    !is.null(.cache_get(object, "DeltaZ_unique_draws"))
  if (has_model_z && is.null(newdata$Z)) {
    stop("Linear prior prediction requires newdata$Z when the model was fit with Z.")
  }

  z_for_delta <- newdata$Z
  if (is.null(z_for_delta)) {
    z_for_delta <- matrix(0, nrow = 1L, ncol = 1L)
  }
  delta_newdata <- list(Z = z_for_delta)
  delta_z <- .predict_structural_common(
    object, delta_newdata, "DeltaZ", burn, r_verbose,
    force_tree_eval = force_tree_eval,
    hetercov_comps = hetercov_comps, ...
  )

  .predict_linear_prior_response(
    object, newdata, delta_z, kept_draws_indices, nsim, r_verbose,
    hetercov_comps = hetercov_comps
  )
}

.predictive_dispatch_negbin <- function(object, newdata, mode, type,
                                        burn = 0, nsim = 10,
                                        r_verbose = TRUE,
                                        force_tree_eval = FALSE, ...) {
  mode <- .resolve_predict_mode(mode)
  if (!(mode %in% c("prior", "posterior"))) {
    stop("NegBin predictive dispatcher supports mode in {'prior','posterior'}.")
  }
  if (!identical(type, "response")) {
    stop("NegBin predictive dispatcher currently supports type='response' only.")
  }

  ndraws_total <- .validate_predict_inputs(
    object, newdata,
    type = if (mode == "prior") "prior_negbin_response" else "posterior_negbin_response",
    burn = burn, nsim = nsim,
    valid_types = c("prior_negbin_response", "posterior_negbin_response")
  )
  kept_draws_indices <- .kept_draw_indices(ndraws_total, burn)

  if (mode == "posterior") {
    return(.predict_negbin_posterior_response(
      object, newdata, kept_draws_indices, r_verbose
    ))
  }

  if (is.null(newdata$X) || !is.list(newdata$X)) {
    stop("NegBin prior prediction requires newdata$X as a list of design matrices.")
  }

  hetercov_comps <- NULL
  if (inherits(object, "bayesm.HART.HeterCov")) {
    if (is.null(newdata$Z))
      stop("Heter-cov negbin prior prediction requires newdata$Z.")
    map_seen <- .match_cached_unique_rows(object, newdata$Z)
    seen_all <- !is.null(map_seen) && map_seen$all_seen
    has_sigma_cache <- !is.null(.cache_get(object, "SigmaZ_unique_draws"))

    need_tree_eval <- .has_tree_payload(object) && (
      force_tree_eval || !(has_sigma_cache && seen_all)
    )
    if (!.has_tree_payload(object) && !(has_sigma_cache && seen_all)) {
      stop("NegBin prior response prediction requires stored tree objects for unseen Z; refit with store_trees=TRUE.")
    }
    if (need_tree_eval) {
      nvar <- dim(object$betadraw)[2]
      npred <- nrow(newdata$Z)
      hetercov_comps <- .hetercov_components(object, newdata$Z, npred, nvar,
                                             ndraws_total, r_verbose)
    }
  }

  has_model_z <- !is.null(object$Deltadraw) ||
    !is.null(object$bart_models) ||
    !is.null(.cache_get(object, "DeltaZ_unique_draws"))
  if (has_model_z && is.null(newdata$Z)) {
    stop("NegBin prior prediction requires newdata$Z when the model was fit with Z.")
  }

  z_for_delta <- newdata$Z
  if (is.null(z_for_delta)) {
    z_for_delta <- matrix(0, nrow = 1L, ncol = 1L)
  }
  delta_newdata <- list(Z = z_for_delta)
  delta_z <- .predict_structural_common(
    object, delta_newdata, "DeltaZ", burn, r_verbose,
    force_tree_eval = force_tree_eval,
    hetercov_comps = hetercov_comps, ...
  )

  .predict_negbin_prior_response(
    object, newdata, delta_z, kept_draws_indices, nsim, r_verbose,
    hetercov_comps = hetercov_comps
  )
}
