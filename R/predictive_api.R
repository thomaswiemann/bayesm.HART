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
    if (is.null(newdata$p)) stop("newdata$p must be provided.")
    return(.hart_predict_posterior(
      object, newdata$nlgtdata, kept_draws_indices,
      .hart_quantity(object, "mnl", p = newdata$p),
      "posterior probabilities", "newdata$nlgtdata", r_verbose))
  }

  hetercov_comps <- .hart_prepare_components(
    object, newdata$Z, ndraws_total, r_verbose, force_tree_eval,
    "prior_probs prediction")

  delta_z <- .predict_structural_common(
    object, newdata, "DeltaZ", burn, r_verbose,
    force_tree_eval = force_tree_eval,
    hetercov_comps = hetercov_comps, ...
  )
  if (is.null(newdata$p)) stop("newdata$p must be provided.")
  .hart_predict_prior(object, newdata, delta_z, kept_draws_indices, nsim,
                      .hart_quantity(object, "mnl", p = newdata$p),
                      reduce = "quantity", label = "prior probabilities",
                      what = "newdata$X", r_verbose = r_verbose,
                      hetercov_comps = hetercov_comps,
                      force_tree_eval = force_tree_eval)
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
    return(.hart_predict_posterior(
      object, newdata$regdata, kept_draws_indices,
      .hart_quantity(object, "linear"),
      "posterior responses", "newdata$regdata", r_verbose))
  }

  if (is.null(newdata$X) || !is.list(newdata$X)) {
    stop("Linear prior prediction requires newdata$X as a list of design matrices.")
  }

  hetercov_comps <- .hart_prepare_components(
    object, newdata$Z, ndraws_total, r_verbose, force_tree_eval,
    "Linear prior response prediction")

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

  .hart_predict_prior(
    object, newdata, delta_z, kept_draws_indices, nsim,
    .hart_quantity(object, "linear",
                   tau_ref = rowMeans(object$taudraw[kept_draws_indices, ,
                                                     drop = FALSE])),
    reduce = "beta", label = "linear prior responses", what = "newdata$X",
    r_verbose = r_verbose, hetercov_comps = hetercov_comps,
    force_tree_eval = force_tree_eval)
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
    return(.hart_predict_posterior(
      object, newdata$regdata, kept_draws_indices,
      .hart_quantity(object, "negbin"),
      "posterior responses", "newdata$regdata", r_verbose))
  }

  if (is.null(newdata$X) || !is.list(newdata$X)) {
    stop("NegBin prior prediction requires newdata$X as a list of design matrices.")
  }

  hetercov_comps <- .hart_prepare_components(
    object, newdata$Z, ndraws_total, r_verbose, force_tree_eval,
    "NegBin prior response prediction")

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

  .hart_predict_prior(
    object, newdata, delta_z, kept_draws_indices, nsim,
    .hart_quantity(object, "negbin"), reduce = "beta",
    label = "negbin prior responses", what = "newdata$X",
    r_verbose = r_verbose, hetercov_comps = hetercov_comps,
    force_tree_eval = force_tree_eval)
}
