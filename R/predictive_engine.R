# ==============================================================================
# Predictive engine
# ==============================================================================
# One posterior loop and one prior loop for all three likelihoods.  A
# likelihood contributes its per-unit quantity and, for the prior, how the
# nsim simulations are reduced.

# The draw index is the last dimension, so its elements are contiguous.
.quantity_slot <- function(d, len) ((d - 1L) * len + 1L):(d * len)

.check_unit_design <- function(X, i, nvar, what) {
  if (is.null(X) || !is.matrix(X))
    stop(sprintf("%s[[%d]]$X must be a matrix.", what, i))
  if (ncol(X) != nvar)
    stop(sprintf("%s[[%d]]$X has %d columns; expected %d.",
                 what, i, ncol(X), nvar))
  X
}#.CHECK_UNIT_DESIGN

# Per-unit predictive quantity: choice probabilities for the logit, a response
# draw for the regression models.  `ctx` carries the kept-draw index `s`, the
# unit `i`, its position `draw_idx`, and the mode.
.hart_quantity <- function(object, likelihood, p = NULL, tau_ref = NULL) {
  switch(likelihood,
    mnl = function(X_i, beta, ctx)
      calculate_mnl_probs_from_beta(X_i, beta, p),
    linear = function(X_i, beta, ctx) {
      v <- if (identical(ctx$mode, "posterior")) object$taudraw[ctx$s, ctx$i]
           else                                  tau_ref[ctx$draw_idx]
      rnorm(nrow(X_i), mean = as.vector(X_i %*% beta), sd = sqrt(v))
    },
    negbin = function(X_i, beta, ctx)
      rnbinom(nrow(X_i), size = object$alphadraw[ctx$s],
              mu = as.vector(exp(X_i %*% beta))),
    stop(sprintf("No predictive quantity for likelihood '%s'.", likelihood)))
}#.HART_QUANTITY

# Posterior predictive for the estimation units, each at its own betadraw.
.hart_predict_posterior <- function(object, units, kept, quantity, label,
                                    what, r_verbose) {
  n    <- dim(object$betadraw)[1]
  nvar <- dim(object$betadraw)[2]
  if (is.null(units) || !is.list(units))
    stop(sprintf("%s must be provided as a list.", what))
  out <- vector("list", n)
  names(out) <- paste0("Unit_", seq_len(n))
  for (i in seq_len(n)) {
    if (r_verbose)
      cat("Calculating", label, "for Unit", i, "of", n, "...\n")
    X_i <- .check_unit_design(units[[i]]$X, i, nvar, what)
    arr <- NULL
    for (d in seq_along(kept)) {
      q <- quantity(X_i, object$betadraw[i, , kept[d], drop = TRUE],
                    list(s = kept[d], i = i, draw_idx = d,
                         mode = "posterior"))
      if (is.null(arr))
        arr <- array(0, dim = c(dim(q) %||% length(q), length(kept)))
      arr[.quantity_slot(d, length(q))] <- q
    }
    out[[i]] <- arr
  }
  if (r_verbose) cat(label, "complete.\n")
  out
}#.HART_PREDICT_POSTERIOR

# Prior predictive at newdata: beta = Delta(Z_i) + mu + Sigma(Z_i)^{1/2} z,
# nsim simulations per prediction unit and kept draw.  `reduce` is "quantity"
# (average the per-simulation quantity) or "beta" (average the simulated
# coefficients, then evaluate once).
.hart_predict_prior <- function(object, newdata, delta_z_array, kept, nsim,
                                quantity, reduce, label, what, r_verbose,
                                hetercov_comps = NULL,
                                force_tree_eval = FALSE) {
  nvar <- dim(object$betadraw)[2]
  if (is.null(newdata$X) || !is.list(newdata$X))
    stop(sprintf("newdata$X must be provided as a list for %s.", label))
  npred <- length(newdata$X)
  nd    <- length(kept)
  if (dim(delta_z_array)[3] != nd)
    stop(paste("delta_z_array draw dimension does not match",
               "kept_draws_indices length."))
  if (nrow(delta_z_array) == 1L && npred > 1L)
    delta_z_array <- delta_z_array[rep(1L, npred), , , drop = FALSE]
  if (nrow(delta_z_array) != npred)
    stop("delta_z_array unit dimension must match length(newdata$X).")
  if (!is.null(newdata$Z) && nrow(newdata$Z) != npred)
    stop("nrow(newdata$Z) must equal length(newdata$X).")

  mu_kept <- .hart_mu(object, kept)
  Sigma   <- .hart_sigma(object, newdata$Z, kept,
                         force_tree_eval = force_tree_eval,
                         r_verbose = r_verbose, comps = hetercov_comps)
  n_sig <- dim(Sigma)[1]

  out <- vector("list", npred)
  names(out) <- paste0("PredUnit_", seq_len(npred))
  for (i in seq_len(npred)) {
    if (r_verbose)
      cat("Simulating", label, "for Prediction Unit", i, "of", npred, "...\n")
    X_i <- .check_unit_design(newdata$X[[i]], i, nvar, "newdata$X")
    arr <- NULL
    for (d in seq_len(nd)) {
      ctx  <- list(s = kept[d], i = i, draw_idx = d, mode = "prior")
      eta  <- .hart_prior_eta(mu_kept[d, ], Sigma[min(i, n_sig), , , d], nsim)
      beta <- eta + matrix(delta_z_array[i, , d], nsim, nvar, byrow = TRUE)
      q <- if (identical(reduce, "beta")) {
        quantity(X_i, colMeans(beta), ctx)
      } else {
        acc <- quantity(X_i, beta[1, ], ctx)
        if (nsim > 1L)
          for (k in 2:nsim) acc <- acc + quantity(X_i, beta[k, ], ctx)
        acc / nsim
      }
      if (is.null(arr)) arr <- array(0, dim = c(dim(q) %||% length(q), nd))
      arr[.quantity_slot(d, length(q))] <- q
    }
    out[[i]] <- arr
  }
  if (r_verbose) cat(label, "complete.\n")
  out
}#.HART_PREDICT_PRIOR
