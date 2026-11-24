# ==============================================================================
# Schema accessors
# ==============================================================================
# The one place that knows where a fit keeps mu, Delta(Z) and Sigma(Z).  The
# predict helpers, the predictive engines and the summary methods read the
# draws through these rather than reaching into fields of their own.

# Mean model: "none" (constant mu), "linear" (Deltadraw), "scalar" (one
# ensemble per whitened coordinate), "vector" (one vector-leaf ensemble).
.hart_mean_kind <- function(object) {
  if (!is.null(object$bart_models$treedraws)) "vector"
  else if (!is.null(object$bart_models))      "scalar"
  else if (!is.null(object$Deltadraw))        "linear"
  else                                        "none"
}#.HART_MEAN_KIND

# Covariance model: "heteroskedastic" (covariance ensembles) or "constant".
.hart_cov_kind <- function(object) {
  if (inherits(object, "bayesm.HART.HeterCov")) "heteroskedastic"
  else                                          "constant"
}#.HART_COV_KIND

# mu draws, nkept x nvar.
.hart_mu <- function(object, kept) {
  if (!is.null(object$mu_draw)) return(object$mu_draw[kept, , drop = FALSE])
  .assert_single_component_nmix(object, "mu draws")
  if (is.null(object$nmix) || is.null(object$nmix$compdraw))
    stop("mu draws require mu_draw or nmix$compdraw.")
  nvar <- dim(object$betadraw)[2]
  t(vapply(object$nmix$compdraw[kept], function(x) {
    if (is.null(x[[1]]$mu)) stop("mu missing in component draw.")
    as.numeric(x[[1]]$mu)
  }, numeric(nvar)))
}#.HART_MU

# Sigma(Z) draws, npred x nvar x nvar x nkept.  A constant Sigma is the same
# matrix on every row of Z.
.hart_sigma <- function(object, Z, kept, force_tree_eval = FALSE,
                        r_verbose = FALSE, comps = NULL) {
  nvar  <- dim(object$betadraw)[2]
  npred <- if (is.null(Z)) 1L else nrow(Z)
  if (.hart_cov_kind(object) == "constant") {
    .assert_single_component_nmix(object, "SigmaZ prediction")
    if (is.null(object$nmix) || is.null(object$nmix$compdraw))
      stop(paste("SigmaZ requires mu_draw with covariance trees,",
                 "or nmix$compdraw."))
    out <- array(0, dim = c(npred, nvar, nvar, length(kept)))
    for (s in seq_along(kept)) {
      rooti <- object$nmix$compdraw[[kept[s]]][[1]]$rooti
      if (is.null(rooti))
        stop(sprintf("Missing nmix$compdraw[[%d]][[1]]$rooti", kept[s]))
      A <- solve(rooti)                 # Sigma = A^T A, since R R^T = Sigma^-1
      out[, , , s] <- rep(t(A) %*% A, each = npred)
    }
    return(out)
  }
  if (is.null(Z)) stop("Heter-cov predictions require newdata$Z.")

  cache <- .cache_get(object, "SigmaZ_unique_draws")
  if (!force_tree_eval && !is.null(cache)) {
    map <- .match_cached_unique_rows(object, Z)
    if (!is.null(map) && map$all_seen)
      return(cache[map$idx, , , kept, drop = FALSE])
    if (!is.null(map) && !map$all_seen && !.has_tree_payload(object))
      stop(paste("Out-of-sample Z prediction requires stored tree objects;",
                 "refit with store_trees=TRUE."))
  }
  nall  <- dim(object$betadraw)[3]
  comps <- comps %||% .hart_components(object, Z, npred, nvar, nall, r_verbose)
  .sigma_from_components(comps, npred, nvar, kept)
}#.HART_SIGMA

# Delta(Z) draws, npred x nvar x nkept.
.hart_delta <- function(object, Z, kept, force_tree_eval = FALSE,
                        r_verbose = FALSE, comps = NULL, ...) {
  nvar <- dim(object$betadraw)[2]
  nall <- dim(object$betadraw)[3]
  kind <- .hart_mean_kind(object)
  has_cache <- !is.null(.cache_get(object, "DeltaZ_unique_draws"))
  if (kind == "none" && !has_cache)
    return(array(0, dim = c(1L, nvar, length(kept))))
  if (is.null(Z)) stop("Model was fit with Z. 'newdata$Z' must be provided.")
  npred <- nrow(Z)

  if (!force_tree_eval && has_cache) {
    map <- .match_cached_unique_rows(object, Z)
    if (!is.null(map) && map$all_seen)
      return(.cache_get(object, "DeltaZ_unique_draws")[map$idx, , kept,
                                                       drop = FALSE])
    if (!is.null(map) && !map$all_seen && !.has_tree_payload(object))
      stop(paste("Out-of-sample Z prediction requires stored tree objects;",
                 "refit with store_trees=TRUE."))
  }

  pred <- switch(kind,
    linear = .delta_z_linear(object, Z, npred, nvar, ncol(Z), nall),
    scalar = .delta_z_scalar(object, Z, npred, nvar, nall, r_verbose, ...),
    vector = .delta_z_vector(object, Z, nvar, comps),
    stop("Delta(Z) prediction requires Deltadraw, bart_models, or a cache."))
  pred[, , kept, drop = FALSE]
}#.HART_DELTA

# --- per-kind evaluation ------------------------------------------------------

.delta_z_linear <- function(object, Z, npred, nvar, nz, ndraws_total) {
  if (ncol(object$Deltadraw) != nz * nvar)
    stop(sprintf("Dimension mismatch: ncol(Deltadraw)=%d, expected nz*nvar=%d",
                 ncol(object$Deltadraw), nz * nvar))
  pred <- array(0, dim = c(npred, nvar, ndraws_total))
  for (i in seq_len(ndraws_total)) {
    Delta_draw <- matrix(object$Deltadraw[i, ], nrow = nz, byrow = TRUE)
    pred[, , i] <- Z %*% Delta_draw
  }
  pred
}#.DELTA_Z_LINEAR

# Normalized HART: Delta(Z*) = Sigma_r^{1/2} f*(Z*), with Sigma_r^{1/2} from
# nmix$compdraw[[r]][[1]]$rooti.
.delta_z_scalar <- function(object, Z, npred, nvar, ndraws_total,
                            r_verbose, ...) {
  .assert_single_component_nmix(object, "DeltaZ BART prediction")
  if (length(object$bart_models) != nvar)
    stop(sprintf("Number of BART models (%d) does not match nvar (%d)",
                 length(object$bart_models), nvar))
  raw_pred <- array(0, dim = c(npred, nvar, ndraws_total))
  passed_args <- list(...)
  pw_arg_names <- names(formals(.pw_trees))
  valid_passed_args <- passed_args[names(passed_args) %in% pw_arg_names]
  base_args <- list(x.test = Z, mu = 0, transposed = FALSE, dodraws = TRUE)
  final_args_template <- utils::modifyList(base_args, valid_passed_args)
  final_args_template <-
    final_args_template[names(final_args_template) %in% pw_arg_names]
  for (j in seq_len(nvar)) {
    if (is.null(object$bart_models[[j]]$treedraws))
      stop(paste("Missing treedraws for BART model of coefficient", j))
    if (r_verbose) cat("Predicting coefficient", j, "with BART model\n")
    args_j <- final_args_template
    args_j$treedraws <- object$bart_models[[j]]$treedraws
    raw_pred[, j, ] <- t(do.call(.pw_trees, args_j))
  }
  if (is.null(object$nmix) || is.null(object$nmix$compdraw))
    stop("Missing nmix$compdraw needed for BART scaling.")
  pred <- array(0, dim = c(npred, nvar, ndraws_total))
  for (s in seq_len(ndraws_total)) {
    rooti <- object$nmix$compdraw[[s]][[1]]$rooti
    if (is.null(rooti))
      stop(paste0("Missing nmix$compdraw[[", s, "]][[1]]$rooti"))
    root_inv <- tryCatch(
      solve(rooti),
      error = function(e)
        stop(sprintf(paste("DeltaZ BART prediction failed at draw %d while",
                           "inverting nmix rooti: %s"), s, e$message))
    )
    pred[, , s] <- raw_pred[, , s] %*% root_inv
  }
  pred
}#.DELTA_Z_SCALAR

# Direct HART: Delta(Z*) = f(Z*) from the vector-leaf ensemble.
.delta_z_vector <- function(object, Z, nvar, comps = NULL) {
  if (!is.null(comps)) return(comps$delta_arr)
  td <- object$bart_models$treedraws
  if (!identical(as.integer(td$leaf_dim), as.integer(nvar)))
    stop(sprintf("Vector-leaf ensemble needs leaf_dim %d.", nvar))
  .pw_vec_trees(Z, td, nvar)
}#.DELTA_Z_VECTOR

# --- covariance ensembles -----------------------------------------------------

# Evaluate the direct-model ensembles at Z: delta_arr = f(Z) from the
# vector-leaf trees (npred x nvar x ndraws), d_arr from the variance trees,
# and phi_arr (npred x nvar x nvar x ndraws, NULL when diagonal) from the
# off-diagonal trees.  One pass shared by the DeltaZ, SigmaZ and prior
# predictive paths.
.hart_components <- function(object, Z, npred, nvar, ndraws_total,
                             r_verbose = FALSE) {
  td <- object$bart_models$treedraws
  if (is.null(td) || !identical(as.integer(td$leaf_dim), as.integer(nvar)))
    stop(sprintf("Direct HART object needs bart_models with leaf_dim %d.",
                 nvar))
  if (is.null(object$var_models) || length(object$var_models) != nvar)
    stop(sprintf(paste("Heter-cov object missing var_models or wrong length",
                       "(expected %d)."), nvar))
  use_full <- !is.null(object$phi_models)
  if (use_full && length(object$phi_models) != nvar)
    stop(sprintf("phi_models length (%d) != nvar (%d)",
                 length(object$phi_models), nvar))

  if (r_verbose) cat("Predicting f(Z) (mean trees)\n")
  delta_arr <- .pw_vec_trees(Z, td, nvar)
  d_arr     <- array(0, dim = c(npred, nvar, ndraws_total))
  for (j in seq_len(nvar)) {
    if (r_verbose) cat("Predicting d_", j, " (var tree)\n", sep = "")
    d_arr[, j, ] <- t(.pw_trees(
      x.test = Z, treedraws = object$var_models[[j]]$treedraws,
      product = TRUE))
  }

  phi_arr <- NULL
  if (use_full) {
    phi_arr <- array(0, dim = c(npred, nvar, nvar, ndraws_total))
    for (j in 2:nvar) {
      pj <- object$phi_models[[j]]
      if (is.null(pj) || length(pj) != j - 1)
        stop(sprintf("phi_models[[%d]] missing or wrong length (expected %d).",
                     j, j - 1))
      for (k in seq_len(j - 1)) {
        if (r_verbose)
          cat("Predicting phi_", j, k, " (off-diag tree)\n", sep = "")
        phi_arr[, j, k, ] <- t(.pw_trees(
          x.test = Z, treedraws = pj[[k]]$treedraws))
      }
    }
  }
  list(delta_arr = delta_arr, d_arr = d_arr, phi_arr = phi_arr,
       use_full = use_full)
}#.HART_COMPONENTS

# Sigma(Z*) = L(Z*)^{-1} D(Z*) L(Z*)^{-T} from the modified-Cholesky pieces.
.sigma_from_components <- function(comps, npred, nvar, kept) {
  out <- array(0, dim = c(npred, nvar, nvar, length(kept)))
  for (out_idx in seq_along(kept)) {
    s <- kept[out_idx]
    for (i in seq_len(npred)) {
      d_vec <- pmax(comps$d_arr[i, , s], .Machine$double.eps)
      L <- diag(nvar)
      if (comps$use_full && nvar > 1L)
        for (j in 2:nvar)
          L[j, seq_len(j - 1L)] <- -comps$phi_arr[i, j, seq_len(j - 1L), s]
      if (any(!is.finite(L)))
        stop(sprintf(paste("SigmaZ prediction failed at draw %d, unit %d:",
                           "non-finite Cholesky factor entries."), s, i))
      A <- tryCatch(
        solve(L, diag(sqrt(d_vec), nvar)),
        error = function(e)
          stop(sprintf(paste("SigmaZ prediction failed at draw %d, unit %d",
                             "while solving L^{-1}D^{1/2}: %s"),
                       s, i, e$message))
      )
      out[i, , , out_idx] <- A %*% t(A)
    }
  }
  out
}#.SIGMA_FROM_COMPONENTS

# Shared tree evaluation for the predictive paths: NULL when the unique-Z
# cache covers Z or the fit has no covariance ensembles, the evaluated
# ensembles otherwise.  `context` names the caller in the error message.
.hart_prepare_components <- function(object, Z, ndraws_total, r_verbose,
                                     force_tree_eval, context) {
  if (.hart_cov_kind(object) != "heteroskedastic") return(NULL)
  if (is.null(Z)) stop("Heter-cov predictions require newdata$Z.")
  map      <- .match_cached_unique_rows(object, Z)
  covered  <- !is.null(.cache_get(object, "SigmaZ_unique_draws")) &&
              !is.null(map) && map$all_seen
  has_tree <- .has_tree_payload(object)
  if (!has_tree && !covered)
    stop(sprintf(paste("%s requires stored tree objects for unseen Z;",
                       "refit with store_trees=TRUE."), context))
  if (!has_tree || !(force_tree_eval || !covered)) return(NULL)
  .hart_components(object, Z, nrow(Z), dim(object$betadraw)[2], ndraws_total,
                   r_verbose)
}#.HART_PREPARE_COMPONENTS

# nsim prior deviations eta = mu + Sigma^{1/2} z at one unit and kept draw,
# as an nsim x nvar matrix.
.hart_prior_eta <- function(mu_s, Sigma_is, nsim) {
  nvar <- length(mu_s)
  U <- tryCatch(chol(Sigma_is), error = function(e) NULL)
  if (is.null(U))
    stop("Prior prediction failed: Sigma(Z) is not positive definite.")
  matrix(rnorm(nsim * nvar), nsim, nvar) %*% U +
    matrix(mu_s, nsim, nvar, byrow = TRUE)
}#.HART_PRIOR_ETA
