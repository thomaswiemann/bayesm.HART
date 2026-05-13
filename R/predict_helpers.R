# ==============================================================================
# Shared Predict Helpers for Hierarchical Models
# ==============================================================================
# These helpers are model-agnostic and used by predict methods for
# rhierMnlRwMixture, rhierLinearMixture, and rhierNegbinRw.

# Helper 1: Minimal Input Validation
.validate_predict_inputs <- function(object, newdata, type, burn, nsim,
                                      valid_types = c("DeltaZ", "DeltaZ+mu")) {
  # Minimal validation: Type and basic object structure
  if (!(type %in% valid_types)) {
    stop(paste("Invalid type specified. Choose one of:",
               paste(valid_types, collapse = ", ")))
  }
  if (is.null(object$betadraw) || !is.array(object$betadraw) ||
      length(dim(object$betadraw)) != 3) {
    stop("Invalid object: object$betadraw missing or not a 3D array.")
  }

  ndraws_total <- dim(object$betadraw)[3]

  # Minimal validation: burn and nsim (assuming numeric)
  if (burn < 0 || burn >= ndraws_total) {
    stop(paste0("burn must be >= 0 and < ndraws_total (", ndraws_total, ")."))
  }
  if (type == "prior_probs" && nsim <= 0) {
    stop("nsim must be > 0 for type='prior_probs'.")
  }

  return(ndraws_total)
}#VALIDATE_PREDICT_INPUTS

# Helper 2: Calculate DeltaZ component
#
# Returns a (npred x nvar x ndraws_kept) array of representative-consumer
# predictions Delta(Z*) at user-supplied newdata$Z.  Three model paths:
#
#   (a) Linear hierarchical (object$Deltadraw):  Delta(Z*) = Z* %*% Delta_r
#   (b) Sum-of-trees, fixed Sigma (object$bart_models, object$nmix):
#         Delta(Z*) = Sigma_r^{1/2} sum-of-trees(Z*),   with Sigma_r^{1/2}
#         derived from nmix$compdraw[[r]][[1]]$rooti (= Sigma_r^{-1/2}).
#   (c) Heteroscedastic Sigma(Z*) (object inherits "bayesm.HART.HeterCov"):
#         Delta(Z*) = Sigma(Z*)^{1/2} sum-of-trees(Z*) where
#         Sigma(Z*)^{1/2} = L(Z*)^{-1} D(Z*)^{1/2} is rebuilt per draw from
#         var_models (product-of-trees -> d) and phi_models (sum-of-trees -> phi).
.calculate_delta_z <- function(object, newdata, burn, r_verbose,
                               hetercov_comps = NULL, ...) {
  nvar         <- dim(object$betadraw)[2]
  ndraws_total <- dim(object$betadraw)[3]

  is_heter <- inherits(object, "bayesm.HART.HeterCov")
  has_bart <- !is.null(object$bart_models)
  has_lin  <- !is.null(object$Deltadraw)
  model_has_Z <- has_lin || has_bart

  if (!model_has_Z) {
    pred <- array(0, dim = c(1L, nvar, ndraws_total))
  } else {
    if (is.null(newdata$Z))
      stop("Model was fit with Z. 'newdata$Z' must be provided.")
    npred <- nrow(newdata$Z)
    nz    <- ncol(newdata$Z)

    if (has_lin) {
      pred <- .delta_z_linear(object, newdata$Z, npred, nvar, nz, ndraws_total)
    } else if (is_heter) {
      pred <- .delta_z_hetercov(object, newdata$Z, npred, nvar,
                                ndraws_total, r_verbose,
                                hetercov_comps = hetercov_comps)
    } else {
      pred <- .delta_z_bart(object, newdata$Z, npred, nvar, ndraws_total,
                            r_verbose, ...)
    }
  }

  if (burn >= ndraws_total) stop("Burn-in period too large.")
  if (burn > 0) {
    kept_draws <- (burn + 1):ndraws_total
    pred <- pred[, , kept_draws, drop = FALSE]
  }
  return(pred)
}#CALCULATE_DELTA_Z

# --- per-path helpers (split out for clarity) -------------------------------

.delta_z_linear <- function(object, Z, npred, nvar, nz, ndraws_total) {
  if (ncol(object$Deltadraw) != nz * nvar)
    stop(sprintf("Dimension mismatch: ncol(Deltadraw)=%d, expected nz*nvar=%d",
                 ncol(object$Deltadraw), nz * nvar))
  pred <- array(0, dim = c(npred, nvar, ndraws_total))
  for (i in seq_len(ndraws_total)) {
    Delta_draw <- matrix(object$Deltadraw[i, ], nrow = nz, byrow = TRUE)
    pred[, , i] <- Z %*% Delta_draw
  }
  return(pred)
}

.delta_z_bart <- function(object, Z, npred, nvar, ndraws_total,
                          r_verbose, ...) {
  if (length(object$bart_models) != nvar)
    stop(sprintf("Number of BART models (%d) does not match nvar (%d)",
                 length(object$bart_models), nvar))

  raw_pred <- array(0, dim = c(npred, nvar, ndraws_total))

  # --- Prepare arguments for pwbart, handling ... ---
  passed_args <- list(...)
  pwbart_arg_names <- names(formals(pwbart))
  valid_passed_args <- passed_args[names(passed_args) %in% pwbart_arg_names]
  base_args <- list(x.test = Z, mu = 0, transposed = FALSE, dodraws = TRUE)
  final_args_template <- utils::modifyList(base_args, valid_passed_args)
  final_args_template <-
    final_args_template[names(final_args_template) %in% pwbart_arg_names]

  for (j in seq_len(nvar)) {
    if (is.null(object$bart_models[[j]]$treedraws))
      stop(paste("Missing treedraws for BART model of coefficient", j))
    if (r_verbose) cat("Predicting coefficient", j, "with BART model\n")
    args_j <- final_args_template
    args_j$treedraws <- object$bart_models[[j]]$treedraws
    raw_pred[, j, ] <- t(do.call(pwbart, args_j))
  }

  if (is.null(object$nmix) || is.null(object$nmix$compdraw))
    stop("Missing nmix$compdraw needed for BART scaling.")

  pred <- array(0, dim = c(npred, nvar, ndraws_total))
  for (s in seq_len(ndraws_total)) {
    rooti <- object$nmix$compdraw[[s]][[1]]$rooti
    if (is.null(rooti))
      stop(paste0("Missing nmix$compdraw[[", s, "]][[1]]$rooti"))
    root_inv <- tryCatch(solve(rooti), error = function(e)
      stop(paste("rooti matrix is singular for draw", s)))
    pred[, , s] <- raw_pred[, , s] %*% root_inv
  }
  return(pred)
}

# Heter-cov path uses precomputed components from .hetercov_components(); see
# below for the data layout.  Sigma(Z*)^{1/2} = L^{-1} D^{1/2} acts on a vector
# in O(D^2): scale by sqrt(d), then apply unit-lower-triangular inverse
# (L y = x solved by forward substitution; L = I - lower(phi)).
.delta_z_hetercov <- function(object, Z, npred, nvar, ndraws_total, r_verbose,
                              hetercov_comps = NULL) {
  comps <- if (is.null(hetercov_comps))
    .hetercov_components(object, Z, npred, nvar, ndraws_total, r_verbose)
  else
    hetercov_comps
  pred     <- array(0, dim = c(npred, nvar, ndraws_total))
  use_full <- comps$use_full
  for (s in seq_len(ndraws_total)) {
    for (i in seq_len(npred)) {
      x <- comps$delta_arr[i, , s] * sqrt(comps$d_arr[i, , s])  # D^{1/2} delta
      if (use_full) {
        for (j in 2:nvar) {
          x[j] <- x[j] + sum(comps$phi_arr[i, j, seq_len(j - 1), s] *
                             x[seq_len(j - 1)])
        }
      }
      pred[i, , s] <- x
    }
  }
  return(pred)
}

# Evaluate every BART / varBART / phi-BART ensemble in `object` at newdata
# locations Z.  Returns a list of (npred x nvar x ndraws_total) arrays plus a
# (npred x nvar x nvar x ndraws_total) phi_arr (or NULL if diagonal-only).
# Centralized so that .delta_z_hetercov() and the heter-cov prior_probs path
# can share a single tree-evaluation pass.
.hetercov_components <- function(object, Z, npred, nvar, ndraws_total,
                                 r_verbose) {
  if (is.null(object$bart_models) || length(object$bart_models) != nvar)
    stop(sprintf("Heter-cov object missing bart_models or wrong length (expected %d).",
                 nvar))
  if (is.null(object$var_models) || length(object$var_models) != nvar)
    stop(sprintf("Heter-cov object missing var_models or wrong length (expected %d).",
                 nvar))
  use_full <- !is.null(object$phi_models)
  if (use_full && length(object$phi_models) != nvar)
    stop(sprintf("phi_models length (%d) != nvar (%d)",
                 length(object$phi_models), nvar))

  delta_arr <- array(0, dim = c(npred, nvar, ndraws_total))
  d_arr     <- array(0, dim = c(npred, nvar, ndraws_total))
  for (j in seq_len(nvar)) {
    if (r_verbose) cat("Predicting delta_", j, " (mean tree)  / d_", j,
                       " (var tree)\n", sep = "")
    delta_arr[, j, ] <- t(pwbart(
      x.test = Z, treedraws = object$bart_models[[j]]$treedraws,
      mu = 0, transposed = FALSE, dodraws = TRUE))
    d_arr[, j, ] <- t(pwvarbart(
      x.test = Z, treedraws = object$var_models[[j]]$treedraws,
      transposed = FALSE, dodraws = TRUE))
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
        if (r_verbose) cat("Predicting phi_", j, k, " (off-diag tree)\n", sep = "")
        phi_arr[, j, k, ] <- t(pwbart(
          x.test = Z, treedraws = pj[[k]]$treedraws,
          mu = 0, transposed = FALSE, dodraws = TRUE))
      }
    }
  }
  list(delta_arr = delta_arr, d_arr = d_arr, phi_arr = phi_arr,
       use_full = use_full)
}


# Helper 3: Add mu component
#
# Two paths:
#   (a) BART / linear with mixture: pull mu from object$nmix$compdraw[[r]][[1]]$mu
#   (b) Heter-cov:                   pull mu from object$mu_draw[r, ]
.add_mu_component <- function(delta_z_array, object, kept_draws_indices) {
  npred      <- dim(delta_z_array)[1]
  nvar       <- dim(delta_z_array)[2]
  ndraws_out <- dim(delta_z_array)[3]

  if (inherits(object, "bayesm.HART.HeterCov")) {
    if (is.null(object$mu_draw))
      stop("Heter-cov object missing mu_draw.")
    mudraw <- t(object$mu_draw[kept_draws_indices, , drop = FALSE])  # nvar x ndraws_out
  } else {
    if (is.null(object$nmix) || is.null(object$nmix$compdraw))
      stop("DeltaZ+mu requires nmix$compdraw.")
    ncomp <- if (!is.null(object$nmix$probdraw)) ncol(object$nmix$probdraw)
             else                                length(object$nmix$compdraw[[1]])
    if (ncomp > 1)
      warning("DeltaZ+mu prediction currently only uses mu from the first mixture component.")
    mudraw <- tryCatch({
      sapply(object$nmix$compdraw[kept_draws_indices], function(x) {
        if (is.null(x[[1]]$mu)) stop("mu missing in component draw")
        x[[1]]$mu
      })
    }, error = function(e)
      stop(paste("Error extracting mu from nmix$compdraw:", e$message)))
  }

  if (!is.matrix(mudraw) || nrow(mudraw) != nvar || ncol(mudraw) != ndraws_out)
    stop("Extracted mudraw has incorrect dimensions.")

  mudraw_array <- aperm(array(mudraw, dim = c(nvar, ndraws_out, npred)),
                        c(3, 1, 2))
  return(delta_z_array + mudraw_array)
}#ADD_MU_COMPONENT
