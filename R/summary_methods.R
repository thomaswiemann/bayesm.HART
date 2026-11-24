# ==============================================================================
# Summary Methods for Hierarchical Models
# ==============================================================================
# Model-agnostic summary and diagnostics.  Works with any fitted object from
# rhierMnlRwMixture, rhierLinearMixture, or rhierNegbinRw (with or without
# HART / heteroscedastic extensions).

# ==============================================================================
# Internal helpers
# ==============================================================================

# Weighted variance (population weights; denominator = sum(w), not sum(w)-1)
.wvar <- function(x, w) {
  sw <- sum(w)
  mu <- sum(x * w) / sw
  sum(w * (x - mu)^2) / sw
}

.compute_r2 <- function(object, Z, burn, coefs, r_verbose) {
  nvar   <- dim(object$betadraw)[2]
  ndraws <- dim(object$betadraw)[3]
  if (is.null(coefs)) coefs <- seq_len(nvar)

  is_heter <- .hart_cov_kind(object) == "heteroskedastic"

  # With Z = NULL the unique-row cache is used, weighted by z_counts.
  cache <- attr(object, "hart_cache", exact = TRUE)
  used_cache <- FALSE

  if (is.null(Z) && !is.null(cache$DeltaZ_unique_draws)) {
    DeltaZ <- cache$DeltaZ_unique_draws   # [n_unique x nvar x ndraws]
    w      <- cache$z_counts
    used_cache <- TRUE
  } else if (!is.null(Z)) {
    DeltaZ <- predict(object, newdata = list(Z = Z), type = "DeltaZ",
                       burn = burn, r_verbose = r_verbose)
    w <- NULL  # uniform weights (standard var)
  } else {
    stop("R-squared requires either Z or a fitted BART model with cached DeltaZ draws.")
  }

  npred <- dim(DeltaZ)[1]
  # The cache holds all R/keep draws.
  if (used_cache) {
    keep_idx <- seq.int(min(burn + 1, ndraws), ndraws)
    DeltaZ <- DeltaZ[, , keep_idx, drop = FALSE]
  }
  ndraws_out <- dim(DeltaZ)[3]

  vfun <- if (!is.null(w)) function(x) .wvar(x, w) else var

  # Sigma_kk draws
  if (is_heter) {
    if (used_cache && !is.null(cache$SigmaZ_unique_draws)) {
      SigmaZ <- cache$SigmaZ_unique_draws[, , , keep_idx, drop = FALSE]
    } else if (!is.null(Z)) {
      SigmaZ <- predict(object, newdata = list(Z = Z), type = "SigmaZ",
                          burn = burn, r_verbose = r_verbose)
    } else {
      stop("Heteroscedastic R-squared requires Z or cached SigmaZ draws.")
    }
  } else {
    compdraw <- object$nmix$compdraw
    cd_idx <- seq.int(min(burn + 1, length(compdraw)), length(compdraw))
    Sigma_diag <- t(vapply(cd_idx, function(ii) {
      rooti <- compdraw[[ii]][[1]]$rooti
      diag(crossprod(solve(rooti)))  # Sigma = A^T A, A = rooti^{-1}
    }, numeric(nvar)))
  }

  # R2 per coefficient
  if (is_heter) {
    R2_draws <- vector("list", length(coefs))
    R2_mean  <- numeric(length(coefs))
    names(R2_draws) <- names(R2_mean) <- paste0("k", coefs)
    for (ci in seq_along(coefs)) {
      k <- coefs[ci]
      r2_mat <- matrix(0, nrow = ndraws_out, ncol = npred)
      for (s in seq_len(ndraws_out)) {
        V_exp <- vfun(DeltaZ[, k, s])
        V_res <- SigmaZ[, k, k, s]
        r2_mat[s, ] <- V_exp / (V_exp + V_res)
      }
      R2_draws[[ci]] <- r2_mat
      R2_mean[ci]    <- if (!is.null(w)) {
        mean(apply(r2_mat, 1, function(row) sum(row * w) / sum(w)))
      } else {
        mean(colMeans(r2_mat))
      }
    }
  } else {
    R2_draws <- matrix(0, nrow = ndraws_out, ncol = length(coefs))
    colnames(R2_draws) <- paste0("k", coefs)
    R2_mean <- numeric(length(coefs))
    names(R2_mean) <- paste0("k", coefs)
    for (ci in seq_along(coefs)) {
      k <- coefs[ci]
      for (s in seq_len(ndraws_out)) {
        V_exp <- vfun(DeltaZ[, k, s])
        R2_draws[s, ci] <- V_exp / (V_exp + Sigma_diag[s, k])
      }
      R2_mean[ci] <- mean(R2_draws[, ci])
    }
  }

  list(R2_draws = R2_draws, R2_mean = R2_mean, heteroscedastic = is_heter)
}

# ==============================================================================
# Shared summary implementation
# ==============================================================================

.summary_hier <- function(object, Z = NULL, burn = 0, coefs = NULL,
                           r_verbose = FALSE, model_type = "MNL") {
  nreg   <- dim(object$betadraw)[1]
  nvar   <- dim(object$betadraw)[2]
  ndraws <- dim(object$betadraw)[3]

  out <- list()
  out$model_type     <- model_type
  out$nreg           <- nreg
  out$nvar           <- nvar
  out$ndraws         <- ndraws
  out$is_bart        <- !is.null(object$bart_models)
  out$is_heter       <- .hart_cov_kind(object) == "heteroskedastic"

  # Acceptance rates (MNL, NegBin)
  if (!is.null(object$acceptrbeta))
    out$accept_beta <- object$acceptrbeta
  if (!is.null(object$acceptralpha))
    out$accept_alpha <- object$acceptralpha

  # Posterior mean coefficients
  beta_mean <- apply(object$betadraw[, , max(1, burn + 1):ndraws, drop = FALSE],
                     c(1, 2), mean)
  out$beta_mean <- beta_mean

  # Variable importance
  if (out$is_bart && !is.null(object$varcount)) {
    out$varcount <- object$varcount
    out$varprob  <- object$varprob
  }

  # R2, computed when Z is supplied or the DeltaZ cache exists.
  cache <- attr(object, "hart_cache", exact = TRUE)
  has_cache <- !is.null(cache$DeltaZ_unique_draws)
  if (!is.null(Z) || has_cache) {
    out$r2 <- .compute_r2(object, Z, burn, coefs, r_verbose)
  }

  class(out) <- "summary.rhierModel"
  out
}

# ==============================================================================
# Print method
# ==============================================================================

#' @export
print.summary.rhierModel <- function(x, digits = 3, ...) {
  cat("Hierarchical", x$model_type, "Model Summary\n")
  cat(strrep("-", 50), "\n")
  cat("  Units:", x$nreg, " | Coefficients:", x$nvar,
      " | Posterior draws:", x$ndraws, "\n")
  cat("  Prior:", if (x$is_heter) "Heteroscedastic HART"
               else if (x$is_bart) "HART" else "Linear", "\n")

  if (!is.null(x$accept_beta))
    cat("  Beta acceptance rate:", round(x$accept_beta, 1), "%\n")
  if (!is.null(x$accept_alpha))
    cat("  Alpha acceptance rate:", round(x$accept_alpha, 1), "%\n")

  if (!is.null(x$r2)) {
    cat("\nVariance Explained (R^2):\n")
    r2_pct <- round(x$r2$R2_mean * 100, digits)
    for (i in seq_along(r2_pct)) {
      label <- if (x$r2$heteroscedastic) " (per-respondent mean)" else ""
      cat(sprintf("  Coefficient %s: %.1f%%%s\n",
                  names(r2_pct)[i], r2_pct[i], label))
    }
  }

  if (!is.null(x$varprob)) {
    cat("\nBART Variable Importance (posterior inclusion probs):\n")
    vp <- colMeans(x$varprob)
    ord <- order(vp, decreasing = TRUE)
    top_n <- min(10, length(vp))
    for (j in ord[seq_len(top_n)]) {
      cat(sprintf("  Z[%d]: %.3f\n", j, vp[j]))
    }
    if (length(vp) > top_n) cat("  ...\n")
  }

  invisible(x)
}

# ==============================================================================
# Model-specific summary dispatchers
# ==============================================================================

#' Summary Method for Hierarchical Models
#'
#' Computes model diagnostics including acceptance rates, posterior mean betas,
#' BART variable importance, and Bayesian R^2.
#'

#'
#' @param object A fitted hierarchical model object.
#' @param Z Optional matrix of unit-level characteristics for R^2 computation.
#'   If \code{NULL} and cached DeltaZ draws exist, in-sample R^2 is still computed.
#' @param burn Number of initial draws to discard (thinned units).
#' @param coefs Integer vector of coefficient indices for R^2. Default: all.
#' @param r_verbose Print progress? Default FALSE.
#' @param ... Ignored.
#' @return An object of class `summary.rhierModel` with print method.
#' @export
summary.rhierMnlRwMixture <- function(object, Z = NULL, burn = 0,
                                       coefs = NULL, r_verbose = FALSE, ...) {
  .summary_hier(object, Z, burn, coefs, r_verbose, model_type = "MNL")
}

#' @rdname summary.rhierMnlRwMixture
#' @export
summary.rhierLinearMixture <- function(object, Z = NULL, burn = 0,
                                        coefs = NULL, r_verbose = FALSE, ...) {
  .summary_hier(object, Z, burn, coefs, r_verbose, model_type = "Linear")
}

#' @rdname summary.rhierMnlRwMixture
#' @export
summary.rhierNegbinRw <- function(object, Z = NULL, burn = 0,
                                    coefs = NULL, r_verbose = FALSE, ...) {
  .summary_hier(object, Z, burn, coefs, r_verbose, model_type = "NegBin")
}
