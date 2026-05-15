# ==============================================================================
# NegBin Predictive Engines
# ==============================================================================
# Internal helpers for negative binomial posterior/prior predictive responses.

.predict_negbin_posterior_response <- function(object, newdata,
                                               kept_draws_indices, r_verbose) {
  if (is.null(newdata$regdata) || !is.list(newdata$regdata)) {
    stop("newdata$regdata must be provided as a list for negbin posterior predictions.")
  }

  ndraws_out <- length(kept_draws_indices)
  nreg <- dim(object$betadraw)[1]
  nvar <- dim(object$betadraw)[2]
  if (length(newdata$regdata) != nreg) {
    stop(sprintf("newdata$regdata must have length %d (one per estimation unit).", nreg))
  }

  out <- vector("list", nreg)
  names(out) <- paste0("Unit_", seq_len(nreg))

  for (i in seq_len(nreg)) {
    X_i <- newdata$regdata[[i]]$X
    if (is.null(X_i) || !is.matrix(X_i)) {
      stop(sprintf("newdata$regdata[[%d]]$X must be a matrix.", i))
    }
    if (ncol(X_i) != nvar) {
      stop(sprintf("newdata$regdata[[%d]]$X has %d columns; expected %d.",
                   i, ncol(X_i), nvar))
    }
    if (r_verbose) {
      cat("Simulating negbin posterior responses for Unit", i, "of", nreg, "...\n")
    }

    nobs_i <- nrow(X_i)
    y_draws_i <- matrix(0.0, nrow = nobs_i, ncol = ndraws_out)
    draw_count <- 0L
    for (s_orig in kept_draws_indices) {
      draw_count <- draw_count + 1L
      beta_is <- object$betadraw[i, , s_orig, drop = TRUE]
      alpha_s <- object$alphadraw[s_orig]
      mu_is <- as.vector(exp(X_i %*% beta_is))
      y_draws_i[, draw_count] <- rnbinom(nobs_i, size = alpha_s, mu = mu_is)
    }
    out[[i]] <- y_draws_i
  }

  if (r_verbose) cat("NegBin posterior response simulation complete.\n")
  out
}

.predict_negbin_prior_response <- function(object, newdata, delta_z_array,
                                           kept_draws_indices, nsim, r_verbose,
                                           hetercov_comps = NULL) {
  if (is.null(newdata$X) || !is.list(newdata$X)) {
    stop("newdata$X must be provided as a list for negbin prior predictions.")
  }

  is_heter <- inherits(object, "bayesm.HART.HeterCov")
  if (!is_heter && !requireNamespace("bayesm", quietly = TRUE)) {
    stop("The 'bayesm' package is required for negbin prior predictions. Please install it.")
  }

  npred <- length(newdata$X)
  nvar <- dim(object$betadraw)[2]
  ndraws_total <- dim(object$betadraw)[3]
  ndraws_out <- length(kept_draws_indices)
  if (dim(delta_z_array)[3] != ndraws_out) {
    stop("delta_z_array draw dimension does not match kept_draws_indices length.")
  }
  if (nrow(delta_z_array) == 1L && npred > 1L) {
    delta_z_array <- delta_z_array[rep(1L, npred), , , drop = FALSE]
  }
  if (nrow(delta_z_array) != npred) {
    stop("delta_z_array unit dimension must match length(newdata$X).")
  }

  sigma_kept <- NULL
  mu_kept <- NULL
  if (is_heter) {
    if (is.null(newdata$Z))
      stop("Heter-cov negbin prior prediction requires newdata$Z.")
    if (nrow(newdata$Z) != npred)
      stop("nrow(newdata$Z) must equal length(newdata$X).")
    if (is.null(object$mu_draw))
      stop("Heter-cov object missing mu_draw.")
    mu_kept <- object$mu_draw[kept_draws_indices, , drop = FALSE]

    sigma_cache <- .cache_get(object, "SigmaZ_unique_draws")
    if (!is.null(sigma_cache)) {
      map <- .match_cached_unique_rows(object, newdata$Z)
      if (!is.null(map) && map$all_seen) {
        sigma_kept <- sigma_cache[map$idx, , , kept_draws_indices, drop = FALSE]
      }
    }
    if (is.null(sigma_kept) && is.null(hetercov_comps)) {
      hetercov_comps <- .hetercov_components(object, newdata$Z, npred, nvar,
                                             ndraws_total, r_verbose)
    }
  } else {
    .assert_single_component_nmix(object, "negbin prior prediction")
  }

  out <- vector("list", npred)
  names(out) <- paste0("PredUnit_", seq_len(npred))

  for (i in seq_len(npred)) {
    X_i <- newdata$X[[i]]
    if (is.null(X_i) || !is.matrix(X_i)) {
      stop(sprintf("newdata$X[[%d]] must be a matrix.", i))
    }
    if (ncol(X_i) != nvar) {
      stop(sprintf("newdata$X[[%d]] has %d columns; expected %d.",
                   i, ncol(X_i), nvar))
    }
    if (r_verbose) {
      cat("Simulating negbin prior responses for Prediction Unit", i, "of", npred, "...\n")
    }

    nobs_i <- nrow(X_i)
    y_draws_i <- matrix(0.0, nrow = nobs_i, ncol = ndraws_out)

    for (draw_idx in seq_len(ndraws_out)) {
      s_orig <- kept_draws_indices[draw_idx]
      delta_is <- delta_z_array[i, , draw_idx, drop = TRUE]

      if (is_heter) {
        mu_s <- mu_kept[draw_idx, ]
        z <- rnorm(nvar)
        if (!is.null(sigma_kept)) {
          Sigma_is <- sigma_kept[i, , , draw_idx]
          U <- tryCatch(chol(Sigma_is), error = function(e) NULL)
          if (is.null(U)) {
            stop(sprintf("negbin prior heter-cov failed at draw %d, unit %d: SigmaZ chol failed.",
                         s_orig, i))
          }
          eta <- as.vector(z %*% U + mu_s)
        } else {
          d_is  <- hetercov_comps$d_arr[i, , s_orig]
          x <- z * sqrt(d_is)
          if (hetercov_comps$use_full) {
            phi_is <- hetercov_comps$phi_arr[i, , , s_orig]
            for (j in 2:nvar) {
              x[j] <- x[j] + sum(phi_is[j, seq_len(j - 1)] * x[seq_len(j - 1)])
            }
          }
          eta <- x + mu_s
        }
      } else {
        comps_s <- object$nmix$compdraw[[s_orig]]
        pvec_s  <- object$nmix$probdraw[s_orig, ]
        eta <- as.vector(bayesm::rmixture(n = 1L, pvec = pvec_s, comps = comps_s)$x)
      }

      if (nsim > 1L) {
        eta_acc <- eta
        for (k in 2:nsim) {
          if (is_heter) {
            mu_s <- mu_kept[draw_idx, ]
            z <- rnorm(nvar)
            if (!is.null(sigma_kept)) {
              Sigma_is <- sigma_kept[i, , , draw_idx]
              U <- tryCatch(chol(Sigma_is), error = function(e) NULL)
              if (is.null(U)) {
                stop(sprintf("negbin prior heter-cov failed at draw %d, unit %d: SigmaZ chol failed.",
                             s_orig, i))
              }
              eta_k <- as.vector(z %*% U + mu_s)
            } else {
              d_is  <- hetercov_comps$d_arr[i, , s_orig]
              x <- z * sqrt(d_is)
              if (hetercov_comps$use_full) {
                phi_is <- hetercov_comps$phi_arr[i, , , s_orig]
                for (j in 2:nvar) {
                  x[j] <- x[j] + sum(phi_is[j, seq_len(j - 1)] * x[seq_len(j - 1)])
                }
              }
              eta_k <- x + mu_s
            }
          } else {
            comps_s <- object$nmix$compdraw[[s_orig]]
            pvec_s  <- object$nmix$probdraw[s_orig, ]
            eta_k <- as.vector(bayesm::rmixture(n = 1L, pvec = pvec_s, comps = comps_s)$x)
          }
          eta_acc <- eta_acc + eta_k
        }
        eta <- eta_acc / nsim
      }

      beta_draw <- delta_is + eta
      mu_y <- as.vector(exp(X_i %*% beta_draw))
      alpha_s <- object$alphadraw[s_orig]
      y_draws_i[, draw_idx] <- rnbinom(nobs_i, size = alpha_s, mu = mu_y)
    }

    out[[i]] <- y_draws_i
  }

  if (r_verbose) cat("NegBin prior response simulation complete.\n")
  out
}
