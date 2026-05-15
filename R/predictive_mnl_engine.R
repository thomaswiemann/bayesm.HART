# ==============================================================================
# MNL Predictive Engines
# ==============================================================================
# Internal helpers for MNL posterior/prior predictive probabilities.
# These are separated from the model wrapper to keep predictive logic modular.

# Predict posterior choice probabilities for original estimation units.
.predict_posterior_probs <- function(object, newdata, kept_draws_indices,
                                     r_verbose) {
  dims_beta <- dim(object$betadraw)
  nlgt <- dims_beta[1]
  ndraws_out <- length(kept_draws_indices)

  if (is.null(newdata$p)) stop("newdata$p must be provided.")
  p <- newdata$p
  if (is.null(newdata$nlgtdata)) stop("newdata$nlgtdata must be provided.")

  prob_pred_list <- vector("list", nlgt)
  names(prob_pred_list) <- paste0("Unit_", seq_len(nlgt))
  betadraw <- object$betadraw

  for (i in seq_len(nlgt)) {
    if (r_verbose) {
      cat("Calculating posterior probabilities for Unit", i, "of", nlgt, "...\n")
    }
    X_i <- newdata$nlgtdata[[i]]$X
    T_i <- nrow(X_i) / p
    prob_array_i <- array(0.0, dim = c(T_i, p, ndraws_out))

    draw_count <- 0L
    for (s_orig in kept_draws_indices) {
      draw_count <- draw_count + 1L
      beta_is <- betadraw[i, , s_orig, drop = TRUE]
      prob_array_i[, , draw_count] <- calculate_mnl_probs_from_beta(X_i, beta_is, p)
    }
    prob_pred_list[[i]] <- prob_array_i
  }
  if (r_verbose) cat("Posterior probability calculation complete.\n")
  prob_pred_list
}

# Predict prior choice probabilities for prediction units.
.predict_prior_probs <- function(object, newdata, delta_z_array,
                                 kept_draws_indices_nmix, nsim, r_verbose,
                                 hetercov_comps = NULL) {
  is_heter <- inherits(object, "bayesm.HART.HeterCov")
  if (!is_heter && !requireNamespace("bayesm", quietly = TRUE)) {
    stop("The 'bayesm' package is required for prior choice_probs prediction. Please install it.")
  }

  dims_beta    <- dim(object$betadraw)
  nvar         <- dims_beta[2]
  ndraws_total <- dims_beta[3]
  npred        <- dim(delta_z_array)[1]
  ndraws_out   <- dim(delta_z_array)[3]

  if (is.null(newdata$p)) stop("newdata$p must be provided.")
  p <- newdata$p
  if (is.null(newdata$X)) stop("newdata$X must be provided as a list.")

  if (is_heter) {
    if (is.null(newdata$Z))
      stop("Heter-cov prior choice_probs requires newdata$Z (same Z as DeltaZ).")
    if (npred != nrow(newdata$Z))
      stop("npred (from delta_z_array) does not match nrow(newdata$Z).")
    if (is.null(object$mu_draw))
      stop("Heter-cov object missing mu_draw.")
    mu_kept <- object$mu_draw[kept_draws_indices_nmix, , drop = FALSE]

    sigma_kept <- NULL
    sigma_cache <- .cache_get(object, "SigmaZ_unique_draws")
    if (!is.null(sigma_cache)) {
      map <- .match_cached_unique_rows(object, newdata$Z)
      if (!is.null(map) && map$all_seen) {
        sigma_kept <- sigma_cache[map$idx, , , kept_draws_indices_nmix, drop = FALSE]
      }
    }
    if (is.null(sigma_kept) && is.null(hetercov_comps)) {
      hetercov_comps <- .hetercov_components(object, newdata$Z, npred, nvar,
                                             ndraws_total, r_verbose)
    }
  } else {
    .assert_single_component_nmix(object, "prior choice_probs prediction")
  }

  prob_pred_list <- vector("list", npred)
  names(prob_pred_list) <- paste0("PredUnit_", seq_len(npred))

  for (i in seq_len(npred)) {
    if (r_verbose) {
      cat("Calculating prior probabilities for Prediction Unit", i, "of", npred, "...\n")
    }
    X_i <- newdata$X[[i]]
    T_i <- nrow(X_i) / p
    prob_array_i <- array(0.0, dim = c(T_i, p, ndraws_out))

    for (draw_idx in seq_len(ndraws_out)) {
      DeltaZ_is <- delta_z_array[i, , draw_idx, drop = TRUE]
      s_orig    <- kept_draws_indices_nmix[draw_idx]

      if (is_heter) {
        mu_s <- mu_kept[draw_idx, ]
        z_mat <- matrix(rnorm(nsim * nvar), nrow = nsim, ncol = nvar)
        if (!is.null(sigma_kept)) {
          Sigma_is <- sigma_kept[i, , , draw_idx]
          U <- tryCatch(chol(Sigma_is), error = function(e) NULL)
          if (is.null(U)) {
            stop(sprintf("prior choice_probs heter-cov failed at draw %d, unit %d: SigmaZ chol failed.",
                         s_orig, i))
          }
          eta_mat <- z_mat %*% U + matrix(mu_s, nrow = nsim, ncol = nvar, byrow = TRUE)
        } else {
          d_is  <- hetercov_comps$d_arr[i, , s_orig]
          x_mat <- z_mat * matrix(sqrt(d_is), nrow = nsim, ncol = nvar, byrow = TRUE)
          if (hetercov_comps$use_full) {
            phi_is <- hetercov_comps$phi_arr[i, , , s_orig]
            for (k in seq_len(nsim)) {
              x <- x_mat[k, ]
              for (j in 2:nvar) {
                x[j] <- x[j] + sum(phi_is[j, seq_len(j - 1)] * x[seq_len(j - 1)])
              }
              x_mat[k, ] <- x
            }
          }
          eta_mat <- x_mat + matrix(mu_s, nrow = nsim, ncol = nvar, byrow = TRUE)
        }

        prob_is_sum <- matrix(0.0, nrow = T_i, ncol = p)
        for (k in seq_len(nsim)) {
          beta_is_k <- DeltaZ_is + eta_mat[k, ]
          probs_mat_is <- calculate_mnl_probs_from_beta(X_i, beta_is_k, p)
          prob_is_sum <- prob_is_sum + probs_mat_is
        }
      } else {
        comps_s <- object$nmix$compdraw[[s_orig]]
        pvec_s  <- object$nmix$probdraw[s_orig, ]
        prob_is_sum <- matrix(0.0, nrow = T_i, ncol = p)
        for (k in seq_len(nsim)) {
          eta_is_k <- as.vector(bayesm::rmixture(n = 1, pvec = pvec_s, comps = comps_s)$x)
          beta_is_k <- DeltaZ_is + eta_is_k
          probs_mat_is <- calculate_mnl_probs_from_beta(X_i, beta_is_k, p)
          prob_is_sum <- prob_is_sum + probs_mat_is
        }
      }

      prob_array_i[, , draw_idx] <- prob_is_sum / nsim
    }
    prob_pred_list[[i]] <- prob_array_i
  }
  if (r_verbose) cat("Prior probability calculation complete.\n")
  prob_pred_list
}
