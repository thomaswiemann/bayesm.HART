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
.calculate_delta_z <- function(object, newdata, burn, r_verbose, ...) {
  nvar <- dim(object$betadraw)[2]
  ndraws_total <- dim(object$betadraw)[3]

  # Check if Z covariates were used in the original model fit
  model_has_Z <- !is.null(object$Deltadraw) || !is.null(object$bart_models)

  pred <- NULL # Initialize pred
  npred <- 1 # Default if no Z

  if (model_has_Z) {
    # Model with Z covariates - Assume newdata$Z is provided and is a matrix
    if (is.null(newdata$Z)) {
      stop("Model was fit with Z. 'newdata$Z' must be provided.")
    }
    npred <- nrow(newdata$Z)
    nz <- ncol(newdata$Z)

    pred <- raw_pred <- array(0, dim = c(npred, nvar, ndraws_total))

    # Calculate DeltaZ component (systematic part)
    if (!is.null(object$Deltadraw)) {
      # Linear model prediction case
      if (ncol(object$Deltadraw) != nz * nvar) {
        stop(paste("Dimension mismatch: ncol(Deltadraw) is",
                   ncol(object$Deltadraw), "but expected nz*nvar =", nz * nvar))
      }
      for (i in 1:ndraws_total) {
        Delta_draw <- matrix(object$Deltadraw[i, ], nrow = nz, byrow = TRUE)
        pred[, , i] <- newdata$Z %*% Delta_draw
      }#FOR
    } else if (!is.null(object$bart_models)) {
      # BART model prediction case
      if (length(object$bart_models) != nvar) {
        stop(paste("Number of BART models (", length(object$bart_models),
                   ") does not match nvar (", nvar, ")"))
      }
      # --- Prepare arguments for pwbart, handling ... ---
      passed_args <- list(...)
      pwbart_arg_names <- names(formals(bayesm.HART:::pwbart))
      valid_passed_args <- passed_args[names(passed_args) %in% pwbart_arg_names]
      base_args <- list(
        x.test = newdata$Z,
        mu = 0,
        transposed = FALSE,
        dodraws = TRUE
      )
      final_args_template <- utils::modifyList(base_args, valid_passed_args)
      final_args_template <-
        final_args_template[names(final_args_template) %in% pwbart_arg_names]
      # --- End argument preparation ---

      for (j in 1:nvar) {
        if (is.null(object$bart_models[[j]]$treedraws)) {
          stop(paste("Missing treedraws for BART model of coefficient", j))
        }
        if (r_verbose) {
          cat("Predicting coefficient", j, "with BART model\n")
        }
        current_iter_args <- final_args_template
        current_iter_args$treedraws <- object$bart_models[[j]]$treedraws
        bart_pred <- do.call(bayesm.HART:::pwbart, current_iter_args)
        raw_pred[, j, ] <- t(bart_pred)
      }#FOR

      # Apply scaling using rooti (specific to this package's BART impl.)
      if (is.null(object$nmix) || is.null(object$nmix$compdraw)) {
        stop("Missing nmix$compdraw needed for BART scaling.")
      }
      for (s in 1:ndraws_total) {
        if (is.null(object$nmix$compdraw[[s]][[1]]$rooti)) {
          stop(paste("Missing nmix$compdraw[[ ", s, " ]][[1]]$rooti", sep = ""))
        }
        root_i_s <- object$nmix$compdraw[[s]][[1]]$rooti
        if (inherits(try(solve(root_i_s), silent = TRUE), "try-error")) {
          stop(paste("rooti matrix is singular for draw", s))
        }
        root_s_inv <- solve(root_i_s)
        pred[, , s] <- raw_pred[, , s] %*% root_s_inv
      }#FOR
    }#ELSEIF
  } else {
    # Model without Z covariates
    npred <- 1 # Predict only the baseline heterogeneity component
    pred <- array(0, dim = c(npred, nvar, ndraws_total))
  }#ELSE (model_has_Z)

  # Apply burn-in
  if (burn >= ndraws_total) stop("Burn-in period too large.") # Should be caught earlier
  if (burn > 0) {
    kept_draws <- (burn + 1):ndraws_total
    pred <- pred[, , kept_draws, drop = FALSE]
  }#IF

  return(pred)
}#CALCULATE_DELTA_Z

# Helper 3: Add mu component
.add_mu_component <- function(delta_z_array, object, kept_draws_indices) {
  if (is.null(object$nmix) || is.null(object$nmix$compdraw)) {
    stop("DeltaZ+mu requires a BART model with nmix$compdraw.")
  }

  # Check ncomp and issue warning if > 1
  ncomp <- 1 # Default
  if (!is.null(object$nmix$probdraw)) {
    ncomp <- ncol(object$nmix$probdraw)
  } else if (!is.null(object$nmix$compdraw)) {
    # Fallback: check length of first draw's components if probdraw missing
    ncomp <- length(object$nmix$compdraw[[1]])
  }
  if (ncomp > 1) {
    warning("DeltaZ+mu prediction currently only uses mu from the first mixture component.")
  }

  # Extract mu from the first component for kept draws
  mudraw <- tryCatch({
    sapply(object$nmix$compdraw[kept_draws_indices], function(x) {
      if (is.null(x[[1]]$mu)) stop("mu missing in component draw")
      x[[1]]$mu
    })
  }, error = function(e) {
    stop(paste("Error extracting mu from nmix$compdraw:", e$message))
  })

  # Get dimensions from inputs
  npred <- dim(delta_z_array)[1]
  nvar <- dim(delta_z_array)[2]
  ndraws_out <- dim(delta_z_array)[3]

  if (!is.matrix(mudraw) || nrow(mudraw) != nvar || ncol(mudraw) != ndraws_out) {
    stop("Extracted mudraw has incorrect dimensions.")
  }

  # Reshape mudraw and add to delta_z_array
  mudraw_array <- array(mudraw, dim = c(nvar, ndraws_out, npred))
  mudraw_array <- aperm(mudraw_array, c(3, 1, 2))

  return(delta_z_array + mudraw_array)
}#ADD_MU_COMPONENT
