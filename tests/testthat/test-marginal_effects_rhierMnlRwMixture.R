# ==============================================================================
# Unit Tests for marginal_effects.rhierMnlRwMixture
# ==============================================================================

library(testthat)

# ==============================================================================
# Helper Function to Generate Data and Fit Model Quickly
# ==============================================================================

# Similar to the one in test-predict, ensures nz > 0
run_quick_model_for_mfx_test <- function(R = 50, keep = 1, useBART = FALSE) {
  set.seed(654) # Ensure reproducibility for tests
  sim_data <- sim_hier_mnl(
    nlgt = 25, nT = 8, p = 3, nz = 5, nXa = 1, nXd = 0, const = TRUE,
    # Use linear function to avoid potential Friedman simulation issues for now
    beta_func_type = "linear",
    beta_func_args = list() # Let sim_hier_mnl use default Delta
  )
  dims <- sim_data$true_values$dimensions
  ncoef <- dims$ncoef
  nz <- dims$nz

  Prior <- list(ncomp = 1)
  if (useBART) {
    Prior$bart <- list(num_trees = 10, numcut = 10) # Fast BART settings
  } else {
    Prior$Deltabar <- matrix(0, nz * ncoef, 1)
    Prior$Ad <- diag(nz * ncoef) * 0.01
  }

  Mcmc <- list(R = R, keep = keep, nprint = 0) # Suppress printing

  # Fit model (Capture non-essential output)
  fit <- NULL
  capture.output({
    fit <- rhierMnlRwMixture(Data = sim_data, Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  # Basic check that fit looks reasonable
  if(is.null(fit$betadraw)) stop("Model fitting failed for mfx test setup")

  # Add nz to the fit object list for testing marginal_effects nz detection
  fit$nz <- nz

  return(list(fit = fit, sim_data = sim_data))
}

# ==============================================================================
# Test Setup: Generate Fit and Data ONCE for all tests in this file
# ==============================================================================

# Run with BART for broader coverage, use fewer draws for speed
n_draws_test <- 30 
quick_fit_bart <- run_quick_model_for_mfx_test(R = n_draws_test, useBART = TRUE)
fit_object <- quick_fit_bart$fit
orig_data <- quick_fit_bart$sim_data

# Define common parameters based on the *actual* quick fit
dims_actual <- orig_data$true_values$dimensions
ncoef <- dims_actual$ncoef
nz <- dims_actual$nz # Number of Z covariates
nlgt <- dims_actual$nlgt
ndraws_total <- n_draws_test # R from the quick run

base_Z <- orig_data$Z # Base Z matrix for context

# Example z_values matrices for tests
n_grid_points_test <- 5 # Number of points for 1D grid

# Vary Z column 1 (-1 to 1), others NA
z_col1_vec <- seq(-1, 1, length.out = n_grid_points_test)
z_values_1d <- matrix(NA_real_, nrow = n_grid_points_test, ncol = nz)
z_values_1d[, 1] <- z_col1_vec
# Add colnames if base_Z has them
if (!is.null(colnames(base_Z))) {
    colnames(z_values_1d) <- colnames(base_Z)
}

# 1 row, all NA (for MEM-like test)
z_values_all_na <- matrix(NA_real_, nrow = 1, ncol = nz)
if (!is.null(colnames(base_Z))) {
    colnames(z_values_all_na) <- colnames(base_Z)
}

# Matrix violating the NA constraint (mix in col 2)
z_values_violating_na <- matrix(NA_real_, nrow = 3, ncol = nz)
z_values_violating_na[, 1] <- c(1, 2, 3) # Non-NA column
z_values_violating_na[1:2, 2] <- c(4, 5) # Mixed NA/non-NA column
if (!is.null(colnames(base_Z))) {
    colnames(z_values_violating_na) <- colnames(base_Z)
}

# Other parameters
burn_val <- 5 # Burn must be < ndraws_total
ndraws_use <- ndraws_total - burn_val


# ==============================================================================
# Test Cases
# ==============================================================================

test_that("marginal_effects: Basic functionality and output structure (BART fit)", {

  mfx <- marginal_effects(
    object = fit_object,
    z_values = z_values_1d, # Use new args
    Z = base_Z[1:10, ],
    burn = burn_val,
    verbose = FALSE
  )

  expect_s3_class(mfx, "marginal_effects")
  expect_type(mfx, "list")
  expect_named(mfx, c("z_values", "avg_betabar_draws", "call", "burn")) # Updated names

  expect_equal(mfx$z_values, z_values_1d) # Check z_values matrix
  expect_type(mfx$avg_betabar_draws, "list")
  expect_length(mfx$avg_betabar_draws, nrow(z_values_1d)) # Use nrow
  # expect_equal(mfx$target_covariate_index, target_idx) # Removed
  expect_equal(mfx$burn, burn_val)
  expect_true(!is.null(mfx$call))
})

test_that("marginal_effects: avg_betabar_draws elements have correct dimensions (BART fit)", {
  
  mfx <- marginal_effects(fit_object, z_values = z_values_1d, Z = base_Z, burn = burn_val, verbose = FALSE) # Updated call

  # Check dimensions of each matrix in the list
  expect_true(length(mfx$avg_betabar_draws) > 0, info = "avg_betabar_draws list should not be empty")
  # Loop based on nrow of z_values
  for (i in 1:nrow(z_values_1d)) {
      expect_true(is.matrix(mfx$avg_betabar_draws[[i]]), 
                  info = paste("Element", i, "should be a matrix"))
      expect_equal(nrow(mfx$avg_betabar_draws[[i]]), ncoef, 
                  info = paste("nrow mismatch for element", i))
      expect_equal(ncol(mfx$avg_betabar_draws[[i]]), ndraws_use, 
                  info = paste("ncol mismatch for element", i))
  }
})

# ==============================================================================
# Edge Case & Functionality Tests
# ==============================================================================

test_that("marginal_effects: Handles z_values edge cases", { # Renamed and simplified

  # z_values with 1 row
  z_values_1row <- z_values_1d[1, , drop = FALSE]
  mfx_1row <- marginal_effects(fit_object, z_values = z_values_1row, Z = base_Z, burn = burn_val, verbose = FALSE)
  expect_equal(nrow(mfx_1row$z_values), 1)
  expect_length(mfx_1row$avg_betabar_draws, 1)
  expect_true(is.matrix(mfx_1row$avg_betabar_draws[[1]]))
  expect_equal(ncol(mfx_1row$avg_betabar_draws[[1]]), ndraws_use)

  # z_values with all NA columns - NOW EXPECTS ERROR
  expect_error(
      marginal_effects(fit_object, z_values = z_values_all_na, Z = base_Z, burn = burn_val, verbose = FALSE),
      "'z_values' must contain at least one non-NA column to define grid points.",
      fixed = TRUE
  )
})

test_that("marginal_effects: Handles burn edge cases", {

  # Zero burn (burn=0)
  mfx_zero_burn <- marginal_effects(fit_object, z_values=z_values_1d, Z=base_Z, burn = 0, verbose = FALSE)
  expect_equal(mfx_zero_burn$burn, 0)
  expect_true(is.matrix(mfx_zero_burn$avg_betabar_draws[[1]]))
  expect_equal(ncol(mfx_zero_burn$avg_betabar_draws[[1]]), ndraws_total)

  # Maximum burn (leaving 1 draw)
  # Ensure ndraws_total > 1 for this test
  if (ndraws_total > 1) {
      mfx_max_burn <- marginal_effects(fit_object, z_values=z_values_1d, Z=base_Z, burn = ndraws_total - 1, verbose = FALSE)
      expect_equal(mfx_max_burn$burn, ndraws_total - 1)
      expect_true(is.matrix(mfx_max_burn$avg_betabar_draws[[1]]))
      expect_equal(ncol(mfx_max_burn$avg_betabar_draws[[1]]), 1)
  } else {
      skip("Skipping max burn test because ndraws_total <= 1")
  }
})

# ==============================================================================
# Setup for Summary Method Tests
# ==============================================================================

# Generate a marginal_effects object using the new function signature
mfx_object <- marginal_effects(
    object = fit_object,
    z_values = z_values_1d, # Use test z_values
    Z = base_Z[1:10, , drop = FALSE],             # Use test Z
    burn = burn_val,
    verbose = FALSE
  )

# Create the summary object needed for plot tests
sum_mfx_1d <- summary(mfx_object)

# Generate a second object (e.g., from a linear fit) - requires updating this too
quick_fit_linear <- run_quick_model_for_mfx_test(R = n_draws_test, useBART = FALSE)
mfx_linear <- marginal_effects(
    object = quick_fit_linear$fit,
    z_values = z_values_1d, # Use same test z_values
    Z = quick_fit_linear$sim_data$Z[1:10, , drop = FALSE],
    burn = burn_val,
    verbose = FALSE
  )
mfx_summary_linear <- summary(mfx_linear) # Will be used later if plot tests are re-enabled


# ==============================================================================
# Test Cases for summary.marginal_effects
# ==============================================================================

test_that("summary.marginal_effects: Basic structure and class", {
  mfx_summary <- summary(mfx_object)

  expect_s3_class(mfx_summary, "summary.marginal_effects")
  expect_type(mfx_summary, "list")
  expect_named(mfx_summary, c("summary_df", "z_values", "probs", "call")) # Updated names

  expect_equal(mfx_summary$z_values, mfx_object$z_values) # Check z_values
  # expect_equal(mfx_summary$target_covariate_index, mfx_object$target_covariate_index) # Removed
  expect_true(!is.null(mfx_summary$call))
})

test_that("summary.marginal_effects: summary_df structure and content", {
  default_probs = c(0.025, 0.05, 0.5, 0.95, 0.975)
  expected_q_names <- paste0("q", gsub("\\.", "", default_probs * 100))

  # Determine expected z columns based on non-NA columns in z_values_1d
  non_na_cols_idx <- which(colSums(is.na(z_values_1d)) < nrow(z_values_1d))
  z_cols_expected_names <- colnames(z_values_1d)[non_na_cols_idx]
  if (is.null(z_cols_expected_names)) { # Fallback if no colnames were set
      z_cols_expected_names <- paste0("z_val_", non_na_cols_idx)
  }

  expected_colnames <- c("coefficient_index", z_cols_expected_names, "mean", expected_q_names)
  
  mfx_summary <- summary(mfx_object) # Use default probs
  df <- mfx_summary$summary_df
  
  expect_s3_class(df, "data.frame")
  expect_named(df, expected_colnames, ignore.order = FALSE)
  expect_equal(nrow(df), nrow(z_values_1d) * ncoef) # Use nrow(z_values_1d)
  
  # Check column types (basic check)
  expect_type(df$coefficient_index, "integer")
  # expect_type(df$zk_value, "double") # Removed
  if (length(z_cols_expected_names) > 0) {
      expect_type(df[[z_cols_expected_names[1]]], "double") # Check first z col type
  }
  expect_type(df$mean, "double")
  expect_type(df[[expected_q_names[1]]], "double") # Check first quantile col type
})

test_that("summary.marginal_effects: Custom probs argument works", {
  custom_probs <- c(0.1, 0.9)
  expected_q_names <- paste0("q", custom_probs * 100) # q10, q90

  # Determine expected z columns based on non-NA columns in z_values_1d
  non_na_cols_idx <- which(colSums(is.na(z_values_1d)) < nrow(z_values_1d))
  z_cols_expected_names <- colnames(z_values_1d)[non_na_cols_idx]
  if (is.null(z_cols_expected_names)) { # Fallback
      z_cols_expected_names <- paste0("z_val_", non_na_cols_idx)
  }

  expected_colnames <- c("coefficient_index", z_cols_expected_names, "mean", expected_q_names)
  
  mfx_summary_custom <- summary(mfx_object, probs = custom_probs)
  df <- mfx_summary_custom$summary_df
  
  expect_equal(mfx_summary_custom$probs, custom_probs)
  expect_named(df, expected_colnames, ignore.order = FALSE)
  expect_type(df$q10, "double")
  expect_type(df$q90, "double")
})

test_that("summary.marginal_effects: Input validation", {
  
  # Invalid probs
  expect_error(
    summary(mfx_object, probs = c(0.1, 1.1)),
    "'probs' must be a numeric vector with values between 0 and 1."
  )
  expect_error(
    summary(mfx_object, probs = "a"),
    "'probs' must be a numeric vector with values between 0 and 1."
  )
  
  # Missing element in input object (check for z_values)
  bad_mfx_object <- mfx_object
  bad_mfx_object$z_values <- NULL # Remove z_values
  expect_error(
    summary(bad_mfx_object),
    "Invalid 'marginal_effects' object. Missing elements: z_values"
  )

  # Missing avg_betabar_draws (kept from before)
  bad_mfx_object_2 <- mfx_object
  bad_mfx_object_2$avg_betabar_draws <- NULL
  expect_error(
    summary(bad_mfx_object_2),
    "Invalid 'marginal_effects' object. Missing elements: avg_betabar_draws"
  )
})

# ==============================================================================
# Plot Method Tests (`plot.summary.marginal_effects`)
# ==============================================================================


test_that("plot.summary.marginal_effects: Basic plotting works (default x-axis)", {
  
  skip_if_not_installed("ggplot2")
  
  # Use summary object created earlier
  p <- plot(sum_mfx_1d, coef_index = 1)
  
  expect_s3_class(p, "ggplot")
  
  # Check default x-axis variable name (should be the first z_val_ col)
  expected_x_var <- "z_val_1"
  expect_equal(rlang::quo_name(p$mapping$x), expected_x_var)
  
  # Check default x-axis label
  expect_equal(p$labels$x, expected_x_var)
  
  # Check other basic labels
  expect_true(grepl("Coefficient 1", p$labels$y))
  expect_true(grepl(paste("vs", expected_x_var), p$labels$title))
})

test_that("plot.summary.marginal_effects: Plotting with specified x-axis works", {
  
  skip_if_not_installed("ggplot2")
  
  # Specify the x-axis variable explicitly
  plot_var <- "z_val_1"
  p <- plot(sum_mfx_1d, coef_index = 1, plot_axis_var = plot_var)
  
  expect_s3_class(p, "ggplot")
  
  # Check specified x-axis variable name
  expect_equal(rlang::quo_name(p$mapping$x), plot_var)
  
  # Check x-axis label reflects the specified variable
  expect_equal(p$labels$x, plot_var)
})

test_that("plot.summary.marginal_effects: Plotting multiple models works", {
  
  skip_if_not_installed("ggplot2")
  
  # Create a second summary object (can just copy/rename for structural test)
  sum_mfx_1d_alt <- sum_mfx_1d # Simple copy for structure test
  
  p <- plot(sum_mfx_1d, linear = sum_mfx_1d_alt, coef_index = 1)
  
  expect_s3_class(p, "ggplot")
  
  # Check that data includes both models
  expect_true("model_name" %in% names(p$data))
  expect_equal(length(unique(p$data$model_name)), 2)
  
  # Check default x-axis variable name again
  expected_x_var <- "z_val_1"
  expect_equal(rlang::quo_name(p$mapping$x), expected_x_var)
})


# test_that("plot.summary.marginal_effects: Input validation - x object", {
#   expect_error(plot(list()), "is not of class \"summary.marginal_effects\"")
#   bad_sum_mfx <- mfx_object
#   bad_sum_mfx$summary_df <- NULL
#   expect_error(plot(bad_sum_mfx), "missing the 'summary_df' data frame")
# })
#
# test_that("plot.summary.marginal_effects: Input validation - coef_index", {
#   expect_error(plot(mfx_object, coef_index = 0), "must be a single integer between")
#   expect_error(plot(mfx_object, coef_index = ncoef + 1), "must be a single integer between")
#   expect_error(plot(mfx_object, coef_index = c(1, 2)), "must be a single integer between")
# })
#
# test_that("plot.summary.marginal_effects: Input validation - ci_level", {
#   expect_error(plot(mfx_object, ci_level = -0.1), "must be a single number between 0 and 1")
#   expect_error(plot(mfx_object, ci_level = 1.1), "must be a single number between 0 and 1")
#   expect_error(plot(mfx_object, ci_level = c(0.9, 0.95)), "must be a single number between 0 and 1")
#
#   # Check missing required quantiles
#   sum_mfx_missing_q <- mfx_object
#   sum_mfx_missing_q$summary_df <- sum_mfx_missing_q$summary_df[, !(names(sum_mfx_missing_q$summary_df) %in% c("q2.5", "q97.5"))]
#   expect_error(plot(sum_mfx_missing_q, ci_level = 0.95), "missing required quantile columns")
# })
#
# test_that("plot.summary.marginal_effects: Handles multiple models and naming", {
#   skip_if_not_installed("ggplot2")
#
#   # Copy for a second model
#   mfx_object_alt <- mfx_object
#
#   # Default naming
#   p_multi <- plot(mfx_object, mfx_object_alt, coef_index = 1)
#   expect_s3_class(p_multi, "ggplot")
#   expect_true("model_name" %in% names(p_multi$data))
#   expect_equal(length(unique(p_multi$data$model_name)), 2)
#   expect_true(all(c("mfx_object", "mfx_object_alt") %in% levels(p_multi$data$model_name)))
#
#   # Explicit naming via ...
#   p_named <- plot(mfx_object, ModelB = mfx_object_alt, coef_index = 1)
#   expect_true(all(c("mfx_object", "ModelB") %in% levels(p_named$data$model_name)))
#
#   # Explicit naming via x_name
#   p_xnamed <- plot(mfx_object, ModelB = mfx_object_alt, coef_index = 1, x_name = "ModelA")
#   expect_true(all(c("ModelA", "ModelB") %in% levels(p_xnamed$data$model_name)))
#
#   # Warns on unnamed objects
#   expect_warning(plot(mfx_object, mfx_object_alt, x_name = "ModelA", coef_index = 1), "Unnamed objects found")
# })
#
# test_that("plot.summary.marginal_effects: Appearance customization works", {
#   skip_if_not_installed("ggplot2")
#   p <- plot(mfx_object, coef_index = 1,
#             xlab = "X Label", ylab = "Y Label", title = "Title",
#             show_legend = FALSE)
#
#   expect_equal(p$labels$x, "X Label")
#   expect_equal(p$labels$y, "Y Label")
#   expect_equal(p$labels$title, "Title")
#   expect_equal(p$theme$legend.position, "none")
#
#   # Check custom colors/linetypes (requires second model)
#   mfx_object_alt <- mfx_object
#   my_colors <- c("mfx_object" = "#FF0000", "mfx_object_alt" = "#0000FF")
#   my_linetypes <- c("mfx_object" = "dashed", "mfx_object_alt" = "dotted")
#   my_fills <- c("mfx_object" = "#FFCCCC", "mfx_object_alt" = "#CCCCFF")
#
#   p_custom <- plot(mfx_object, mfx_object_alt, coef_index = 1,
#                    color_values = my_colors,
#                    linetype_values = my_linetypes,
#                    fill_values = my_fills)
#
#   # Check if scale_manual was added (indirect check)
#   expect_true("ScaleManual" %in% class(p_custom$scales$get_scales("colour")))
#   expect_true("ScaleManual" %in% class(p_custom$scales$get_scales("linetype")))
#   expect_true("ScaleManual" %in% class(p_custom$scales$get_scales("fill")))
# }) 