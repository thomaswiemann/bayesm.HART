# ==============================================================================
# Unit Tests for marginal_effects.rhierMnlRwMixture
# ==============================================================================

library(testthat)

# ==============================================================================
# Short fit shared by the tests
# ==============================================================================

run_quick_model_for_mfx_test <- function(R = 50, keep = 1, useBART = FALSE) {
  set.seed(654)
  sim_data <- sim_hier_mnl(
    nlgt = 25, nT = 8, p = 3, nz = 5, nXa = 1, nXd = 0, const = TRUE,
    beta_func_type = "linear",
    beta_func_args = list()
  )
  dims <- sim_data$true_values$dimensions
  ncoef <- dims$ncoef
  nz <- dims$nz

  Prior <- list(ncomp = 1)
  if (useBART) {
    Prior$hart_model <- "normalized"
    Prior$bart <- list(num_trees = 10, numcut = 10)
  } else {
    Prior$Deltabar <- matrix(0, nz * ncoef, 1)
    Prior$Ad <- diag(nz * ncoef) * 0.01
  }

  Mcmc <- list(R = R, keep = keep, nprint = 0)

  fit <- NULL
  capture.output({
    fit <- rhierMnlRwMixture(Data = sim_data, Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  if(is.null(fit$betadraw)) stop("Model fitting failed for mfx test setup")

  fit$nz <- nz

  return(list(fit = fit, sim_data = sim_data))
}

# ==============================================================================
# Shared fit and grids
# ==============================================================================

n_draws_test <- 30 
quick_fit_bart <- run_quick_model_for_mfx_test(R = n_draws_test, useBART = TRUE)
fit_object <- quick_fit_bart$fit
orig_data <- quick_fit_bart$sim_data

dims_actual <- orig_data$true_values$dimensions
ncoef <- dims_actual$ncoef
nz <- dims_actual$nz
nlgt <- dims_actual$nlgt
ndraws_total <- n_draws_test

base_Z <- orig_data$Z

n_grid_points_test <- 5

z_col1_vec <- seq(-1, 1, length.out = n_grid_points_test)
z_values_1d <- matrix(NA_real_, nrow = n_grid_points_test, ncol = nz)
z_values_1d[, 1] <- z_col1_vec
if (!is.null(colnames(base_Z))) {
    colnames(z_values_1d) <- colnames(base_Z)
}

# One row, all NA.
z_values_all_na <- matrix(NA_real_, nrow = 1, ncol = nz)
if (!is.null(colnames(base_Z))) {
    colnames(z_values_all_na) <- colnames(base_Z)
}

z_values_violating_na <- matrix(NA_real_, nrow = 3, ncol = nz)
z_values_violating_na[, 1] <- c(1, 2, 3)
z_values_violating_na[1:2, 2] <- c(4, 5)
if (!is.null(colnames(base_Z))) {
    colnames(z_values_violating_na) <- colnames(base_Z)
}

burn_val <- 5
ndraws_use <- ndraws_total - burn_val

# ==============================================================================
# Test Cases
# ==============================================================================

test_that("marginal_effects: Basic functionality and output structure (BART fit)", {

  mfx <- marginal_effects(
    object = fit_object,
    z_values = z_values_1d,
    Z = base_Z[1:10, ],
    burn = burn_val,
    r_verbose = FALSE
  )

  expect_s3_class(mfx, "marginal_effects")
  expect_type(mfx, "list")
  expect_named(mfx, c("z_values", "avg_betabar_draws", "call", "burn"))

  expect_equal(mfx$z_values, z_values_1d)
  expect_type(mfx$avg_betabar_draws, "list")
  expect_length(mfx$avg_betabar_draws, nrow(z_values_1d))
  expect_equal(mfx$burn, burn_val)
  expect_true(!is.null(mfx$call))
})

test_that("marginal_effects: avg_betabar_draws elements have correct dimensions (BART fit)", {
  
  mfx <- marginal_effects(fit_object, z_values = z_values_1d, Z = base_Z, burn = burn_val, r_verbose = FALSE)

  expect_true(length(mfx$avg_betabar_draws) > 0, info = "avg_betabar_draws list should not be empty")
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

test_that("marginal_effects: Handles z_values edge cases", {

  z_values_1row <- z_values_1d[1, , drop = FALSE]
  mfx_1row <- marginal_effects(fit_object, z_values = z_values_1row, Z = base_Z, burn = burn_val, r_verbose = FALSE)
  expect_equal(nrow(mfx_1row$z_values), 1)
  expect_length(mfx_1row$avg_betabar_draws, 1)
  expect_true(is.matrix(mfx_1row$avg_betabar_draws[[1]]))
  expect_equal(ncol(mfx_1row$avg_betabar_draws[[1]]), ndraws_use)

  expect_error(
      marginal_effects(fit_object, z_values = z_values_all_na, Z = base_Z, burn = burn_val, r_verbose = FALSE),
      "'z_values' must contain at least one non-NA column to define grid points.",
      fixed = TRUE
  )
})

test_that("marginal_effects: Handles burn edge cases", {

  mfx_zero_burn <- marginal_effects(fit_object, z_values=z_values_1d, Z=base_Z, burn = 0, r_verbose = FALSE)
  expect_equal(mfx_zero_burn$burn, 0)
  expect_true(is.matrix(mfx_zero_burn$avg_betabar_draws[[1]]))
  expect_equal(ncol(mfx_zero_burn$avg_betabar_draws[[1]]), ndraws_total)

  if (ndraws_total > 1) {
      mfx_max_burn <- marginal_effects(fit_object, z_values=z_values_1d, Z=base_Z, burn = ndraws_total - 1, r_verbose = FALSE)
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

mfx_object <- marginal_effects(
    object = fit_object,
    z_values = z_values_1d,
    Z = base_Z[1:10, , drop = FALSE],
    burn = burn_val,
    r_verbose = FALSE
  )

sum_mfx_1d <- summary(mfx_object)

# ==============================================================================
# Test Cases for summary.marginal_effects
# ==============================================================================

test_that("summary.marginal_effects: Basic structure and class", {
  mfx_summary <- summary(mfx_object)

  expect_s3_class(mfx_summary, "summary.marginal_effects")
  expect_type(mfx_summary, "list")
  expect_named(mfx_summary, c("summary_df", "z_values", "probs", "call"))

  expect_equal(mfx_summary$z_values, mfx_object$z_values)
  expect_true(!is.null(mfx_summary$call))
})

test_that("summary.marginal_effects: summary_df structure and content", {
  default_probs = c(0.025, 0.05, 0.5, 0.95, 0.975)
  expected_q_names <- paste0("q", gsub("\\.", "", default_probs * 100))

  non_na_cols_idx <- which(colSums(is.na(z_values_1d)) < nrow(z_values_1d))
  z_cols_expected_names <- colnames(z_values_1d)[non_na_cols_idx]
  if (is.null(z_cols_expected_names)) {
      z_cols_expected_names <- paste0("z_val_", non_na_cols_idx)
  }

  expected_colnames <- c("coefficient_index", z_cols_expected_names, "mean", expected_q_names)
  
  mfx_summary <- summary(mfx_object)
  df <- mfx_summary$summary_df
  
  expect_s3_class(df, "data.frame")
  expect_named(df, expected_colnames, ignore.order = FALSE)
  expect_equal(nrow(df), nrow(z_values_1d) * ncoef)
  
  expect_type(df$coefficient_index, "integer")
  if (length(z_cols_expected_names) > 0) {
      expect_type(df[[z_cols_expected_names[1]]], "double")
  }
  expect_type(df$mean, "double")
  expect_type(df[[expected_q_names[1]]], "double")
})

test_that("summary.marginal_effects: Custom probs argument works", {
  custom_probs <- c(0.1, 0.9)
  expected_q_names <- paste0("q", custom_probs * 100)

  non_na_cols_idx <- which(colSums(is.na(z_values_1d)) < nrow(z_values_1d))
  z_cols_expected_names <- colnames(z_values_1d)[non_na_cols_idx]
  if (is.null(z_cols_expected_names)) {
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
  
  expect_error(
    summary(mfx_object, probs = c(0.1, 1.1)),
    "'probs' must be a numeric vector with values between 0 and 1."
  )
  expect_error(
    summary(mfx_object, probs = "a"),
    "'probs' must be a numeric vector with values between 0 and 1."
  )
  
  bad_mfx_object <- mfx_object
  bad_mfx_object$z_values <- NULL
  expect_error(
    summary(bad_mfx_object),
    "Invalid 'marginal_effects' object. Missing elements: z_values"
  )

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
  
  p <- plot(sum_mfx_1d, coef_index = 1)
  
  expect_s3_class(p, "ggplot")
  
  expected_x_var <- "z_val_1"
  expect_equal(rlang::quo_name(p$mapping$x), expected_x_var)
  
  expect_equal(p$labels$x, expected_x_var)
  
  expect_true(grepl("Coefficient 1", p$labels$y))
  expect_true(grepl(paste("vs", expected_x_var), p$labels$title))
})

test_that("plot.summary.marginal_effects: Plotting with specified x-axis works", {
  
  skip_if_not_installed("ggplot2")
  
  plot_var <- "z_val_1"
  p <- plot(sum_mfx_1d, coef_index = 1, plot_axis_var = plot_var)
  
  expect_s3_class(p, "ggplot")
  
  expect_equal(rlang::quo_name(p$mapping$x), plot_var)
  
  expect_equal(p$labels$x, plot_var)
})

test_that("plot.summary.marginal_effects: Plotting multiple models works", {
  
  skip_if_not_installed("ggplot2")
  
  sum_mfx_1d_alt <- sum_mfx_1d
  
  p <- plot(sum_mfx_1d, linear = sum_mfx_1d_alt, coef_index = 1)
  
  expect_s3_class(p, "ggplot")
  
  expect_true("model_name" %in% names(p$data))
  expect_equal(length(unique(p$data$model_name)), 2)
  
  expected_x_var <- "z_val_1"
  expect_equal(rlang::quo_name(p$mapping$x), expected_x_var)
})
