library(testthat)
library(bayesm.HART)

# ==============================================================================
# Test sim_hier_mnl Function
# ==============================================================================

# ==============================================================================
# Test Case 1: Basic Linear, ncomp = 1 (Step 5)
# ==============================================================================
test_that("Basic linear simulation works correctly", {

  # Parameters for a small test case
  nlgt_test <- 10
  nT_test <- 5
  p_test <- 3
  nz_test <- 2
  nXa_test <- 1
  nXd_test <- 0 # No choice-invariant vars
  const_test <- TRUE
  seed_test <- 42

  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  # Simulate data
  sim_data <- sim_hier_mnl(
    nlgt = nlgt_test,
    nT = nT_test,
    p = p_test,
    nz = nz_test,
    nXa = nXa_test,
    nXd = nXd_test,
    const = const_test,
    beta_func_type = "linear",
    ncomp = 1,
    seed = seed_test,
    standardize_Z = FALSE # Easier to check linear if not standardized
  )

  # --- Check Overall Structure ---
  expect_true(is.list(sim_data))
  expect_named(sim_data, c("p", "lgtdata", "Z", "true_values"))
  expect_equal(sim_data$p, p_test)

  # --- Check Z Dimensions ---
  expect_true(is.matrix(sim_data$Z))
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_equal(ncol(sim_data$Z), nz_test)

  # --- Check lgtdata Structure ---
  expect_true(is.list(sim_data$lgtdata))
  expect_length(sim_data$lgtdata, nlgt_test)
  # Check structure of the first element
  expect_named(sim_data$lgtdata[[1]], c("y", "X", "beta", "betabar"))
  expect_length(sim_data$lgtdata[[1]]$y, nT_test)
  expect_true(is.matrix(sim_data$lgtdata[[1]]$X))
  # Rows = nT * p (alternatives stacked)
  expect_equal(nrow(sim_data$lgtdata[[1]]$X), nT_test * p_test)
  # Cols = ncoef
  expect_equal(ncol(sim_data$lgtdata[[1]]$X), ncoef_test)
  expect_length(sim_data$lgtdata[[1]]$beta, ncoef_test)
  expect_length(sim_data$lgtdata[[1]]$betabar, ncoef_test)

  # --- Check true_values Structure & Dimensions ---
  expect_true(is.list(sim_data$true_values))
  expect_named(sim_data$true_values,
               c("beta_true", "betabar_true", "true_params", "dimensions"))

  expect_true(is.matrix(sim_data$true_values$beta_true))
  expect_equal(dim(sim_data$true_values$beta_true), c(nlgt_test, ncoef_test))

  expect_true(is.matrix(sim_data$true_values$betabar_true))
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  # --- Check true_params & dimensions match input ---
  true_params <- sim_data$true_values$true_params
  dims <- sim_data$true_values$dimensions

  expect_equal(true_params$beta_func_type, "linear")
  expect_named(true_params$beta_func_args, "Delta") # Should have Delta
  expect_equal(dim(true_params$beta_func_args$Delta), c(nz_test, ncoef_test))
  expect_length(true_params$mixture_comps, 1)
  expect_named(true_params$mixture_comps[[1]], c("mu", "rooti"))
  expect_length(true_params$pvec, 1)
  expect_equal(true_params$pvec, 1)

  expect_equal(dims$nlgt, nlgt_test)
  expect_equal(dims$nT, nT_test)
  expect_equal(dims$p, p_test)
  expect_equal(dims$nz, nz_test)
  expect_equal(dims$nXa, nXa_test)
  expect_equal(dims$nXd, nXd_test)
  expect_equal(dims$ncoef, ncoef_test)
  expect_equal(dims$const, const_test)
  expect_equal(dims$ncomp, 1)

  # --- Check relationship beta = betabar + eps (for first person) ---
  # Note: This requires knowing the generated mixture component
  # We know ncomp=1, so eps ~ N(mu, Sigma = (rooti'rooti)^-1)
  mu1 <- true_params$mixture_comps[[1]]$mu
  # Calculating Sigma requires inverting rooti etc., maybe too complex for basic test
  # Instead, check if beta - betabar is roughly centered around mu
  # For ncomp=1 and default sigma_inv_diag=1, mu should be zeros
  expect_equal(mu1, rep(0, ncoef_test))
  # Calculate residuals
  eps_est <- sim_data$true_values$beta_true - sim_data$true_values$betabar_true
  # Check means are somewhat close to mu (zeros)
  # Allow larger tolerance due to small sample size (nlgt=10)
  expect_true(all(abs(colMeans(eps_est) - mu1) < 0.8)) # Increased tol

}) #TEST_THAT basic linear 

# ==============================================================================
# Test Case 2: Step Function (Step 6)
# ==============================================================================
test_that("Step function simulation works correctly", {

  # Parameters
  nlgt_test <- 15
  nT_test <- 2
  p_test <- 2
  nz_test <- 3
  nXa_test <- 1
  nXd_test <- 0
  const_test <- FALSE # Simpler coef structure
  seed_test <- 123
  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  # Step function arguments
  step_args <- list(
    cutoff = 0.1,        # Custom cutoff
    beta_1 = rep(2, ncoef_test), # Value above cutoff
    beta_2 = rep(-2, ncoef_test),# Value below cutoff
    Z_index = 1          # Use first Z variable
  )

  sim_data <- sim_hier_mnl(
    nlgt = nlgt_test,
    nT = nT_test,
    p = p_test,
    nz = nz_test,
    nXa = nXa_test,
    nXd = nXd_test,
    const = const_test,
    beta_func_type = "step",
    beta_func_args = step_args,
    ncomp = 1,
    seed = seed_test,
    standardize_Z = FALSE # Keep Z scale predictable for cutoff
  )

  # --- Basic Checks ---
  expect_true(is.list(sim_data))
  expect_equal(sim_data$p, p_test)
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_equal(ncol(sim_data$Z), nz_test)
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  # --- Check Parameters ---
  true_params <- sim_data$true_values$true_params
  expect_equal(true_params$beta_func_type, "step")
  expect_equal(true_params$beta_func_args, step_args)

  # --- Verify Step Logic ---
  Z_col <- sim_data$Z[, step_args$Z_index]
  betabar <- sim_data$true_values$betabar_true
  for (i in 1:nlgt_test) {
    expected_beta <- if (Z_col[i] > step_args$cutoff) {
      step_args$beta_1
    } else {
      step_args$beta_2
    }
    expect_equal(betabar[i, ], expected_beta)
  }

}) #TEST_THAT step function

# ==============================================================================
# Test Case 3: Friedman Function (Step 6)
# ==============================================================================
test_that("Friedman function simulation works correctly", {

  # Parameters
  nlgt_test <- 15
  nT_test <- 3
  p_test <- 3
  nz_test <- 5 # Friedman needs at least 5
  nXa_test <- 2
  nXd_test <- 1
  const_test <- TRUE
  seed_test <- 444

  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  # Friedman function arguments - specify target coef and expected defaults
  friedman_args <- list(
    coef_index = 3,
    coefs = c(10, 20, 10, 5), # Default Friedman coefficients
    shift = 0,                # Default Friedman shift
    scale = 8                 # Default Friedman scale
  )
  if (friedman_args$coef_index > ncoef_test) {
      skip("Test skipped: coef_index is larger than generated ncoef")
  }

  sim_data <- sim_hier_mnl(
    nlgt = nlgt_test,
    nT = nT_test,
    p = p_test,
    nz = nz_test,
    nXa = nXa_test,
    nXd = nXd_test,
    const = const_test,
    beta_func_type = "friedman",
    beta_func_args = friedman_args,
    ncomp = 1,
    seed = seed_test,
    standardize_Z = TRUE # Friedman often used with standardized inputs
  )

  # --- Basic Checks ---
  expect_true(is.list(sim_data))
  expect_equal(sim_data$p, p_test)
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_equal(ncol(sim_data$Z), nz_test)
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  # --- Check Parameters ---
  true_params <- sim_data$true_values$true_params
  expect_equal(true_params$beta_func_type, "friedman")
  expect_equal(true_params$beta_func_args, friedman_args)

  # --- Verify Friedman Logic ---
  Z_mat <- sim_data$Z
  betabar <- sim_data$true_values$betabar_true
  coef_idx <- friedman_args$coef_index

  # Check the target coefficient matches internal calculation
  expected_friedman_vals <- apply(Z_mat, 1, bayesm.HART:::.beta_Z_friedman)
  expect_equal(betabar[, coef_idx], expected_friedman_vals)

  # Check that *other* coefficients are zero for betabar
  other_indices <- setdiff(1:ncoef_test, coef_idx)
  if (length(other_indices) > 0) {
      expect_true(all(betabar[, other_indices] == 0))
  }

}) #TEST_THAT friedman function

# ==============================================================================
# Test Case 4: Mixture Components (ncomp > 1) (Step 7)
# ==============================================================================
test_that("Mixture component simulation (ncomp > 1) works", {

  # Parameters
  nlgt_test <- 20
  nT_test <- 4
  p_test <- 3
  nz_test <- 2
  nXa_test <- 1
  nXd_test <- 1
  const_test <- TRUE
  seed_test <- 555
  ncomp_test <- 3 # Test with 3 components

  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  sim_data <- sim_hier_mnl(
    nlgt = nlgt_test,
    nT = nT_test,
    p = p_test,
    nz = nz_test,
    nXa = nXa_test,
    nXd = nXd_test,
    const = const_test,
    beta_func_type = "linear", # Keep f(Z) simple
    ncomp = ncomp_test,
    seed = seed_test,
    standardize_Z = TRUE
  )

  # --- Basic Checks ---
  expect_true(is.list(sim_data))
  expect_equal(sim_data$p, p_test)
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_equal(dim(sim_data$true_values$beta_true), c(nlgt_test, ncoef_test))
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  # --- Check Parameters ---
  true_params <- sim_data$true_values$true_params
  dims <- sim_data$true_values$dimensions

  expect_equal(dims$ncomp, ncomp_test)
  expect_length(true_params$mixture_comps, ncomp_test)
  # Check structure of first generated component
  expect_named(true_params$mixture_comps[[1]], c("mu", "rooti"))
  expect_length(true_params$mixture_comps[[1]]$mu, ncoef_test)
  expect_equal(dim(true_params$mixture_comps[[1]]$rooti), c(ncoef_test, ncoef_test))
  # Check pvec
  expect_length(true_params$pvec, ncomp_test)
  expect_equal(sum(true_params$pvec), 1.0)

  # --- Check Residuals (beta - betabar) ---
  eps_est <- sim_data$true_values$beta_true - sim_data$true_values$betabar_true
  expect_equal(dim(eps_est), c(nlgt_test, ncoef_test))
  # Expect residuals not to be exactly zero when ncomp > 1
  # (unless by extreme chance or mu=0 and sigma=0)
  # Use a tolerance check
  expect_false(all(abs(eps_est) < 1e-9))

}) #TEST_THAT mixture components

# ==============================================================================
# Test Case 5: Custom Beta Function (Step 8)
# ==============================================================================
test_that("Custom beta function simulation works", {

  # Parameters
  nlgt_test <- 10
  nT_test <- 2
  p_test <- 2
  nz_test <- 4
  nXa_test <- 2 # ncoef = 2
  nXd_test <- 0
  const_test <- FALSE
  seed_test <- 789
  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  # Define a simple custom function: betabar_i = [ sum(Zi), prod(Zi) ]
  custom_beta_func <- function(Zi) {
    if(length(Zi) != nz_test) stop("Incorrect Zi length passed to custom func")
    out_vec <- c(sum(Zi), prod(Zi[1:min(3, length(Zi))])) # Use max 3 for prod
    # Ensure length matches ncoef
    if(length(out_vec) > ncoef_test) out_vec <- out_vec[1:ncoef_test]
    if(length(out_vec) < ncoef_test) {
        out_vec <- c(out_vec, rep(0, ncoef_test - length(out_vec)))
    }
    return(out_vec)
  }

  custom_args = list(func = custom_beta_func)

  sim_data <- sim_hier_mnl(
    nlgt = nlgt_test,
    nT = nT_test,
    p = p_test,
    nz = nz_test,
    nXa = nXa_test,
    nXd = nXd_test,
    const = const_test,
    beta_func_type = "custom",
    beta_func_args = custom_args,
    ncomp = 1,
    seed = seed_test,
    standardize_Z = FALSE # Keep Z scale simpler for checking
  )

  # --- Basic Checks ---
  expect_true(is.list(sim_data))
  expect_equal(sim_data$p, p_test)
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  # --- Check Parameters ---
  true_params <- sim_data$true_values$true_params
  expect_equal(true_params$beta_func_type, "custom")
  # Can't directly compare functions, check it exists
  expect_true(is.function(true_params$beta_func_args$func))

  # --- Verify Custom Logic ---
  Z_mat <- sim_data$Z
  betabar <- sim_data$true_values$betabar_true

  # Calculate expected betabar using the custom function
  expected_betabar <- t(apply(Z_mat, 1, custom_beta_func))

  expect_equal(betabar, expected_betabar)

}) #TEST_THAT custom function

# ==============================================================================
# Test Case 6: No Z variables (nz = 0)
# ==============================================================================
test_that("Simulation works with nz = 0", {

  # Parameters
  nlgt_test <- 8
  nT_test <- 3
  p_test <- 4
  nz_test <- 0 # Explicitly testing nz = 0
  nXa_test <- 2
  nXd_test <- 1
  const_test <- TRUE
  seed_test <- 111
  ncomp_test <- 2 # Use multiple components for heterogeneity

  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  # Simulate data
  sim_data <- sim_hier_mnl(
    nlgt = nlgt_test,
    nT = nT_test,
    p = p_test,
    nz = nz_test, # Set to zero
    nXa = nXa_test,
    nXd = nXd_test,
    const = const_test,
    beta_func_type = "linear", # Should be ignored
    ncomp = ncomp_test,
    seed = seed_test
  )

  # --- Check Overall Structure ---
  expect_true(is.list(sim_data))
  expect_named(sim_data, c("p", "lgtdata", "Z", "true_values"))
  expect_equal(sim_data$p, p_test)

  # --- Check Z is NULL ---
  expect_null(sim_data$Z)

  # --- Check lgtdata Structure ---
  expect_true(is.list(sim_data$lgtdata))
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_named(sim_data$lgtdata[[1]], c("y", "X", "beta", "betabar"))

  # --- Check true_values Structure & Dimensions ---
  expect_true(is.list(sim_data$true_values))
  expect_named(sim_data$true_values,
               c("beta_true", "betabar_true", "true_params", "dimensions"))

  # --- Check betabar_true is all zeros ---
  expect_true(is.matrix(sim_data$true_values$betabar_true))
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))
  expect_true(all(sim_data$true_values$betabar_true == 0))

  # --- Check beta_true is not all zeros (due to mixture) ---
  expect_true(is.matrix(sim_data$true_values$beta_true))
  expect_equal(dim(sim_data$true_values$beta_true), c(nlgt_test, ncoef_test))
  expect_false(all(abs(sim_data$true_values$beta_true) < 1e-9))

  # --- Check beta_true = betabar_true + eps ---
  # This means beta_true should equal the residuals
  eps_est <- sim_data$true_values$beta_true - sim_data$true_values$betabar_true
  expect_equal(sim_data$true_values$beta_true, eps_est)

  # --- Check true_params & dimensions match input ---
  true_params <- sim_data$true_values$true_params
  dims <- sim_data$true_values$dimensions

  expect_equal(true_params$beta_func_type, "linear") # Still stored
  expect_equal(true_params$beta_func_args, list()) # Args ignored, stored as empty
  expect_length(true_params$mixture_comps, ncomp_test)
  expect_length(true_params$pvec, ncomp_test)

  expect_equal(dims$nz, 0)
  expect_equal(dims$ncomp, ncomp_test)

}) #TEST_THAT nz = 0 

# ==============================================================================
# Test Case 7: User-provided mixture components
# ==============================================================================
test_that("User-provided mixture components work", {

  # Parameters
  nlgt_test <- 5
  nT_test <- 2
  p_test <- 3 
  nz_test <- 1
  nXa_test <- 1 
  nXd_test <- 0
  const_test <- TRUE
  seed_test <- 666
  ncomp_test <- 2 
  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test # Should be 3
  expect_equal(ncoef_test, 3) # Verify ncoef calc

  # Define the components to pass
  user_comps <- list(
    list(mu = c(1, 0, 0), rooti = diag(1.0, ncoef_test)),
    list(mu = c(-1, 1, 0), rooti = diag(2.0, ncoef_test))
  )
  
  # Simulate data, passing the components
  sim_data <- NULL
  expect_no_error({
      sim_data <- sim_hier_mnl(
        nlgt = nlgt_test,
        nT = nT_test,
        p = p_test,
        nz = nz_test,
        nXa = nXa_test,
        nXd = nXd_test,
        const = const_test,
        ncomp = ncomp_test,
        mixture_comps = user_comps, # Pass the user list
        seed = seed_test
      )
  })
  
  # --- Check if simulation ran --- 
  expect_true(!is.null(sim_data), "Simulation failed to run with user components.")
  
  # --- Check true_params match input --- 
  true_params <- sim_data$true_values$true_params
  expect_equal(true_params$mixture_comps, user_comps)
  expect_equal(sim_data$true_values$dimensions$ncomp, ncomp_test)

}) #TEST_THAT user mixture components 