# Tests for predict.rhierMnlRwMixture

# Helper function to run model quickly for tests
run_model_for_pred_test <- function(sim_data, R = 50, keep = 1, useBART = FALSE) {
  # Check dimensions structure (now checking within true_values)
  if (is.null(sim_data$true_values)) stop("sim_data$true_values is NULL")
  if (is.null(sim_data$true_values$dimensions)) stop("sim_data$true_values$dimensions is NULL")
  dims <- sim_data$true_values$dimensions
  if (is.null(dims$nz)) stop("sim_data$true_values$dimensions$nz is NULL")
  if (!is.numeric(dims$nz) || length(dims$nz) != 1) {
      stop(paste("sim_data$true_values$dimensions$nz is not a single number. Value:", 
                 paste(dims$nz, collapse=", ")))
  }
  if (is.null(dims$nlgt)) stop("sim_data$true_values$dimensions$nlgt is NULL")
  if (is.null(dims$ncoef)) stop("sim_data$true_values$dimensions$ncoef is NULL")
  
  nlgt <- dims$nlgt
  ncoef <- dims$ncoef
  nz <- dims$nz
  
  Prior <- list(ncomp = 1) # Single component for simplicity in tests
  if (nz > 0) { 
      Prior$Deltabar <- matrix(0, nz*ncoef, 1)
      Prior$Ad <- diag(nz*ncoef) * 0.01 # Weak prior
  }
  # Use bart settings only if requested
  if (useBART && nz > 0) { 
    Prior$bart <- list(num_trees = 20) # Small number of trees
  } else {
    Prior$bart <- NULL # Ensure BART is off otherwise
  }
  
  Mcmc <- list(R = R, keep = keep, nprint = 0) # Suppress printing
  
  # Fit model
  fit <- rhierMnlRwMixture(Data = sim_data, Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  return(fit)
}

# Helper to run models quickly - suppress messages
run_quick_hmnl <- function(data, prior, mcmc) {
  capture.output({
    fit <- rhierMnlRwMixture(Data = data, Prior = prior, Mcmc = mcmc)
  })
  return(fit)
}

test_that("predict works for linear model (nz > 0)", {
  # 1. Simulate data (linear beta function)
  set.seed(123)
  sim_data_linear <- sim_hier_mnl(
    nlgt = 100, nT = 10, p = 3, nz = 2, nXa = 1, nXd = 0, const = TRUE,
    beta_func_type = "linear"
  )
  dims <- sim_data_linear$true_values$dimensions 
  
  # 2. Fit linear model (Capture output)
  n_draws <- 1000
  fit_linear <- NULL
  capture.output({
      fit_linear <- run_model_for_pred_test(sim_data_linear, R = n_draws, useBART = FALSE)
  })
  expect_false(is.null(fit_linear), info="Model fitting should succeed")
  
  # 3. Predict with type="DeltaZ"
  pred_deltaz <- predict(fit_linear, newdata = list(Z = sim_data_linear$Z), type = "DeltaZ", r_verbose=FALSE)
  
  # Check dimensions
  expect_equal(dim(pred_deltaz), c(dims$nlgt, dims$ncoef, n_draws))
  
  # Check plausibility (mean prediction vs true betabar)
  mean_pred_deltaz <- apply(pred_deltaz, c(1, 2), mean)
  corr_deltaz <- cor(as.vector(mean_pred_deltaz), as.vector(sim_data_linear$true_values$betabar_true))
  # Correlation should be reasonably high if model captures the linear trend
  expect_true(corr_deltaz > 0.5, 
              label = paste("Correlation between mean DeltaZ pred and true betabar_true:", round(corr_deltaz, 3)))

  # 4. Predict with type="DeltaZ+mu"
  pred_deltaz_mu <- predict(fit_linear, newdata = list(Z = sim_data_linear$Z), type = "DeltaZ+mu", r_verbose=FALSE)

  # Check dimensions
  expect_equal(dim(pred_deltaz_mu), c(dims$nlgt, dims$ncoef, n_draws))
  
  # Check plausibility (mean prediction vs true beta)
  mean_pred_deltaz_mu <- apply(pred_deltaz_mu, c(1, 2), mean)
  corr_beta <- cor(as.vector(mean_pred_deltaz_mu), as.vector(sim_data_linear$true_values$beta_true))
  # Correlation should be reasonably high 
  expect_true(corr_beta > 0.5, 
              label = paste("Correlation between mean DeltaZ+mu pred and true beta_true:", round(corr_beta, 3)))

})

test_that("predict works for BART model (nz > 0)", {
  # 1. Simulate data (Friedman beta function requires nz >= 5)
  set.seed(456)
  sim_data_bart <- sim_hier_mnl(
    nlgt = 100, nT = 10, p = 3, nz = 5, nXa = 1, nXd = 0, const = TRUE,
    beta_func_type = "friedman",
    beta_func_args = list(coef_index = 1) # Apply Friedman to 1st coef
  )
  dims <- sim_data_bart$true_values$dimensions
  
  # 2. Fit BART model (Capture output)
  n_draws <- 100 
  fit_bart <- NULL
  capture.output({
      fit_bart <- run_model_for_pred_test(sim_data_bart, R = n_draws, useBART = TRUE)
  })
  expect_false(is.null(fit_bart), info="BART model fitting should succeed")

  # Check if BART models were actually produced
  expect_true(!is.null(fit_bart$bart_models), info = "BART models should exist in the output")
  expect_equal(length(fit_bart$bart_models), dims$ncoef, info = "Number of BART models should match ncoef")

  # 3. Predict with type="DeltaZ"
  pred_deltaz <- predict(fit_bart, newdata = list(Z = sim_data_bart$Z), type = "DeltaZ", r_verbose=FALSE)
  
  # Check dimensions
  expect_equal(dim(pred_deltaz), c(dims$nlgt, dims$ncoef, n_draws))
  
  # Check plausibility (mean prediction vs true betabar)
  # Note: Correlation might be lower for BART with few draws
  mean_pred_deltaz <- apply(pred_deltaz, c(1, 2), mean)
  # Check only the coefficient affected by Friedman for stronger correlation
  coef_idx <- sim_data_bart$true_values$true_params$beta_func_args$coef_index
  corr_deltaz <- cor(mean_pred_deltaz[, coef_idx], 
                     sim_data_bart$true_values$betabar_true[, coef_idx])
  expect_true(corr_deltaz > 0.3, # Lower threshold for BART with short run
              label = paste("Correlation (coef ", coef_idx, ") between mean DeltaZ pred and true betabar_true:", round(corr_deltaz, 3)))

  # 4. Predict with type="DeltaZ+mu"
  pred_deltaz_mu <- predict(fit_bart, newdata = list(Z = sim_data_bart$Z), type = "DeltaZ+mu", r_verbose=FALSE)

  # Check dimensions
  expect_equal(dim(pred_deltaz_mu), c(dims$nlgt, dims$ncoef, n_draws))
  
  # Check plausibility (mean prediction vs true beta)
  mean_pred_deltaz_mu <- apply(pred_deltaz_mu, c(1, 2), mean)
  corr_beta <- cor(mean_pred_deltaz_mu[, coef_idx], 
                     sim_data_bart$true_values$beta_true[, coef_idx])
  expect_true(corr_beta > 0.3, # Lower threshold for BART with short run
              label = paste("Correlation (coef ", coef_idx, ") between mean DeltaZ+mu pred and true beta_true:", round(corr_beta, 3)))

})

test_that("predict handles edge cases (burn-in, different npred)", {
  # --- Setup: Use the linear model fit from the first test --- 
  # (Regenerate quickly here for self-containment & capture output)
  set.seed(123)
  sim_data_linear <- sim_hier_mnl(
    nlgt = 20, nT = 10, p = 3, nz = 2, nXa = 1, nXd = 0, const = TRUE,
    beta_func_type = "linear"
  )
  dims <- sim_data_linear$true_values$dimensions 
  n_draws <- 50 
  fit_linear <- NULL
  capture.output({
      fit_linear <- run_model_for_pred_test(sim_data_linear, R = n_draws, useBART = FALSE)
  })
  expect_false(is.null(fit_linear), info="Model fitting for edge cases should succeed")
  
  # --- Test burn-in parameter --- 
  burn_in <- 10
  expect_true(burn_in < n_draws, info = "Ensure burn_in is less than n_draws for test validity")
  
  pred_burned <- predict(fit_linear, newdata = list(Z = sim_data_linear$Z), 
                         type = "DeltaZ", burn = burn_in, r_verbose=FALSE)
  
  # Check dimension after burn-in
  expect_equal(dim(pred_burned), c(dims$nlgt, dims$ncoef, n_draws - burn_in),
               info = "Third dimension should be n_draws - burn_in")

  # Check error for too large burn-in
  expect_error(predict(fit_linear, newdata = list(Z = sim_data_linear$Z), burn = n_draws), 
               regexp = "burn must be >= 0 and < ndraws_total")
  expect_error(predict(fit_linear, newdata = list(Z = sim_data_linear$Z), burn = n_draws + 5),
               regexp = "burn must be >= 0 and < ndraws_total")

  # --- Test prediction with different number of rows in newdata$Z --- 
  npred_new <- 5
  expect_true(npred_new < dims$nlgt, info = "Ensure npred_new is less than nlgt for test validity")
  Z_new <- sim_data_linear$Z[1:npred_new, , drop = FALSE]
  
  pred_new_n <- predict(fit_linear, newdata = list(Z = Z_new), type = "DeltaZ", r_verbose=FALSE)
  
  # Check dimension for new number of predictions
  expect_equal(dim(pred_new_n), c(npred_new, dims$ncoef, n_draws),
               info = "First dimension should match nrow(newdata$Z)")
  
  # Also check with burn-in
  pred_new_n_burned <- predict(fit_linear, newdata = list(Z = Z_new), 
                               type = "DeltaZ", burn = burn_in, r_verbose=FALSE)
  expect_equal(dim(pred_new_n_burned), c(npred_new, dims$ncoef, n_draws - burn_in),
                info = "Dimensions should reflect new npred and burn-in")

})

test_that("predict works for model with nz=0", {
  # 1. Simulate data with nz=0
  set.seed(789)
  sim_data_no_z <- sim_hier_mnl(
    nlgt = 15, nT = 10, p = 3, nz = 0, nXa = 1, nXd = 1, const = TRUE
  )
  dims <- sim_data_no_z$true_values$dimensions 
  
  # 2. Fit model (Capture output)
  n_draws <- 55
  fit_no_z <- NULL
  capture.output({
      fit_no_z <- run_model_for_pred_test(sim_data_no_z, R = n_draws, useBART = FALSE)
  })
  expect_false(is.null(fit_no_z), info="Model fitting (nz=0) should succeed")
  
  # Check that Deltadraw and bart_models are NULL
  expect_null(fit_no_z$Deltadraw, info = "Deltadraw should be NULL for nz=0")
  expect_null(fit_no_z$bart_models, info = "bart_models should be NULL for nz=0")

  # 3. Predict with type="DeltaZ"
  # newdata should not be required
  pred_deltaz_no_z <- predict(fit_no_z, type = "DeltaZ", r_verbose=FALSE)
  
  # Check dimensions - should be 1 x ncoef x ndraws
  expect_equal(dim(pred_deltaz_no_z), c(1, dims$ncoef, n_draws))
  
  # Check values - should be all zeros
  expect_true(all(pred_deltaz_no_z == 0), info = "Predictions for DeltaZ with nz=0 should be zero")

  # 4. Predict with type="DeltaZ+mu"
  pred_deltaz_mu_no_z <- predict(fit_no_z, type = "DeltaZ+mu", r_verbose=FALSE)

  # Check dimensions - should be 1 x ncoef x ndraws
  expect_equal(dim(pred_deltaz_mu_no_z), c(1, dims$ncoef, n_draws))
  
  # Check values - should be equal to the mu draws
  # Extract mu draws from the fit object (assuming ncomp=1 from helper)
  mu_draws <- sapply(fit_no_z$nmix$compdraw, function(comp) comp[[1]]$mu)
  # Reshape mu_draws to match prediction dimensions [1, ncoef, ndraws]
  expected_mu_array <- array(mu_draws, dim = c(dims$ncoef, n_draws, 1))
  expected_mu_array <- aperm(expected_mu_array, c(3, 1, 2))
  
  expect_equal(pred_deltaz_mu_no_z, expected_mu_array, 
               info = "Predictions for DeltaZ+mu with nz=0 should match mu draws")

})

test_that("predict method with type = 'posterior_probs' works", {

  # --- Setup Mock Data ---
  set.seed(1001)
  nlgt_test <- 2
  nvar_test <- 3
  ndraws_test <- 10
  p_test <- 4
  T1_test <- 1
  T2_test <- 2

  # Mock rhierMnlRwMixture object
  mock_object <- list()
  mock_object$betadraw <- array(rnorm(nlgt_test * nvar_test * ndraws_test), 
                                dim = c(nlgt_test, nvar_test, ndraws_test))
  class(mock_object) <- "rhierMnlRwMixture" # Add class for method dispatch

  # Mock newdata for T_i=1 case
  mock_newdata_t1 <- list(
    p = p_test,
    nlgtdata = list(
      list(X = matrix(rnorm(T1_test * p_test * nvar_test), 
                      nrow = T1_test * p_test, ncol = nvar_test)),
      list(X = matrix(rnorm(T1_test * p_test * nvar_test), 
                      nrow = T1_test * p_test, ncol = nvar_test))
    )
  )
  
  # Mock newdata for T_i>1 case (using nlgt=1 for simplicity here)
  nlgt_test_t2 <- 1
  ndraws_test_t2 <- 5
  mock_object_t2 <- list()
  mock_object_t2$betadraw <- array(rnorm(nlgt_test_t2 * nvar_test * ndraws_test_t2), 
                                    dim = c(nlgt_test_t2, nvar_test, ndraws_test_t2))
  class(mock_object_t2) <- "rhierMnlRwMixture"
  mock_newdata_t2 <- list(
    p = p_test,
    nlgtdata = list(
       list(X = matrix(rnorm(T2_test * p_test * nvar_test), 
                       nrow = T2_test * p_test, ncol = nvar_test))
    )
  )

  # --- Test Case 1: Basic Functionality (T_i = 1) ---
  pred_t1 <- predict(mock_object, newdata = mock_newdata_t1, type = "posterior_probs", r_verbose=FALSE)
  
  # Check output structure
  expect_type(pred_t1, "list")
  expect_length(pred_t1, nlgt_test)
  expect_true(is.array(pred_t1[[1]]), info = "output[[1]] should be an array (T=1)")
  expect_equal(dim(pred_t1[[1]]), c(T1_test, p_test, ndraws_test))
  expect_equal(dim(pred_t1[[2]]), c(T1_test, p_test, ndraws_test))
  
  # Check calculation for one unit/draw
  # Manually calculate using the helper (assuming it's loaded)
  beta_1_1 <- mock_object$betadraw[1, , 1]
  X_1 <- mock_newdata_t1$nlgtdata[[1]]$X
  expected_prob_1_1 <- bayesm.HART:::calculate_mnl_probs_from_beta(X_1, beta_1_1, p_test)
  
  # Use drop=FALSE when subsetting the 3D array to keep dimensions
  expect_equal(pred_t1[[1]][,,1, drop = FALSE], array(expected_prob_1_1, dim=c(T1_test, p_test, 1)), 
               info = "Probability calculation mismatch for unit 1, draw 1 (T=1)")
               
  # Check probabilities sum to 1
  expect_equal(sum(pred_t1[[1]][1,,1]), 1.0, tolerance = 1e-6, 
               info = "Probabilities do not sum to 1 for unit 1, draw 1 (T=1)")
  expect_equal(sum(pred_t1[[2]][1,,5]), 1.0, tolerance = 1e-6, 
               info = "Probabilities do not sum to 1 for unit 2, draw 5 (T=1)")

  # --- Test Case 2: Multi-Observation Functionality (T_i > 1) ---
  pred_t2 <- predict(mock_object_t2, newdata = mock_newdata_t2, type = "posterior_probs", r_verbose=FALSE)
  
  # Check output structure
  expect_type(pred_t2, "list")
  expect_length(pred_t2, nlgt_test_t2)
  expect_true(is.array(pred_t2[[1]]), info = "output[[1]] should be an array (T>1)")
  expect_equal(dim(pred_t2[[1]]), c(T2_test, p_test, ndraws_test_t2))
  
  # Check calculation for one unit/obs/draw
  beta_1_3_t2 <- mock_object_t2$betadraw[1, , 3]
  X_1_t2 <- mock_newdata_t2$nlgtdata[[1]]$X
  expected_prob_1_3_t2 <- bayesm.HART:::calculate_mnl_probs_from_beta(X_1_t2, beta_1_3_t2, p_test)
  
  # Use drop=FALSE when subsetting the 3D array
  expect_equal(pred_t2[[1]][,,3, drop = FALSE], array(expected_prob_1_3_t2, dim=c(T2_test, p_test, 1)), 
               info = "Probability calculation mismatch for unit 1, draw 3 (T>1)")
               
  # Check probabilities sum to 1 for one observation
  # Subset with drop=FALSE to get matrix, then apply rowSums
  probs_draw3_t2 <- pred_t2[[1]][,,3, drop = FALSE][,,1] # Extract the matrix for draw 3
  expect_equal(rowSums(probs_draw3_t2), rep(1.0, T2_test), tolerance = 1e-6, 
               info = "Probabilities do not sum to 1 for unit 1, draw 3 (T>1)")
  probs_draw5_t2 <- pred_t2[[1]][,,5, drop = FALSE][,,1] # Extract the matrix for draw 5
  expect_equal(rowSums(probs_draw5_t2), rep(1.0, T2_test), tolerance = 1e-6, 
               info = "Probabilities do not sum to 1 for unit 1, draw 5 (T>1)")

  # --- Test Case 3: Ignored newdata$Z ---
  mock_newdata_t1_with_Z <- mock_newdata_t1
  mock_newdata_t1_with_Z$Z <- matrix(0, nrow=nlgt_test, ncol=1) # Add dummy Z
  
  expect_no_error({ 
      pred_t1_Z <- predict(mock_object, newdata = mock_newdata_t1_with_Z, type = "posterior_probs", r_verbose=FALSE) 
  }, message="Predicting with type='posterior_probs' and present Z failed.")
  
  # Compare results with and without Z
  expect_identical(pred_t1, pred_t1_Z, 
                   info="Result differs when newdata$Z is present for type='posterior_probs'")

}) 

test_that("predict method with type = 'prior_probs' works", {

  # --- Setup --- 
  set.seed(1101)
  nlgt_fit <- 10 # Number of units used for fitting
  npred_test <- 10
  nT1_test <- 1
  nT2_test <- 2
  # --- ADJUSTED PARAMETERS ---
  p_test <- 2       # Change: Use binary choice for simplicity
  ncoef_test <- 2   # Target number of coefficients
  nXa_test <- 1     # Derived: Need 1 Xa to get ncoef=2 with p=2, const=T
  nXd_test <- 0     # Be explicit
  const_test <- TRUE # Be explicit
  # --- END ADJUSTMENTS ---
  nz_test <- 1   # Use nz=1 for Z tests
  R_test <- 10
  R_burn_test <- 20 # For burn-in test
  keep_test <- 1
  ndraws_out <- R_test / keep_test # = 10
  ndraws_burn_out <- (R_burn_test - 10) / keep_test # = 10 (burn first 10)
  ncomp_test <- 1
  num_trees_test <- 5

  # Fit Minimal Models (suppress messages)
  # 1. Linear with Z
  sim_lin_Z <- sim_hier_mnl(nlgt = nlgt_fit, nT = nT2_test, p = p_test, nz = nz_test, 
                            nXa = nXa_test, nXd = nXd_test, const = const_test,
                            seed = 1102)
  fit_lin_Z <- run_quick_hmnl(sim_lin_Z, list(ncomp = ncomp_test), 
                              list(R = R_test, keep = keep_test, nprint=0))
  fit_lin_Z_R20 <- run_quick_hmnl(sim_lin_Z, list(ncomp = ncomp_test), 
                                  list(R = R_burn_test, keep = keep_test, nprint=0))
  ncoef_lin_Z <- dim(fit_lin_Z$betadraw)[2] # Should be ncoef_test
  expect_equal(ncoef_lin_Z, ncoef_test, info="ncoef mismatch for fit_lin_Z")

  # 2. BART with Z
  sim_bart_Z <- sim_lin_Z # Can use same data
  fit_bart_Z <- run_quick_hmnl(sim_bart_Z, list(ncomp = ncomp_test, bart = list(num_trees = num_trees_test)), 
                               list(R = R_test, keep = keep_test, nprint=0))
  ncoef_bart_Z <- dim(fit_bart_Z$betadraw)[2] # Should be ncoef_test
  expect_equal(ncoef_bart_Z, ncoef_test, info="ncoef mismatch for fit_bart_Z")

  # 3. Linear without Z
  sim_lin_noZ <- sim_hier_mnl(nlgt = nlgt_fit, nT = nT2_test, p = p_test, nz = 0, 
                              nXa = nXa_test, nXd = nXd_test, const = const_test,
                              seed = 1103)
  fit_lin_noZ <- run_quick_hmnl(sim_lin_noZ, list(ncomp = ncomp_test), 
                                list(R = R_test, keep = keep_test, nprint=0))
  ncoef_lin_noZ <- dim(fit_lin_noZ$betadraw)[2] # Should be ncoef_test
  expect_equal(ncoef_lin_noZ, ncoef_test, info="ncoef mismatch for fit_lin_noZ")

  # 4. BART without Z (should run as linear internally)
  sim_bart_noZ <- sim_lin_noZ
  fit_bart_noZ <- run_quick_hmnl(sim_bart_noZ, list(ncomp = ncomp_test, bart = list(num_trees = num_trees_test)), 
                                 list(R = R_test, keep = keep_test, nprint=0))
  ncoef_bart_noZ <- dim(fit_bart_noZ$betadraw)[2] # Should be ncoef_test
  expect_equal(ncoef_bart_noZ, ncoef_test, info="ncoef mismatch for fit_bart_noZ")

  # Create newdata structures
  # newdata for models WITH Z
  newdata_Z <- list(
    p = p_test,
    Z = matrix(rnorm(npred_test * nz_test), nrow = npred_test, ncol = nz_test),
    X = lapply(1:npred_test, function(i) {
          matrix(rnorm(nT2_test * p_test * ncoef_lin_Z), # Use ncoef here 
                 nrow = nT2_test * p_test, ncol = ncoef_lin_Z)
        })
  )
  
  # newdata for models WITHOUT Z (only 1 prediction unit expected)
  newdata_noZ <- list(
    p = p_test,
    X = list(matrix(rnorm(nT2_test * p_test * ncoef_lin_noZ), # Use ncoef here
                    nrow = nT2_test * p_test, ncol = ncoef_lin_noZ))
    # Z should be NULL or absent
  )
  
  # newdata for T=1 test
  newdata_Z_T1 <- list(
    p = p_test,
    Z = matrix(rnorm(npred_test * nz_test), nrow = npred_test, ncol = nz_test),
    X = lapply(1:npred_test, function(i) {
          matrix(rnorm(nT1_test * p_test * ncoef_lin_Z), # Use ncoef here
                 nrow = nT1_test * p_test, ncol = ncoef_lin_Z)
        })
  )
  
  # newdata for npred=1 test
  newdata_Z_npred1 <- list(
    p = p_test,
    Z = matrix(rnorm(1 * nz_test), nrow = 1, ncol = nz_test),
    X = list(matrix(rnorm(nT2_test * p_test * ncoef_lin_Z), # Use ncoef here
                    nrow = nT2_test * p_test, ncol = ncoef_lin_Z))
  )

  # --- Core Functionality Tests --- 
  
  # Test 2.1: Basic Execution (Linear, with Z)
  pred_lin_Z <- predict(fit_lin_Z, newdata = newdata_Z, type = "prior_probs", nsim = 1, r_verbose=FALSE)
  expect_type(pred_lin_Z, "list")
  expect_length(pred_lin_Z, npred_test)
  expect_true(is.array(pred_lin_Z[[1]]), info = "pred_lin_Z[[1]] is not array")
  expect_equal(dim(pred_lin_Z[[1]]), c(nT2_test, p_test, ndraws_out))
  expect_equal(dim(pred_lin_Z[[npred_test]]), c(nT2_test, p_test, ndraws_out))
  
  # Test 2.2: Basic Execution (BART, with Z)
  pred_bart_Z <- predict(fit_bart_Z, newdata = newdata_Z, type = "prior_probs", nsim = 1, r_verbose=FALSE)
  expect_type(pred_bart_Z, "list")
  expect_length(pred_bart_Z, npred_test)
  expect_true(is.array(pred_bart_Z[[1]]), info = "pred_bart_Z[[1]] is not array")
  expect_equal(dim(pred_bart_Z[[1]]), c(nT2_test, p_test, ndraws_out))
  
  # Test 2.3: Basic Execution (Linear, without Z)
  pred_lin_noZ <- predict(fit_lin_noZ, newdata = newdata_noZ, type = "prior_probs", nsim = 1, r_verbose=FALSE)
  expect_type(pred_lin_noZ, "list")
  expect_length(pred_lin_noZ, 1)
  expect_true(is.array(pred_lin_noZ[[1]]), info = "pred_lin_noZ[[1]] is not array")
  expect_equal(dim(pred_lin_noZ[[1]]), c(nT2_test, p_test, ndraws_out))
  
  # Test 2.4: Basic Execution (BART, without Z)
  pred_bart_noZ <- predict(fit_bart_noZ, newdata = newdata_noZ, type = "prior_probs", nsim = 1, r_verbose=FALSE)
  expect_type(pred_bart_noZ, "list")
  expect_length(pred_bart_noZ, 1)
  expect_true(is.array(pred_bart_noZ[[1]]), info = "pred_bart_noZ[[1]] is not array")
  expect_equal(dim(pred_bart_noZ[[1]]), c(nT2_test, p_test, ndraws_out))

  # --- nsim Argument Tests --- 
  
  # Test 3.1: nsim > 1 produces different results
  set.seed(1104)
  pred_nsim1 <- predict(fit_lin_Z, newdata = newdata_Z, type = "prior_probs", nsim = 1, r_verbose=FALSE)
  set.seed(1104) # Reset seed
  pred_nsim5 <- predict(fit_lin_Z, newdata = newdata_Z, type = "prior_probs", nsim = 5, r_verbose=FALSE)
  
  expect_type(pred_nsim5, "list")
  expect_length(pred_nsim5, npred_test)
  expect_equal(dim(pred_nsim1[[1]]), dim(pred_nsim5[[1]]), info = "Dimensions differ for nsim=1 vs nsim=5")
  # Results should differ due to averaging different eta draws
  expect_false(identical(pred_nsim1[[1]][1,1,1], pred_nsim5[[1]][1,1,1]), 
               info = "Results unexpectedly identical for nsim=1 vs nsim=5")

  # --- Edge Case Tests --- 

  # Test 5.1: npred = 1 (with Z)
  pred_npred1 <- predict(fit_lin_Z, newdata = newdata_Z_npred1, type = "prior_probs", nsim = 1, r_verbose=FALSE)
  expect_type(pred_npred1, "list")
  expect_length(pred_npred1, 1)
  expect_true(is.array(pred_npred1[[1]]), info = "pred_npred1[[1]] is not array")
  expect_equal(dim(pred_npred1[[1]]), c(nT2_test, p_test, ndraws_out))

  # Test 5.2: T_i = 1 (single observation per unit)
  pred_T1 <- predict(fit_lin_Z, newdata = newdata_Z_T1, type = "prior_probs", nsim = 1, r_verbose=FALSE)
  expect_type(pred_T1, "list")
  expect_length(pred_T1, npred_test)
  expect_true(is.array(pred_T1[[1]]), info = "pred_T1[[1]] is not array")
  expect_equal(dim(pred_T1[[1]]), c(nT1_test, p_test, ndraws_out))
  # Check sum-to-one (already a row vector, so sum directly)
  expect_equal(sum(pred_T1[[1]][1,,1]), 1.0, tolerance = 1e-6, 
               info = "T=1: Probabilities don't sum to 1")

}) 
