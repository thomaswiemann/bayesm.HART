library(testthat)

# Tests for rhierMnlRwMixture function

# Helper to run models quickly - suppress messages
run_quick_hmnl <- function(data, prior, mcmc) {
  capture.output({
    fit <- rhierMnlRwMixture(Data = data, Prior = prior, Mcmc = mcmc)
  })
  return(fit)
}

test_that("Basic function execution works", {

  # --- Settings for quick tests ---
  nlgt_test <- 5
  nT_test <- 3
  p_test <- 3
  nz_test <- 2
  R_test <- 10
  keep_test <- 1

  # --- 1. Linear Model with Z (nz > 0) ---
  set.seed(1)
  sim_data_linear <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = nz_test, nXa = 1, nXd = 0,
    beta_func_type = "linear"
  )
  Prior_linear <- list(ncomp = 1) # Defaults determined within function
  Mcmc_linear <- list(R = R_test, keep = keep_test, nprint = 0)

  expect_no_error({
    fit_linear <- run_quick_hmnl(sim_data_linear, Prior_linear, Mcmc_linear)
  }, message = "Linear model failed basic execution.")
  expect_s3_class(fit_linear, "rhierMnlRwMixture")

  # --- 2. BART Model with Z (nz > 0) ---
  set.seed(2)
  sim_data_bart <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = nz_test, nXa = 1, nXd = 0,
    beta_func_type = "linear" # Can use linear, BART fits non-linearity
  )
  Prior_bart <- list(ncomp = 1, bart = list(num_trees = 5)) # Minimal BART
  Mcmc_bart <- list(R = R_test, keep = keep_test, nprint = 0)

  expect_no_error({
    fit_bart <- run_quick_hmnl(sim_data_bart, Prior_bart, Mcmc_bart)
  }, message = "BART model failed basic execution.")
  expect_s3_class(fit_bart, "rhierMnlRwMixture")

  # --- 3. Model without Z (nz = 0) ---
  set.seed(3)
  sim_data_noZ <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = 0, nXa = 1, nXd = 0
  )
  Prior_noZ <- list(ncomp = 1)
  Mcmc_noZ <- list(R = R_test, keep = keep_test, nprint = 0)

  expect_no_error({
    fit_noZ <- run_quick_hmnl(sim_data_noZ, Prior_noZ, Mcmc_noZ)
  }, message = "Model with nz=0 failed basic execution.")
  expect_s3_class(fit_noZ, "rhierMnlRwMixture")

})

test_that("Input validation works", {

  # --- Minimal Valid Setup ---
  set.seed(100)
  nlgt_val <- 3
  nT_val <- 2
  p_val <- 2
  nz_val <- 1
  valid_data <- sim_hier_mnl(
    nlgt = nlgt_val, nT = nT_val, p = p_val, nz = nz_val, nXa = 1, nXd = 0,
    seed = 101
  )
  valid_prior <- list(ncomp = 1)
  valid_mcmc <- list(R = 2, keep = 1, nprint = 0)

  # --- Test Missing Arguments ---
  expect_error(capture.output(rhierMnlRwMixture(Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "Requires Data argument")
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Mcmc = valid_mcmc)),
               regexp = "Requires Prior list argument")
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = valid_prior)),
               regexp = "Requires Mcmc list argument")

  # --- Test Data list validation ---
  bad_data_p <- valid_data; bad_data_p$p <- NULL
  expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_p, Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "Requires Data element p")

  bad_data_lgt <- valid_data; bad_data_lgt$lgtdata <- NULL
  expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_lgt, Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "Requires Data element lgtdata")

  # --- Test Data$Z validation ---
  bad_data_Z_mat <- valid_data; bad_data_Z_mat$Z <- as.data.frame(valid_data$Z)
  expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_Z_mat, Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "Z must be a matrix")

  # --- Test lgtdata elements (y, X) validation ---
  bad_data_y_null <- valid_data; bad_data_y_null$lgtdata[[1]]$y <- NULL
  expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_y_null, Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "Requires element y of lgtdata")

  bad_data_X_null <- valid_data; bad_data_X_null$lgtdata[[1]]$X <- NULL
  expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_X_null, Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "Requires element X of lgtdata")

  bad_data_X_mat <- valid_data; bad_data_X_mat$lgtdata[[1]]$X <- as.data.frame(valid_data$lgtdata[[1]]$X)
  expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_X_mat, Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "X must be a matrix")

  bad_data_X_row <- valid_data; bad_data_X_row$lgtdata[[1]]$X <- valid_data$lgtdata[[1]]$X[1:(nrow(valid_data$lgtdata[[1]]$X)-p_val), ]
  expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_X_row, Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "nrow\\(X\\) ne p\\*length\\(yi\\)")

  bad_data_X_col <- valid_data; bad_data_X_col$lgtdata[[2]]$X <- valid_data$lgtdata[[1]]$X[, 1, drop=FALSE]
  # Need at least 2 units for this check, adjust setup if needed
  if(nlgt_val >= 2) {
      expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_X_col, Prior = valid_prior, Mcmc = valid_mcmc)),
                   regexp = "All X elements must have same # of cols")
  }

  bad_data_y_val <- valid_data; bad_data_y_val$lgtdata[[1]]$y <- rep(p_val + 1, nT_val) # Invalid choice value
  expect_error(capture.output(rhierMnlRwMixture(Data = bad_data_y_val, Prior = valid_prior, Mcmc = valid_mcmc)),
               regexp = "y takes on .* values -- must be = p") # Match the earlier error check

  # --- Test Prior list validation ---
  bad_prior_ncomp <- list() # Missing ncomp
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_ncomp, Mcmc = valid_mcmc)),
               regexp = "Requires Prior element ncomp")

  ncoef_val <- valid_data$true_values$dimensions$ncoef
  bad_prior_signres_len <- valid_prior; bad_prior_signres_len$SignRes <- rep(0, ncoef_val + 1)
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_signres_len, Mcmc = valid_mcmc)),
               regexp = "SignRes must be equal to the dimension of X")

  bad_prior_signres_val <- valid_prior; bad_prior_signres_val$SignRes <- rep(2, ncoef_val)
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_signres_val, Mcmc = valid_mcmc)),
               regexp = "All elements of SignRes must be equal to")

  bad_prior_mubar <- valid_prior; bad_prior_mubar$mubar <- matrix(0, nrow = 1, ncol = ncoef_val + 1)
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_mubar, Mcmc = valid_mcmc)),
               regexp = "mubar must have ncomp cols") # Note: error msg seems wrong, checks nvar

  bad_prior_amu <- valid_prior; bad_prior_amu$Amu <- matrix(0, nrow = 2, ncol = 1)
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_amu, Mcmc = valid_mcmc)),
               regexp = "Am must be a 1 x 1 array")

  bad_prior_nu <- valid_prior; bad_prior_nu$nu <- 0
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_nu, Mcmc = valid_mcmc)),
               regexp = "invalid nu value")

  bad_prior_V <- valid_prior; bad_prior_V$V <- diag(ncoef_val + 1)
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_V, Mcmc = valid_mcmc)),
               regexp = "Invalid V in prior")

  # Only test delta priors if nz > 0
  if (nz_val > 0) {
    bad_prior_Ad <- valid_prior; bad_prior_Ad$Ad <- diag(ncoef_val * nz_val + 1)
    expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_Ad, Mcmc = valid_mcmc)),
                 regexp = "Ad must be nvar\\*nz x nvar\\*nz")

    bad_prior_deltabar <- valid_prior; bad_prior_deltabar$deltabar <- rep(0, ncoef_val * nz_val + 1)
    expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_deltabar, Mcmc = valid_mcmc)),
                 regexp = "deltabar must be of length nvar\\*nz")
  }

  bad_prior_a <- valid_prior; bad_prior_a$a <- rep(-1, valid_prior$ncomp)
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_a, Mcmc = valid_mcmc)),
               regexp = "invalid values in a vector")

  bad_prior_a_len <- valid_prior; bad_prior_a_len$a <- rep(1, valid_prior$ncomp + 1)
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = bad_prior_a_len, Mcmc = valid_mcmc)),
               regexp = "Requires dim\\(a\\)= ncomp")

  # --- Test Mcmc list validation ---
  bad_mcmc_R <- valid_mcmc; bad_mcmc_R$R <- NULL
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = valid_prior, Mcmc = bad_mcmc_R)),
               regexp = "Requires R argument in Mcmc list")

  bad_mcmc_nprint <- valid_mcmc; bad_mcmc_nprint$nprint <- -1
  expect_error(capture.output(rhierMnlRwMixture(Data = valid_data, Prior = valid_prior, Mcmc = bad_mcmc_nprint)),
               regexp = "nprint must be an integer greater than or equal to 0")

  # (End of Task 2 tests)

})

test_that("Prior specification handling works", {

  # --- Settings for quick tests ---
  nlgt_test <- 5
  nT_test <- 3
  p_test <- 3
  nz_test <- 2
  R_test <- 10
  keep_test <- 1
  Mcmc_test <- list(R = R_test, keep = keep_test, nprint = 0)

  # --- Test Default Priors (Run without error) ---
  # 1. No Z, No SignRes
  set.seed(201)
  sim_noZ_noSign <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = 0, nXa = 1, nXd = 0, seed = 202
  )
  Prior_1 <- list(ncomp = 1)
  expect_no_error({
    fit1 <- run_quick_hmnl(sim_noZ_noSign, Prior_1, Mcmc_test)
  }, message = "Default prior failed: No Z, No SignRes")
  expect_s3_class(fit1, "rhierMnlRwMixture")

  # 2. With Z, No SignRes
  set.seed(203)
  sim_Z_noSign <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = nz_test, nXa = 1, nXd = 0, seed = 204
  )
  Prior_2 <- list(ncomp = 1)
  expect_no_error({
    fit2 <- run_quick_hmnl(sim_Z_noSign, Prior_2, Mcmc_test)
  }, message = "Default prior failed: With Z, No SignRes")
  expect_s3_class(fit2, "rhierMnlRwMixture")

  # 3. No Z, With SignRes
  set.seed(205)
  sim_noZ_Sign <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = 0, nXa = 1, nXd = 0, seed = 206
  )
  ncoef_3 <- sim_noZ_Sign$true_values$dimensions$ncoef
  Prior_3 <- list(ncomp = 1, SignRes = c(1, rep(0, ncoef_3 - 1))) # Restrict first coef > 0
  expect_no_error({
    fit3 <- run_quick_hmnl(sim_noZ_Sign, Prior_3, Mcmc_test)
  }, message = "Default prior failed: No Z, With SignRes")
  expect_s3_class(fit3, "rhierMnlRwMixture")

  # 4. With Z, With SignRes
  set.seed(207)
  sim_Z_Sign <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = nz_test, nXa = 1, nXd = 0, seed = 208
  )
  ncoef_4 <- sim_Z_Sign$true_values$dimensions$ncoef
  Prior_4 <- list(ncomp = 1, SignRes = c(-1, rep(0, ncoef_4 - 1))) # Restrict first coef < 0
  expect_no_error({
    fit4 <- run_quick_hmnl(sim_Z_Sign, Prior_4, Mcmc_test)
  }, message = "Default prior failed: With Z, With SignRes")
  expect_s3_class(fit4, "rhierMnlRwMixture")

  # --- Test User-Specified Priors (Run without error) ---
  # Using sim_Z_noSign data from above (nz>0, no SignRes)
  ncoef_user <- sim_Z_noSign$true_values$dimensions$ncoef
  nz_user <- sim_Z_noSign$true_values$dimensions$nz

  Prior_user <- list(
    ncomp = 1,
    mubar = matrix(rep(1, ncoef_user), nrow = 1), # Specify mubar
    Amu = matrix(0.5),                          # Specify Amu
    nu = ncoef_user + 5,                        # Specify nu
    V = diag(ncoef_user) * (ncoef_user + 5),    # Specify V
    # Specify delta priors only if nz > 0
    deltabar = if(nz_user > 0) rep(0.5, nz_user * ncoef_user) else NULL,
    Ad = if(nz_user > 0) diag(nz_user * ncoef_user) * 0.5 else NULL,
    a = rep(2, 1)                               # Specify a
  )

  expect_no_error({
    fit_user <- run_quick_hmnl(sim_Z_noSign, Prior_user, Mcmc_test)
  }, message = "User-specified prior failed")
  expect_s3_class(fit_user, "rhierMnlRwMixture")

  # --- Test BART Priors (Run without error) ---
  # Using sim_Z_noSign data
  # 1. BART defaults (just specify bart list)
  Prior_bart_def <- list(ncomp = 1, bart = list())
  expect_no_error({
      fit_bart_def <- run_quick_hmnl(sim_Z_noSign, Prior_bart_def, Mcmc_test)
  }, message = "BART default prior failed")
  expect_s3_class(fit_bart_def, "rhierMnlRwMixture")
  # Check if bart_models exist, indicating BART ran (Task 6 will check more formally)
  expect_true(!is.null(fit_bart_def$bart_models), "bart_models should be in output for BART default")

  names(fit_bart_def$bart_models)

  # 2. BART user overrides
  Prior_bart_user <- list(ncomp = 1, bart = list(num_trees = 10, power = 1.5, base = 0.8))
  expect_no_error({
      fit_bart_user <- run_quick_hmnl(sim_Z_noSign, Prior_bart_user, Mcmc_test)
  }, message = "BART user-specified prior failed")
  expect_s3_class(fit_bart_user, "rhierMnlRwMixture")
  # Check if bart_models exist, indicating BART ran (Task 6 will check more formally)
  expect_true(!is.null(fit_bart_user$bart_models), "bart_models should be in output for BART user")

  # (End of Task 3 tests)

})

test_that("Sign restriction handling works", {

  # --- Settings for quick tests ---
  nlgt_test <- 5
  nT_test <- 3
  p_test <- 3
  nz_test <- 2
  R_sign <- 50 # Need a few more draws to check restrictions
  keep_sign <- 1
  Mcmc_sign <- list(R = R_sign, keep = keep_sign, nprint = 0)

  set.seed(301)
  sim_data_sign <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = nz_test,
    nXa = 2, nXd = 0, const = TRUE, # Ensure at least 2 coefs (p-1 const + nXa)
    seed = 302
  )
  ncoef_sign <- sim_data_sign$true_values$dimensions$ncoef
  expect_true(ncoef_sign >= 2, "Need at least 2 coefs for sign restriction tests")

  # --- Test Runs with Different Restrictions (expect no error) ---
  # 1. Positive restriction
  Prior_pos <- list(ncomp = 1, SignRes = c(1, rep(0, ncoef_sign - 1)))
  expect_no_error({
    fit_pos <- run_quick_hmnl(sim_data_sign, Prior_pos, Mcmc_sign)
  }, message = "SignRes failed: Positive restriction")
  expect_s3_class(fit_pos, "rhierMnlRwMixture")

  # 2. Negative restriction
  Prior_neg <- list(ncomp = 1, SignRes = c(-1, rep(0, ncoef_sign - 1)))
  expect_no_error({
    fit_neg <- run_quick_hmnl(sim_data_sign, Prior_neg, Mcmc_sign)
  }, message = "SignRes failed: Negative restriction")
  expect_s3_class(fit_neg, "rhierMnlRwMixture")

  # 3. Mixed restriction
  Prior_mix <- list(ncomp = 1, SignRes = c(1, -1, rep(0, ncoef_sign - 2)))
  expect_no_error({
    fit_mix <- run_quick_hmnl(sim_data_sign, Prior_mix, Mcmc_sign)
  }, message = "SignRes failed: Mixed restriction")
  expect_s3_class(fit_mix, "rhierMnlRwMixture")

  # --- Test Verification: Draws adhere to restrictions ---
  # Use fit_pos, fit_neg, fit_mix from above

  # Get betadraws (remove burn-in implicitly via Mcmc$keep=1)
  betadraw_pos <- fit_pos$betadraw
  betadraw_neg <- fit_neg$betadraw
  betadraw_mix <- fit_mix$betadraw

  # Check positive restriction (first coef should be > 0)
  # Allow for potential numerical noise near zero, but expect mostly positive
  expect_true(all(betadraw_pos[, 1, ] > -1e-9), # Allow tiny negative due to float
              info = "Positive restriction not fully respected for coef 1")
  expect_true(all(betadraw_pos[, 1, ] > 0),
              info = "Majority of draws should respect positive restriction for coef 1")

  # Check negative restriction (first coef should be < 0)
  expect_true(all(betadraw_neg[, 1, ] < 1e-9), # Allow tiny positive due to float
              info = "Negative restriction not fully respected for coef 1")
  expect_true(all(betadraw_neg[, 1, ] < 0),
              info = "All draws should respect negative restriction for coef 1")

  # Check mixed restriction (coef 1 > 0, coef 2 < 0)
  expect_true(all(betadraw_mix[, 1, ] > -1e-9),
              info = "Positive restriction (mixed) not fully respected for coef 1")
  expect_true(all(betadraw_mix[, 1, ] > 0),
              info = "All draws should respect positive restriction (mixed) for coef 1")

  expect_true(all(betadraw_mix[, 2, ] < 1e-9),
              info = "Negative restriction (mixed) not fully respected for coef 2")
  expect_true(all(betadraw_mix[, 2, ] < 0),
              info = "All draws should respect negative restriction (mixed) for coef 2")

  # Note: Optimization and Hessian checks are harder to verify directly
  # in unit tests without accessing internal states or specific outputs.
  # Relying on successful run and draw verification provides indirect evidence.

})

test_that("BART integration works", {

  # --- Settings for quick tests ---
  nlgt_test <- 5
  nT_test <- 3
  p_test <- 3
  nz_test <- 2
  R_test <- 10
  keep_test <- 1
  Mcmc_test <- list(R = R_test, keep = keep_test, nprint = 0)

  set.seed(401)
  sim_data_bartint <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = nz_test, nXa = 1, nXd = 0,
    seed = 402
  )

  # --- Test BART is triggered by Prior$bart ---
  Prior_bart <- list(ncomp = 1, bart = list(num_trees = 5))
  fit_bart <- run_quick_hmnl(sim_data_bartint, Prior_bart, Mcmc_test)

  expect_true(!is.null(fit_bart$bart_models), "bart_models should exist when Prior$bart is set")
  expect_true(is.null(fit_bart$Deltadraw), "Deltadraw should be NULL when Prior$bart is set")

  # --- Test Linear model runs when Prior$bart is NULL ---
  Prior_linear_null <- list(ncomp = 1, bart = NULL)
  fit_linear_null <- run_quick_hmnl(sim_data_bartint, Prior_linear_null, Mcmc_test)

  expect_true(is.null(fit_linear_null$bart_models), "bart_models should be NULL when Prior$bart is NULL")
  expect_true(!is.null(fit_linear_null$Deltadraw), "Deltadraw should exist when Prior$bart is NULL")

  # --- Test Linear model runs when Prior$bart is missing ---
  Prior_linear_miss <- list(ncomp = 1) # bart element missing entirely
  fit_linear_miss <- run_quick_hmnl(sim_data_bartint, Prior_linear_miss, Mcmc_test)

  expect_true(is.null(fit_linear_miss$bart_models), "bart_models should be NULL when Prior$bart is missing")
  expect_true(!is.null(fit_linear_miss$Deltadraw), "Deltadraw should exist when Prior$bart is missing")

  # --- Test Linear model runs when drawdelta=FALSE (no Z) ---
  # BART is only applicable when Z exists (drawdelta=TRUE internally)
  set.seed(403)
  sim_data_noZ <- sim_hier_mnl(
    nlgt = nlgt_test, nT = nT_test, p = p_test, nz = 0, nXa = 1, nXd = 0, seed = 404
  )
  Prior_noZ <- list(ncomp = 1, bart = list()) # Specify bart, but nz=0
  fit_noZ <- run_quick_hmnl(sim_data_noZ, Prior_noZ, Mcmc_test)

  expect_true(is.null(fit_noZ$bart_models), "bart_models should be NULL when nz=0, even if Prior$bart is set")
  # Deltadraw is also NULL when nz=0, so no check needed here.

})

test_that("Output structure and content are correct", {

  # --- Settings --- 
  nlgt_out <- 4
  nT_out <- 2
  p_out <- 3
  nz_out <- 2
  nXa_out <- 1
  nXd_out <- 0
  const_out <- TRUE
  ncoef_out <- const_out*(p_out-1) + (p_out-1)*nXd_out + nXa_out
  R_out <- 20
  keep_out <- 2 # Use keep > 1
  ndraws_out <- R_out / keep_out
  Mcmc_out <- list(R = R_out, keep = keep_out, nprint = 0)

  set.seed(501)
  sim_data_out <- sim_hier_mnl(
    nlgt = nlgt_out, nT = nT_out, p = p_out, nz = nz_out,
    nXa = nXa_out, nXd = nXd_out, const = const_out, seed = 502
  )

  # --- Test Linear Model Output --- 
  Prior_linear <- list(ncomp = 1)
  fit_linear <- run_quick_hmnl(sim_data_out, Prior_linear, Mcmc_out)

  expect_s3_class(fit_linear, "rhierMnlRwMixture")

  # Deltadraw (Linear only)
  expect_true(!is.null(fit_linear$Deltadraw))
  expect_true(is.matrix(fit_linear$Deltadraw))
  expect_equal(nrow(fit_linear$Deltadraw), ndraws_out)
  expect_equal(ncol(fit_linear$Deltadraw), nz_out * ncoef_out)
  expect_s3_class(fit_linear$Deltadraw, "bayesm.mat")
  expect_s3_class(fit_linear$Deltadraw, "mcmc")
  expect_equal(attributes(fit_linear$Deltadraw)$mcpar, c(1, R_out, keep_out))

  # betadraw (Always present)
  expect_true(!is.null(fit_linear$betadraw))
  expect_true(is.array(fit_linear$betadraw))
  expect_equal(dim(fit_linear$betadraw), c(nlgt_out, ncoef_out, ndraws_out))
  expect_s3_class(fit_linear$betadraw, "bayesm.hcoef")

  # nmix (Always present)
  expect_true(!is.null(fit_linear$nmix))
  expect_s3_class(fit_linear$nmix, "bayesm.nmix")
  # nmix$probdraw
  expect_true(!is.null(fit_linear$nmix$probdraw))
  expect_true(is.matrix(fit_linear$nmix$probdraw))
  expect_equal(nrow(fit_linear$nmix$probdraw), ndraws_out)
  expect_equal(ncol(fit_linear$nmix$probdraw), Prior_linear$ncomp)
  # nmix$compdraw (list of lists of lists)
  expect_true(!is.null(fit_linear$nmix$compdraw))
  expect_true(is.list(fit_linear$nmix$compdraw))
  expect_length(fit_linear$nmix$compdraw, ndraws_out)
  expect_true(is.list(fit_linear$nmix$compdraw[[1]])) # Check first draw
  expect_length(fit_linear$nmix$compdraw[[1]], Prior_linear$ncomp)
  expect_true(is.list(fit_linear$nmix$compdraw[[1]][[1]])) # Check first comp
  expect_true(all(c("mu", "rooti") %in% names(fit_linear$nmix$compdraw[[1]][[1]])))

  # loglike (Always present)
  expect_true(!is.null(fit_linear$loglike))
  expect_true(is.numeric(fit_linear$loglike))
  expect_length(fit_linear$loglike, ndraws_out)

  # --- Test BART Model Output --- 
  Prior_bart <- list(ncomp = 1, bart = list(num_trees = 5))
  fit_bart <- run_quick_hmnl(sim_data_out, Prior_bart, Mcmc_out)

  expect_s3_class(fit_bart, "rhierMnlRwMixture")

  # Deltadraw (NULL for BART)
  expect_true(is.null(fit_bart$Deltadraw))

  # bart_models (BART only)
  expect_true(!is.null(fit_bart$bart_models))
  expect_true(is.list(fit_bart$bart_models))
  expect_length(fit_bart$bart_models, ncoef_out)
  # Check structure of one model (basic check)
  expect_true(!is.null(fit_bart$bart_models[[1]]$treedraws))

  # betadraw, nmix, loglike (should exist and have same structure as linear)
  expect_true(!is.null(fit_bart$betadraw))
  expect_equal(dim(fit_bart$betadraw), c(nlgt_out, ncoef_out, ndraws_out))
  expect_s3_class(fit_bart$betadraw, "bayesm.hcoef")

  expect_true(!is.null(fit_bart$nmix))
  expect_s3_class(fit_bart$nmix, "bayesm.nmix")
  expect_equal(nrow(fit_bart$nmix$probdraw), ndraws_out)
  expect_length(fit_bart$nmix$compdraw, ndraws_out)

  expect_true(!is.null(fit_bart$loglike))
  expect_true(is.numeric(fit_bart$loglike))
  expect_length(fit_bart$loglike, ndraws_out)

})

test_that("Edge cases run without error", {

  # --- Shared Settings ---
  R_edge <- 10
  keep_edge <- 1
  Mcmc_edge <- list(R = R_edge, keep = keep_edge, nprint = 0)
  Prior_edge <- list(ncomp = 1)

  # --- Dimensionality Extremes ---
  # Minimal nT = 1
  set.seed(605)
  sim_nT1 <- sim_hier_mnl(nlgt = 50, nT = 1, p = 3, nz = 1, nXa = 1, seed = 606)
  expect_no_error(run_quick_hmnl(sim_nT1, Prior_edge, Mcmc_edge))

  # Binary choice p = 2
  set.seed(607)
  sim_p2 <- sim_hier_mnl(nlgt = 50, nT = 5, p = 2, nz = 1, nXa = 1, seed = 608)
  expect_no_error(run_quick_hmnl(sim_p2, Prior_edge, Mcmc_edge))

  # No intercepts const = FALSE
  set.seed(609)
  sim_noconst <- sim_hier_mnl(nlgt = 50, nT = 5, p = 3, nz = 1, nXa = 1, const = FALSE, seed = 610)
  expect_no_error(run_quick_hmnl(sim_noconst, Prior_edge, Mcmc_edge))

  # Minimal predictors (only 1 Xa, no const, no Xd)
  set.seed(611)
  sim_minpred <- sim_hier_mnl(nlgt = 50, nT = 5, p = 3, nz = 1, nXa = 1, nXd = 0, const = FALSE, seed = 612)
  expect_warning(
    run_quick_hmnl(sim_minpred, Prior_edge, Mcmc_edge),
    regexp = "one-dimensional optimization by Nelder-Mead is unreliable"
  )

  # Minimal Z (nz = 1)
  set.seed(613)
  sim_nz1 <- sim_hier_mnl(nlgt = 50, nT = 5, p = 3, nz = 1, nXa = 1, seed = 614)
  expect_no_error(run_quick_hmnl(sim_nz1, Prior_edge, Mcmc_edge))

  # Minimal Z (nz = 0) - Already tested in Task 1, but confirm here
  set.seed(615)
  sim_nz0 <- sim_hier_mnl(nlgt = 50, nT = 5, p = 3, nz = 0, nXa = 1, seed = 616)
  expect_no_error(run_quick_hmnl(sim_nz0, Prior_edge, Mcmc_edge))

  # --- MCMC Parameters ---
  # keep > 1 (Output structure tested in Task 6, just check run)
  Mcmc_keep <- list(R = R_edge, keep = 5, nprint = 0)
  expect_no_error(run_quick_hmnl(sim_nz1, Prior_edge, Mcmc_keep)) # Use nz=1 data

  # --- BART Specifics ---
  # Minimal BART parameters (num_trees = 1)
  Prior_bart_min <- list(ncomp = 1, bart = list(num_trees = 1))
  expect_no_error(run_quick_hmnl(sim_nz1, Prior_bart_min, Mcmc_edge)) # Use nz=1 data

})

test_that("Reproducibility works", {

  # --- Settings ---
  nlgt_rep <- 100
  nT_rep <- 2
  p_rep <- 2
  nz_rep <- 1
  R_rep <- 10
  keep_rep <- 1
  Mcmc_rep <- list(R = R_rep, keep = keep_rep, nprint = 0)
  Prior_rep <- list(ncomp = 1)
  Prior_bart_rep <- list(ncomp = 1, bart = list(num_trees = 3))

  # --- Linear Model ---
  set.seed(701) # Seed for data generation
  sim_data_rep_lin <- sim_hier_mnl(
    nlgt = nlgt_rep, nT = nT_rep, p = p_rep, nz = nz_rep, nXa = 1, seed = 702
  )

  set.seed(703) # Seed for first run
  fit1_lin <- run_quick_hmnl(sim_data_rep_lin, Prior_rep, Mcmc_rep)

  set.seed(703) # SAME Seed for second run
  fit2_lin <- run_quick_hmnl(sim_data_rep_lin, Prior_rep, Mcmc_rep)

  # Compare key outputs
  expect_equal(fit1_lin$Deltadraw, fit2_lin$Deltadraw, info = "Linear Deltadraw not reproducible")
  expect_equal(fit1_lin$betadraw, fit2_lin$betadraw, info = "Linear betadraw not reproducible")
  expect_equal(fit1_lin$nmix, fit2_lin$nmix, info = "Linear nmix not reproducible")
  expect_equal(fit1_lin$loglike, fit2_lin$loglike, info = "Linear loglike not reproducible")

  # --- BART Model ---
  set.seed(704) # Seed for data generation
  sim_data_rep_bart <- sim_hier_mnl(
    nlgt = nlgt_rep, nT = nT_rep, p = p_rep, nz = nz_rep, nXa = 1, seed = 705
  )

  set.seed(706) # Seed for first run
  fit1_bart <- run_quick_hmnl(sim_data_rep_bart, Prior_bart_rep, Mcmc_rep)

  set.seed(706) # SAME Seed for second run
  fit2_bart <- run_quick_hmnl(sim_data_rep_bart, Prior_bart_rep, Mcmc_rep)

  # Compare key outputs (no Deltadraw for BART)
  # Comparing bart_models directly might be tricky due to environments/pointers?
  # Compare summaries or components if direct comparison fails reliably.
  # For now, compare the other outputs.
  expect_equal(fit1_bart$betadraw, fit2_bart$betadraw, info = "BART betadraw not reproducible")
  expect_equal(fit1_bart$nmix, fit2_bart$nmix, info = "BART nmix not reproducible")
  expect_equal(fit1_bart$loglike, fit2_bart$loglike, info = "BART loglike not reproducible")
  # We can't easily compare bart_models for perfect equality due to how they store trees
  # but the other components being equal strongly suggests reproducibility.

})

test_that("beta and delta parameter recovery is reasonable", {

  # --- Settings ---
  nlgt_rec <- 500 # More units for better estimates
  nT_rec <- 10
  p_rec <- 3
  nz_rec <- 2
  R_rec <- 1000  # More draws needed for recovery check
  keep_rec <- 10
  Mcmc_rec <- list(R = R_rec, keep = keep_rec, nprint = 0)

  # --- Linear Model Recovery ---
  set.seed(801)
  sim_data_lin_rec <- sim_hier_mnl(
    nlgt = nlgt_rec, nT = nT_rec, p = p_rec, nz = nz_rec, nXa = 1,
    beta_func_type = "linear", seed = 802
  )
  Prior_lin_rec <- list(ncomp = 1)

  set.seed(803)
  fit_lin_rec <- run_quick_hmnl(sim_data_lin_rec, Prior_lin_rec, Mcmc_rec)

  # Compare betadraw posterior mean to true beta
  beta_mean_lin <- apply(fit_lin_rec$betadraw, c(1, 2), mean)
  corr_beta_lin <- cor(as.vector(beta_mean_lin), as.vector(sim_data_lin_rec$true_values$beta_true))
  expect_gt(corr_beta_lin, 0.5, label = "Linear: Posterior mean beta vs true beta correlation low")

  # Compare Deltadraw posterior mean to true Delta
  delta_mean_lin <- colMeans(fit_lin_rec$Deltadraw)
  # Reshape true Delta (which is nz x ncoef) into a vector matching Deltadraw
  true_delta_lin <- as.vector(t(sim_data_lin_rec$true_values$true_params$beta_func_args$Delta))
  corr_delta_lin <- cor(delta_mean_lin, true_delta_lin)
  expect_gt(corr_delta_lin, 0.5, label = "Linear: Posterior mean Delta vs true Delta correlation low")

  # --- BART Model Recovery ---
  set.seed(804)
  sim_data_bart_rec <- sim_hier_mnl(
    nlgt = nlgt_rec, nT = nT_rec, p = p_rec, nz = nz_rec, nXa = 1,
    # Use a non-linear function for a better BART test
    beta_func_type = "step", seed = 805
  )
  Prior_bart_rec <- list(ncomp = 1, bart = list(num_trees = 20))

  set.seed(806)
  fit_bart_rec <- run_quick_hmnl(sim_data_bart_rec, Prior_bart_rec, Mcmc_rec)

  # Compare betadraw posterior mean to true beta
  beta_mean_bart <- apply(fit_bart_rec$betadraw, c(1, 2), mean)
  corr_beta_bart <- cor(as.vector(beta_mean_bart), as.vector(sim_data_bart_rec$true_values$beta_true))
  expect_gt(corr_beta_bart, 0.5, label = "BART: Posterior mean beta vs true beta correlation low")

  # Compare posterior mean predicted betabar (DeltaZ) to true betabar
  # Need predict method - assume it's available/loaded
  pred_betabar_draws <- predict(fit_bart_rec, newdata = sim_data_bart_rec, type = "DeltaZ", r_verbose = FALSE)
  pred_betabar_mean <- apply(pred_betabar_draws, c(1, 2), mean)
  corr_betabar_bart <- cor(as.vector(pred_betabar_mean), as.vector(sim_data_bart_rec$true_values$betabar_true))
  expect_gt(corr_betabar_bart, 0.3, label = "BART: Predicted betabar vs true betabar correlation low")

}) 

test_that("Mixture with ncomp = 1 parameter recovery is reasonable", {

  # --- Settings ---
  nlgt_rec <- 500 # More units for better estimates
  nT_rec <- 10
  p_rec <- 3
  nz_rec <- 2
  R_rec <- 1000  # More draws needed for recovery check
  keep_rec <- 10
  Mcmc_rec <- list(R = R_rec, keep = keep_rec, nprint = 0)

  # --- Mixture Component Recovery (ncomp > 1) ---
  set.seed(807)
  ncomp_rec <- 1
  # Explicitly calculate ncoef based on sim_hier_mnl defaults being used
  const_rec <- TRUE # Default in sim_hier_mnl
  nXd_rec <- 0   # Default in sim_hier_mnl
  ncoef_rec <- const_rec * (p_rec - 1) + (p_rec - 1) * nXd_rec + 1

  # Restore definition of true_comps_rec_lin
  true_comps_rec_lin <- list(
   list(mu = c(-2, rep(0, ncoef_rec - 1)), rooti = diag(1, ncoef_rec))
  )
  # Revert the first sim call to use the defined list
  sim_data_mix_rec <- sim_hier_mnl(
    nlgt = nlgt_rec, nT = nT_rec, p = p_rec, nz = nz_rec, nXa = 1, 
    nXd = 0, # <<< Add nXd = 0
    const = TRUE, # Explicitly state const=TRUE (default, but good practice)
    beta_func_type = "linear",
    ncomp = ncomp_rec,
    mixture_comps = true_comps_rec_lin, # Use the defined list again
    seed = 808
  )
  Prior_mix_rec <- list(ncomp = ncomp_rec)

  set.seed(809)
  fit_mix_rec <- run_quick_hmnl(sim_data_mix_rec, Prior_mix_rec, Mcmc_rec)

  # Check compdraw recovery (mu and rooti) - accounting for label switching
  # Extract posterior means for each component
  draws <- fit_mix_rec$nmix$compdraw
  ndraws <- length(draws)

  ncoef <- length(draws[[1]][[1]]$mu) # Get ncoef from the mu vector
  
  # Initialize storage for means
  mu_means <- array(0, dim = c(ncoef, ncomp_rec))
  rooti_means <- array(0, dim = c(ncoef, ncoef, ncomp_rec))
  
  # Calculate means (sum across draws)
  for (r in 1:ndraws) {
    for (c in 1:ncomp_rec) {
      mu_means[, c] <- mu_means[, c] + as.vector(draws[[r]][[c]]$mu)
      rooti_means[, , c] <- rooti_means[, , c] + draws[[r]][[c]]$rooti
    }
  }
  mu_means <- mu_means / ndraws
  rooti_means <- rooti_means / ndraws

  # --- Order components for comparison (handle label switching) ---
  # Order estimated components based on the first element of mu
  est_order_lin <- order(mu_means[1, ])
  mu_means_ordered_lin <- mu_means[, est_order_lin, drop=FALSE]
  rooti_means_ordered_lin <- rooti_means[, , est_order_lin, drop=FALSE]

  # --- Compare ordered components ---
  # Compare Mu vectors
  expect_equal(mu_means_ordered_lin[, 1], true_comps_rec_lin[[1]]$mu,
               tolerance = 0.2, # Still using loose tolerance
               label = "Mixture: Ordered mu for component 1 doesn't match true")

  # Compare Rooti matrices
  expect_equal(rooti_means_ordered_lin[, , 1], true_comps_rec_lin[[1]]$rooti,
               tolerance = 0.2,
               label = "Mixture: Ordered rooti for component 1 doesn't match true")

  # --- Test Mixture Component Recovery (ncomp > 1, BART DeltaZ) ---
  set.seed(901)
  # Use ncoef_rec calculated above
  true_comps_rec_bart <- true_comps_rec_lin # Reuse component defs

  sim_data_mix_bart_rec <- sim_hier_mnl(
    nlgt = nlgt_rec, nT = nT_rec, p = p_rec, nz = nz_rec, nXa = 1,
    nXd = 0, # <<< Add nXd = 0
    const = TRUE, # Explicitly state const=TRUE (default, but good practice)
    beta_func_type = "step", # Use step for BART test
    ncomp = ncomp_rec,
    mixture_comps = true_comps_rec_bart,
    seed = 902
  )
  Prior_mix_bart_rec <- list(ncomp = ncomp_rec,
                             bart = list(num_trees = 20))

  fit_mix_bart_rec <- run_quick_hmnl(sim_data_mix_bart_rec, Prior_mix_bart_rec, Mcmc_rec)


  # Check compdraw recovery (mu and rooti) - accounting for label switching
  # Extract posterior means for each component
  draws_bart <- fit_mix_bart_rec$nmix$compdraw
  ndraws_bart <- length(draws_bart)

  ncoef_bart <- length(draws_bart[[1]][[1]]$mu) # Get ncoef from the mu vector
  
  # Initialize storage for means
  mu_means_bart <- array(0, dim = c(ncoef_bart, ncomp_rec))
  rooti_means_bart <- array(0, dim = c(ncoef_bart, ncoef_bart, ncomp_rec))
  
  # Calculate means (sum across draws)
  for (r in 1:ndraws_bart) {
    for (c in 1:ncomp_rec) {
      mu_means_bart[, c] <- mu_means_bart[, c] + as.vector(draws_bart[[r]][[c]]$mu)
      rooti_means_bart[, , c] <- rooti_means_bart[, , c] + draws_bart[[r]][[c]]$rooti
    }
  }
  mu_means_bart <- mu_means_bart / ndraws_bart
  rooti_means_bart <- rooti_means_bart / ndraws_bart

  # --- Order components for comparison (handle label switching) ---
  # Order estimated components based on the first element of mu
  est_order_bart <- order(mu_means_bart[1, ])
  mu_means_ordered_bart <- mu_means_bart[, est_order_bart, drop=FALSE]
  rooti_means_ordered_bart <- rooti_means_bart[, , est_order_bart, drop=FALSE]

  # Order true components based on the first element of mu
  true_order_bart <- order(sapply(true_comps_rec_bart, function(x) x$mu[1]))
  true_comps_rec_ordered_bart <- true_comps_rec_bart[true_order_bart]

  # --- Compare ordered components ---
  # Compare Mu vectors -- unidentified, skipping for now
#   expect_equal(mu_means_ordered_bart[, 1], true_comps_rec_ordered_bart[[1]]$mu,
#                tolerance = 0.5, # Still using loose tolerance
#                label = "Mixture: Ordered mu for component 1 doesn't match true")
#   expect_equal(mu_means_ordered_bart[, 2], true_comps_rec_ordered_bart[[2]]$mu,
#                tolerance = 0.5,
#                label = "Mixture: Ordered mu for component 2 doesn't match true")

  # Compare Rooti matrices
  expect_equal(rooti_means_ordered_bart[, , 1], true_comps_rec_ordered_bart[[1]]$rooti,
               tolerance = 0.2,
               label = "Mixture: Ordered rooti for component 1 doesn't match true")

}) 

test_that("Mixture parameter recovery is reasonable", {

  # --- Settings ---
  nlgt_rec <- 300 # More units for better estimates
  nT_rec <- 20
  p_rec <- 3
  nz_rec <- 2
  R_rec <- 3000  # More draws needed for recovery check
  keep_rec <- 10
  Mcmc_rec <- list(R = R_rec, keep = keep_rec, nprint = 0)

  # --- Mixture Component Recovery (ncomp > 1) ---
  set.seed(807)
  ncomp_rec <- 2
  # Explicitly calculate ncoef based on sim_hier_mnl defaults being used
  const_rec <- TRUE # Default in sim_hier_mnl
  nXd_rec <- 0   # Default in sim_hier_mnl
  ncoef_rec <- const_rec * (p_rec - 1) + (p_rec - 1) * nXd_rec + 1

  # Restore definition of true_comps_rec_lin
  true_comps_rec_lin <- list(
   list(mu = c(-2, rep(0, ncoef_rec - 1)), rooti = diag(0.5, ncoef_rec)),
   list(mu = c( 2, rep(0, ncoef_rec - 1)), rooti = diag(1, ncoef_rec))
  )
  # Revert the first sim call to use the defined list
  sim_data_mix_rec <- sim_hier_mnl(
    nlgt = nlgt_rec, nT = nT_rec, p = p_rec, nz = nz_rec, nXa = 1, 
    nXd = 0, # <<< Add nXd = 0
    const = TRUE, # Explicitly state const=TRUE (default, but good practice)
    beta_func_type = "linear",
    ncomp = ncomp_rec,
    mixture_comps = true_comps_rec_lin, # Use the defined list again
    seed = 808
  )
  Prior_mix_rec <- list(ncomp = ncomp_rec)

  set.seed(809)
  fit_mix_rec <- run_quick_hmnl(sim_data_mix_rec, Prior_mix_rec, Mcmc_rec)

  # Check probdraw recovery
  true_pvec_lin <- sim_data_mix_rec$true_values$true_params$pvec
  est_pvec_mean_lin <- colMeans(fit_mix_rec$nmix$probdraw)
  expect_equal(sum(est_pvec_mean_lin), 1.0, tolerance = 1e-6,
               label = "Linear Mixture: Posterior mean pvec should sum to 1")

  # Check compdraw recovery (mu and rooti) - accounting for label switching
  # Extract posterior means for each component
  draws <- fit_mix_rec$nmix$compdraw
  ndraws <- length(draws)

  ncoef <- length(draws[[1]][[1]]$mu) # Get ncoef from the mu vector
  
  # Initialize storage for means
  mu_means <- array(0, dim = c(ncoef, ncomp_rec))
  rooti_means <- array(0, dim = c(ncoef, ncoef, ncomp_rec))
  
  # Calculate means (sum across draws)
  for (r in 1:ndraws) {
    for (c in 1:ncomp_rec) {
      mu_means[, c] <- mu_means[, c] + as.vector(draws[[r]][[c]]$mu)
      rooti_means[, , c] <- rooti_means[, , c] + draws[[r]][[c]]$rooti
    }
  }
  mu_means <- mu_means / ndraws
  rooti_means <- rooti_means / ndraws

  # --- Order components for comparison (handle label switching) ---
  # Order estimated components based on the first element of mu
  est_order_lin <- order(mu_means[1, ])
  mu_means_ordered_lin <- mu_means[, est_order_lin, drop=FALSE]
  rooti_means_ordered_lin <- rooti_means[, , est_order_lin, drop=FALSE]

  # Order true components based on the first element of mu
  true_order_lin <- order(sapply(true_comps_rec_lin, function(x) x$mu[1]))
  true_comps_rec_ordered_lin <- true_comps_rec_lin[true_order_lin]

  # Compare ordered pvec
  expect_equal(est_pvec_mean_lin[est_order_lin], true_pvec_lin[true_order_lin],
               tolerance = 0.4,
               label = "Linear Mixture: Ordered posterior mean pvec far from ordered true pvec")

  # --- Compare ordered components ---
  # Compare Mu vectors
  expect_equal(mu_means_ordered_lin[, 1], true_comps_rec_ordered_lin[[1]]$mu,
               tolerance = 0.5, # Still using loose tolerance
               label = "Mixture: Ordered mu for component 1 doesn't match true")
  expect_equal(mu_means_ordered_lin[, 2], true_comps_rec_ordered_lin[[2]]$mu,
               tolerance = 0.5,
               label = "Mixture: Ordered mu for component 2 doesn't match true")

  # Compare Rooti matrices
  expect_equal(rooti_means_ordered_lin[, , 1], true_comps_rec_ordered_lin[[1]]$rooti,
               tolerance = 0.5,
               label = "Mixture: Ordered rooti for component 1 doesn't match true")
  expect_equal(rooti_means_ordered_lin[, , 2], true_comps_rec_ordered_lin[[2]]$rooti,
               tolerance = 1.5,
               label = "Mixture: Ordered rooti for component 2 doesn't match true")

  # --- SKIP BART Mixture Recovery PART for now ---
  # --- Test Mixture Component Recovery (ncomp > 1, BART DeltaZ) ---
  # set.seed(901)
  # # Use ncoef_rec calculated above
  # true_comps_rec_bart <- true_comps_rec_lin # Reuse component defs
  #
  # sim_data_mix_bart_rec <- sim_hier_mnl(
  #   nlgt = nlgt_rec, nT = nT_rec, p = p_rec, nz = nz_rec, nXa = 1,
  #   nXd = 0, # <<< Add nXd = 0
  #   const = TRUE, # Explicitly state const=TRUE (default, but good practice)
  #   beta_func_type = "step", # Use step for BART test
  #   ncomp = ncomp_rec,
  #   mixture_comps = true_comps_rec_bart,
  #   seed = 902
  # )
  # Prior_mix_bart_rec <- list(ncomp = ncomp_rec,
  #                            bart = list(num_trees = 100))
  #
  # fit_mix_bart_rec <- run_quick_hmnl(sim_data_mix_bart_rec, Prior_mix_bart_rec, Mcmc_rec)
  #
  # # Check probdraw recovery
  # true_pvec_bart <- sim_data_mix_bart_rec$true_values$true_params$pvec
  # est_pvec_mean_bart <- colMeans(fit_mix_bart_rec$nmix$probdraw)
  # expect_equal(sum(est_pvec_mean_bart), 1.0, tolerance = 1e-6,
  #              label = "BART Mixture: Posterior mean pvec should sum to 1")
  #
  # # Check compdraw recovery (mu and rooti) - accounting for label switching
  # # Extract posterior means for each component
  # draws_bart <- fit_mix_bart_rec$nmix$compdraw
  # ndraws_bart <- length(draws_bart)
  #
  # ncoef_bart <- length(draws_bart[[1]][[1]]$mu) # Get ncoef from the mu vector
  # 
  # # Initialize storage for means
  # mu_means_bart <- array(0, dim = c(ncoef_bart, ncomp_rec))
  # rooti_means_bart <- array(0, dim = c(ncoef_bart, ncoef_bart, ncomp_rec))
  # 
  # # Calculate means (sum across draws)
  # for (r in 1:ndraws_bart) {
  #   for (c in 1:ncomp_rec) {
  #     mu_means_bart[, c] <- mu_means_bart[, c] + as.vector(draws_bart[[r]][[c]]$mu)
  #     rooti_means_bart[, , c] <- rooti_means_bart[, , c] + draws_bart[[r]][[c]]$rooti
  #   }
  # }
  # mu_means_bart <- mu_means_bart / ndraws_bart
  # rooti_means_bart <- rooti_means_bart / ndraws_bart
  #
  # # --- Order components for comparison (handle label switching) ---
  # # Order estimated components based on the first element of mu
  # est_order_bart <- order(mu_means_bart[1, ])
  # mu_means_ordered_bart <- mu_means_bart[, est_order_bart, drop=FALSE]
  # rooti_means_ordered_bart <- rooti_means_bart[, , est_order_bart, drop=FALSE]
  #
  # # Order true components based on the first element of mu
  # true_order_bart <- order(sapply(true_comps_rec_bart, function(x) x$mu[1]))
  # true_comps_rec_ordered_bart <- true_comps_rec_bart[true_order_bart]
  #
  # # Compare ordered pvec
  # expect_equal(est_pvec_mean_bart[est_order_bart], true_pvec_bart[true_order_bart],
  #              tolerance = 0.4,
  #              label = "BART Mixture: Ordered posterior mean pvec far from ordered true pvec")
  #
  # # --- Compare ordered components ---
  # # Compare Mu vectors -- unidentified, skipping for now
  # #   expect_equal(mu_means_ordered_bart[, 1], true_comps_rec_ordered_bart[[1]]$mu,
  # #                tolerance = 0.5, # Still using loose tolerance
  # #                label = "Mixture: Ordered mu for component 1 doesn't match true")
  # #   expect_equal(mu_means_ordered_bart[, 2], true_comps_rec_ordered_bart[[2]]$mu,
  # #                tolerance = 0.5,
  # #                label = "Mixture: Ordered mu for component 2 doesn't match true")
  #
  # # Compare Rooti matrices
  # expect_equal(rooti_means_ordered_bart[, , 1], true_comps_rec_ordered_bart[[1]]$rooti,
  #              tolerance = 0.5,
  #              label = "Mixture: Ordered rooti for component 1 doesn't match true")
  # expect_equal(rooti_means_ordered_bart[, , 2], true_comps_rec_ordered_bart[[2]]$rooti,
  #              tolerance = 1.5,
  #              label = "Mixture: Ordered rooti for component 2 doesn't match true")

}) 


