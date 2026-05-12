# Tests for rhierNegbinRw and sim_hier_negbin

# ==============================================================================
# Test Data Generation Helpers
# ==============================================================================

generate_negbin_test_data <- function(nreg = 50, nobs = 30, nvar = 1, nz = 3,
                                       het = "linear", alpha = 5.0, seed = 42) {
  sim <- sim_hier_negbin(
    nreg = nreg, nobs = nobs, nvar = nvar, nz = nz,
    const = TRUE, het_observed = het,
    target_var_betabar = 0.5, target_var_eps = 0.25,
    alpha = alpha, seed = seed
  )
  return(sim)
}

# ==============================================================================
# 1. sim_hier_negbin Tests
# ==============================================================================

test_that("sim_hier_negbin returns correct structure", {
  sim <- sim_hier_negbin(nreg = 30, nobs = 20, nvar = 1, nz = 3,
                          het_observed = "linear", seed = 1)
  expect_type(sim, "list")
  expect_named(sim, c("regdata", "hessdata", "Z", "true_values"))
  expect_length(sim$regdata, 30)
  expect_length(sim$hessdata, 30)
  expect_equal(nrow(sim$Z), 30)
  expect_equal(ncol(sim$Z), 3)

  # Check regdata structure
  rd1 <- sim$regdata[[1]]
  expect_true(!is.null(rd1$y))
  expect_true(!is.null(rd1$X))
  expect_equal(ncol(rd1$X), 2)  # nvar=1 + intercept
  expect_true(all(rd1$y >= 0))  # counts are non-negative
  expect_true(all(rd1$y == floor(rd1$y)))  # counts are integers

  # Check hessdata structure
  hd1 <- sim$hessdata[[1]]
  expect_true(!is.null(hd1$hess))
  expect_equal(dim(hd1$hess), c(2, 2))
})

test_that("sim_hier_negbin with het_observed='none' works", {
  sim <- sim_hier_negbin(nreg = 20, nobs = 10, nvar = 1, nz = 0,
                          het_observed = "none", seed = 2)
  expect_null(sim$Z)
  expect_equal(nrow(sim$true_values$beta_true), 20)
})

test_that("sim_hier_negbin with het_observed='step' works", {
  sim <- sim_hier_negbin(nreg = 40, nobs = 20, nvar = 1, nz = 3,
                          het_observed = "step", seed = 3)
  expect_equal(nrow(sim$true_values$beta_true), 40)
})

test_that("sim_hier_negbin is reproducible with same seed", {
  sim1 <- sim_hier_negbin(nreg = 20, nobs = 10, nvar = 1, nz = 2,
                           het_observed = "linear", seed = 99)
  sim2 <- sim_hier_negbin(nreg = 20, nobs = 10, nvar = 1, nz = 2,
                           het_observed = "linear", seed = 99)
  expect_equal(sim1$true_values$beta_true, sim2$true_values$beta_true)
  expect_equal(sim1$regdata[[1]]$y, sim2$regdata[[1]]$y)
})

# ==============================================================================
# 2. rhierNegbinRw Basic Execution (non-BART)
# ==============================================================================

test_that("rhierNegbinRw runs with Z (non-BART)", {
  sim <- generate_negbin_test_data(nreg = 20, nobs = 20, nvar = 1, nz = 2)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 50, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierNegbinRw")
  expect_true(!is.null(out$betadraw))
  expect_true(!is.null(out$alphadraw))
  expect_true(!is.null(out$loglike))
  expect_true(!is.null(out$Deltadraw))
  expect_true(!is.null(out$nmix))
})

test_that("rhierNegbinRw runs without Z", {
  sim <- sim_hier_negbin(nreg = 15, nobs = 15, nvar = 1, nz = 0,
                          het_observed = "none", seed = 10)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierNegbinRw")
  expect_equal(dim(out$betadraw)[1], 15)  # nreg
  expect_equal(dim(out$betadraw)[2], 2)   # nvar + intercept
  expect_equal(dim(out$betadraw)[3], 30)  # R/keep
})

# ==============================================================================
# 3. Output Structure
# ==============================================================================

test_that("rhierNegbinRw output dimensions are correct (non-BART)", {
  nreg <- 20; nvar_x <- 1; nz <- 2; R <- 40; keep <- 2
  ncoef <- nvar_x + 1
  sim <- generate_negbin_test_data(nreg = nreg, nobs = 15, nvar = nvar_x, nz = nz)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = R, keep = keep, nprint = 0),
    r_verbose = FALSE
  )

  expect_equal(dim(out$betadraw), c(nreg, ncoef, R / keep))
  expect_equal(length(out$alphadraw), R / keep)
  expect_equal(length(out$loglike), R / keep)
  expect_equal(dim(out$Deltadraw), c(R / keep, ncoef * nz))
  expect_true(!is.null(out$nmix))
  expect_true(!is.null(out$nmix$probdraw))
  expect_equal(nrow(out$nmix$probdraw), R / keep)
})

# ==============================================================================
# 4. Input Validation
# ==============================================================================

test_that("rhierNegbinRw validates missing Data", {
  expect_error(rhierNegbinRw(Prior = list(ncomp = 1), Mcmc = list(R = 10)),
               "Requires Data argument")
})

test_that("rhierNegbinRw validates missing regdata", {
  expect_error(rhierNegbinRw(Data = list(), Prior = list(ncomp = 1), Mcmc = list(R = 10)),
               "Requires Data element regdata")
})

test_that("rhierNegbinRw validates missing R", {
  sim <- generate_negbin_test_data(nreg = 10, nobs = 10, nvar = 1, nz = 2)
  expect_error(rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(keep = 1)
  ), "Requires element R of Mcmc")
})

test_that("rhierNegbinRw validates Z dimensions", {
  sim <- generate_negbin_test_data(nreg = 10, nobs = 10, nvar = 1, nz = 2)
  bad_Z <- sim$Z[1:5, ]
  expect_error(rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = bad_Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 10, nprint = 0)
  ), "ne number units")
})

# ==============================================================================
# 5. Acceptance Rate
# ==============================================================================

test_that("rhierNegbinRw returns acceptance rates", {
  sim <- generate_negbin_test_data(nreg = 15, nobs = 15, nvar = 1, nz = 2)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_true(out$acceptrbeta >= 0 && out$acceptrbeta <= 100)
  expect_true(out$acceptralpha >= 0 && out$acceptralpha <= 100)
})

# ==============================================================================
# 6. fixalpha
# ==============================================================================

test_that("rhierNegbinRw with fixalpha keeps alpha constant", {
  sim <- generate_negbin_test_data(nreg = 15, nobs = 15, nvar = 1, nz = 2)
  alpha_fixed <- 3.0
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0, alpha = alpha_fixed, fixalpha = TRUE),
    r_verbose = FALSE
  )
  expect_true(all(out$alphadraw == alpha_fixed))
})

# ==============================================================================
# 7. BART Integration Tests
# ==============================================================================

test_that("rhierNegbinRw runs with BART prior", {
  sim <- generate_negbin_test_data(nreg = 25, nobs = 20, nvar = 1, nz = 3)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10)),
    Mcmc = list(R = 50, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierNegbinRw")
  expect_true(!is.null(out$bart_models))
  expect_true(!is.null(out$nmix))
  expect_true(!is.null(out$betadraw))
  expect_true(!is.null(out$alphadraw))
  # BART path should NOT have Deltadraw
  expect_null(out$Deltadraw)
  expect_true(!is.null(out$nmix))
})

test_that("rhierNegbinRw BART output dimensions are correct", {
  nreg <- 20; nvar_x <- 1; nz <- 3; R <- 40; keep <- 2
  ncoef <- nvar_x + 1
  sim <- generate_negbin_test_data(nreg = nreg, nobs = 15, nvar = nvar_x, nz = nz)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10)),
    Mcmc = list(R = R, keep = keep, nprint = 0),
    r_verbose = FALSE
  )

  expect_equal(dim(out$betadraw), c(nreg, ncoef, R / keep))
  expect_equal(length(out$alphadraw), R / keep)
  expect_equal(length(out$bart_models), ncoef)
  expect_equal(dim(out$varcount), c(nz, ncoef, R / keep))
  expect_equal(dim(out$varprob), c(nz, ncoef, R / keep))
})

test_that("rhierNegbinRw BART with sparse/DART works", {
  sim <- generate_negbin_test_data(nreg = 15, nobs = 15, nvar = 1, nz = 3)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10, sparse = TRUE, burn = 5)),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierNegbinRw")
  expect_true(!is.null(out$bart_models))
})

test_that("rhierNegbinRw BART with fixalpha works", {
  sim <- generate_negbin_test_data(nreg = 15, nobs = 15, nvar = 1, nz = 3)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10)),
    Mcmc = list(R = 30, keep = 1, nprint = 0, alpha = 5.0, fixalpha = TRUE),
    r_verbose = FALSE
  )
  expect_true(all(out$alphadraw == 5.0))
})

# ==============================================================================
# 8. predict.rhierNegbinRw Tests
# ==============================================================================

test_that("predict.rhierNegbinRw DeltaZ works (non-BART)", {
  sim <- generate_negbin_test_data(nreg = 20, nobs = 20, nvar = 1, nz = 2)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  pred <- predict(out, newdata = list(Z = sim$Z), type = "DeltaZ", burn = 0)
  expect_true(is.array(pred))
  expect_equal(length(dim(pred)), 3)
  expect_equal(dim(pred)[1], 20)  # npred = nreg
  expect_equal(dim(pred)[2], 2)   # ncoef = nvar + intercept
  expect_equal(dim(pred)[3], 30)  # ndraws = R/keep
})

test_that("predict.rhierNegbinRw burn-in discards draws", {
  sim <- generate_negbin_test_data(nreg = 15, nobs = 15, nvar = 1, nz = 2)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 40, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  pred_full <- predict(out, newdata = list(Z = sim$Z), type = "DeltaZ", burn = 0)
  pred_burn <- predict(out, newdata = list(Z = sim$Z), type = "DeltaZ", burn = 10)
  expect_equal(dim(pred_full)[3], 40)
  expect_equal(dim(pred_burn)[3], 30)
})

test_that("predict.rhierNegbinRw DeltaZ+mu works (BART)", {
  sim <- generate_negbin_test_data(nreg = 20, nobs = 20, nvar = 1, nz = 3)
  set.seed(222)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10)),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  pred <- predict(out, newdata = list(Z = sim$Z), type = "DeltaZ+mu",
                  burn = 5, r_verbose = FALSE)
  expect_true(is.array(pred))
  expect_equal(dim(pred)[1], 20)  # npred
  expect_equal(dim(pred)[2], 2)   # ncoef
  expect_equal(dim(pred)[3], 25)  # 30 - 5 burn
})

test_that("predict.rhierNegbinRw errors on invalid type", {
  sim <- generate_negbin_test_data(nreg = 10, nobs = 10, nvar = 1, nz = 2)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 10, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_error(predict(out, newdata = list(Z = sim$Z), type = "prior_probs"),
               "Invalid type")
})

# ==============================================================================
# 9. marginal_effects.rhierNegbinRw Tests
# ==============================================================================

test_that("marginal_effects.rhierNegbinRw works (non-BART)", {
  sim <- generate_negbin_test_data(nreg = 20, nobs = 20, nvar = 1, nz = 2)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  z_grid <- matrix(NA, nrow = 3, ncol = 2)
  z_grid[, 1] <- c(-1, 0, 1)

  mfx <- marginal_effects(out, z_values = z_grid, Z = sim$Z,
                            burn = 5, verbose = FALSE)
  expect_s3_class(mfx, "marginal_effects")
  expect_length(mfx$avg_betabar_draws, 3)
  expect_equal(nrow(mfx$avg_betabar_draws[[1]]), 2)  # ncoef
  expect_equal(ncol(mfx$avg_betabar_draws[[1]]), 25)  # 30 - 5
})

test_that("summary.marginal_effects works on negbin model results", {
  sim <- generate_negbin_test_data(nreg = 15, nobs = 15, nvar = 1, nz = 2)
  out <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 20, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  z_grid <- matrix(NA, nrow = 3, ncol = 2)
  z_grid[, 1] <- c(-1, 0, 1)

  mfx <- marginal_effects(out, z_values = z_grid, Z = sim$Z,
                            burn = 5, verbose = FALSE)
  smry <- summary(mfx, probs = c(0.025, 0.5, 0.975))
  expect_s3_class(smry, "summary.marginal_effects")
  expect_true("summary_df" %in% names(smry))
  expect_true("mean" %in% names(smry$summary_df))
  expect_true("q25" %in% names(smry$summary_df))
})
