# Tests for rhierLinearMixture (mixture-of-normals linear hierarchical model)
# Mirrors structure of test-rhierLinearModel.R with API adjustments:
#   - Prior must include ncomp
#   - Output has nmix instead of Vbetadraw
#   - Class is "rhierLinearMixture"

# ==============================================================================
# Test Data Generation Helpers
# ==============================================================================

generate_mixture_linear_data <- function(nreg = 50, nobs = 20, nvar = 2, nz = 3,
                                         het = "linear", seed = 42) {
  sim <- sim_hier_linear(
    nreg = nreg, nobs = nobs, nvar = nvar, nz = nz,
    const = TRUE, het_observed = het,
    target_var_betabar = 1.0, target_var_eps = 0.5,
    sigma_sq = 1.0, seed = seed
  )
  return(sim)
}

# ==============================================================================
# 1. Basic Execution
# ==============================================================================

test_that("rhierLinearMixture runs with Z and ncomp=1", {
  sim <- generate_mixture_linear_data(nreg = 30, nobs = 15, nvar = 2, nz = 2)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 50, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierLinearMixture")
  expect_true(!is.null(out$betadraw))
  expect_true(!is.null(out$taudraw))
  expect_true(!is.null(out$Deltadraw))
  expect_true(!is.null(out$nmix))
  expect_true(!is.null(out$loglike))
  expect_equal(length(out$loglike), 50)
})

test_that("rhierLinearMixture runs without Z", {
  sim <- sim_hier_linear(nreg = 20, nobs = 10, nvar = 1, nz = 0,
                          het_observed = "none", seed = 10)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierLinearMixture")
  expect_equal(dim(out$betadraw)[1], 20)  # nreg
  expect_equal(dim(out$betadraw)[2], 2)   # nvar + intercept
  expect_equal(dim(out$betadraw)[3], 30)  # R/keep
  # no-Z path should not have Deltadraw
  expect_null(out$Deltadraw)
  expect_true(!is.null(out$nmix))
})

test_that("rhierLinearMixture rejects ncomp > 1", {
  sim <- generate_mixture_linear_data(nreg = 40, nobs = 15, nvar = 1, nz = 2)
  expect_error(
    rhierLinearMixture(
      Data = list(regdata = sim$regdata, Z = sim$Z),
      Prior = list(ncomp = 2),
      Mcmc = list(R = 50, keep = 1, nprint = 0),
      r_verbose = FALSE
    ),
    regexp = "Only ncomp = 1 is currently supported"
  )
})

# ==============================================================================
# 2. Output Structure
# ==============================================================================

test_that("rhierLinearMixture output dimensions are correct", {
  nreg <- 25; nvar_x <- 2; nz <- 3; R <- 40; keep <- 2
  ncoef <- nvar_x + 1  # with intercept
  sim <- generate_mixture_linear_data(nreg = nreg, nobs = 10, nvar = nvar_x, nz = nz)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = R, keep = keep, nprint = 0),
    r_verbose = FALSE
  )

  expect_equal(dim(out$betadraw), c(nreg, ncoef, R / keep))
  expect_equal(dim(out$taudraw), c(R / keep, nreg))
  expect_equal(dim(out$Deltadraw), c(R / keep, nz * ncoef))
  expect_equal(length(out$loglike), R / keep)
  # nmix structure
  expect_true(!is.null(out$nmix$probdraw))
  expect_true(!is.null(out$nmix$compdraw))
  expect_equal(nrow(out$nmix$probdraw), R / keep)
})

test_that("rhierLinearMixture has correct S3 class attributes", {
  sim <- generate_mixture_linear_data(nreg = 20, nobs = 10, nvar = 1, nz = 2)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 20, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierLinearMixture")
  expect_true("bayesm.hcoef" %in% class(out$betadraw))
  expect_true("bayesm.nmix" %in% class(out$nmix))
  expect_true("bayesm.mat" %in% class(out$taudraw))
  expect_true("mcmc" %in% class(out$taudraw))
})

# ==============================================================================
# 3. Input Validation
# ==============================================================================

test_that("rhierLinearMixture validates missing Data", {
  expect_error(rhierLinearMixture(Prior = list(ncomp = 1), Mcmc = list(R = 10)),
               "Requires Data argument")
})

test_that("rhierLinearMixture validates missing regdata", {
  expect_error(rhierLinearMixture(Data = list(), Prior = list(ncomp = 1), Mcmc = list(R = 10)),
               "Requires Data element regdata")
})

test_that("rhierLinearMixture requires ncomp in Prior", {
  sim <- generate_mixture_linear_data(nreg = 10, nobs = 5, nvar = 1, nz = 2)
  expect_error(rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(),
    Mcmc = list(R = 10, nprint = 0)
  ), "ncomp")
})

test_that("rhierLinearMixture validates missing R", {
  sim <- generate_mixture_linear_data(nreg = 10, nobs = 5, nvar = 1, nz = 2)
  expect_error(rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(keep = 1)
  ), "Requires R argument")
})

test_that("rhierLinearMixture validates Z dimensions", {
  sim <- generate_mixture_linear_data(nreg = 10, nobs = 5, nvar = 1, nz = 2)
  bad_Z <- sim$Z[1:5, ]  # Wrong number of rows
  expect_error(rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = bad_Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 10, nprint = 0)
  ), "ne number regressions")
})

test_that("rhierLinearMixture validates Z must be matrix", {
  sim <- generate_mixture_linear_data(nreg = 10, nobs = 5, nvar = 1, nz = 2)
  expect_error(rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = as.data.frame(sim$Z)),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 10, nprint = 0)
  ), "Z must be a matrix")
})

# ==============================================================================
# 4. Prior Specification
# ==============================================================================

test_that("rhierLinearMixture accepts custom priors", {
  sim <- generate_mixture_linear_data(nreg = 20, nobs = 10, nvar = 1, nz = 2)
  ncoef <- 2  # 1 + intercept
  nz <- 2
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(
      ncomp = 1,
      nu = 10, V = 10 * diag(ncoef),
      nu.e = 5, ssq = rep(2.0, 20),
      deltabar = rep(0, nz * ncoef),
      Ad = 0.1 * diag(nz * ncoef)
    ),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierLinearMixture")
})

# ==============================================================================
# 5. Reproducibility
# ==============================================================================

test_that("rhierLinearMixture is reproducible with same seed", {
  sim <- generate_mixture_linear_data(nreg = 15, nobs = 8, nvar = 1, nz = 2, seed = 77)

  set.seed(123)
  out1 <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 20, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  set.seed(123)
  out2 <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 20, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  expect_equal(out1$betadraw, out2$betadraw)
  expect_equal(out1$taudraw, out2$taudraw)
  expect_equal(out1$Deltadraw, out2$Deltadraw)
})

# ==============================================================================
# 6. Parameter Recovery
# ==============================================================================

test_that("rhierLinearMixture recovers beta (linear DGP, moderate R)", {
  sim <- sim_hier_linear(
    nreg = 100, nobs = 30, nvar = 2, nz = 3,
    het_observed = "linear",
    target_var_betabar = 2.0, target_var_eps = 0.5,
    sigma_sq = 0.5, seed = 42
  )

  set.seed(999)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 2000, keep = 2, nprint = 0),
    r_verbose = FALSE
  )

  # Posterior mean of beta for each unit (drop burn-in)
  ndraws <- dim(out$betadraw)[3]
  burn <- floor(ndraws / 2)
  beta_post_mean <- apply(out$betadraw[, , (burn + 1):ndraws], c(1, 2), mean)

  # Check correlation between posterior mean and true beta
  ncoef <- ncol(sim$true_values$beta_true)
  for (k in 1:ncoef) {
    cor_k <- cor(beta_post_mean[, k], sim$true_values$beta_true[, k])
    expect_gt(cor_k, 0.5,
              label = paste("Correlation for coefficient", k, "=", round(cor_k, 3)))
  }
})

test_that("rhierLinearMixture precision matrix convention is correct (M3.2 bug check)", {
  # This test uses nvar = 2 (ncoef = 3) with an asymmetric Cholesky root to
  # ensure that the precision matrix is computed as rootpi * t(rootpi) and
  # not t(rootpi) * rootpi.
  set.seed(42)
  nreg <- 1
  nobs <- 10
  ncoef <- 3
  
  X <- cbind(1, matrix(rnorm(nobs * 2), nobs, 2))
  beta_true <- c(1.5, -0.5, 0.8)
  y <- X %*% beta_true + rnorm(nobs, 0, 0.5)
  
  # Asymmetric V ensures R'R != RR'
  Sigma_k <- matrix(c(1, 0.8, 0.2, 
                      0.8, 2, -0.5,
                      0.2, -0.5, 1.5), 3, 3)
  
  Prior <- list(
    ncomp = 1,
    mubar = matrix(c(0, 0, 0), 3, 1),
    Amu = matrix(100, 1, 1),
    nu = 100,
    V = Sigma_k * 100,
    nu.e = 3,
    ssq = rep(1.0, 1)
  )
  
  Data <- list(regdata = list(list(y=y, X=X)))
  
  # Run the sampler
  set.seed(123)
  out <- rhierLinearMixture(Data=Data, Prior=Prior, 
                            Mcmc=list(R=1000, keep=1, nprint=0), r_verbose=FALSE)
  
  # Get posterior mean
  beta_post_mean <- apply(out$betadraw[1,,501:1000], 1, mean)
  
  # Check against known correct values from the patched implementation
  expect_equal(beta_post_mean[1], 1.04265, tolerance = 0.05)
  expect_equal(beta_post_mean[2], -0.04129, tolerance = 0.05)
  expect_equal(beta_post_mean[3],  0.82022, tolerance = 0.05)
})

# ==============================================================================
# 7. BART Integration
# ==============================================================================

test_that("rhierLinearMixture runs with BART prior", {
  sim <- generate_mixture_linear_data(nreg = 30, nobs = 15, nvar = 2, nz = 3)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10)),
    Mcmc = list(R = 50, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierLinearMixture")
  expect_true(!is.null(out$bart_models))
  expect_true(!is.null(out$nmix))
  expect_true(!is.null(out$betadraw))
  expect_true(!is.null(out$taudraw))
  # BART path should NOT have Deltadraw
  expect_null(out$Deltadraw)
})

test_that("rhierLinearMixture BART output dimensions are correct", {
  nreg <- 25; nvar_x <- 2; nz <- 3; R <- 40; keep <- 2
  ncoef <- nvar_x + 1
  sim <- generate_mixture_linear_data(nreg = nreg, nobs = 10, nvar = nvar_x, nz = nz)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10)),
    Mcmc = list(R = R, keep = keep, nprint = 0),
    r_verbose = FALSE
  )

  expect_equal(dim(out$betadraw), c(nreg, ncoef, R / keep))
  expect_equal(dim(out$taudraw), c(R / keep, nreg))
  expect_equal(length(out$bart_models), ncoef)
  expect_equal(dim(out$varcount), c(nz, ncoef, R / keep))
  expect_equal(dim(out$varprob), c(nz, ncoef, R / keep))
})

test_that("rhierLinearMixture BART with sparse/DART works", {
  sim <- generate_mixture_linear_data(nreg = 20, nobs = 10, nvar = 1, nz = 3)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10, sparse = TRUE, burn = 5)),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierLinearMixture")
  expect_true(!is.null(out$bart_models))
})

test_that("rhierLinearMixture BART with minimal trees works", {
  sim <- generate_mixture_linear_data(nreg = 15, nobs = 8, nvar = 1, nz = 2)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 2)),
    Mcmc = list(R = 20, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_s3_class(out, "rhierLinearMixture")
  expect_equal(length(out$bart_models), 2)  # ncoef = 1 + intercept
})

test_that("rhierLinearMixture BART recovers beta (step DGP)", {
  sim <- sim_hier_linear(
    nreg = 100, nobs = 30, nvar = 2, nz = 3,
    het_observed = "step",
    target_var_betabar = 2.0, target_var_eps = 0.5,
    sigma_sq = 0.5, seed = 42
  )

  set.seed(888)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 50)),
    Mcmc = list(R = 2000, keep = 2, nprint = 0),
    r_verbose = FALSE
  )

  ndraws <- dim(out$betadraw)[3]
  burn <- floor(ndraws / 2)
  beta_post_mean <- apply(out$betadraw[, , (burn + 1):ndraws], c(1, 2), mean)

  ncoef <- ncol(sim$true_values$beta_true)
  for (k in 1:ncoef) {
    cor_k <- cor(beta_post_mean[, k], sim$true_values$beta_true[, k])
    expect_gt(cor_k, 0.5,
              label = paste("BART recovery: coefficient", k, "cor =", round(cor_k, 3)))
  }
})

# ==============================================================================
# 8. predict.rhierLinearMixture Tests
# ==============================================================================

test_that("predict.rhierLinearMixture DeltaZ works (non-BART)", {
  sim <- generate_mixture_linear_data(nreg = 30, nobs = 15, nvar = 2, nz = 3)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  pred <- predict(out, newdata = list(Z = sim$Z), type = "DeltaZ", burn = 0)
  expect_true(is.array(pred))
  expect_equal(length(dim(pred)), 3)
  expect_equal(dim(pred)[1], 30)  # npred = nreg
  expect_equal(dim(pred)[2], 3)   # ncoef = nvar + intercept
  expect_equal(dim(pred)[3], 30)  # ndraws = R/keep
})

test_that("predict.rhierLinearMixture DeltaZ+mu works (non-BART)", {
  sim <- generate_mixture_linear_data(nreg = 25, nobs = 15, nvar = 1, nz = 2)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  pred <- predict(out, newdata = list(Z = sim$Z), type = "DeltaZ+mu",
                  burn = 5, r_verbose = FALSE)
  expect_true(is.array(pred))
  expect_equal(dim(pred)[1], 25)
  expect_equal(dim(pred)[2], 2)   # ncoef
  expect_equal(dim(pred)[3], 25)  # 30 - 5
})

test_that("predict.rhierLinearMixture burn-in discards draws", {
  sim <- generate_mixture_linear_data(nreg = 20, nobs = 10, nvar = 1, nz = 2)
  out <- rhierLinearMixture(
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

test_that("predict.rhierLinearMixture DeltaZ+mu works (BART)", {
  sim <- generate_mixture_linear_data(nreg = 25, nobs = 15, nvar = 2, nz = 3)
  set.seed(111)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 10)),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  pred <- predict(out, newdata = list(Z = sim$Z), type = "DeltaZ+mu",
                  burn = 5, r_verbose = FALSE)
  expect_true(is.array(pred))
  expect_equal(dim(pred)[1], 25)  # npred
  expect_equal(dim(pred)[2], 3)   # ncoef
  expect_equal(dim(pred)[3], 25)  # 30 - 5 burn
})

test_that("predict.rhierLinearMixture errors on invalid type", {
  sim <- generate_mixture_linear_data(nreg = 15, nobs = 10, nvar = 1, nz = 2)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 10, keep = 1, nprint = 0),
    r_verbose = FALSE
  )
  expect_error(predict(out, newdata = list(Z = sim$Z), type = "posterior_probs"),
               "Invalid type")
  expect_error(
    predict(out, newdata = list(Z = sim$Z), type = "SigmaZ"),
    "only available for heteroscedastic covariance models"
  )
})

# ==============================================================================
# 9. marginal_effects.rhierLinearMixture Tests
# ==============================================================================

test_that("marginal_effects.rhierLinearMixture works (non-BART)", {
  sim <- generate_mixture_linear_data(nreg = 25, nobs = 15, nvar = 2, nz = 3)
  out <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1),
    Mcmc = list(R = 30, keep = 1, nprint = 0),
    r_verbose = FALSE
  )

  z_grid <- matrix(NA, nrow = 3, ncol = 3)
  z_grid[, 1] <- c(-1, 0, 1)

  mfx <- marginal_effects(out, z_values = z_grid, Z = sim$Z,
                            burn = 5, verbose = FALSE)
  expect_s3_class(mfx, "marginal_effects")
  expect_length(mfx$avg_betabar_draws, 3)
  expect_equal(nrow(mfx$avg_betabar_draws[[1]]), 3)  # ncoef
  expect_equal(ncol(mfx$avg_betabar_draws[[1]]), 25)  # 30 - 5
})

test_that("summary.marginal_effects works on rhierLinearMixture results", {
  sim <- generate_mixture_linear_data(nreg = 20, nobs = 10, nvar = 1, nz = 2)
  out <- rhierLinearMixture(
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
