# ==============================================================================
# marginal_effects() smoke tests under heter-cov for all three model families.
#
# These tests verify that:
#   1. The generic dispatches correctly on heter-cov subclasses
#      (rhierMnlRwMixtureHeterCov / rhierLinearMixtureHeterCov /
#       rhierNegbinRwHeterCov inherit from their respective base classes).
#   2. .marginal_effects_impl() picks `predict_type = "DeltaZ+mu"` for objects
#      that carry mu in `mu_draw` (heter-cov) just as it does for those that
#      carry it in `nmix$compdraw` (linear / mix-BART).
#   3. summary() returns a well-formed `summary.marginal_effects` object.
#   4. Effects vary across grid points (sanity: not constant in z).
# ==============================================================================

library(testthat)

# --- shared helper ----------------------------------------------------------
.run_mfx_check <- function(fit, Z, ncoef) {
  nz <- ncol(Z)
  z_grid <- matrix(NA_real_, nrow = 4L, ncol = nz)
  z_grid[, 1] <- seq(min(Z[, 1]), max(Z[, 1]), length.out = 4L)

  mfx <- marginal_effects(fit, z_values = z_grid, Z = Z,
                          burn = 5L, verbose = FALSE)
  expect_s3_class(mfx, "marginal_effects")
  expect_equal(nrow(mfx$z_values), 4L)
  expect_length(mfx$avg_betabar_draws, 4L)
  expect_equal(nrow(mfx$avg_betabar_draws[[1]]), ncoef)

  mfx_summary <- summary(mfx, probs = c(0.05, 0.5, 0.95))
  expect_s3_class(mfx_summary, "summary.marginal_effects")
  expect_true(all(c("coefficient_index", "z_val_1",
                    "mean", "q5", "q50", "q95") %in% names(mfx_summary$summary_df)))
  expect_equal(nrow(mfx_summary$summary_df), 4L * ncoef)

  # Effects should vary across grid points (heter-cov delta_Z is non-trivial
  # in Z so the average across units shouldn't collapse to a single value).
  for (k in seq_len(ncoef)) {
    means_k <- subset(mfx_summary$summary_df,
                      mfx_summary$summary_df$coefficient_index == k)$mean
    expect_true(diff(range(means_k)) > 1e-8,
                info = sprintf("coef %d means are constant in z", k))
  }
  invisible(mfx_summary)
}

# --- MNL --------------------------------------------------------------------
test_that("marginal_effects() works for rhierMnlRwMixtureHeterCov", {
  set.seed(20260513L)
  nlgt <- 25L; nvar <- 2L; nz <- 3L; p <- 3L; T_i <- 6L
  Z <- cbind(runif(nlgt, -1, 1), runif(nlgt, -1, 1), runif(nlgt, -1, 1))
  Z <- scale(Z, center = TRUE, scale = FALSE)
  theta_true <- matrix(rnorm(nlgt * nvar, sd = 0.5), nlgt, nvar)
  lgtdata <- lapply(seq_len(nlgt), function(i) {
    X    <- matrix(rnorm(T_i * p * nvar), T_i * p, nvar)
    eta  <- X %*% theta_true[i, ]
    prob <- matrix(exp(eta), T_i, p, byrow = TRUE)
    prob <- prob / rowSums(prob)
    y    <- apply(prob, 1, function(pp) sample.int(p, 1, prob = pp))
    list(y = y, X = X)
  })
  Data <- list(p = p, lgtdata = lgtdata, Z = Z)
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 30L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierMnlRwMixture(Data = Data, Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })
  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
  .run_mfx_check(fit, Z = Z, ncoef = nvar)
})

# --- Linear -----------------------------------------------------------------
test_that("marginal_effects() works for rhierLinearMixtureHeterCov", {
  sim <- sim_hier_linear(
    nreg = 25L, nobs = 12L, nvar = 1L, nz = 3L,
    const = TRUE, het_observed = "linear",
    target_var_betabar = 1.0, target_var_eps = 0.5,
    sigma_sq = 0.5, seed = 20260513L)
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 30L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierLinearMixture(
      Data = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })
  expect_s3_class(fit, "rhierLinearMixtureHeterCov")
  ncoef <- ncol(sim$regdata[[1]]$X)
  .run_mfx_check(fit, Z = sim$Z, ncoef = ncoef)
})

# --- Negbin -----------------------------------------------------------------
test_that("marginal_effects() works for rhierNegbinRwHeterCov", {
  sim <- sim_hier_negbin(
    nreg = 25L, nobs = 12L, nvar = 1L, nz = 3L,
    const = TRUE, het_observed = "linear",
    target_var_betabar = 0.5, target_var_eps = 0.25,
    alpha = 5.0, seed = 20260513L)
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 30L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierNegbinRw(
      Data = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })
  expect_s3_class(fit, "rhierNegbinRwHeterCov")
  ncoef <- ncol(sim$regdata[[1]]$X)
  .run_mfx_check(fit, Z = sim$Z, ncoef = ncoef)
})
