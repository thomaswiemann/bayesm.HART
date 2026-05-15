# ==============================================================================
# predict() extrapolation regression tests under heter-cov.
#
# Verifies that the heter-cov predict.* methods evaluate Sigma(Z*) and
# delta_Z(Z*) at user-supplied Z* (including Z* values OUTSIDE the support of
# the training Z) without crashing or producing non-finite output.
#
# These tests guard against:
#   * dimension mismatches when Z* has different nrow than the training Z,
#   * NaN propagation through Sigma(Z*) for unseen Z values,
#   * tree extrapolation past the training boundaries.
# ==============================================================================

library(testthat)

# Build an extrapolation Z*: 5 in-support rows, 5 well outside the training
# Z range (chosen so that BART must extrapolate in every Z column).
.make_extrapolation_Z <- function(Z) {
  nz <- ncol(Z)
  z_in  <- Z[seq_len(min(5L, nrow(Z))), , drop = FALSE]
  rng   <- apply(Z, 2, range)
  span  <- rng[2, ] - rng[1, ]
  # 5 extrapolation rows: 2.5 * span beyond each boundary, alternating side.
  z_out <- matrix(0, 5L, nz)
  for (j in seq_len(nz)) {
    z_out[, j] <- c(rng[2, j] + 2.5 * span[j],
                    rng[1, j] - 2.5 * span[j],
                    rng[2, j] + 1.5 * span[j],
                    rng[1, j] - 1.5 * span[j],
                    rng[2, j] + 3.0 * span[j])
  }
  rbind(z_in, z_out)
}

.check_pred <- function(pred, expected_nrow, expected_ncoef, expected_ndraws) {
  expect_equal(dim(pred), c(expected_nrow, expected_ncoef, expected_ndraws))
  expect_true(all(is.finite(pred)),
              info = sprintf("predictions contain non-finite values (n=%d)",
                             sum(!is.finite(pred))))
}

.check_sigma_pred <- function(pred, expected_nrow, expected_ncoef, expected_ndraws) {
  expect_equal(dim(pred), c(expected_nrow, expected_ncoef, expected_ncoef, expected_ndraws))
  expect_true(all(is.finite(pred)),
              info = sprintf("SigmaZ contains non-finite values (n=%d)",
                             sum(!is.finite(pred))))
  for (i in seq_len(expected_nrow)) {
    for (s in seq_len(expected_ndraws)) {
      expect_equal(pred[i, , , s], t(pred[i, , , s]), tolerance = 1e-8,
                   ignore_attr = TRUE)
      expect_true(all(diag(pred[i, , , s]) > 0))
    }
  }
}

# --- MNL --------------------------------------------------------------------
test_that("predict.rhierMnlRwMixture extrapolates over Z* outside training support", {
  set.seed(20260513L)
  nlgt <- 30L; nvar <- 2L; nz <- 3L; p <- 3L; T_i <- 6L
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
  Z_star <- .make_extrapolation_Z(Z)
  # Focus on Sigma(Z*) / delta_Z(Z*) extrapolation; predictive choice_probs
  # need an additional `X` list keyed to the new units and are covered by the
  # dedicated test-rhierMnlRwMixture-hetercov.R predict.* dispatch tests.
  for (tp in c("DeltaZ", "DeltaZ+mu")) {
    pred <- predict(fit, newdata = list(Z = Z_star), type = tp,
                    burn = 5L, r_verbose = FALSE)
    .check_pred(pred, nrow(Z_star), nvar, 30L - 5L)
  }
  pred_sigma <- predict(fit, newdata = list(Z = Z_star), type = "SigmaZ",
                        burn = 5L, r_verbose = FALSE)
  .check_sigma_pred(pred_sigma, nrow(Z_star), nvar, 30L - 5L)
})

# --- Linear -----------------------------------------------------------------
test_that("predict.rhierLinearMixture extrapolates over Z* outside training support", {
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
  Z_star <- .make_extrapolation_Z(sim$Z)
  ncoef  <- ncol(sim$regdata[[1]]$X)
  for (tp in c("DeltaZ", "DeltaZ+mu")) {
    pred <- predict(fit, newdata = list(Z = Z_star), type = tp,
                    burn = 5L, r_verbose = FALSE)
    .check_pred(pred, nrow(Z_star), ncoef, 30L - 5L)
  }
  pred_sigma <- predict(fit, newdata = list(Z = Z_star), type = "SigmaZ",
                        burn = 5L, r_verbose = FALSE)
  .check_sigma_pred(pred_sigma, nrow(Z_star), ncoef, 30L - 5L)
})

# --- Negbin -----------------------------------------------------------------
test_that("predict.rhierNegbinRw extrapolates over Z* outside training support", {
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
  Z_star <- .make_extrapolation_Z(sim$Z)
  ncoef  <- ncol(sim$regdata[[1]]$X)
  for (tp in c("DeltaZ", "DeltaZ+mu")) {
    pred <- predict(fit, newdata = list(Z = Z_star), type = tp,
                    burn = 5L, r_verbose = FALSE)
    .check_pred(pred, nrow(Z_star), ncoef, 30L - 5L)
  }
  pred_sigma <- predict(fit, newdata = list(Z = Z_star), type = "SigmaZ",
                        burn = 5L, r_verbose = FALSE)
  .check_sigma_pred(pred_sigma, nrow(Z_star), ncoef, 30L - 5L)
})
