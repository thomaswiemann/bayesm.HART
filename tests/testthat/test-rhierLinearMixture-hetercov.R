# Tests for the heteroscedastic-covariance Linear path
# (Prior$vartree, optionally Prior$phitree).
#
# Smoke tests focused on:
#   * input validation (vartree without bart, ncomp != 1, no Z),
#   * end-to-end execution via the user wrapper,
#   * class assignment ("bayesm.HART.HeterCov" marker) and return-slot shape,
#   * predict.* dispatch on the marker class for "DeltaZ" and "DeltaZ+mu",
#   * the covariance-tree banner under r_verbose = TRUE.
#
# Recovery accuracy is exercised by dev/smoke-heter-linear.R.

library(testthat)

# --- shared minimal-cost simulation -----------------------------------------
make_small_linear_sim <- function(seed = 20260513L,
                                  nreg = 25L, nobs = 15L,
                                  nvar = 1L, nz = 2L,
                                  sigma_sq = 0.5) {
  sim_hier_linear(
    nreg = nreg, nobs = nobs, nvar = nvar, nz = nz,
    const = TRUE, het_observed = "linear",
    target_var_betabar = 1.0, target_var_eps = 0.5,
    sigma_sq = sigma_sq, seed = seed
  )
}

# --- input validation --------------------------------------------------------
test_that("vartree requires hart_model, ncomp == 1, and Z (drawdelta = TRUE)", {
  sim <- make_small_linear_sim()

  # vartree without hart_model -> error
  expect_error(
    capture.output(
      rhierLinearMixture(
        Data  = list(regdata = sim$regdata, Z = sim$Z),
        Prior = list(ncomp = 1L, vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "require Prior\\$hart_model"
  )

  # ncomp > 1 with vartree -> error
  expect_error(
    capture.output(
      rhierLinearMixture(
        Data  = list(regdata = sim$regdata, Z = sim$Z),
        Prior = list(ncomp = 2L, hart_model = "direct",
                     bart    = list(num_trees = 5L),
                     vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "Only ncomp = 1 is currently supported"
  )

  # vartree without Z -> error (drawdelta auto-flips to FALSE)
  expect_error(
    capture.output(
      rhierLinearMixture(
        Data  = list(regdata = sim$regdata),       # no Z
        Prior = list(ncomp = 1L, hart_model = "direct",
                     bart    = list(num_trees = 5L),
                     vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "requires Data\\$Z"
  )
})

# --- vartree only with nvar > 1 (auto-promotes to full Cholesky) ------------
test_that("vartree only with nvar > 1 auto-promotes to full Cholesky", {
  # Default sim has nvar = 1 + intercept = ncoef = 2, so auto-promotion fires.
  sim   <- make_small_linear_sim()
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))   # NO phitree
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierLinearMixture(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierLinearMixtureHeterCov")
  expect_s3_class(fit, "bayesm.HART.HeterCov")
  expect_s3_class(fit, "rhierLinearMixture")

  expected_slots <- c("bart_models", "var_models", "phi_models", "mu_draw",
                      "mean_scale_draw",
                      "varcount", "varprob", "var_varcount", "var_varprob",
                      "betadraw", "taudraw", "loglike", "hart")
  expect_named(fit, expected_slots, ignore.order = TRUE)
  expect_null(fit$nmix)
  expect_null(fit$Deltadraw)

  ncoef <- ncol(sim$regdata[[1]]$X)
  # phi_models is jagged.
  expect_length(fit$phi_models, ncoef)
  expect_length(fit$phi_models[[2]], 1L)

  expect_equal(dim(fit$mu_draw),  c(20L, ncoef))
  expect_equal(dim(fit$betadraw), c(length(sim$regdata), ncoef, 20L))
  expect_equal(dim(fit$taudraw),  c(20L, length(sim$regdata)))
  expect_equal(fit$bart_models$treedraws$leaf_dim, ncoef)
  expect_length(fit$var_models,  ncoef)

  expect_s3_class(fit$mu_draw, "mcmc")
  expect_s3_class(fit$taudraw, "mcmc")
})

test_that("vartree only with ncoef == 1 stays diagonal (scalar case)", {
  # const = FALSE + nvar = 1  =>  ncoef = 1.  No off-diagonals to model.
  sim <- sim_hier_linear(
    nreg = 25L, nobs = 15L, nvar = 1L, nz = 2L,
    const = FALSE, het_observed = "linear",
    target_var_betabar = 1.0, target_var_eps = 0.5,
    sigma_sq = 0.5, seed = 20260513L)
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierLinearMixture(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierLinearMixtureHeterCov")
  expect_null(fit$phi_models)
  expect_equal(fit$bart_models$treedraws$leaf_dim, 1L)
  expect_length(fit$var_models,  1L)
  expect_equal(dim(fit$mu_draw), c(20L, 1L))
})

# --- vartree + phitree ------------------------------------------------------
test_that("vartree + phitree returns jagged phi_models", {
  # nvar = 2: one off-diagonal tree phi_21 per coefficient hierarchy.
  sim   <- make_small_linear_sim(nvar = 2L, nz = 2L,
                                 nreg = 25L, nobs = 15L,
                                 seed = 20260514L)
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L),
                phitree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierLinearMixture(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  expect_s3_class(fit, "bayesm.HART.HeterCov")
  ncoef <- ncol(sim$regdata[[1]]$X)
  expect_length(fit$phi_models, ncoef)

  for (j in seq_len(ncoef)) {
    if (j == 1L) {
      expect_true(is.null(fit$phi_models[[1]]) ||
                  length(fit$phi_models[[1]]) == 0L)
    } else {
      expect_length(fit$phi_models[[j]], j - 1L)
      for (k in seq_len(j - 1L)) {
        expect_true(!is.null(fit$phi_models[[j]][[k]]$treedraws))
      }
    }
  }
})

# --- predict.* dispatch -------------------------------------------------------
test_that("predict.rhierLinearMixture dispatches via marker class for heter-cov", {
  sim   <- make_small_linear_sim(nvar = 2L, nz = 2L,
                                 nreg = 25L, nobs = 15L,
                                 seed = 20260515L)
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L),
                phitree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierLinearMixture(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  ncoef <- ncol(sim$regdata[[1]]$X)
  npred <- 4L
  Zstar <- matrix(runif(npred * ncol(sim$Z), -1, 1), npred, ncol(sim$Z))
  newdata <- list(Z = Zstar)

  pred_dz <- predict(fit, newdata = newdata, type = "DeltaZ", burn = 5L,
                     r_verbose = FALSE)
  expect_equal(dim(pred_dz), c(npred, ncoef, 15L))
  expect_true(all(is.finite(pred_dz)))

  pred_dzmu <- predict(fit, newdata = newdata, type = "DeltaZ+mu", burn = 5L,
                       r_verbose = FALSE)
  expect_equal(dim(pred_dzmu), dim(pred_dz))

  diff_arr <- pred_dzmu - pred_dz
  mu_kept  <- t(fit$mu_draw[6:20, , drop = FALSE])
  for (i in seq_len(npred)) {
    expect_equal(diff_arr[i, , ], mu_kept, tolerance = 1e-8,
                 ignore_attr = TRUE)
  }

  pred_sigma <- predict(fit, newdata = newdata, type = "SigmaZ", burn = 5L,
                        r_verbose = FALSE)
  expect_equal(dim(pred_sigma), c(npred, ncoef, ncoef, 15L))
  expect_true(all(is.finite(pred_sigma)))
  for (i in seq_len(npred)) {
    for (s in 1:15) {
      expect_equal(pred_sigma[i, , , s], t(pred_sigma[i, , , s]),
                   tolerance = 1e-8, ignore_attr = TRUE)
      expect_true(all(diag(pred_sigma[i, , , s]) > 0))
    }
  }
})

# --- caches -------------------------------------------------------------------
test_that("DeltaZ and SigmaZ caches match tree evaluation", {
  sim   <- make_small_linear_sim(nvar = 2L, nz = 2L, nreg = 25L, nobs = 15L,
                                 seed = 20260516L)
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 10L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierLinearMixture(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })
  newdata <- list(Z = sim$Z)
  dz <- predict(fit, newdata = newdata, type = "DeltaZ", burn = 0L,
                r_verbose = FALSE)
  expect_equal(dz, bayesm.HART:::.pw_vec_trees(sim$Z, fit$bart_models$treedraws,
                                           dim(fit$betadraw)[2]),
               tolerance = 1e-7, ignore_attr = TRUE)
  sz <- predict(fit, newdata = newdata, type = "SigmaZ", burn = 0L,
                r_verbose = FALSE)
  sz_trees <- bayesm.HART:::.calculate_sigma_z(
    fit, newdata, burn = 0, r_verbose = FALSE, force_tree_eval = TRUE)
  expect_equal(sz, sz_trees, tolerance = 1e-7, ignore_attr = TRUE)
})

# --- verbose banner -----------------------------------------------------------
test_that("verbose output lists the covariance-tree priors", {
  sim   <- make_small_linear_sim()
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 6L, keep = 1L, nprint = 0L)

  out <- capture.output({
    fit <- rhierLinearMixture(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = TRUE)
  })

  expect_s3_class(fit, "rhierLinearMixtureHeterCov")
  expect_true(any(grepl("Covariance trees", out)))
  expect_true(any(grepl("vartree: num_trees = 5", out)))
})

# --- OLS warm-start verification ---------------------------------------------
test_that("Beta_init warm-start gives finite betadraws even with tiny R", {
  sim   <- make_small_linear_sim()
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 2L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierLinearMixture(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })
  expect_true(all(is.finite(fit$betadraw)))
  expect_true(all(is.finite(fit$mu_draw)))
  expect_true(all(fit$taudraw > 0))
})
