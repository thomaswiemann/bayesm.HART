# Tests for the heteroscedastic-covariance Negbin path
# (Prior$vartree, optionally Prior$phitree).
#
# Smoke tests focused on:
#   * input validation (vartree without bart, ncomp != 1, no Z),
#   * end-to-end execution via the user wrapper,
#   * class assignment ("bayesm.HART.HeterCov" marker) and return-slot shape,
#   * predict.* dispatch on the marker class for "DeltaZ" and "DeltaZ+mu",
#   * lambda auto-calibration message under r_verbose = TRUE.
#
# Recovery accuracy is exercised by dev/smoke-heter-negbin.R.

library(testthat)

# --- shared minimal-cost simulation -----------------------------------------
make_small_negbin_sim <- function(seed = 20260513L,
                                  nreg = 25L, nobs = 15L,
                                  nvar = 1L, nz = 2L) {
  sim_hier_negbin(
    nreg = nreg, nobs = nobs, nvar = nvar, nz = nz,
    const = TRUE, het_observed = "linear",
    target_var_betabar = 0.5, target_var_eps = 0.25,
    alpha = 5.0, seed = seed
  )
}

# --- input validation --------------------------------------------------------
test_that("vartree requires bart, ncomp == 1, and Z (drawdelta = TRUE)", {
  sim <- make_small_negbin_sim()

  # vartree without bart -> error
  expect_error(
    capture.output(
      rhierNegbinRw(
        Data  = list(regdata = sim$regdata, Z = sim$Z),
        Prior = list(ncomp = 1L, vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "requires Prior\\$bart"
  )

  # ncomp > 1 with vartree -> error
  expect_error(
    capture.output(
      rhierNegbinRw(
        Data  = list(regdata = sim$regdata, Z = sim$Z),
        Prior = list(ncomp = 2L,
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
      rhierNegbinRw(
        Data  = list(regdata = sim$regdata),       # no Z
        Prior = list(ncomp = 1L,
                     bart    = list(num_trees = 5L),
                     vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "requires Z"
  )
})

# --- vartree only with nvar > 1 (auto-promotes to full Cholesky) ------------
test_that("vartree only with nvar > 1 auto-promotes to full Cholesky", {
  # Default sim has nvar = 1 + intercept = ncoef = 2, so auto-promotion fires.
  sim   <- make_small_negbin_sim()
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))   # NO phitree
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierNegbinRw(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierNegbinRwHeterCov")
  expect_s3_class(fit, "bayesm.HART.HeterCov")
  expect_s3_class(fit, "rhierNegbinRw")

  expected_slots <- c("bart_models", "var_models", "phi_models", "mu_draw",
                      "varcount", "varprob", "var_varcount", "var_varprob",
                      "betadraw", "alphadraw", "loglike",
                      "acceptrbeta", "acceptralpha")
  expect_named(fit, expected_slots, ignore.order = TRUE)
  expect_null(fit$nmix)
  expect_null(fit$Deltadraw)

  ncoef <- ncol(sim$regdata[[1]]$X)
  # Auto-promotion: phi_models is non-null and jagged (length-1 at j=2 here).
  expect_length(fit$phi_models, ncoef)
  expect_length(fit$phi_models[[2]], 1L)

  expect_equal(dim(fit$mu_draw),  c(20L, ncoef))
  expect_equal(dim(fit$betadraw), c(length(sim$regdata), ncoef, 20L))
  expect_length(fit$alphadraw, 20L)
  expect_length(fit$bart_models, ncoef)
  expect_length(fit$var_models,  ncoef)
  expect_s3_class(fit$mu_draw,    "mcmc")
  expect_s3_class(fit$alphadraw,  "mcmc")
})

test_that("vartree only with ncoef == 1 stays diagonal (scalar case)", {
  # const = FALSE + nvar = 1  =>  ncoef = 1.  No off-diagonals to model.
  sim <- sim_hier_negbin(
    nreg = 25L, nobs = 15L, nvar = 1L, nz = 2L,
    const = FALSE, het_observed = "linear",
    target_var_betabar = 0.5, target_var_eps = 0.25,
    alpha = 5.0, seed = 20260513L)
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierNegbinRw(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierNegbinRwHeterCov")
  expect_null(fit$phi_models)
  expect_length(fit$bart_models, 1L)
  expect_length(fit$var_models,  1L)
  expect_equal(dim(fit$mu_draw), c(20L, 1L))
})

# --- Phase 2: vartree + phitree ---------------------------------------------
test_that("Phase 2 (vartree + phitree) returns jagged phi_models", {
  # nvar = 2 so we have one off-diagonal phi_{2,1}(.) tree
  sim   <- make_small_negbin_sim(nvar = 2L, nz = 2L,
                                 nreg = 25L, nobs = 15L,
                                 seed = 20260514L)
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L),
                phitree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierNegbinRw(
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
test_that("predict.rhierNegbinRw dispatches via marker class for heter-cov", {
  sim   <- make_small_negbin_sim(nvar = 2L, nz = 2L,
                                 nreg = 25L, nobs = 15L,
                                 seed = 20260515L)
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L),
                phitree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierNegbinRw(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
  })

  ncoef <- ncol(sim$regdata[[1]]$X)
  npred <- 4L
  Zstar <- matrix(runif(npred * ncol(sim$Z), -1, 1), npred, ncol(sim$Z))
  newdata <- list(Z = Zstar)

  # DeltaZ: heter-cov path (rebuilds Sigma(Z*)^{1/2} per draw)
  pred_dz <- predict(fit, newdata = newdata, type = "DeltaZ", burn = 5L,
                     r_verbose = FALSE)
  expect_equal(dim(pred_dz), c(npred, ncoef, 15L))
  expect_true(all(is.finite(pred_dz)))

  # DeltaZ+mu adds the per-draw mu (from object$mu_draw under marker class)
  pred_dzmu <- predict(fit, newdata = newdata, type = "DeltaZ+mu", burn = 5L,
                       r_verbose = FALSE)
  expect_equal(dim(pred_dzmu), dim(pred_dz))

  # Difference equals mu_kept broadcast across npred and draws
  diff_arr <- pred_dzmu - pred_dz
  mu_kept  <- t(fit$mu_draw[6:20, , drop = FALSE])  # ncoef x 15
  for (i in seq_len(npred)) {
    expect_equal(diff_arr[i, , ], mu_kept, tolerance = 1e-8,
                 ignore_attr = TRUE)
  }

  # SigmaZ: covariance draws [npred, ncoef, ncoef, ndraws_kept]
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

# --- lambda auto-calibration --------------------------------------------------
test_that("auto-calibrated lambda is announced when Prior$vartree$lambda is NULL", {
  sim   <- make_small_negbin_sim()
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))   # no lambda -> auto
  Mcmc  <- list(R = 6L, keep = 1L, nprint = 0L)

  out <- capture.output({
    fit <- rhierNegbinRw(
      Data  = list(regdata = sim$regdata, Z = sim$Z),
      Prior = Prior, Mcmc = Mcmc, r_verbose = TRUE)
  })

  expect_s3_class(fit, "rhierNegbinRwHeterCov")
  expect_true(any(grepl("auto-calibrated", out)),
              info = "expected lambda auto-calibration line in verbose output")
  expect_true(any(grepl("Heteroscedastic Sigma", out)),
              info = "expected heter-cov banner in verbose output")
})
