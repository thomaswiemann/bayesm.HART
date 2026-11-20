library(testthat)
library(bayesm.HART)

# ==============================================================================
# Test sim_hier_mnl Function
# ==============================================================================

# ==============================================================================
# Basic Linear, ncomp = 1
# ==============================================================================
test_that("Basic linear simulation works correctly", {

  nlgt_test <- 10
  nT_test <- 5
  p_test <- 3
  nz_test <- 2
  nXa_test <- 1
  nXd_test <- 0
  const_test <- TRUE
  seed_test <- 42

  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

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
    standardize_Z = FALSE
  )

  expect_true(is.list(sim_data))
  expect_named(sim_data, c("p", "lgtdata", "Z", "true_values"))
  expect_equal(sim_data$p, p_test)

  expect_true(is.matrix(sim_data$Z))
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_equal(ncol(sim_data$Z), nz_test)

  expect_true(is.list(sim_data$lgtdata))
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_named(sim_data$lgtdata[[1]], c("y", "X", "beta", "betabar"))
  expect_length(sim_data$lgtdata[[1]]$y, nT_test)
  expect_true(is.matrix(sim_data$lgtdata[[1]]$X))
  # Rows = nT * p (alternatives stacked)
  expect_equal(nrow(sim_data$lgtdata[[1]]$X), nT_test * p_test)
  # Cols = ncoef
  expect_equal(ncol(sim_data$lgtdata[[1]]$X), ncoef_test)
  expect_length(sim_data$lgtdata[[1]]$beta, ncoef_test)
  expect_length(sim_data$lgtdata[[1]]$betabar, ncoef_test)

  expect_true(is.list(sim_data$true_values))
  expect_named(sim_data$true_values,
               c("beta_true", "betabar_true", "true_params", "dimensions"))

  expect_true(is.matrix(sim_data$true_values$beta_true))
  expect_equal(dim(sim_data$true_values$beta_true), c(nlgt_test, ncoef_test))

  expect_true(is.matrix(sim_data$true_values$betabar_true))
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  true_params <- sim_data$true_values$true_params
  dims <- sim_data$true_values$dimensions

  expect_equal(true_params$beta_func_type, "linear")
  expect_named(true_params$beta_func_args, "Delta")
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

  # ncomp = 1: eps ~ N(mu, Sigma) with Sigma = (rooti' rooti)^-1.
  mu1 <- true_params$mixture_comps[[1]]$mu
  # beta - betabar is centered at mu (zero under the default sigma_inv_diag).
  expect_equal(mu1, rep(0, ncoef_test))
  eps_est <- sim_data$true_values$beta_true - sim_data$true_values$betabar_true
  # Loose tolerance: nlgt = 10.
  expect_true(all(abs(colMeans(eps_est) - mu1) < 0.8))

})

# ==============================================================================
# Step Function
# ==============================================================================
test_that("Step function simulation works correctly", {

  nlgt_test <- 15
  nT_test <- 2
  p_test <- 2
  nz_test <- 3
  nXa_test <- 1
  nXd_test <- 0
  const_test <- FALSE
  seed_test <- 123
  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  step_args <- list(
    cutoff = 0.1,
    beta_1 = rep(2, ncoef_test),
    beta_2 = rep(-2, ncoef_test),
    Z_index = 1
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
    standardize_Z = FALSE
  )

  expect_true(is.list(sim_data))
  expect_equal(sim_data$p, p_test)
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_equal(ncol(sim_data$Z), nz_test)
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  true_params <- sim_data$true_values$true_params
  expect_equal(true_params$beta_func_type, "step")
  expect_equal(true_params$beta_func_args, step_args)

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

})

# ==============================================================================
# Friedman Function
# ==============================================================================
test_that("Friedman function simulation works correctly", {

  nlgt_test <- 15
  nT_test <- 3
  p_test <- 3
  nz_test <- 5
  nXa_test <- 2
  nXd_test <- 1
  const_test <- TRUE
  seed_test <- 444

  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  friedman_args <- list(
    coef_index = 3,
    coefs = c(10, 20, 10, 5),
    shift = 0,
    scale = 8
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
    standardize_Z = TRUE
  )

  expect_true(is.list(sim_data))
  expect_equal(sim_data$p, p_test)
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_equal(ncol(sim_data$Z), nz_test)
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  true_params <- sim_data$true_values$true_params
  expect_equal(true_params$beta_func_type, "friedman")
  expect_equal(true_params$beta_func_args, friedman_args)

  Z_mat <- sim_data$Z
  betabar <- sim_data$true_values$betabar_true
  coef_idx <- friedman_args$coef_index

  expected_friedman_vals <- apply(Z_mat, 1, .beta_Z_friedman)
  expect_equal(betabar[, coef_idx], expected_friedman_vals)

  other_indices <- setdiff(1:ncoef_test, coef_idx)
  if (length(other_indices) > 0) {
      expect_true(all(betabar[, other_indices] == 0))
  }

})

# ==============================================================================
# Mixture Components (ncomp > 1)
# ==============================================================================
test_that("Mixture component simulation (ncomp > 1) works", {

  nlgt_test <- 20
  nT_test <- 4
  p_test <- 3
  nz_test <- 2
  nXa_test <- 1
  nXd_test <- 1
  const_test <- TRUE
  seed_test <- 555
  ncomp_test <- 3

  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  sim_data <- sim_hier_mnl(
    nlgt = nlgt_test,
    nT = nT_test,
    p = p_test,
    nz = nz_test,
    nXa = nXa_test,
    nXd = nXd_test,
    const = const_test,
    beta_func_type = "linear",
    ncomp = ncomp_test,
    seed = seed_test,
    standardize_Z = TRUE
  )

  expect_true(is.list(sim_data))
  expect_equal(sim_data$p, p_test)
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_equal(dim(sim_data$true_values$beta_true), c(nlgt_test, ncoef_test))
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  true_params <- sim_data$true_values$true_params
  dims <- sim_data$true_values$dimensions

  expect_equal(dims$ncomp, ncomp_test)
  expect_length(true_params$mixture_comps, ncomp_test)
  expect_named(true_params$mixture_comps[[1]], c("mu", "rooti"))
  expect_length(true_params$mixture_comps[[1]]$mu, ncoef_test)
  expect_equal(dim(true_params$mixture_comps[[1]]$rooti), c(ncoef_test, ncoef_test))
  expect_length(true_params$pvec, ncomp_test)
  expect_equal(sum(true_params$pvec), 1.0)

  eps_est <- sim_data$true_values$beta_true - sim_data$true_values$betabar_true
  expect_equal(dim(eps_est), c(nlgt_test, ncoef_test))
  expect_false(all(abs(eps_est) < 1e-9))

})

# ==============================================================================
# Custom Beta Function
# ==============================================================================
test_that("Custom beta function simulation works", {

  nlgt_test <- 10
  nT_test <- 2
  p_test <- 2
  nz_test <- 4
  nXa_test <- 2
  nXd_test <- 0
  const_test <- FALSE
  seed_test <- 789
  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  custom_beta_func <- function(Zi) {
    if(length(Zi) != nz_test) stop("Incorrect Zi length passed to custom func")
    out_vec <- c(sum(Zi), prod(Zi[1:min(3, length(Zi))]))
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
    standardize_Z = FALSE
  )

  expect_true(is.list(sim_data))
  expect_equal(sim_data$p, p_test)
  expect_equal(nrow(sim_data$Z), nlgt_test)
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))

  true_params <- sim_data$true_values$true_params
  expect_equal(true_params$beta_func_type, "custom")
  # Functions are not compared; only presence.
  expect_true(is.function(true_params$beta_func_args$func))

  Z_mat <- sim_data$Z
  betabar <- sim_data$true_values$betabar_true

  expected_betabar <- t(apply(Z_mat, 1, custom_beta_func))

  expect_equal(betabar, expected_betabar)

})

# ==============================================================================
# No Z variables (nz = 0)
# ==============================================================================
test_that("Simulation works with nz = 0", {

  nlgt_test <- 8
  nT_test <- 3
  p_test <- 4
  nz_test <- 0
  nXa_test <- 2
  nXd_test <- 1
  const_test <- TRUE
  seed_test <- 111
  ncomp_test <- 2

  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test

  sim_data <- sim_hier_mnl(
    nlgt = nlgt_test,
    nT = nT_test,
    p = p_test,
    nz = nz_test,
    nXa = nXa_test,
    nXd = nXd_test,
    const = const_test,
    beta_func_type = "linear",
    ncomp = ncomp_test,
    seed = seed_test
  )

  expect_true(is.list(sim_data))
  expect_named(sim_data, c("p", "lgtdata", "Z", "true_values"))
  expect_equal(sim_data$p, p_test)

  expect_null(sim_data$Z)

  expect_true(is.list(sim_data$lgtdata))
  expect_length(sim_data$lgtdata, nlgt_test)
  expect_named(sim_data$lgtdata[[1]], c("y", "X", "beta", "betabar"))

  expect_true(is.list(sim_data$true_values))
  expect_named(sim_data$true_values,
               c("beta_true", "betabar_true", "true_params", "dimensions"))

  expect_true(is.matrix(sim_data$true_values$betabar_true))
  expect_equal(dim(sim_data$true_values$betabar_true), c(nlgt_test, ncoef_test))
  expect_true(all(sim_data$true_values$betabar_true == 0))

  expect_true(is.matrix(sim_data$true_values$beta_true))
  expect_equal(dim(sim_data$true_values$beta_true), c(nlgt_test, ncoef_test))
  expect_false(all(abs(sim_data$true_values$beta_true) < 1e-9))

  eps_est <- sim_data$true_values$beta_true - sim_data$true_values$betabar_true
  expect_equal(sim_data$true_values$beta_true, eps_est)

  true_params <- sim_data$true_values$true_params
  dims <- sim_data$true_values$dimensions

  expect_equal(true_params$beta_func_type, "linear")
  expect_equal(true_params$beta_func_args, list())
  expect_length(true_params$mixture_comps, ncomp_test)
  expect_length(true_params$pvec, ncomp_test)

  expect_equal(dims$nz, 0)
  expect_equal(dims$ncomp, ncomp_test)

})

# ==============================================================================
# User-provided mixture components
# ==============================================================================
test_that("User-provided mixture components work", {

  nlgt_test <- 5
  nT_test <- 2
  p_test <- 3 
  nz_test <- 1
  nXa_test <- 1 
  nXd_test <- 0
  const_test <- TRUE
  seed_test <- 666
  ncomp_test <- 2 
  ncoef_test <- const_test*(p_test - 1) + (p_test - 1)*nXd_test + nXa_test
  expect_equal(ncoef_test, 3)

  user_comps <- list(
    list(mu = c(1, 0, 0), rooti = diag(1.0, ncoef_test)),
    list(mu = c(-1, 1, 0), rooti = diag(2.0, ncoef_test))
  )
  
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
        mixture_comps = user_comps,
        seed = seed_test
      )
  })
  
  expect_true(!is.null(sim_data), "Simulation failed to run with user components.")
  
  true_params <- sim_data$true_values$true_params
  expect_equal(true_params$mixture_comps, user_comps)
  expect_equal(sim_data$true_values$dimensions$ncomp, ncomp_test)

})
