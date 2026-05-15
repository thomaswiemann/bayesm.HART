library(testthat)

test_that("DeltaZ+mu contract holds exactly for heter-cov helper path", {
  set.seed(20260514)
  npred <- 3L
  nvar <- 2L
  ndraws <- 4L

  delta_z <- array(rnorm(npred * nvar * ndraws), dim = c(npred, nvar, ndraws))
  mu_draw <- matrix(rnorm(ndraws * nvar), nrow = ndraws, ncol = nvar)

  object <- list(mu_draw = mu_draw)
  class(object) <- "bayesm.HART.HeterCov"

  out <- bayesm.HART:::.add_mu_component(delta_z, object, kept_draws_indices = seq_len(ndraws))
  diff <- out - delta_z
  mu_kept <- t(mu_draw)

  for (i in seq_len(npred)) {
    expect_equal(diff[i, , ], mu_kept, tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("SigmaZ helper returns symmetric positive-diagonal covariances", {
  set.seed(20260514)
  npred <- 2L
  nvar <- 3L
  ndraws <- 5L

  object <- list(betadraw = array(0, dim = c(1L, nvar, ndraws)))
  class(object) <- "bayesm.HART.HeterCov"
  newdata <- list(Z = matrix(rnorm(npred * 2L), nrow = npred, ncol = 2L))

  comps <- list(
    d_arr = array(runif(npred * nvar * ndraws, min = 0.2, max = 2.0),
                  dim = c(npred, nvar, ndraws)),
    phi_arr = array(0, dim = c(npred, nvar, nvar, ndraws)),
    use_full = FALSE
  )

  sigma_arr <- bayesm.HART:::.calculate_sigma_z(
    object, newdata, burn = 1L, r_verbose = FALSE, hetercov_comps = comps
  )

  expect_equal(dim(sigma_arr), c(npred, nvar, nvar, ndraws - 1L))
  for (i in seq_len(npred)) {
    for (s in seq_len(ndraws - 1L)) {
      expect_equal(sigma_arr[i, , , s], t(sigma_arr[i, , , s]),
                   tolerance = 1e-12, ignore_attr = TRUE)
      expect_true(all(diag(sigma_arr[i, , , s]) > 0))
    }
  }
})

test_that("single-component guard rejects ncomp > 1 in prediction helpers", {
  object <- list(nmix = list(probdraw = matrix(c(0.5, 0.5), nrow = 1L)))
  expect_error(
    bayesm.HART:::.assert_single_component_nmix(object, "test context"),
    regexp = "supports only ncomp = 1"
  )
})

test_that("posterior_probs helper preserves simplex rows", {
  set.seed(20260514)
  nlgt <- 2L
  p <- 2L
  nvar <- 2L
  ndraws <- 3L
  T_i <- 3L

  object <- list(betadraw = array(rnorm(nlgt * nvar * ndraws), dim = c(nlgt, nvar, ndraws)))
  x1 <- matrix(rnorm(T_i * p * nvar), nrow = T_i * p, ncol = nvar)
  x2 <- matrix(rnorm(T_i * p * nvar), nrow = T_i * p, ncol = nvar)
  newdata <- list(
    p = p,
    nlgtdata = list(list(X = x1), list(X = x2))
  )

  pp <- bayesm.HART:::.predict_posterior_probs(
    object, newdata, kept_draws_indices = seq_len(ndraws), r_verbose = FALSE
  )

  expect_length(pp, nlgt)
  for (i in seq_len(nlgt)) {
    sums <- apply(pp[[i]], c(1, 3), sum)
    expect_equal(as.numeric(sums), rep(1, length(sums)), tolerance = 1e-10)
  }
})

