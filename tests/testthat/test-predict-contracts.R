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

test_that("MNL predictive dispatcher matches predict() prior/posterior modes", {
  set.seed(20260515)
  sim <- sim_hier_mnl(
    nlgt = 20, nT = 4, p = 3, nz = 2, nXa = 1, nXd = 0, const = TRUE,
    beta_func_type = "linear"
  )
  prior <- list(ncomp = 1L, bart = list(num_trees = 5L))
  mcmc <- list(R = 12L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierMnlRwMixture(Data = sim, Prior = prior, Mcmc = mcmc, r_verbose = FALSE)
  })

  newdata_post <- list(p = sim$p, nlgtdata = sim$lgtdata)
  out_pred <- predict(fit, newdata = newdata_post, mode = "posterior",
                      type = "choice_probs", burn = 1L, r_verbose = FALSE)
  out_dispatch <- bayesm.HART:::.predictive_dispatch_mnl(
    fit, newdata_post, mode = "posterior", type = "choice_probs",
    burn = 1L, nsim = 2L, r_verbose = FALSE
  )
  expect_equal(out_pred, out_dispatch, tolerance = 1e-10, ignore_attr = TRUE)

  X_template <- sim$lgtdata[[1]]$X
  newdata_prior <- list(
    Z = sim$Z, p = sim$p,
    X = replicate(nrow(sim$Z), X_template, simplify = FALSE)
  )
  set.seed(20260516)
  out_pred_prior <- predict(fit, newdata = newdata_prior, mode = "prior",
                            type = "choice_probs", burn = 1L, nsim = 2L, r_verbose = FALSE)
  set.seed(20260516)
  out_dispatch_prior <- bayesm.HART:::.predictive_dispatch_mnl(
    fit, newdata_prior, mode = "prior", type = "choice_probs",
    burn = 1L, nsim = 2L, r_verbose = FALSE
  )
  expect_equal(out_pred_prior, out_dispatch_prior, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("MNL unified mode/type contract works via predict()", {
  set.seed(20260525)
  sim <- sim_hier_mnl(
    nlgt = 16, nT = 3, p = 3, nz = 2, nXa = 1, nXd = 0, const = TRUE,
    beta_func_type = "linear"
  )
  prior <- list(ncomp = 1L, bart = list(num_trees = 5L))
  mcmc <- list(R = 12L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierMnlRwMixture(Data = sim, Prior = prior, Mcmc = mcmc, r_verbose = FALSE)
  })

  newdata_post <- list(p = sim$p, nlgtdata = sim$lgtdata)
  out_mode <- predict(fit, newdata = newdata_post, mode = "posterior",
                      type = "choice_probs", burn = 1L, r_verbose = FALSE)
  out_direct <- bayesm.HART:::.predictive_dispatch(
    fit, newdata_post, mode = "posterior", type = "choice_probs",
    burn = 1L, nsim = 2L, r_verbose = FALSE
  )
  expect_equal(out_mode, out_direct, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("Linear predictive dispatcher matches predict() in posterior mode", {
  set.seed(20260517)
  sim <- sim_hier_linear(
    nreg = 12, nobs = 6, nvar = 2, nz = 2,
    const = TRUE, het_observed = "linear"
  )
  prior <- list(ncomp = 1L)
  mcmc <- list(R = 14L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierLinearMixture(
      Data = list(regdata = sim$regdata, Z = sim$Z),
      Prior = prior, Mcmc = mcmc, r_verbose = FALSE
    )
  })

  newdata_post <- list(regdata = lapply(sim$regdata, function(u) list(X = u$X)))
  set.seed(20260518)
  out_pred <- predict(fit, newdata = newdata_post, mode = "posterior",
                      burn = 2L, r_verbose = FALSE)
  set.seed(20260518)
  out_dispatch <- bayesm.HART:::.predictive_dispatch(
    fit, newdata_post, mode = "posterior", type = "response",
    burn = 2L, nsim = 1L, r_verbose = FALSE
  )
  expect_equal(out_pred, out_dispatch, tolerance = 1e-10, ignore_attr = TRUE)
  expect_length(out_dispatch, length(sim$regdata))
  expect_equal(ncol(out_dispatch[[1]]), dim(fit$betadraw)[3] - 2L)
})

test_that("Linear prior predictive response path is reproducible and shaped", {
  set.seed(20260519)
  sim <- sim_hier_linear(
    nreg = 10, nobs = 5, nvar = 1, nz = 2,
    const = TRUE, het_observed = "linear"
  )
  prior <- list(ncomp = 1L)
  mcmc <- list(R = 12L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierLinearMixture(
      Data = list(regdata = sim$regdata, Z = sim$Z),
      Prior = prior, Mcmc = mcmc, r_verbose = FALSE
    )
  })

  npred <- 4L
  new_Z <- sim$Z[seq_len(npred), , drop = FALSE]
  X_template <- sim$regdata[[1]]$X
  newdata_prior <- list(
    Z = new_Z,
    X = replicate(npred, X_template, simplify = FALSE)
  )

  set.seed(20260520)
  out1 <- predict(fit, newdata = newdata_prior, mode = "prior",
                  burn = 1L, nsim = 3L, r_verbose = FALSE)
  set.seed(20260520)
  out2 <- bayesm.HART:::.predictive_dispatch(
    fit, newdata_prior, mode = "prior", type = "response",
    burn = 1L, nsim = 3L, r_verbose = FALSE
  )
  expect_equal(out1, out2, tolerance = 1e-10, ignore_attr = TRUE)
  expect_length(out1, npred)
  expect_equal(dim(out1[[1]]), c(nrow(X_template), dim(fit$betadraw)[3] - 1L))
})

test_that("NegBin predictive dispatcher matches predict() in posterior mode", {
  set.seed(20260521)
  sim <- sim_hier_negbin(
    nreg = 10, nobs = 8, nvar = 1, nz = 2,
    const = TRUE, het_observed = "linear"
  )
  prior <- list(ncomp = 1L)
  mcmc <- list(R = 14L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierNegbinRw(
      Data = list(regdata = sim$regdata, hessdata = sim$hessdata, Z = sim$Z),
      Prior = prior, Mcmc = mcmc, r_verbose = FALSE
    )
  })

  newdata_post <- list(regdata = lapply(sim$regdata, function(u) list(X = u$X)))
  set.seed(20260522)
  out_pred <- predict(fit, newdata = newdata_post, mode = "posterior",
                      burn = 2L, r_verbose = FALSE)
  set.seed(20260522)
  out_dispatch <- bayesm.HART:::.predictive_dispatch(
    fit, newdata_post, mode = "posterior", type = "response",
    burn = 2L, nsim = 1L, r_verbose = FALSE
  )
  expect_equal(out_pred, out_dispatch, tolerance = 1e-10, ignore_attr = TRUE)
  expect_length(out_dispatch, length(sim$regdata))
  expect_equal(ncol(out_dispatch[[1]]), dim(fit$betadraw)[3] - 2L)
})

test_that("NegBin prior predictive response path is reproducible and shaped", {
  set.seed(20260523)
  sim <- sim_hier_negbin(
    nreg = 8, nobs = 6, nvar = 1, nz = 2,
    const = TRUE, het_observed = "linear"
  )
  prior <- list(ncomp = 1L)
  mcmc <- list(R = 12L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierNegbinRw(
      Data = list(regdata = sim$regdata, hessdata = sim$hessdata, Z = sim$Z),
      Prior = prior, Mcmc = mcmc, r_verbose = FALSE
    )
  })

  npred <- 4L
  new_Z <- sim$Z[seq_len(npred), , drop = FALSE]
  X_template <- sim$regdata[[1]]$X
  newdata_prior <- list(
    Z = new_Z,
    X = replicate(npred, X_template, simplify = FALSE)
  )

  set.seed(20260524)
  out1 <- predict(fit, newdata = newdata_prior, mode = "prior",
                  burn = 1L, nsim = 3L, r_verbose = FALSE)
  set.seed(20260524)
  out2 <- bayesm.HART:::.predictive_dispatch(
    fit, newdata_prior, mode = "prior", type = "response",
    burn = 1L, nsim = 3L, r_verbose = FALSE
  )
  expect_equal(out1, out2, tolerance = 1e-10, ignore_attr = TRUE)
  expect_length(out1, npred)
  expect_equal(dim(out1[[1]]), c(nrow(X_template), dim(fit$betadraw)[3] - 1L))
})
