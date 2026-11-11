# Tests for R/hart_priors.R: the explicit HART model selector, validation,
# hyperparameter defaults, and the configuration object the wrappers consume.

library(testthat)

cfg_of <- function(Prior, nvar = 2L, nz = 2L, drawdelta = TRUE) {
  bayesm.HART:::.hart_config(Prior, nvar = nvar, nz = nz,
                             drawdelta = drawdelta, coef_names = c("a", "b"))
}#CFG_OF

test_that(".hart_config defaults to the bayesm paths", {
  cfg <- cfg_of(list(ncomp = 1))
  expect_equal(cfg$model, "none")
  expect_null(cfg$bart)
  expect_null(cfg$vartree)
  expect_null(cfg$phitree)
  expect_false(cfg$full_cov)
  expect_true(cfg$store_trees)
  # "none" with no Z is the bayesm constant-mean path, not an error.
  expect_equal(cfg_of(list(ncomp = 1), drawdelta = FALSE)$model, "none")
})#TEST_THAT

test_that(".hart_config rejects unknown selectors", {
  expect_error(cfg_of(list(ncomp = 1, hart_model = "hart")),
               regexp = "must be one of")
  expect_error(cfg_of(list(ncomp = 1, hart_model = c("none", "direct"))),
               regexp = "must be one of")
})#TEST_THAT

test_that("tree lists without a HART selector are errors", {
  expect_error(cfg_of(list(ncomp = 1, bart = list(num_trees = 5))),
               regexp = "Prior\\$bart requires Prior\\$hart_model")
  expect_error(cfg_of(list(ncomp = 1, vartree = list())),
               regexp = "require Prior\\$hart_model = 'direct'")
  expect_error(cfg_of(list(ncomp = 1, phitree = list())),
               regexp = "require Prior\\$hart_model = 'direct'")
  # bart = NULL is the same as no bart.
  expect_equal(cfg_of(list(ncomp = 1, bart = NULL))$model, "none")
})#TEST_THAT

test_that("HART models require Z", {
  expect_error(cfg_of(list(ncomp = 1, hart_model = "normalized"), drawdelta = FALSE),
               regexp = "requires Data\\$Z")
  expect_error(cfg_of(list(ncomp = 1, hart_model = "direct"), drawdelta = FALSE),
               regexp = "requires Data\\$Z")
})#TEST_THAT

test_that("normalized model: bart optional, covariance lists rejected", {
  cfg <- cfg_of(list(ncomp = 1, hart_model = "normalized"))
  expect_equal(cfg$model, "normalized")
  expect_equal(cfg$bart$num_trees, 200)
  expect_equal(cfg$bart$tau, 1)
  expect_equal(cfg$bart$nmin, 5L)
  expect_true(cfg$bart$store_trees)
  expect_null(cfg$vartree)
  expect_error(cfg_of(list(ncomp = 1, hart_model = "normalized",
                           vartree = list())),
               regexp = "require Prior\\$hart_model = 'direct'")
  expect_error(cfg_of(list(ncomp = 1, hart_model = "normalized",
                           phitree = list())),
               regexp = "require Prior\\$hart_model = 'direct'")
})#TEST_THAT

test_that("direct model: all tree lists optional, phi iff nvar > 1", {
  cfg <- cfg_of(list(ncomp = 1, hart_model = "direct"))
  expect_equal(cfg$model, "direct")
  expect_true(cfg$full_cov)
  expect_equal(cfg$vartree, list())       # parsed later by the wrapper
  expect_equal(cfg$phitree, list())
  cfg1 <- cfg_of(list(ncomp = 1, hart_model = "direct"), nvar = 1L)
  expect_false(cfg1$full_cov)
  expect_null(cfg1$phitree)
  expect_error(cfg_of(list(ncomp = 1, hart_model = "direct", phitree = list()),
                      nvar = 1L),
               regexp = "requires nvar > 1")
})#TEST_THAT

test_that("user hyperparameters and store_trees are passed through", {
  cfg <- cfg_of(list(ncomp = 1, hart_model = "direct", store_trees = FALSE,
                     bart = list(num_trees = 10, power = 1.5, nmin = 2L),
                     vartree = list(num_trees = 7L)))
  expect_false(cfg$store_trees)
  expect_false(cfg$bart$store_trees)
  expect_equal(cfg$bart$num_trees, 10)
  expect_equal(cfg$bart$power, 1.5)
  expect_equal(cfg$bart$nmin, 2L)
  expect_equal(cfg$bart$tau, 1)
  expect_equal(cfg$vartree, list(num_trees = 7L))
  expect_equal(cfg$coef_names, c("a", "b"))
  expect_equal(cfg$nvar, 2L)
  expect_equal(cfg$nz, 2L)
})#TEST_THAT

test_that(".hart_Amu builds the independent-prior precision", {
  A <- bayesm.HART:::.hart_Amu(NULL, 3L, default = 0.01)
  expect_equal(A, 0.01 * diag(3))
  expect_equal(bayesm.HART:::.hart_Amu(2, 2L, default = 0.01), 2 * diag(2))
  M <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  expect_equal(bayesm.HART:::.hart_Amu(M, 2L, default = 0.01), M)
  expect_error(bayesm.HART:::.hart_Amu(-1, 2L, default = 0.01),
               regexp = "positive scalar")
  expect_error(bayesm.HART:::.hart_Amu(diag(3), 2L, default = 0.01),
               regexp = "nvar x nvar")
  expect_error(bayesm.HART:::.hart_Amu(matrix(c(1, 2, 0, 1), 2, 2), 2L,
                                       default = 0.01),
               regexp = "symmetric positive definite")
  expect_error(bayesm.HART:::.hart_Amu(matrix(c(1, 2, 2, 1), 2, 2), 2L,
                                       default = 0.01),
               regexp = "symmetric positive definite")
})#TEST_THAT

test_that(".parse_var_params and .parse_phi_params supply defaults", {
  oldbetas <- matrix(c(1, 2, 3, 4, 6, 8), ncol = 2)
  v <- bayesm.HART:::.parse_var_params(list(), oldbetas, nz = 3L)
  expect_equal(v$num_trees, 40L)
  expect_equal(v$nu, 10)
  expect_equal(v$lambda, mean(apply(oldbetas, 2, var)))
  expect_equal(v$nmin, 5L)
  expect_equal(v$rho, 3L)
  v_user <- bayesm.HART:::.parse_var_params(list(lambda = 0.3, nmin = 1L),
                                            oldbetas, nz = 3L)
  expect_equal(v_user$lambda, 0.3)
  expect_equal(v_user$nmin, 1L)

  p <- bayesm.HART:::.parse_phi_params(list(), nz = 3L, nu = 5, nvar = 2L)
  expect_equal(p$num_trees, 40L)
  expect_equal(p$nmin, 2L)
  expect_equal(p$ess_min, 5)
  expect_true(p$tau > 0)
})#TEST_THAT

test_that("fitted objects carry hart metadata only for HART models", {
  set.seed(20261110)
  sim <- sim_hier_mnl(nlgt = 12, nT = 4, p = 3, nz = 2, nXa = 1, nXd = 0,
                      const = TRUE, beta_func_type = "linear")
  Mcmc <- list(R = 4L, keep = 1L, nprint = 0L)
  capture.output({
    fit_none <- rhierMnlRwMixture(sim, Prior = list(ncomp = 1), Mcmc,
                                  r_verbose = FALSE)
    fit_norm <- rhierMnlRwMixture(
      sim, Prior = list(ncomp = 1, hart_model = "normalized",
                        bart = list(num_trees = 3L)),
      Mcmc, r_verbose = FALSE)
    fit_dir <- rhierMnlRwMixture(
      sim, Prior = list(ncomp = 1, hart_model = "direct",
                        bart = list(num_trees = 3L),
                        vartree = list(num_trees = 3L)),
      Mcmc, r_verbose = FALSE)
  })
  expect_null(fit_none$hart)
  expect_equal(fit_norm$hart$model, "normalized")
  expect_equal(fit_norm$hart$schema, 1L)
  expect_equal(fit_norm$hart$nvar, dim(fit_norm$betadraw)[2])
  expect_equal(fit_norm$hart$nz, 2L)
  expect_equal(fit_norm$hart$priors$bart$num_trees, 3L)
  expect_null(fit_norm$hart$priors$vartree)
  expect_equal(fit_dir$hart$model, "direct")
  expect_equal(fit_dir$hart$priors$vartree$num_trees, 3L)
  expect_true(is.list(fit_dir$hart$priors$phitree))
  expect_s3_class(fit_dir, "bayesm.HART.HeterCov")
})#TEST_THAT
