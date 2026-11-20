# Covariate-dependent heterogeneity in the simulators (Sigma_func).

library(testthat)

sig2 <- function(Zi) {
  v <- c(0.5 * exp(Zi[1]), 0.5)
  r <- 0.6 * tanh(2 * Zi[2]) * sqrt(v[1] * v[2])
  matrix(c(v[1], r, r, v[2]), 2, 2)
}#SIG2

test_that("Sigma_func drives eps_i in all three simulators", {
  lin <- sim_hier_linear(nreg = 2000, nobs = 2, nvar = 1, nz = 2,
                         het_observed = "linear", Sigma_func = sig2, seed = 1)
  expect_equal(dim(lin$true_values$Sigma_true), c(2000L, 2L, 2L))
  expect_equal(lin$true_values$Sigma_true[7, , ], sig2(lin$Z[7, ]))
  e <- lin$true_values$eps_true
  S <- lin$true_values$Sigma_true
  expect_equal(colMeans(e^2), colMeans(cbind(S[, 1, 1], S[, 2, 2])),
               tolerance = 0.1)

  nb <- sim_hier_negbin(nreg = 300, nobs = 3, nvar = 1, nz = 2,
                        het_observed = "linear", Sigma_func = sig2, seed = 2)
  expect_equal(dim(nb$true_values$Sigma_true), c(300L, 2L, 2L))

  mnl <- sim_hier_mnl(nlgt = 300, nT = 2, p = 3, nz = 2, nXa = 2, nXd = 0,
                      const = FALSE, beta_func_type = "linear",
                      Sigma_func = sig2, seed = 3)
  expect_equal(dim(mnl$true_values$Sigma_true), c(300L, 2L, 2L))
  expect_equal(mnl$true_values$Sigma_true[5, , ], sig2(mnl$Z[5, ]))
})#TEST_THAT

test_that("Sigma_func requires Z and one component", {
  expect_error(sim_hier_linear(nreg = 10, nobs = 2, nvar = 1, nz = 0,
                               het_observed = "none", Sigma_func = sig2),
               regexp = "requires Z")
  expect_error(sim_hier_mnl(nlgt = 10, nT = 2, p = 3, nz = 2, nXa = 2,
                            nXd = 0, const = FALSE, ncomp = 2,
                            Sigma_func = sig2),
               regexp = "ncomp = 1")
  expect_error(sim_hier_linear(nreg = 10, nobs = 2, nvar = 1, nz = 2,
                               het_observed = "linear",
                               Sigma_func = function(Zi) diag(3)),
               regexp = "symmetric ncoef x ncoef")
})#TEST_THAT
