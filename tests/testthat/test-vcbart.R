# Tests for the varying-coefficient tree machinery (src/vcbartfuns.cpp,
# src/vcbd.cpp, src/vcbart.cpp): y_i = a_i f(Z_i) + e_i, e_i ~ N(0, var_i),
# leaf prior N(0, tau^2).  The trees are scalar, so .pw_trees() reads them.

library(testthat)

# Exact log marginal of the units in one leaf: y ~ N(0, tau^2 a a' + diag(var)).
log_marginal_vc <- function(y, a, v, tau) {
  V <- tau^2 * tcrossprod(a) + diag(v, length(v))
  -0.5 * length(y) * log(2 * pi) -
    0.5 * as.numeric(determinant(V, logarithm = TRUE)$modulus) -
    0.5 * drop(t(y) %*% solve(V, y))
}#LOG_MARGINAL_VC

test_that("vclh matches the exact leaf marginal up to unit constants", {
  set.seed(1)
  n <- 5; a <- rnorm(n); v <- runif(n, 0.5, 2); y <- rnorm(n); tau <- 0.7
  P <- sum(a^2 / v); b <- sum(a * y / v)
  # Constants dropped by the score: per-unit normalising terms and y_i^2 / v_i.
  unit_const <- sum(-0.5 * log(2 * pi * v) - 0.5 * y^2 / v)
  expect_equal(vclh_R(P, b, tau), log_marginal_vc(y, a, v, tau) - unit_const,
               tolerance = 1e-10)
  # A leaf with no information (all regressors zero) has marginal one.
  expect_equal(vclh_R(0, 0, tau), 0)
  # Weighted-mean kernel at w = a^2 / v, ytilde = y / a.
  w <- a^2 / v
  expect_equal(c(P, b), c(sum(w), sum(w * (y / a))))
})#TEST_THAT

test_that("stump leaves are drawn from the leaf conditional", {
  set.seed(2)
  n <- 30; Zt <- matrix(runif(n, -1, 1), 1)
  a <- rnorm(n); v <- runif(n, 0.5, 2); tau <- 0.8
  y <- a * 0.4 + rnorm(n, sd = sqrt(v))
  # alpha = 0 removes all prior mass from splits, so tree 1 stays a stump and
  # every unit shares the single leaf value.
  out <- vcbart_chain_R(Zt, y, a, v, m = 1L, tau = tau, alpha = 0, beta = 1,
                        nmin = 1L, numcut = 10L, niter = 20000L)
  expect_true(all(out$keys == ""))
  leaf <- out$fits[1, ]
  expect_true(all(abs(sweep(out$fits, 2, leaf)) < 1e-12))
  prec <- 1 / tau^2 + sum(a^2 / v)
  expect_lt(abs(mean(leaf) - sum(a * y / v) / prec), 4 * sqrt(1 / prec / 20000))
  expect_equal(var(leaf), 1 / prec, tolerance = 0.1)
})#TEST_THAT

test_that("a zero-information leaf draws from the prior", {
  set.seed(3)
  n <- 20; Zt <- matrix(runif(n), 1)
  a <- rep(0, n); v <- rep(1, n); y <- rnorm(n); tau <- 0.6
  out <- vcbart_chain_R(Zt, y, a, v, m = 1L, tau = tau, alpha = 0, beta = 1,
                        nmin = 1L, numcut = 5L, niter = 20000L)
  leaf <- out$fits[1, ]
  expect_true(all(is.finite(leaf)))
  expect_lt(abs(mean(leaf)), 4 * tau / sqrt(20000))
  expect_equal(var(leaf), tau^2, tolerance = 0.1)
})#TEST_THAT

vc_problem <- function(seed = 4L, n = 50L, p = 2L) {
  set.seed(seed)
  Zt  <- matrix(runif(p * n, -1, 1), p, n)
  phi <- ifelse(Zt[1, ] > 0, 0.8, -0.4)
  a   <- rnorm(n); v <- runif(n, 0.3, 1)
  list(Zt = Zt, y = a * phi + rnorm(n, sd = sqrt(v)), a = a, v = v, tau = 0.5)
}#VC_PROBLEM

test_that("vcbart keeps split counts consistent and reads back via trees", {
  pb <- vc_problem()
  out <- with(pb, vcbart_chain_R(Zt, y, a, v, m = 4L, tau = tau, alpha = 0.95,
                                 beta = 1, nmin = 3L, numcut = 15L,
                                 niter = 40L))
  expect_equal(sum(out$nv), out$nleaves[40L] - 4L)
  expect_true(all(is.finite(out$fits)))
  pred <- .pw_trees(t(pb$Zt), out$treedraws, mu = 0, transposed = FALSE,
                    dodraws = TRUE)                   # ndraws x npred
  expect_equal(t(pred), out$fits, tolerance = 1e-7)
})#TEST_THAT

test_that("DART bookkeeping carries over to vcbart", {
  pb <- vc_problem(seed = 5L)
  out <- with(pb, vcbart_chain_R(Zt, y, a, v, m = 3L, tau = tau, alpha = 0.95,
                                 beta = 2, nmin = 5L, numcut = 20L, niter = 40L,
                                 sparse = TRUE, burn = 2L))
  expect_equal(sum(out$pv), 1, tolerance = 1e-10)
  expect_true(all(out$pv > 0))
  expect_equal(sum(out$nv), out$nleaves[40L] - 3L)
})#TEST_THAT
