# Block conditionals of the direct HART sampler (src/direct_hart_backend.cpp)
# against closed forms.  Every ensemble is one stump (base = 0, one tree), so
# f(Z) = m, v_j(Z) = s_j, phi_jk(Z) = lambda_jk and each block is a conjugate
# draw whose parameters follow from the pinned state.

library(testthat)

make_problem <- function(seed = 1L, n = 40L, D = 3L, nz = 2L) {
  set.seed(seed)
  S <- matrix(c(1, 0.5, 0.2, 0, 1, -0.3, 0, 0, 0.8), D, D)
  list(theta = matrix(rnorm(n * D), n, D) %*% S +
                 matrix(c(0.5, -0.5, 1), n, D, byrow = TRUE),
       Z = matrix(runif(n * nz, -1, 1), n, nz),
       n = n, D = D, nz = nz, mubar = rep(0, D), Amu = 0.1 * diag(D))
}#MAKE_PROBLEM

# `a_v` is the per-leaf shape; with one tree, tau = sqrt(trigamma(a_v)).
stump_params <- function(pb, a_tau = 3, tau0 = 1, a_v = 3, v0 = 1,
                         tau_phi = 1, base = 0, num_trees = 1L) {
  list(bart = bayesm.HART:::.parse_bart_params(
         list(num_trees = num_trees, base = base, a_tau = a_tau, tau0 = tau0),
         pb$nz, "direct", pb$D),
       vartree = bayesm.HART:::.parse_var_params(
         list(num_trees = num_trees, base = base,
              tau = sqrt(num_trees * trigamma(a_v)), v0 = v0),
         pb$nz, pb$D),
       phitree = bayesm.HART:::.parse_phi_params(
         list(num_trees = num_trees, base = base, tau = tau_phi), pb$nz))
}#STUMP_PARAMS

# update = c(mu, v, phi, f, tau)
run_blocks <- function(pb, pr, update, warm = 50L, nsweep = 4000L) {
  direct_hart_R(pb$theta, pb$Z, pb$mubar, pb$Amu, pr$bart, pr$vartree,
                pr$phitree, warm = warm, nsweep = nsweep, update = update)
}#RUN_BLOCKS

# update = c(mu, Sigma, f, tau)
run_constant <- function(pb, pr, nu, V, update, warm = 50L, nsweep = 4000L) {
  direct_hart_constant_R(pb$theta, pb$Z, pb$mubar, pb$Amu, nu, V, pr$bart,
                       warm = warm, nsweep = nsweep, update = update)
}#RUN_CONSTANT

unit_lower <- function(phi) {
  L <- diag(nrow(phi)); L[lower.tri(L)] <- -phi[lower.tri(phi)]; L
}#UNIT_LOWER

test_that("the kernels handed to the driver are mu + f and T^T V^{-1} T", {
  pb <- make_problem(); pr <- stump_params(pb, base = 0.95, num_trees = 3L)
  st <- run_blocks(pb, pr, update = rep(TRUE, 5), warm = 0L,
                   nsweep = 20L)$state
  expect_equal(st$prior_mean, st$f + matrix(st$mu, pb$n, pb$D, byrow = TRUE),
               tolerance = 1e-12)
  L <- unit_lower(st$phi)
  expect_equal(st$Q[, , 1], t(L) %*% diag(1 / st$v[1, ]) %*% L,
               tolerance = 1e-12)
  expect_true(all(is.finite(st$tau2)) && all(st$tau2 > 0))
  expect_gte(st$nleaves, 3L)
})#TEST_THAT

test_that("mu | rest is N(P^{-1} b, P^{-1}) with P = A + sum Q_i", {
  pb <- make_problem(); pr <- stump_params(pb)
  out <- run_blocks(pb, pr, update = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  st <- out$state
  P <- pb$Amu; b <- pb$Amu %*% pb$mubar
  for (i in seq_len(pb$n)) {
    Qi <- st$Q[, , i]
    P <- P + Qi; b <- b + Qi %*% (pb$theta[i, ] - st$f[i, ])
  }
  m <- solve(P, b); V <- solve(P)
  se <- sqrt(diag(V) / nrow(out$mu))
  expect_true(all(abs(colMeans(out$mu) - as.vector(m)) < 4 * se))
  expect_equal(cov(out$mu), V, tolerance = 0.15)
})#TEST_THAT

test_that("the stump mean leaf is N_D(P^{-1} sum Q (theta - mu), P^{-1})", {
  pb <- make_problem(); pr <- stump_params(pb)
  out <- run_blocks(pb, pr, update = c(FALSE, FALSE, FALSE, TRUE, FALSE))
  st <- out$state
  P <- diag(1 / st$tau2); b <- rep(0, pb$D)       # C = diag(tau2) / 1
  for (i in seq_len(pb$n)) {
    Qi <- st$Q[, , i]
    P <- P + Qi; b <- b + Qi %*% (pb$theta[i, ] - st$mu)
  }
  m <- solve(P, b); V <- solve(P)
  # One leaf: every unit carries the same vector.
  expect_true(all(abs(sweep(out$f, c(2, 3), out$f[1, , ])) < 1e-12))
  draws <- t(out$f[1, , ])
  se <- sqrt(diag(V) / nrow(draws))
  expect_true(all(abs(colMeans(draws) - as.vector(m)) < 4 * se))
  expect_equal(cov(draws), V, tolerance = 0.15)
})#TEST_THAT

test_that("tau_d^2 | f is IG(a_tau + L/2, b_tau + (H/2) sum m^2)", {
  pb <- make_problem(); pr <- stump_params(pb, a_tau = 4, tau0 = 0.7)
  out <- run_blocks(pb, pr, update = c(FALSE, FALSE, FALSE, FALSE, TRUE),
                    nsweep = 20000L)
  st <- out$state
  expect_equal(st$nleaves, 1L)
  shape <- 4 + 0.5
  rate  <- (4 - 1) * 0.7^2 + 0.5 * st$leaf_sumsq
  for (d in seq_len(pb$D)) {
    # 1 / tau^2 ~ Gamma(shape, rate).
    prec <- 1 / out$tau2[, d]
    expect_lt(abs(mean(prec) - shape / rate[d]),
              4 * sqrt(shape) / rate[d] / sqrt(length(prec)))
    expect_equal(var(prec), shape / rate[d]^2, tolerance = 0.1)
  }
})#TEST_THAT

test_that("the stump variance leaf is IG(a_v + n/2, b_v + sum r^2 / 2)", {
  pb <- make_problem(); pr <- stump_params(pb, a_v = 4, v0 = 0.5)
  out <- run_blocks(pb, pr, update = c(FALSE, TRUE, FALSE, FALSE, FALSE),
                    nsweep = 20000L)
  st <- out$state
  e <- pb$theta - st$f - matrix(st$mu, pb$n, pb$D, byrow = TRUE)
  r <- e %*% t(unit_lower(st$phi))
  shape <- 4 + pb$n / 2
  rate  <- (4 - 1) * 0.5 + colSums(r^2) / 2       # H_v = 1: b_v = (a_v - 1) v0
  for (j in seq_len(pb$D)) {
    prec <- 1 / out$v[1, j, ]
    expect_lt(abs(mean(prec) - shape / rate[j]),
              4 * sqrt(shape) / rate[j] / sqrt(length(prec)))
    expect_equal(var(prec), shape / rate[j]^2, tolerance = 0.1)
  }
})#TEST_THAT

test_that("the stump phi leaves target the Gaussian posterior row by row", {
  pb <- make_problem(); pr <- stump_params(pb, tau_phi = 0.8)
  out <- run_blocks(pb, pr, update = c(FALSE, FALSE, TRUE, FALSE, FALSE),
                    nsweep = 30000L)
  st <- out$state
  e <- pb$theta - st$f - matrix(st$mu, pb$n, pb$D, byrow = TRUE)
  v <- st$v[1, ]
  # Row j: r_j = e_j - sum_{k<j} lambda_jk e_k, r_ij ~ N(0, v_j),
  # lambda ~ N(0, tau_phi^2); the pairs of one row form a Gibbs chain.
  for (j in 2:pb$D) {
    E <- e[, seq_len(j - 1), drop = FALSE]
    P <- diag(1 / 0.8^2, j - 1) + crossprod(E) / v[j]
    b <- crossprod(E, e[, j]) / v[j]
    m <- solve(P, b); V <- solve(P)
    draws <- t(matrix(out$phi[j, seq_len(j - 1), ], nrow = j - 1))
    expect_true(all(abs(colMeans(draws) - as.vector(m)) <
                      8 * sqrt(diag(V) / nrow(draws))))
    expect_equal(cov(draws), V, tolerance = 0.15)
  }
})#TEST_THAT

test_that("constant Sigma | rest is IW(nu + n, V + E^T E)", {
  pb <- make_problem(); pr <- stump_params(pb)
  nu <- pb$D + 3; V <- nu * diag(pb$D)
  out <- run_constant(pb, pr, nu, V, update = c(FALSE, TRUE, FALSE, FALSE))
  st <- out$state
  E <- pb$theta - st$f - matrix(st$mu, pb$n, pb$D, byrow = TRUE)
  S <- V + crossprod(E)
  prec <- apply(out$Sigma, 3, solve)
  expect_equal(matrix(rowMeans(prec), pb$D, pb$D), (nu + pb$n) * solve(S),
               tolerance = 0.05)
  expect_equal(apply(out$Sigma, 1:2, mean), S / (nu + pb$n - pb$D - 1),
               tolerance = 0.05)
  expect_equal(st$Q[, , 1], solve(st$Sigma), tolerance = 1e-10)
  expect_equal(st$rooti %*% t(st$rooti), solve(st$Sigma), tolerance = 1e-10)
  expect_equal(st$prior_mean, st$f + matrix(st$mu, pb$n, pb$D, byrow = TRUE),
               tolerance = 1e-12)
})#TEST_THAT

test_that("constant model: mu | rest is N(P^{-1} b, P^{-1}) with P = A + n Q", {
  pb <- make_problem(); pr <- stump_params(pb)
  nu <- pb$D + 3; V <- nu * diag(pb$D)
  out <- run_constant(pb, pr, nu, V, update = c(TRUE, FALSE, FALSE, FALSE))
  st <- out$state
  Q <- st$Q[, , 1]
  P <- pb$Amu + pb$n * Q
  b <- pb$Amu %*% pb$mubar + Q %*% colSums(pb$theta - st$f)
  expect_equal(colMeans(out$mu), as.vector(solve(P, b)), tolerance = 0.05)
  expect_equal(cov(out$mu), solve(P), tolerance = 0.15)
})#TEST_THAT
