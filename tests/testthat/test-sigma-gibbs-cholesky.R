## Tests for the normalized HART (mu, L) blocks in src/sigma_gibbs.cpp.
##
##   theta_i = mu + L^{-1} (f*_i + eps_i),  eps_i ~ N(0, I),
##   L lower triangular with L^T L = Sigma^{-1},
##   mu ~ N(mubar, Amu^{-1}),  Sigma ~ IW(nu, Psi).
##
## Harness: normalized_kernel_R() returns the conditional kernels at a state;
## normalized_gibbs_R() runs sweeps with either block optionally pinned.

library(testthat)

# Independent joint log-density of (mu, L) given (theta, f*), up to a constant:
# the IW prior on Sigma = (L^T L)^{-1} mapped to L with its Jacobian, the
# Gaussian prior on mu, and the whitened likelihood with its |L|^n factor.
# Jacobians: Sigma = Q^{-1} gives |Q|^{-(D+1)}; Q = L^T L = U U^T with
# U = L^T upper triangular gives 2^D prod_j L_jj^j.
log_joint <- function(L, mu, theta, fstar, Psi, nu, mubar, Amu) {
  D <- nrow(L); n <- nrow(theta)
  Q <- crossprod(L)
  logdetQ <- as.numeric(determinant(Q, logarithm = TRUE)$modulus)
  log_iw  <- (nu + D + 1) / 2 * logdetQ - 0.5 * sum(diag(Psi %*% Q))
  log_jac <- -(D + 1) * logdetQ + sum(seq_len(D) * log(diag(L)))
  log_mu  <- -0.5 * drop(t(mu - mubar) %*% Amu %*% (mu - mubar))
  resid   <- (theta - matrix(mu, n, D, byrow = TRUE)) %*% t(L) - fstar
  log_lik <- n * sum(log(diag(L))) - 0.5 * sum(resid^2)
  log_iw + log_jac + log_mu + log_lik
}#LOG_JOINT

# A small problem with nothing special: nonzero mu, nondiagonal L, nonzero f*,
# non-identity Psi, non-diagonal Amu.
make_problem <- function(seed = 1L, D = 3L, n = 15L) {
  set.seed(seed)
  L <- diag(c(1.3, 0.7, 1.1))
  L[2, 1] <- 0.4; L[3, 1] <- -0.3; L[3, 2] <- 0.2
  mu    <- c(0.5, -1, 0.25)
  theta <- matrix(rnorm(n * D), n, D) %*% solve(L) +
           matrix(mu, n, D, byrow = TRUE)
  fstar <- matrix(rnorm(n * D, sd = 0.4), n, D)
  Psi   <- crossprod(matrix(rnorm(D * D), D, D)) + diag(D)
  Amu   <- crossprod(matrix(rnorm(D * D, sd = 0.3), D, D)) + 0.5 * diag(D)
  list(L = L, mu = mu, theta = theta, fstar = fstar, Psi = Psi, nu = D + 2,
       mubar = c(0.1, 0.2, -0.1), Amu = Amu)
}#MAKE_PROBLEM

run_chain <- function(theta, fstar, Psi, nu, mubar, Amu, n_iter, burn,
                      L_init, mu_init, update_mu = TRUE, update_L = TRUE,
                      seed = 1L) {
  set.seed(seed)
  D <- ncol(L_init)
  L <- L_init; mu <- mu_init
  L_chain  <- array(NA_real_, dim = c(n_iter, D, D))
  mu_chain <- matrix(NA_real_, n_iter, D)
  Sigma_sum <- matrix(0, D, D); kept <- 0L
  for (r in seq_len(n_iter)) {
    out <- normalized_gibbs_R(L, mu, theta, fstar, Psi, nu, mubar, Amu,
                              nsweep = 1L, update_mu = update_mu,
                              update_L = update_L)
    L <- out$L; mu <- as.vector(out$mu)
    L_chain[r, , ] <- L; mu_chain[r, ] <- mu
    if (r > burn) { Sigma_sum <- Sigma_sum + solve(crossprod(L)); kept <- kept + 1L }
  }#FOR
  list(L_chain = L_chain, mu_chain = mu_chain, Sigma_mean = Sigma_sum / kept,
       mu_mean = colMeans(mu_chain[(burn + 1):n_iter, , drop = FALSE]))
}#RUN_CHAIN

## ----------------------------------------------------------------------------
## Full-density checks: each conditional kernel used by the sampler must
## reproduce the differences of the independently coded joint log-density.
## ----------------------------------------------------------------------------
test_that("off-diagonal kernels match the joint density", {
  pb <- make_problem()
  k <- with(pb, normalized_kernel_R(L, theta, fstar, mu, Psi, nu, mubar, Amu))
  set.seed(11)
  for (j in 2:3) for (kk in seq_len(j - 1)) {
    for (rep in 1:3) {
      L1 <- pb$L; L1[j, kk] <- pb$L[j, kk] + rnorm(1)
      dJ <- with(pb, log_joint(L1, mu, theta, fstar, Psi, nu, mubar, Amu) -
                     log_joint(L,  mu, theta, fstar, Psi, nu, mubar, Amu))
      m <- k$offdiag_mean[j, kk]; p <- k$offdiag_prec[j, kk]
      dK <- -0.5 * p * (L1[j, kk] - m)^2 + 0.5 * p * (pb$L[j, kk] - m)^2
      expect_equal(dJ, dK, tolerance = 1e-8,
                   label = sprintf("L[%d,%d] kernel difference", j, kk))
    }
  }
})#TEST_THAT

test_that("diagonal kernels match the joint density", {
  pb <- make_problem()
  k <- with(pb, normalized_kernel_R(L, theta, fstar, mu, Psi, nu, mubar, Amu))
  set.seed(12)
  for (j in 1:3) for (rep in 1:3) {
    L1 <- pb$L; L1[j, j] <- pb$L[j, j] * exp(rnorm(1, sd = 0.5))
    dJ <- with(pb, log_joint(L1, mu, theta, fstar, Psi, nu, mubar, Amu) -
                   log_joint(L,  mu, theta, fstar, Psi, nu, mubar, Amu))
    kern <- function(x) k$diag_a[j] * log(x) - 0.5 * k$diag_b[j] * x^2 - k$diag_c[j] * x
    expect_equal(dJ, kern(L1[j, j]) - kern(pb$L[j, j]), tolerance = 1e-8,
                 label = sprintf("L[%d,%d] kernel difference", j, j))
  }
  # Exponent a_j = nu + n - D + j - 1 (j 1-based).
  expect_equal(as.vector(k$diag_a), pb$nu + nrow(pb$theta) - 3 + (1:3) - 1)
})#TEST_THAT

test_that("mu kernel matches the joint density", {
  pb <- make_problem()
  k <- with(pb, normalized_kernel_R(L, theta, fstar, mu, Psi, nu, mubar, Amu))
  m <- solve(k$mu_prec, k$mu_b)
  set.seed(13)
  for (rep in 1:5) {
    mu1 <- pb$mu + rnorm(3)
    dJ <- with(pb, log_joint(L, mu1, theta, fstar, Psi, nu, mubar, Amu) -
                   log_joint(L, mu,  theta, fstar, Psi, nu, mubar, Amu))
    dK <- -0.5 * drop(t(mu1 - m) %*% k$mu_prec %*% (mu1 - m)) +
           0.5 * drop(t(pb$mu - m) %*% k$mu_prec %*% (pb$mu - m))
    expect_equal(dJ, dK, tolerance = 1e-8)
  }
  # Closed form: P = Amu + n L^T L, b = Amu mubar + L^T sum_i (L theta_i - f_i).
  n <- nrow(pb$theta)
  expect_equal(k$mu_prec, pb$Amu + n * crossprod(pb$L))
  s <- colSums(pb$theta %*% t(pb$L) - pb$fstar)
  expect_equal(as.vector(k$mu_b), as.vector(pb$Amu %*% pb$mubar + t(pb$L) %*% s))
})#TEST_THAT

## ----------------------------------------------------------------------------
## Draws
## ----------------------------------------------------------------------------
test_that("mu draws follow the conditional at a pinned L", {
  pb <- make_problem(seed = 2L, n = 40L)
  k <- with(pb, normalized_kernel_R(L, theta, fstar, mu, Psi, nu, mubar, Amu))
  m <- solve(k$mu_prec, k$mu_b); V <- solve(k$mu_prec)
  set.seed(21)
  draws <- t(replicate(4000, as.vector(
    with(pb, normalized_gibbs_R(L, mu, theta, fstar, Psi, nu, mubar, Amu,
                                nsweep = 1L, update_L = FALSE))$mu)))
  se <- sqrt(diag(V) / 4000)
  expect_true(all(abs(colMeans(draws) - m) < 4 * se))
  expect_equal(cov(draws), V, tolerance = 0.15)
})#TEST_THAT

test_that("f* = 0 at a pinned mu reproduces the IW(nu + n, Psi + S) posterior", {
  set.seed(123)
  D <- 3; n <- 200
  Sigma_true <- diag(c(1, 2, 0.5))
  mu_true <- c(0, 1, -1)
  theta <- matrix(rnorm(n * D), n, D) %*% chol(Sigma_true) +
           matrix(mu_true, n, D, byrow = TRUE)
  fstar <- matrix(0, n, D)
  nu <- D + 2; Psi <- diag(D)
  fit <- run_chain(theta, fstar, Psi, nu, mubar = rep(0, D), Amu = 0.01 * diag(D),
                   n_iter = 4000L, burn = 1000L, L_init = diag(D), mu_init = mu_true,
                   update_mu = FALSE, seed = 42L)
  S <- crossprod(theta - matrix(mu_true, n, D, byrow = TRUE))
  Sigma_iw_mean <- (Psi + S) / (nu + n - D - 1)
  expect_equal(diag(fit$Sigma_mean), diag(Sigma_iw_mean), tolerance = 0.05)
  expect_lt(max(abs(fit$Sigma_mean - Sigma_iw_mean)) / max(abs(Sigma_iw_mean)), 0.10)
})#TEST_THAT

test_that("n = 0 reproduces the Bartlett prior diagonals", {
  D <- 4; nu <- 12; Psi <- diag(D)
  theta <- matrix(0, 0, D); fstar <- matrix(0, 0, D)
  fit <- run_chain(theta, fstar, Psi, nu, mubar = rep(0, D), Amu = diag(D),
                   n_iter = 8000L, burn = 1000L, L_init = diag(D),
                   mu_init = rep(0, D), update_mu = FALSE, seed = 13L)
  for (j in seq_len(D)) {
    # L_jj ~ chi_{nu - D + j}: E = sqrt(2) Gamma((k+1)/2) / Gamma(k/2).
    k <- nu - D + j
    expected_mean <- sqrt(2) * exp(lgamma((k + 1) / 2) - lgamma(k / 2))
    expect_equal(mean(fit$L_chain[1001:8000, j, j]), expected_mean, tolerance = 0.05,
                 label = sprintf("L[%d,%d] mean ~ chi_%d", j, j, k))
  }
  for (j in 2:D) for (k in 1:(j - 1))
    expect_lt(abs(mean(fit$L_chain[1001:8000, j, k])), 0.05)
})#TEST_THAT

test_that("the (mu, Sigma) posterior is permutation invariant", {
  set.seed(2028)
  D <- 3; n <- 250
  Sigma_true <- matrix(c(2, 0.5, 0, 0.5, 1, 0.3, 0, 0.3, 1.5), D, D)
  mu_true <- c(0.5, -0.5, 1)
  theta <- matrix(rnorm(n * D), n, D) %*% chol(Sigma_true) +
           matrix(mu_true, n, D, byrow = TRUE)
  fstar <- matrix(rnorm(n * D, sd = 0.3), n, D)
  Psi <- diag(D); nu <- D + 3; mubar <- rep(0, D); Amu <- 0.1 * diag(D)
  perm <- c(3, 1, 2); P <- diag(D)[perm, ]
  fit_orig <- run_chain(theta, fstar, Psi, nu, mubar, Amu, 3000L, 500L,
                        diag(D), mubar, seed = 100L)
  fit_perm <- run_chain(theta[, perm], fstar[, perm], Psi[perm, perm], nu,
                        mubar[perm], Amu[perm, perm], 3000L, 500L,
                        diag(D), mubar[perm], seed = 100L)
  expect_equal(t(P) %*% fit_perm$Sigma_mean %*% P, fit_orig$Sigma_mean,
               tolerance = 0.10)
  expect_equal(as.vector(t(P) %*% fit_perm$mu_mean), fit_orig$mu_mean,
               tolerance = 0.05)
})#TEST_THAT

test_that("the slice sampler is robust to extremes", {
  set.seed(2029)
  D <- 3; n <- 1000L
  theta <- matrix(rnorm(n * D), n, D) %*% chol(diag(c(0.1, 10, 1)))
  fstar <- matrix(rnorm(n * D, sd = 2), n, D)
  fit <- run_chain(theta, fstar, Psi = diag(c(0.01, 10, 1)), nu = D + 2,
                   mubar = rep(0, D), Amu = 1e-3 * diag(D), n_iter = 200L,
                   burn = 0L, L_init = diag(D), mu_init = rep(0, D), seed = 17L)
  expect_true(all(is.finite(fit$L_chain)))
  for (j in seq_len(D)) expect_true(all(fit$L_chain[, j, j] > 0))
})#TEST_THAT

## ----------------------------------------------------------------------------
## End to end through rhierLinearMixture (smoke test of the wiring)
## ----------------------------------------------------------------------------
test_that("rhierLinearMixture normalized path runs end-to-end", {
  set.seed(2030)
  sim <- sim_hier_linear(nreg = 25, nobs = 15, nvar = 2, nz = 2,
                         het_observed = "linear", seed = 2030)
  fit <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, hart_model = "normalized", bart = list(num_trees = 5)),
    Mcmc  = list(R = 30, keep = 1, nprint = 0), r_verbose = FALSE)
  expect_equal(fit$hart$model, "normalized")
  expect_equal(fit$hart$priors$Amu, 0.01 * diag(ncol(fit$betadraw)))
  expect_true(!is.null(fit$nmix$compdraw))
  comp_last <- fit$nmix$compdraw[[length(fit$nmix$compdraw)]][[1]]
  rooti_last <- comp_last$rooti
  expect_true(all(is.finite(rooti_last)))
  expect_true(all(diag(rooti_last) > 0))
  expect_equal(rooti_last[lower.tri(rooti_last)],
               rep(0, sum(lower.tri(rooti_last))))
  Sigma_last <- solve(rooti_last %*% t(rooti_last))
  expect_true(all(eigen(Sigma_last)$values > 0))
  expect_equal(dim(fit$nmix$probdraw), c(30L, 1L))
})#TEST_THAT
