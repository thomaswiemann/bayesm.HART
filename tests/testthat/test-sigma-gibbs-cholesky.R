## Tests for sigma_mu_block_gibbs (Cholesky-Gibbs (mu, Sigma) update for HART).
##
## Math: discussions/2026-05-14-cholesky-gibbs-sigma.md
## Plan: discussions/2026-05-14-cholesky-gibbs-implementation-plan.md
##
## Helper: `run_chain` runs a sigma_mu_block_gibbs chain at fixed (theta, delta)
## and returns the posterior mean of (Sigma, mu).

run_chain <- function(theta, delta, Psi, nu, mubar0, Amu,
                      n_iter = 2000L, burn = 500L, seed = 1L,
                      L_init = NULL, mu_init = NULL) {
  set.seed(seed)
  D <- ncol(theta)
  if (is.null(L_init)) {
    Sigma_init <- diag(D)
    L_init <- t(chol(solve(Sigma_init)))     # lower-tri, L^T L = Sigma^{-1}
  }
  if (is.null(mu_init)) mu_init <- mubar0

  L <- L_init
  mu <- mu_init
  L_sum <- matrix(0, D, D)
  mu_sum <- numeric(D)
  Sigma_sum <- matrix(0, D, D)
  kept <- 0L
  L_chain <- array(NA_real_, dim = c(n_iter, D, D))
  mu_chain <- matrix(NA_real_, n_iter, D)

  for (r in seq_len(n_iter)) {
    out <- sigma_mu_block_gibbs_R(L, mu, theta, delta, Psi, nu, mubar0, Amu)
    L <- out$L
    mu <- as.vector(out$mu)
    L_chain[r, , ] <- L
    mu_chain[r, ] <- mu
    if (r > burn) {
      L_sum <- L_sum + L
      mu_sum <- mu_sum + mu
      Sigma_sum <- Sigma_sum + solve(crossprod(L))    # (L^T L)^{-1}
      kept <- kept + 1L
    }
  }

  list(
    L_mean     = L_sum / kept,
    mu_mean    = mu_sum / kept,
    Sigma_mean = Sigma_sum / kept,
    L_chain    = L_chain,
    mu_chain   = mu_chain
  )
}

## ----------------------------------------------------------------------------
## Test 1: delta = 0 reduces to IW + Gaussian conjugate posterior
##
## With delta = 0, our kernel is the standard IW posterior on Sigma and the
## standard Gaussian conjugate on mu.  Verify that the Cholesky-Gibbs chain
## reproduces the analytical IW(nu+n, Psi+S) posterior mean for Sigma.
## ----------------------------------------------------------------------------
test_that("delta = 0 chain reproduces IW conjugate posterior", {
  set.seed(123)
  D <- 3
  n <- 200
  Sigma_true <- diag(c(1, 2, 0.5))
  Sigma_half <- chol(Sigma_true)                  # upper-tri C, C^T C = Sigma
  mu_true <- c(0, 1, -1)
  theta <- matrix(rnorm(n * D), n, D) %*% Sigma_half +
           matrix(mu_true, n, D, byrow = TRUE)
  delta <- matrix(0, n, D)

  nu <- D + 2
  Psi <- diag(D)
  mubar0 <- rep(0, D)
  Amu <- 0.01

  fit <- run_chain(theta, delta, Psi, nu, mubar0, Amu,
                   n_iter = 4000L, burn = 1000L, seed = 42L)

  # Analytical IW posterior mean: (Psi + S) / (nu + n - D - 1)
  S <- crossprod(scale(theta, center = colMeans(theta), scale = FALSE))
  Sigma_iw_mean <- (Psi + S) / (nu + n - D - 1)

  # Compare diagonals first (most robust)
  expect_equal(diag(fit$Sigma_mean), diag(Sigma_iw_mean), tolerance = 0.05,
               label = "Sigma diagonals match IW posterior mean")
  # Off-diagonals (relative tolerance, allow more slack)
  expect_lt(max(abs(fit$Sigma_mean - Sigma_iw_mean)) /
              max(abs(Sigma_iw_mean)), 0.10,
            label = "max relative deviation from IW posterior < 10%")
})

## ----------------------------------------------------------------------------
## Test 2: mu posterior matches closed-form gold standard
##
## With L pinned (one-step draw), the mu posterior is N(mu_hat, Sigma/(n+a_mu))
## with mu_hat = (n*theta_bar + a_mu*mubar0 - L^{-1} S_delta) / (n + a_mu).
## ----------------------------------------------------------------------------
test_that("mu posterior matches closed-form (centered & shifted cases)", {
  set.seed(2026)
  D <- 4
  n <- 100

  # Construct fixed L, mu, theta, delta with known structure
  L_fixed <- matrix(0, D, D)
  diag(L_fixed) <- c(1.2, 0.8, 1.0, 1.5)
  L_fixed[2, 1] <- 0.3
  L_fixed[3, 1] <- -0.2; L_fixed[3, 2] <- 0.1
  L_fixed[4, 1] <- 0.05; L_fixed[4, 2] <- -0.1; L_fixed[4, 3] <- 0.2
  Sigma_fixed <- solve(crossprod(L_fixed))     # (L^T L)^{-1}

  mu_true <- c(1, -1, 0.5, 2)
  theta <- matrix(rnorm(n * D), n, D) %*% chol(Sigma_fixed) +
           matrix(mu_true, n, D, byrow = TRUE)

  Psi <- diag(D)
  nu <- D + 2
  mubar0 <- rep(0, D)
  Amu <- 0.5

  # ---- Case A: delta = 0  =>  mu_hat reduces to standard conjugate
  delta_zero <- matrix(0, n, D)
  S_delta_a <- rep(0, D)
  mu_hat_a <- (n * colMeans(theta) + Amu * mubar0 - solve(L_fixed, S_delta_a)) /
              (n + Amu)
  Sigma_post_a <- Sigma_fixed / (n + Amu)

  # Run a chain that reuses the fixed L by re-injecting it each iteration
  set.seed(7L)
  n_draws <- 5000L
  mu_draws_a <- matrix(NA_real_, n_draws, D)
  for (r in seq_len(n_draws)) {
    out <- sigma_mu_block_gibbs_R(L_fixed, mu_true, theta, delta_zero,
                                  Psi, nu, mubar0, Amu)
    mu_draws_a[r, ] <- as.vector(out$mu)
  }
  expect_equal(colMeans(mu_draws_a), mu_hat_a, tolerance = 0.05,
               label = "delta=0: empirical mu mean matches closed-form")
  # Empirical variance close to Sigma_post diagonal
  expect_equal(apply(mu_draws_a, 2, var), diag(Sigma_post_a), tolerance = 0.20,
               label = "delta=0: empirical mu variance matches Sigma/(n+a)")

  # ---- Case B: delta != 0  =>  mu_hat has -L^{-1} S_delta shift
  delta_nonzero <- matrix(rnorm(n * D, sd = 0.5), n, D)
  S_delta_b <- colSums(delta_nonzero)
  mu_hat_b <- (n * colMeans(theta) + Amu * mubar0 -
               solve(L_fixed, S_delta_b)) / (n + Amu)

  set.seed(8L)
  mu_draws_b <- matrix(NA_real_, n_draws, D)
  for (r in seq_len(n_draws)) {
    out <- sigma_mu_block_gibbs_R(L_fixed, mu_true, theta, delta_nonzero,
                                  Psi, nu, mubar0, Amu)
    mu_draws_b[r, ] <- as.vector(out$mu)
  }
  expect_equal(colMeans(mu_draws_b), mu_hat_b, tolerance = 0.05,
               label = "delta!=0: empirical mu mean matches closed-form with shift")

  # ---- Sanity: the shift difference between (B) and (A) is exactly -L^{-1} S_delta / (n+Amu)
  shift_expected <- -solve(L_fixed, S_delta_b) / (n + Amu)
  shift_observed <- colMeans(mu_draws_b) - colMeans(mu_draws_a)
  expect_equal(shift_observed, shift_expected, tolerance = 0.05,
               label = "shift A->B equals -L^{-1} S_delta / (n + Amu)")
})

## ----------------------------------------------------------------------------
## Test 3: Bartlett prior reproduction (n = 0, delta = 0)
##
## With no data and Psi = I, the L_jj marginal should approximately follow
## chi_{nu - D + j} (df = nu - D + j).  Diagonal exponent in the kernel is
## a_j = nu + n - D + j - 1 = nu - D + j - 1 (for n = 0), matching chi.
## ----------------------------------------------------------------------------
test_that("n=0 reproduces Bartlett prior diagonals (chi distributions)", {
  set.seed(2027)
  D <- 4
  nu <- 12        # generous df so the prior is well-behaved
  Psi <- diag(D)
  mubar0 <- rep(0, D)
  Amu <- 1.0

  theta <- matrix(0, 0, D)
  delta <- matrix(0, 0, D)

  L_init <- diag(D)
  fit <- run_chain(theta, delta, Psi, nu, mubar0, Amu,
                   n_iter = 8000L, burn = 1000L, seed = 13L,
                   L_init = L_init, mu_init = mubar0)

  for (j in seq_len(D)) {
    # E[chi_k] = sqrt(2) * Gamma((k+1)/2) / Gamma(k/2)
    k <- nu - D + j     # df under our kernel for j (1-indexed)
    expected_mean <- sqrt(2) * exp(lgamma((k + 1) / 2) - lgamma(k / 2))
    empirical_mean <- mean(fit$L_chain[1001:8000, j, j])
    expect_equal(empirical_mean, expected_mean, tolerance = 0.05,
                 label = sprintf("L[%d,%d] mean ~ chi_%d", j, j, k))
  }

  # Off-diagonals: prior is iid N(0, 1) for Psi = I.
  for (j in 2:D) {
    for (k in 1:(j - 1)) {
      empirical_mean <- mean(fit$L_chain[1001:8000, j, k])
      expect_lt(abs(empirical_mean), 0.05,
                label = sprintf("L[%d,%d] off-diag mean ~ 0", j, k))
    }
  }
})

## ----------------------------------------------------------------------------
## Test 4: Permutation invariance of Sigma posterior
##
## The IW posterior on Sigma is invariant under joint permutation of theta,
## delta, Psi, mubar0.  Our Cholesky-Gibbs sampler should produce permuted
## Sigma posterior means (modulo MC error).
## ----------------------------------------------------------------------------
test_that("Sigma posterior is permutation-invariant", {
  set.seed(2028)
  D <- 3
  n <- 250
  Sigma_true <- matrix(c(2, 0.5, 0,
                         0.5, 1, 0.3,
                         0, 0.3, 1.5), D, D)
  mu_true <- c(0.5, -0.5, 1)
  theta <- matrix(rnorm(n * D), n, D) %*% chol(Sigma_true) +
           matrix(mu_true, n, D, byrow = TRUE)
  delta <- matrix(rnorm(n * D, sd = 0.3), n, D)

  Psi <- diag(D)
  nu <- D + 3
  mubar0 <- rep(0, D)
  Amu <- 0.1

  perm <- c(3, 1, 2)
  P <- diag(D)[perm, ]    # permutation matrix: rows of P pick perm

  # Original chain
  fit_orig <- run_chain(theta, delta, Psi, nu, mubar0, Amu,
                        n_iter = 3000L, burn = 500L, seed = 100L)

  # Permuted chain (theta, delta, Psi, mubar0 all permuted)
  fit_perm <- run_chain(theta[, perm], delta[, perm],
                        Psi[perm, perm], nu, mubar0[perm], Amu,
                        n_iter = 3000L, burn = 500L, seed = 100L)

  # Compare permuted -> un-permuted: P^T Sigma_perm P should equal Sigma_orig
  Sigma_back <- t(P) %*% fit_perm$Sigma_mean %*% P
  expect_equal(Sigma_back, fit_orig$Sigma_mean, tolerance = 0.10,
               label = "Sigma posterior permutation-invariant (within MC tol)")

  mu_back <- as.vector(t(P) %*% fit_perm$mu_mean)
  expect_equal(mu_back, fit_orig$mu_mean, tolerance = 0.05,
               label = "mu posterior permutation-invariant")
})

## ----------------------------------------------------------------------------
## Test 5: Slice sampler robustness on extremes
##
## Hammer the diagonal slice with extreme regimes (large n, small Psi, large
## |c|) and confirm no infinite loops, no negative draws, no NaN.
## ----------------------------------------------------------------------------
test_that("slice sampler is robust to extremes", {
  ## Drive the slice with regimes that stress (a) large n -> large b,
  ## (b) large |c| via cross-correlated delta, (c) varying Sigma scales.
  ## Verify positivity / finiteness; skip the post-hoc Sigma inversion since
  ## extreme L can be near-singular by design.
  set.seed(2029)
  D <- 3
  n <- 1000L
  Sigma_true <- diag(c(0.1, 10, 1))           # 100x scale span
  mu_true <- rep(0, D)
  theta <- matrix(rnorm(n * D), n, D) %*% chol(Sigma_true)
  delta <- matrix(rnorm(n * D, sd = 2), n, D)  # heavy cross terms

  nu <- D + 2
  Psi <- diag(c(0.01, 10, 1))
  mubar0 <- rep(0, D)
  Amu <- 1e-3

  set.seed(17L)
  L <- t(chol(solve(diag(D))))
  mu <- mubar0
  L_chain <- array(NA_real_, dim = c(200L, D, D))
  for (r in seq_len(200L)) {
    out <- sigma_mu_block_gibbs_R(L, mu, theta, delta, Psi, nu, mubar0, Amu)
    L <- out$L
    mu <- as.vector(out$mu)
    L_chain[r, , ] <- L
  }

  expect_true(all(is.finite(L_chain)),
              label = "no NaN/Inf in L draws under extremes")
  for (j in seq_len(D)) {
    expect_true(all(L_chain[, j, j] > 0),
                label = sprintf("L[%d,%d] strictly positive across all draws", j, j))
  }
})

## ----------------------------------------------------------------------------
## Test 6: End-to-end via rhierLinearMixture (smoke + sanity)
##
## Run the full BART-path MCMC on a small synthetic dataset and check
## that the chain produces sensible Sigma draws (positive-definite, finite).
## This is mostly a smoke test that the wiring in mixbart_block.cpp works.
## ----------------------------------------------------------------------------
test_that("rhierLinearMixture BART path runs end-to-end with new sampler", {
  skip_if_not(exists("rhierLinearMixture"))
  skip_if_not(exists("sim_hier_linear"))

  set.seed(2030)
  sim <- sim_hier_linear(nreg = 25, nobs = 15, nvar = 2, nz = 2,
                         seed = 2030)
  fit <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z),
    Prior = list(ncomp = 1, bart = list(num_trees = 5)),
    Mcmc  = list(R = 30, keep = 1, nprint = 0))

  expect_true(!is.null(fit$nmix$compdraw), label = "compdraw populated")
  comp_last <- fit$nmix$compdraw[[length(fit$nmix$compdraw)]][[1]]
  rooti_last <- comp_last$rooti
  Sigma_last <- solve(rooti_last %*% t(rooti_last))

  expect_true(all(is.finite(rooti_last)), label = "rooti finite")
  expect_true(all(diag(rooti_last) > 0), label = "rooti diag positive")
  expect_true(all(eigen(Sigma_last)$values > 0),
              label = "Sigma positive definite")
})
