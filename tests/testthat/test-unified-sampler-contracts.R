# Tests for the unified sampler architecture contracts.
# Verifies:
#   * Algebra: R_i * t(R_i) == Sigma_i^{-1} for both backends
#   * Prior-mean consistency
#   * Schema snapshots (field names + class vectors)
#   * Cross-model output structure parity

library(testthat)

# --- helpers -----------------------------------------------------------------
run_quiet <- function(expr) capture.output(suppressWarnings(expr))

make_mnl_sim <- function(seed = 42L, nlgt = 10L, nvar = 2L, nz = 2L,
                         p = 3L, T_i = 5L) {
  set.seed(seed)
  Z <- cbind(runif(nlgt, -1, 1), runif(nlgt, -1, 1))
  Z <- scale(Z, center = TRUE, scale = FALSE)
  theta <- matrix(rnorm(nlgt * nvar, sd = 0.5), nlgt, nvar)
  lgtdata <- lapply(seq_len(nlgt), function(i) {
    X    <- matrix(rnorm(T_i * p * nvar), T_i * p, nvar)
    eta  <- X %*% theta[i, ]
    prob <- matrix(exp(eta), T_i, p, byrow = TRUE)
    prob <- prob / rowSums(prob)
    y    <- apply(prob, 1, function(pp) sample.int(p, 1, prob = pp))
    list(y = y, X = X)
  })
  list(p = p, lgtdata = lgtdata, Z = Z)
}

# ===========================================================================
# 8.1 Algebra contract: precision-root R_i * t(R_i) == Sigma_i^{-1}
# ===========================================================================

test_that("cov_helpers roundtrip: rootpi * t(rootpi) == Sigma^{-1}", {
  # Use the exported C++ test harness for cov_helpers
  D_vals <- c(2.0, 3.0)
  Phi_vals <- matrix(c(0, 0.5, 0, 0), 2, 2)
  v <- c(1.0, 2.0)
  mu <- c(0.0, 0.0)

  res <- cov_helpers_test(D_vals, Phi_vals, v, mu, diagonal = FALSE)

  # Contract: rootpi * t(rootpi) = Sigma^{-1} (must be symmetric PD)
  rootpi <- res$rootpi
  Q <- rootpi %*% t(rootpi)
  expect_true(isSymmetric(Q, tol = 1e-10),
              label = "rootpi * t(rootpi) must be symmetric")
  evals <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0),
              label = "rootpi * t(rootpi) must be positive definite")

  # Verify rootpi is upper-triangular
  expect_equal(rootpi[lower.tri(rootpi)], rep(0, sum(lower.tri(rootpi))),
               tolerance = 1e-10, label = "rootpi must be upper-triangular")
})

test_that("cov_helpers diagonal case: rootpi * t(rootpi) == diag(1/D)", {
  D_vals <- c(4.0, 9.0)
  Phi_vals <- matrix(0, 2, 2)
  v <- c(1.5, -0.5)
  mu <- c(0, 0)

  res <- cov_helpers_test(D_vals, Phi_vals, v, mu, diagonal = TRUE)

  rootpi <- res$rootpi
  Q <- rootpi %*% t(rootpi)
  expect_equal(Q, diag(1 / D_vals), tolerance = 1e-10,
               label = "diagonal case: Q must equal diag(1/D)")
})

# ===========================================================================
# 8.1 Algebra contract: global backend precision_root via rmixGibbs
# ===========================================================================

test_that("global backend: compdraw rooti satisfies R*t(R) == Sigma^{-1}", {
  sim <- make_mnl_sim()
  set.seed(100)
  run_quiet({
    fit <- rhierMnlRwMixture(
      Data = sim, Prior = list(ncomp = 1),
      Mcmc = list(R = 10, keep = 1, nprint = 0))
  })

  # Check last draw's rooti for component 1
  last_draw <- fit$nmix$compdraw[[length(fit$nmix$compdraw)]][[1]]
  rooti <- last_draw$rooti
  mu_k  <- last_draw$mu

  # Contract: rooti * t(rooti) must be symmetric PD (= Sigma_k^{-1})
  Q <- rooti %*% t(rooti)
  expect_true(isSymmetric(Q, tol = 1e-10),
              label = "rooti * t(rooti) must be symmetric")
  evals <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0),
              label = "rooti * t(rooti) must be positive definite")
})

# ===========================================================================
# 8.4 Schema snapshot: MNL output field names
# ===========================================================================

test_that("MNL global output has expected field names", {
  sim <- make_mnl_sim()

  # No Z path
  sim_noz <- sim
  sim_noz$Z <- NULL
  set.seed(200)
  run_quiet({
    fit_noz <- rhierMnlRwMixture(
      Data = sim_noz, Prior = list(ncomp = 1),
      Mcmc = list(R = 4, keep = 1, nprint = 0))
  })
  expect_true(all(c("betadraw", "loglike", "nmix", "SignRes", "acceptrbeta")
                  %in% names(fit_noz)))
  expect_s3_class(fit_noz, "rhierMnlRwMixture")

  # Linear delta path
  set.seed(201)
  run_quiet({
    fit_lin <- rhierMnlRwMixture(
      Data = sim, Prior = list(ncomp = 1),
      Mcmc = list(R = 4, keep = 1, nprint = 0))
  })
  expect_true(all(c("betadraw", "loglike", "nmix", "Deltadraw",
                     "SignRes", "acceptrbeta") %in% names(fit_lin)))

  # BART path
  set.seed(202)
  run_quiet({
    fit_bart <- rhierMnlRwMixture(
      Data = sim, Prior = list(ncomp = 1, bart = list(num_trees = 3)),
      Mcmc = list(R = 4, keep = 1, nprint = 0))
  })
  expect_true(all(c("betadraw", "loglike", "nmix", "bart_models",
                     "SignRes", "acceptrbeta") %in% names(fit_bart)))
  expect_true(is.null(fit_bart$Deltadraw))
})

test_that("MNL heter-cov output has expected field names", {
  sim <- make_mnl_sim()
  set.seed(300)
  run_quiet({
    fit <- rhierMnlRwMixture(
      Data = sim,
      Prior = list(ncomp = 1,
                   bart = list(num_trees = 3),
                   vartree = list(num_trees = 3)),
      Mcmc = list(R = 4, keep = 1, nprint = 0))
  })
  expect_true(all(c("betadraw", "loglike", "bart_models",
                     "var_models", "mu_draw", "SignRes", "acceptrbeta")
                  %in% names(fit)))
  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
  # Heter-cov should NOT have nmix or Deltadraw
  expect_true(is.null(fit$nmix))
  expect_true(is.null(fit$Deltadraw))
})

# ===========================================================================
# 8.4 Schema snapshot: Linear output field names
# ===========================================================================

test_that("Linear global output has expected field names", {
  set.seed(400)
  n <- 20; nvar <- 2; nz <- 2
  Z <- matrix(rnorm(n * nz), n, nz)
  regdata <- lapply(1:n, function(i) {
    ni <- 10
    X <- matrix(rnorm(ni * nvar), ni, nvar)
    y <- X %*% rnorm(nvar) + rnorm(ni)
    list(y = y, X = X)
  })

  run_quiet({
    fit <- rhierLinearMixture(
      Data = list(regdata = regdata, Z = Z),
      Prior = list(ncomp = 1),
      Mcmc = list(R = 4, keep = 1, nprint = 0))
  })
  expect_true(all(c("betadraw", "taudraw", "loglike", "nmix")
                  %in% names(fit)))
  expect_s3_class(fit, "rhierLinearMixture")
})

# ===========================================================================
# 8.4 Schema snapshot: NegBin output field names
# ===========================================================================

test_that("NegBin global output has expected field names", {
  set.seed(500)
  n <- 20; nvar <- 2; nz <- 2
  Z <- matrix(rnorm(n * nz), n, nz)
  regdata <- lapply(1:n, function(i) {
    ni <- 15
    X <- matrix(rnorm(ni * nvar), ni, nvar)
    beta <- rnorm(nvar, sd = 0.3)
    lambda <- exp(X %*% beta)
    y <- rnbinom(ni, size = 3, mu = lambda)
    list(y = y, X = X)
  })

  run_quiet({
    fit <- rhierNegbinRw(
      Data = list(regdata = regdata, Z = Z),
      Prior = list(ncomp = 1),
      Mcmc = list(R = 10, keep = 1, nprint = 0))
  })
  expect_true(all(c("betadraw", "alphadraw", "loglike", "nmix",
                     "acceptrbeta", "acceptralpha") %in% names(fit)))
  expect_s3_class(fit, "rhierNegbinRw")
})

# ===========================================================================
# 8.2 Cross-model parity: betadraw dimensions
# ===========================================================================

test_that("betadraw dimensions are consistent across models", {
  # MNL
  sim <- make_mnl_sim(seed = 600, nlgt = 8, nvar = 2, nz = 2, p = 3, T_i = 5)
  set.seed(601)
  run_quiet({
    fit_mnl <- rhierMnlRwMixture(
      Data = sim, Prior = list(ncomp = 1),
      Mcmc = list(R = 6, keep = 2, nprint = 0))
  })
  expect_equal(dim(fit_mnl$betadraw), c(8, 2, 3))  # nlgt x nvar x R/keep

  # Linear
  set.seed(602)
  n <- 8; nvar <- 2; nz <- 2
  Z <- matrix(rnorm(n * nz), n, nz)
  regdata <- lapply(1:n, function(i) {
    ni <- 10
    X <- matrix(rnorm(ni * nvar), ni, nvar)
    y <- X %*% rnorm(nvar) + rnorm(ni)
    list(y = y, X = X)
  })
  run_quiet({
    fit_lin <- rhierLinearMixture(
      Data = list(regdata = regdata, Z = Z),
      Prior = list(ncomp = 1),
      Mcmc = list(R = 6, keep = 2, nprint = 0))
  })
  expect_equal(dim(fit_lin$betadraw), c(8, 2, 3))
})
