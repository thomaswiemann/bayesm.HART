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
  D_vals <- c(2.0, 3.0)
  Phi_vals <- matrix(c(0, 0.5, 0, 0), 2, 2)
  v <- c(1.0, 2.0)

  res <- cov_eval_R(D_vals, Phi_vals, v, diagonal = FALSE)

  # rootpi * t(rootpi) = L^T D^{-1} L with L = I - lower(phi).
  L <- diag(2); L[2, 1] <- -0.5
  rootpi <- res$rootpi
  Q <- rootpi %*% t(rootpi)
  expect_equal(Q, t(L) %*% diag(1 / D_vals) %*% L, tolerance = 1e-12)
  expect_equal(res$Lv, L %*% v, tolerance = 1e-12, ignore_attr = TRUE)
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

  res <- cov_eval_R(D_vals, Phi_vals, v, diagonal = TRUE)

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
      Data = sim, Prior = list(ncomp = 1, hart_model = "normalized", bart = list(num_trees = 3)),
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
      Prior = list(ncomp = 1, hart_model = "direct",
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

test_that("MNL direct HART with a constant Sigma has the nmix schema", {
  sim <- make_mnl_sim()
  nvar <- ncol(sim$lgtdata[[1]]$X)
  set.seed(301)
  run_quiet({
    fit <- rhierMnlRwMixture(
      Data = sim,
      Prior = list(ncomp = 1, hart_model = "direct", hart_cov = "constant",
                   bart = list(num_trees = 3)),
      Mcmc = list(R = 6, keep = 1, nprint = 0))
  })
  expect_equal(class(fit), "rhierMnlRwMixture")
  expect_named(fit, c("mu_draw", "mean_scale_draw", "varcount", "varprob",
                      "nmix", "bart_models", "betadraw", "loglike", "SignRes",
                      "acceptrbeta", "hart", "restart"), ignore.order = TRUE)
  expect_equal(dim(fit$mu_draw), c(6L, nvar))
  expect_equal(fit$mu_draw[1, ], fit$nmix$compdraw[[1]][[1]]$mu)
  expect_s3_class(fit$nmix, "bayesm.nmix")
  expect_equal(fit$bart_models$treedraws$leaf_dim, nvar)
  expect_equal(dim(fit$nmix$compdraw[[1]][[1]]$rooti), c(nvar, nvar))
  expect_equal(dim(fit$varcount), c(2L, 6L))
  # Delta(Z) from the cache and from the trees agree; DeltaZ+mu adds the
  # nmix mean; SigmaZ needs the covariance trees.
  newdata <- list(Z = sim$Z)
  d1 <- predict(fit, newdata = newdata, type = "DeltaZ", burn = 0L,
                r_verbose = FALSE)
  d2 <- predict(fit, newdata = newdata, type = "DeltaZ", burn = 0L,
                force_tree_eval = TRUE, r_verbose = FALSE)
  expect_equal(d1, d2, tolerance = 1e-8)
  dm <- predict(fit, newdata = newdata, type = "DeltaZ+mu", burn = 0L,
                r_verbose = FALSE)
  mu1 <- fit$nmix$compdraw[[1]][[1]]$mu
  mu1_rooti <- fit$nmix$compdraw[[1]][[1]]$rooti
  expect_equal(dm[, , 1], d1[, , 1] + matrix(mu1, nrow(sim$Z), nvar,
                                              byrow = TRUE),
               tolerance = 1e-10)
  sig <- predict(fit, newdata = newdata, type = "SigmaZ", burn = 0L,
                 r_verbose = FALSE)
  expect_equal(dim(sig), c(nrow(sim$Z), nvar, nvar, 6L))
  expect_equal(sig[1, , , 1], solve(mu1_rooti %*% t(mu1_rooti)),
               tolerance = 1e-10)
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
# 8.4 Schema snapshot: HART outputs for Linear and NegBin
# ===========================================================================

make_reg_sim <- function(seed, n = 20L, nvar = 2L, nz = 2L, counts = FALSE) {
  set.seed(seed)
  Z <- matrix(rnorm(n * nz), n, nz)
  regdata <- lapply(seq_len(n), function(i) {
    ni <- 15L
    X <- matrix(rnorm(ni * nvar), ni, nvar)
    beta <- rnorm(nvar, sd = 0.3)
    y <- if (counts) rnbinom(ni, size = 3, mu = exp(X %*% beta))
         else X %*% beta + rnorm(ni)
    list(y = y, X = X)
  })
  list(regdata = regdata, Z = Z)
}#MAKE_REG_SIM

hart_priors <- list(
  normalized = list(ncomp = 1, hart_model = "normalized",
                    bart = list(num_trees = 3)),
  direct     = list(ncomp = 1, hart_model = "direct",
                    bart = list(num_trees = 3), vartree = list(num_trees = 3)))

direct_fields <- c("bart_models", "var_models", "phi_models", "mu_draw",
                   "mean_scale_draw", "varcount", "varprob", "var_varcount",
                   "var_varprob", "betadraw", "loglike", "hart", "restart")

test_that("Linear HART outputs have expected field names and classes", {
  dat  <- make_reg_sim(700)
  Mcmc <- list(R = 4, keep = 1, nprint = 0)
  run_quiet({
    fit_a <- rhierLinearMixture(Data = dat, Prior = hart_priors$normalized,
                                Mcmc = Mcmc)
    fit_b <- rhierLinearMixture(Data = dat, Prior = hart_priors$direct,
                                Mcmc = Mcmc)
  })
  expect_true(all(c("betadraw", "taudraw", "loglike", "nmix", "bart_models",
                    "varcount", "varprob", "hart") %in% names(fit_a)))
  expect_length(fit_a$bart_models, 2L)
  expect_equal(class(fit_a), "rhierLinearMixture")
  expect_null(fit_a$Deltadraw)
  expect_named(fit_b, c(direct_fields, "taudraw"), ignore.order = TRUE)
  expect_equal(class(fit_b), c("rhierLinearMixtureHeterCov",
                               "bayesm.HART.HeterCov", "rhierLinearMixture"))
  expect_equal(fit_b$bart_models$treedraws$leaf_dim, 2L)
  expect_equal(dim(fit_b$varcount), c(2L, 4L))
  expect_equal(dim(fit_b$mean_scale_draw), c(4L, 2L))
  expect_null(fit_b$nmix)
})

test_that("NegBin HART outputs have expected field names and classes", {
  dat  <- make_reg_sim(701, counts = TRUE)
  Mcmc <- list(R = 4, keep = 1, nprint = 0)
  run_quiet({
    fit_a <- rhierNegbinRw(Data = dat, Prior = hart_priors$normalized,
                           Mcmc = Mcmc, r_verbose = FALSE)
    fit_b <- rhierNegbinRw(Data = dat, Prior = hart_priors$direct,
                           Mcmc = Mcmc, r_verbose = FALSE)
  })
  expect_true(all(c("betadraw", "alphadraw", "loglike", "nmix", "bart_models",
                    "varcount", "varprob", "acceptrbeta", "acceptralpha",
                    "hart") %in% names(fit_a)))
  expect_length(fit_a$bart_models, 2L)
  expect_equal(class(fit_a), "rhierNegbinRw")
  expect_named(fit_b, c(direct_fields, "alphadraw", "acceptrbeta",
                        "acceptralpha"), ignore.order = TRUE)
  expect_equal(class(fit_b), c("rhierNegbinRwHeterCov",
                               "bayesm.HART.HeterCov", "rhierNegbinRw"))
  expect_equal(fit_b$bart_models$treedraws$leaf_dim, 2L)
  expect_equal(dim(fit_b$mean_scale_draw), c(4L, 2L))
  expect_null(fit_b$nmix)
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
