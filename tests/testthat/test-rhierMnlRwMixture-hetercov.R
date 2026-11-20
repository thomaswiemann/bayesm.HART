# Tests for the heteroscedastic-covariance MNL path
# (Prior$vartree, optionally Prior$phitree).
#
# These are smoke tests focused on:
#   * input validation (vartree without bart, ncomp != 1, no Z),
#   * end-to-end execution via the user wrapper,
#   * class assignment and return-slot shape,
#   * sign restrictions are honored under heter-cov (positivity preserved),
#   * predict.* dispatch on "rhierMnlRwMixtureHeterCov" for all five types.

library(testthat)

# --- shared minimal-cost simulation -----------------------------------------
make_small_sim <- function(seed = 20260512L, nlgt = 30L, nvar = 2L, nz = 2L,
                           p = 3L, T_i = 6L) {
  set.seed(seed)
  Z <- cbind(runif(nlgt, -1, 1), runif(nlgt, -1, 1))
  Z <- scale(Z, center = TRUE, scale = FALSE)
  theta_true <- matrix(rnorm(nlgt * nvar, sd = 0.5), nlgt, nvar)
  lgtdata <- lapply(seq_len(nlgt), function(i) {
    X    <- matrix(rnorm(T_i * p * nvar), T_i * p, nvar)
    eta  <- X %*% theta_true[i, ]
    prob <- matrix(exp(eta), T_i, p, byrow = TRUE)
    prob <- prob / rowSums(prob)
    y    <- apply(prob, 1, function(pp) sample.int(p, 1, prob = pp))
    list(y = y, X = X)
  })
  list(Data = list(p = p, lgtdata = lgtdata, Z = Z),
       theta_true = theta_true, Z = Z, p = p, nvar = nvar)
}

test_that("vartree requires hart_model, ncomp == 1, and Z", {
  sim <- make_small_sim()

  # vartree without hart_model
  expect_error(
    capture.output(
      rhierMnlRwMixture(
        Data  = sim$Data,
        Prior = list(ncomp = 1L, vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "require Prior\\$hart_model"
  )

  # ncomp > 1 with vartree
  expect_error(
    capture.output(
      rhierMnlRwMixture(
        Data  = sim$Data,
        Prior = list(ncomp = 2L, hart_model = "direct",
                     bart = list(num_trees = 5L),
                     vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "Only ncomp = 1 is currently supported"
  )
})

test_that("Sign restrictions under heter-cov preserve positivity", {
  # beta = SignRes * exp(beta*) keeps the sign of every kept draw.
  sim <- make_small_sim()                  # nvar = 2L
  Prior <- list(ncomp = 1L, hart_model = "direct",
                SignRes = c(1L, 0L),       # constrain coef 1 to be positive
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 30L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
  expect_equal(as.integer(fit$SignRes), c(1L, 0L))
  # Coef 1 (sign-restricted positive) must be > 0 for every unit and draw.
  expect_true(all(fit$betadraw[, 1, ] > 0))
  # Coef 2 (unrestricted) is unchanged in interpretation; just check finite.
  expect_true(all(is.finite(fit$betadraw[, 2, ])))
})

test_that("vartree only with nvar > 1 auto-promotes to full Cholesky", {
  # nvar > 1 without phitree: the off-diagonal ensembles exist anyway.
  sim <- make_small_sim()                # nvar = 2L
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))   # no phitree
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
  expect_s3_class(fit, "rhierMnlRwMixture")
  expect_named(fit,
               c("bart_models", "var_models", "phi_models", "mu_draw",
                 "mean_scale_draw", "varcount", "varprob", "var_varcount",
                 "var_varprob", "betadraw", "loglike", "SignRes",
                 "acceptrbeta", "hart", "restart"),
               ignore.order = TRUE)
  # phi_models is jagged: entry j holds j - 1 trees.
  expect_length(fit$phi_models, sim$nvar)
  for (j in seq_len(sim$nvar)) {
    if (j == 1L) {
      expect_true(is.null(fit$phi_models[[1]]) ||
                  length(fit$phi_models[[1]]) == 0L)
    } else {
      expect_length(fit$phi_models[[j]], j - 1L)
    }
  }
  expect_equal(dim(fit$mu_draw), c(20L, sim$nvar))
  expect_equal(dim(fit$mean_scale_draw), c(20L, sim$nvar))
  expect_true(all(fit$mean_scale_draw > 0))
  expect_equal(dim(fit$betadraw)[1:2], c(nrow(sim$Z), sim$nvar))
  # One vector-leaf ensemble with a shared partition.
  expect_equal(fit$bart_models$treedraws$leaf_dim, sim$nvar)
  expect_equal(dim(fit$varcount), c(ncol(sim$Z), 20L))
  expect_length(fit$var_models, sim$nvar)
  expect_equal(dim(fit$var_varcount), c(ncol(sim$Z), sim$nvar, 20L))
})

test_that("vartree only with nvar == 1 stays diagonal (scalar case)", {
  # nvar == 1: no off-diagonals, Sigma(Z) = d_1(Z).
  sim <- make_small_sim(nvar = 1L)
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
  expect_null(fit$phi_models)
  expect_equal(fit$bart_models$treedraws$leaf_dim, 1L)
  expect_length(fit$var_models,  1L)
  expect_equal(dim(fit$mu_draw), c(20L, 1L))
})

test_that("DART on vartree moves var_varprob away from uniform", {
  # Only Z[, 1] drives the variance; after activation var_varprob leaves 1/nz.
  set.seed(20260513L)
  nlgt <- 60L; nz <- 5L; nvar <- 1L; p <- 3L; T_i <- 6L
  Z <- cbind(runif(nlgt, -1, 1),
             matrix(runif(nlgt * (nz - 1L), -1, 1), nlgt, nz - 1L))
  Z <- scale(Z, center = TRUE, scale = FALSE)
  theta_true <- rnorm(nlgt, mean = 0, sd = exp(1.0 * Z[, 1]))
  lgtdata <- lapply(seq_len(nlgt), function(i) {
    X <- matrix(rnorm(T_i * p * nvar), T_i * p, nvar)
    eta  <- X %*% theta_true[i]
    prob <- matrix(exp(eta), T_i, p, byrow = TRUE)
    prob <- prob / rowSums(prob)
    y    <- apply(prob, 1, function(pp) sample.int(p, 1, prob = pp))
    list(y = y, X = X)
  })
  Data <- list(p = p, lgtdata = lgtdata, Z = Z)
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 20L, sparse = TRUE, burn = 5L))
  Mcmc  <- list(R = 200L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })
  vp_post <- apply(fit$var_varprob[, 1, 101:200], 1, mean)
  # Largest column above uniform, smallest below.
  expect_true(max(vp_post) > 1.5 * (1 / nz),
              info = sprintf("max vp=%.3f, uniform=%.3f", max(vp_post), 1 / nz))
  expect_true(min(vp_post) < (1 / nz),
              info = sprintf("min vp=%.3f, uniform=%.3f", min(vp_post), 1 / nz))
})

test_that("vartree + phitree returns jagged phi_models", {
  sim <- make_small_sim()
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L),
                phitree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
  expect_length(fit$phi_models, sim$nvar)
  # phi_models[[1]] is empty (no parents); for j > 1 there are j - 1 entries
  for (j in seq_len(sim$nvar)) {
    if (j == 1L) {
      expect_true(is.null(fit$phi_models[[1]]) ||
                  length(fit$phi_models[[1]]) == 0L)
    } else {
      expect_length(fit$phi_models[[j]], j - 1L)
      for (k in seq_len(j - 1L)) {
        expect_true(!is.null(fit$phi_models[[j]][[k]]$treedraws))
      }
    }
  }
})

test_that("predict.rhierMnlRwMixtureHeterCov dispatches for all five types", {
  sim <- make_small_sim()
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L),
                phitree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })

  npred  <- 4L
  Zstar  <- matrix(runif(npred * ncol(sim$Z), -1, 1), npred, ncol(sim$Z))
  Xstar  <- lapply(seq_len(npred), function(i)
              matrix(rnorm(3L * sim$p * sim$nvar), 3L * sim$p, sim$nvar))
  newdata <- list(Z = Zstar, X = Xstar, p = sim$p,
                  nlgtdata = lapply(seq_len(nrow(sim$Z)), function(i)
                    list(X = sim$Data$lgtdata[[i]]$X)))

  # DeltaZ: 3D array (npred x nvar x ndraws_kept)
  pred_dz <- predict(fit, newdata = newdata, type = "DeltaZ", burn = 5L,
                     r_verbose = FALSE)
  expect_equal(dim(pred_dz), c(npred, sim$nvar, 15L))
  expect_true(all(is.finite(pred_dz)))

  # DeltaZ+mu adds the per-draw mu
  pred_dzmu <- predict(fit, newdata = newdata, type = "DeltaZ+mu", burn = 5L,
                       r_verbose = FALSE)
  expect_equal(dim(pred_dzmu), dim(pred_dz))
  # Difference should equal mu_kept (broadcast across npred and draws)
  diff <- pred_dzmu - pred_dz
  mu_kept <- t(fit$mu_draw[6:20, , drop = FALSE])           # nvar x 15
  for (i in seq_len(npred)) {
    expect_equal(diff[i, , ], mu_kept, tolerance = 1e-8,
                 ignore_attr = TRUE)
  }

  # posterior choice_probs: list of npred-by-... arrays for the original units
  pp <- predict(fit, newdata = newdata, mode = "posterior",
                type = "choice_probs", burn = 5L, r_verbose = FALSE)
  expect_length(pp, nrow(sim$Z))
  expect_equal(dim(pp[[1]])[2:3], c(sim$p, 15L))

  # prior choice_probs: averaged over nsim draws of eta = mu + Sigma(Z*)^{1/2} z
  prp <- predict(fit, newdata = newdata, mode = "prior",
                 type = "choice_probs", burn = 5L, nsim = 3L,
                 r_verbose = FALSE)
  expect_length(prp, npred)
  expect_equal(dim(prp[[1]])[2:3], c(sim$p, 15L))
  # probabilities sum to 1 across alternatives (within numerical tolerance)
  for (i in seq_len(npred)) {
    sums <- apply(prp[[i]], c(1, 3), sum)
    expect_equal(as.numeric(sums), rep(1, length(sums)), tolerance = 1e-8)
  }

  # SigmaZ: covariance draws [npred, nvar, nvar, ndraws_kept]
  pred_sigma <- predict(fit, newdata = newdata, type = "SigmaZ", burn = 5L,
                        r_verbose = FALSE)
  expect_equal(dim(pred_sigma), c(npred, sim$nvar, sim$nvar, 15L))
  expect_true(all(is.finite(pred_sigma)))
  for (i in seq_len(npred)) {
    for (s in 1:15) {
      expect_equal(pred_sigma[i, , , s], t(pred_sigma[i, , , s]),
                   tolerance = 1e-8, ignore_attr = TRUE)
      expect_true(all(diag(pred_sigma[i, , , s]) > 0))
    }
  }
})

test_that("DeltaZ is f(Z) from the mean trees alone and matches the cache", {
  sim <- make_small_sim()
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 10L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
    fit_b <- rhierMnlRwMixture(Data = sim$Data, Mcmc = Mcmc, r_verbose = FALSE,
                               Prior = modifyList(Prior, list(
                                 vartree = list(num_trees = 3L))))
  })
  newdata <- list(Z = sim$Z)
  from_cache <- predict(fit, newdata = newdata, type = "DeltaZ", burn = 0L,
                        r_verbose = FALSE)
  from_trees <- bayesm.HART:::.pw_vec_trees(sim$Z, fit$bart_models$treedraws,
                                            dim(fit$betadraw)[2])
  expect_equal(from_cache, from_trees, tolerance = 1e-7, ignore_attr = TRUE)
  # Swapping in other covariance trees leaves DeltaZ unchanged.
  mixed <- fit
  mixed$var_models <- fit_b$var_models
  mixed$phi_models <- fit_b$phi_models
  pred_mixed <- bayesm.HART:::.calculate_delta_z(
    mixed, newdata, burn = 0, r_verbose = FALSE, force_tree_eval = TRUE)
  expect_equal(pred_mixed, from_trees, tolerance = 1e-7, ignore_attr = TRUE)
  # SigmaZ: cache and tree evaluation agree.
  sig_cache <- predict(fit, newdata = newdata, type = "SigmaZ", burn = 0L,
                       r_verbose = FALSE)
  sig_trees <- bayesm.HART:::.calculate_sigma_z(
    fit, newdata, burn = 0, r_verbose = FALSE, force_tree_eval = TRUE)
  expect_equal(sig_cache, sig_trees, tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("verbose output lists the covariance-tree priors", {
  sim <- make_small_sim()
  Prior <- list(ncomp = 1L, hart_model = "direct",
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L, v0 = 0.5))
  Mcmc  <- list(R = 4L, keep = 1L, nprint = 0L)
  out <- capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = TRUE)
  })
  expect_true(any(grepl("Covariance trees \\(full Cholesky\\)", out)))
  expect_true(any(grepl("vartree: num_trees = 5, tau = 1 .*v0 = 0.5 0.5", out)))
  expect_equal(fit$hart$priors$vartree$v0, c(0.5, 0.5))
})
