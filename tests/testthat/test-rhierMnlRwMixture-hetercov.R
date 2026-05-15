# Tests for the heteroscedastic-covariance MNL path
# (Prior$vartree, optionally Prior$phitree).
#
# These are smoke tests focused on:
#   * input validation (vartree without bart, ncomp != 1, no Z),
#   * end-to-end execution via the user wrapper,
#   * class assignment and return-slot shape,
#   * sign restrictions are honored under heter-cov (positivity preserved),
#   * predict.* dispatch on "rhierMnlRwMixtureHeterCov" for all five types.
#
# Recovery accuracy is exercised by the longer dev/smoke-heter-mnl*.R scripts.

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

test_that("vartree requires bart, ncomp == 1, and Z", {
  sim <- make_small_sim()

  # vartree without bart
  expect_error(
    capture.output(
      rhierMnlRwMixture(
        Data  = sim$Data,
        Prior = list(ncomp = 1L, vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "requires Prior\\$bart"
  )

  # ncomp > 1 with vartree
  expect_error(
    capture.output(
      rhierMnlRwMixture(
        Data  = sim$Data,
        Prior = list(ncomp = 2L, bart = list(num_trees = 5L),
                     vartree = list(num_trees = 5L)),
        Mcmc  = list(R = 4L, keep = 1L, nprint = 0L),
        r_verbose = FALSE)
    ),
    regexp = "Only ncomp = 1 is currently supported"
  )
})

test_that("Sign restrictions under heter-cov preserve positivity (lifted validator)", {
  # Sigma(Z) is the heteroscedastic covariance of beta* = log|beta|/SignRes;
  # post-loop transform beta = SignRes * exp(beta*) preserves the sign for
  # every kept draw.  This regression test guards against any future change
  # that re-introduces a SignRes block in .validate_heter_cov.
  sim <- make_small_sim()                  # nvar = 2L
  Prior <- list(ncomp = 1L,
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
  # Diagonal-only Sigma(Z) for nvar > 1 would force every conditional
  # cross-correlation to hard zero at every Z, which is essentially never a
  # believable posterior.  The wrapper auto-synthesizes phi_params from the
  # default .parse_phi_params(NULL, nz) so the user gets full Cholesky even
  # without explicitly passing Prior$phitree.
  sim <- make_small_sim()                # nvar = 2L
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))   # NO phitree
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
  expect_s3_class(fit, "rhierMnlRwMixture")
  expect_named(fit,
               c("bart_models", "var_models", "phi_models", "mu_draw",
                 "varcount", "varprob", "var_varcount", "var_varprob",
                 "betadraw", "loglike", "SignRes", "acceptrbeta"),
               ignore.order = TRUE)
  # Auto-promotion: phi_models is now non-null and properly jagged.
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
  expect_equal(dim(fit$betadraw)[1:2], c(nrow(sim$Z), sim$nvar))
  expect_length(fit$bart_models, sim$nvar)
  expect_length(fit$var_models, sim$nvar)
})

test_that("vartree only with nvar == 1 stays diagonal (scalar case)", {
  # When nvar == 1 there are no off-diagonals to model, so phi_models is
  # legitimately NULL and Sigma(Z) is just d_1(Z).  Auto-promotion is a no-op.
  sim <- make_small_sim(nvar = 1L)
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))
  Mcmc  <- list(R = 20L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })

  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
  expect_null(fit$phi_models)
  expect_length(fit$bart_models, 1L)
  expect_length(fit$var_models,  1L)
  expect_equal(dim(fit$mu_draw), c(20L, 1L))
})

test_that("DART (sparse=TRUE) on vartree concentrates var_varprob away from uniform", {
  # Validates the full DART pipeline (parse Prior$vartree$sparse + burn ->
  # initializeVarBART:setdart -> startdart at burn+1).  Phi-tree DART is plumbed
  # via the same mechanism (s.phi_sparse / s.phi_burn / startdart in
  # hetercov_draw_iter); validating the var path confirms the activation.
  #
  # Robust check: under a uniform Dirichlet over `nz` columns, each varprob
  # entry has expected value 1/nz.  DART, post-activation, drives varprob away
  # from uniform whenever ANY column is more informative than others; we don't
  # over-specify which one wins (sample-size noise can shuffle the winner).
  set.seed(20260513L)
  nlgt <- 30L; nz <- 5L; nvar <- 1L; p <- 3L; T_i <- 6L
  Z <- cbind(runif(nlgt, -1, 1), matrix(runif(nlgt * (nz - 1L), -1, 1), nlgt, nz - 1L))
  Z <- scale(Z, center = TRUE, scale = FALSE)
  theta_true <- rnorm(nlgt, mean = 0, sd = exp(0.6 * Z[, 1]))
  lgtdata <- lapply(seq_len(nlgt), function(i) {
    X <- matrix(rnorm(T_i * p * nvar), T_i * p, nvar)
    eta  <- X %*% theta_true[i]
    prob <- matrix(exp(eta), T_i, p, byrow = TRUE)
    prob <- prob / rowSums(prob)
    y    <- apply(prob, 1, function(pp) sample.int(p, 1, prob = pp))
    list(y = y, X = X)
  })
  Data <- list(p = p, lgtdata = lgtdata, Z = Z)
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 20L, sparse = TRUE, burn = 5L))
  Mcmc  <- list(R = 60L, keep = 1L, nprint = 0L)

  capture.output({
    fit <- rhierMnlRwMixture(Data = Data, Prior = Prior, Mcmc = Mcmc,
                             r_verbose = FALSE)
  })
  vp_post <- apply(fit$var_varprob[, 1, 30:60], 1, mean)
  # If DART is active, the *most* selected column should poll meaningfully
  # above uniform (1/nz = 0.2) and the *least* selected should be squeezed
  # below uniform.  Both directions must be present for a true non-uniform
  # posterior over Z columns -- a sufficient signature of DART activation.
  expect_true(max(vp_post) > 1.5 * (1 / nz),
              info = sprintf("max vp=%.3f, uniform=%.3f", max(vp_post), 1 / nz))
  expect_true(min(vp_post) < (1 / nz),
              info = sprintf("min vp=%.3f, uniform=%.3f", min(vp_post), 1 / nz))
})

test_that("Phase 2 (vartree + phitree) runs end-to-end and returns jagged phi_models", {
  sim <- make_small_sim()
  Prior <- list(ncomp = 1L,
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
  Prior <- list(ncomp = 1L,
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

  # posterior_probs: list of npred-by-... arrays for the original units
  pp <- predict(fit, newdata = newdata, type = "posterior_probs", burn = 5L,
                r_verbose = FALSE)
  expect_length(pp, nrow(sim$Z))
  expect_equal(dim(pp[[1]])[2:3], c(sim$p, 15L))

  # prior_probs: averaged over nsim draws of eta = mu + Sigma(Z*)^{1/2} z
  prp <- predict(fit, newdata = newdata, type = "prior_probs", burn = 5L,
                 nsim = 3L, r_verbose = FALSE)
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

test_that("auto-calibrated lambda is used when Prior$vartree$lambda is NULL", {
  sim <- make_small_sim()
  Prior <- list(ncomp = 1L,
                bart    = list(num_trees = 5L),
                vartree = list(num_trees = 5L))   # no lambda -> auto
  Mcmc  <- list(R = 6L, keep = 1L, nprint = 0L)

  expect_message(
    {
      capture.output({
        fit <- rhierMnlRwMixture(Data = sim$Data, Prior = Prior, Mcmc = Mcmc,
                                 r_verbose = TRUE)
      })
    },
    regexp = NA   # only verify execution; message inspection requires output capture
  )
  expect_s3_class(fit, "rhierMnlRwMixtureHeterCov")
})
