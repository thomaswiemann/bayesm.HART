library(testthat)

.fit_mnl_for_cache_tests <- function(use_heter = FALSE, store_trees = TRUE) {
  set.seed(20260515)
  sim <- sim_hier_mnl(
    nlgt = 24, nT = 4, p = 3, nz = 2, nXa = 1, nXd = 0, const = TRUE,
    beta_func_type = "linear"
  )
  prior <- list(
    ncomp = 1L,
    bart = list(num_trees = 5L),
    store_trees = store_trees
  )
  if (use_heter) prior$vartree <- list(num_trees = 4L)
  mcmc <- list(R = 12L, keep = 1L, nprint = 0L)
  capture.output({
    fit <- rhierMnlRwMixture(Data = sim, Prior = prior, Mcmc = mcmc, r_verbose = FALSE)
  })
  list(sim = sim, fit = fit)
}

test_that("store_trees=FALSE supports seen-Z DeltaZ cache and guards unseen rows", {
  obj <- .fit_mnl_for_cache_tests(use_heter = FALSE, store_trees = FALSE)
  fit <- obj$fit
  Z <- obj$sim$Z

  expect_null(fit$bart_models)
  pred_seen <- predict(fit, newdata = list(Z = Z), type = "DeltaZ", r_verbose = FALSE)
  expect_true(is.array(pred_seen))
  expect_equal(dim(pred_seen)[1], nrow(Z))

  Z_unseen <- Z
  Z_unseen[1, ] <- Z_unseen[1, ] + 100
  expect_error(
    predict(fit, newdata = list(Z = Z_unseen), type = "DeltaZ", r_verbose = FALSE),
    regexp = "Out-of-sample Z prediction requires stored tree objects"
  )
})

test_that("DeltaZ cache path matches force_tree_eval on seen rows", {
  obj <- .fit_mnl_for_cache_tests(use_heter = FALSE, store_trees = TRUE)
  fit <- obj$fit
  Z <- obj$sim$Z

  pred_cache <- predict(fit, newdata = list(Z = Z), type = "DeltaZ", r_verbose = FALSE)
  pred_tree <- predict(
    fit, newdata = list(Z = Z), type = "DeltaZ",
    force_tree_eval = TRUE, r_verbose = FALSE
  )
  expect_equal(pred_cache, pred_tree, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("SigmaZ cache path matches force_tree_eval on seen rows", {
  obj <- .fit_mnl_for_cache_tests(use_heter = TRUE, store_trees = TRUE)
  fit <- obj$fit
  Z <- obj$sim$Z

  sigma_cache <- predict(fit, newdata = list(Z = Z), type = "SigmaZ",
                         burn = 1L, r_verbose = FALSE)
  sigma_tree <- predict(
    fit, newdata = list(Z = Z), type = "SigmaZ", burn = 1L,
    force_tree_eval = TRUE, r_verbose = FALSE
  )
  expect_equal(sigma_cache, sigma_tree, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("row order is preserved for permuted and mixed seen/unseen Z", {
  obj <- .fit_mnl_for_cache_tests(use_heter = TRUE, store_trees = TRUE)
  fit <- obj$fit
  Z <- obj$sim$Z

  idx <- c(3L, 1L, 5L, 2L)
  Z_perm <- Z[idx, , drop = FALSE]
  sigma_perm <- predict(fit, newdata = list(Z = Z_perm), type = "SigmaZ",
                        burn = 1L, r_verbose = FALSE)
  sigma_full <- predict(fit, newdata = list(Z = Z), type = "SigmaZ",
                        burn = 1L, r_verbose = FALSE)
  expect_equal(sigma_perm, sigma_full[idx, , , , drop = FALSE],
               tolerance = 1e-6, ignore_attr = TRUE)

  Z_mix <- rbind(Z[2, , drop = FALSE], Z[4, , drop = FALSE] + 5, Z[6, , drop = FALSE])
  delta_mix <- predict(fit, newdata = list(Z = Z_mix), type = "DeltaZ", r_verbose = FALSE)
  delta_mix_tree <- predict(
    fit, newdata = list(Z = Z_mix), type = "DeltaZ",
    force_tree_eval = TRUE, r_verbose = FALSE
  )
  expect_equal(delta_mix, delta_mix_tree, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("heter prior choice_probs uses SigmaZ cache on seen Z without trees", {
  obj <- .fit_mnl_for_cache_tests(use_heter = TRUE, store_trees = FALSE)
  fit <- obj$fit
  Z <- obj$sim$Z
  p <- obj$sim$p

  expect_null(fit$bart_models)
  expect_null(fit$var_models)

  X_template <- obj$sim$lgtdata[[1]]$X
  newdata_seen <- list(
    Z = Z,
    p = p,
    X = replicate(nrow(Z), X_template, simplify = FALSE)
  )
  pp_seen <- predict(
    fit, newdata = newdata_seen, mode = "prior", type = "choice_probs",
    nsim = 2L, burn = 1L, r_verbose = FALSE
  )
  expect_type(pp_seen, "list")
  expect_length(pp_seen, nrow(Z))
  expect_true(is.array(pp_seen[[1]]))

  Z_unseen <- Z
  Z_unseen[1, ] <- Z_unseen[1, ] + 7
  newdata_unseen <- list(
    Z = Z_unseen,
    p = p,
    X = replicate(nrow(Z_unseen), X_template, simplify = FALSE)
  )
  expect_error(
    predict(
      fit, newdata = newdata_unseen, mode = "prior", type = "choice_probs",
      nsim = 2L, burn = 1L, r_verbose = FALSE
    ),
    regexp = "stored tree objects for unseen Z"
  )
})
