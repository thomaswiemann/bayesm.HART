# Fresh-fit restart metadata, resolved priors and unit ordering.

library(testthat)

for (case in .guard_cases()) {
  if (case$likelihood != "linear" && case$model != "heteroskedastic") next
  test_that(paste("fresh restart metadata:", case$name), {
    data <- .guard_data(case$likelihood)
    units <- if (case$likelihood == "mnl") "lgtdata" else "regdata"
    data$unit_id <- seq_along(data[[units]]) + 100L
    before <- serialize(data, NULL)
    prior <- .guard_priors[[case$model]]
    fit <- .guard_fit(case$likelihood, data, prior, case$seed)
    r <- fit$restart
    expect_named(r, c("schema", "likelihood", "data", "config", "mcmc",
                      "preparation", "state"))
    expect_identical(r$schema, 1L)
    expect_identical(r$likelihood, case$likelihood)
    expect_false("state" %in% names(fit))
    st <- r$state
    expect_named(st, c("completed_sweeps", "betas",
      if (case$likelihood == "linear") "tau"
      else if (case$likelihood == "negbin") "alpha", "backend"))
    expect_equal(st$completed_sweeps, 30L)
    expect_equal(st$betas, unclass(fit$betadraw)[, , 30L])
    if (case$likelihood == "linear")
      expect_equal(st$tau, as.numeric(fit$taudraw[30L, ]))
    if (case$likelihood == "negbin")
      expect_equal(st$alpha, as.numeric(fit$alphadraw[30L]))
    if (case$model == "none") {
      expect_named(st$backend, c("delta", "prob", "ind"))
      expect_equal(dim(st$backend$delta), c(ncol(data[[units]][[1]]$X), ncol(data$Z)))
      expect_length(st$backend$ind, length(data[[units]]))
      expect_equal(st$backend$prob, 1)
    } else if (case$model == "normalized") {
      expect_named(st$backend, c("mu", "L", "mean"))
      expect_length(st$backend$mean, ncol(data[[units]][[1]]$X))
      expect_true(all(vapply(st$backend$mean, function(x) nzchar(x$trees), logical(1))))
    } else {
      expect_named(st$backend, c("mu", "tau2", "mean", "covariance"))
      expect_true(nzchar(st$backend$mean$trees))
      if (case$model == "constant") {
        expect_named(st$backend$covariance, "Sigma")
        expect_true(isSymmetric(st$backend$covariance$Sigma))
      } else {
        expect_named(st$backend$covariance, c("variance", "phi"))
        expect_length(st$backend$covariance$variance, ncol(data[[units]][[1]]$X))
        expect_equal(lengths(st$backend$covariance$phi),
                     seq_len(ncol(data[[units]][[1]]$X)) - 1L)
      }
    }
    expect_identical(r$data$unit_id, as.character(data$unit_id))
    expect_identical(r$data[[units]],
                     lapply(data[[units]], function(u) list(y = u$y, X = u$X)))
    expect_identical(r$data$Z, data$Z)
    expect_identical(serialize(data, NULL), before)
    expect_equal(r$config$model, prior$hart_model %||% "none")
    expect_equal(r$config$cov, if (case$model %in% c("constant", "heteroskedastic"))
      prior$hart_cov else NULL)
    expect_equal(r$config$nvar, ncol(data[[units]][[1]]$X))
    expect_equal(r$config$nz, ncol(data$Z))
    expect_equal(r$config$ncomp, 1L)
    expect_true(r$config$drawdelta)
    expect_identical(r$config$coef_names, colnames(data[[units]][[1]]$X))
    expect_named(r$config$prior, c("mubar", "Amu", "nu", "V", "Ad", "deltabar",
      if (case$likelihood == "mnl") c("a", "SignRes")
      else if (case$likelihood == "linear") c("a", "nu.e", "ssq")
      else c("a", "a_mix", "b")))
    expect_named(r$mcmc, c("keep", "nprint",
      if (case$likelihood == "mnl") c("s", "w")
      else if (case$likelihood == "negbin") c("s_alpha", "s_beta", "w")))
    expect_equal(r$mcmc$keep, 1L)
    expect_equal(r$mcmc$nprint, 0L)
    for (name in c("bart", "vartree", "phitree")) {
      expect_false(any(c("z_index", "n_unique") %in% names(r$config[[name]])))
    }
    if (case$model != "none") {
      expect_equal(r$config$bart$num_trees, 5L)
      expect_true(r$config$bart$store_trees)
    }
    if (case$likelihood == "mnl") expect_identical(r$data$p, data$p)
    if (case$likelihood == "linear") {
      expect_identical(r$preparation, list())
      expect_equal(r$config$prior$ssq, sapply(data$regdata, function(u) var(u$y)))
    } else {
      expect_named(r$preparation, c("hess",
        if (case$likelihood == "negbin") "alphacroot"))
      expect_length(r$preparation$hess, length(data[[units]]))
      expect_true(all(vapply(r$preparation$hess, function(h)
        is.matrix(h) && all(dim(h) == r$config$nvar) && all(is.finite(h)), logical(1))))
      if (case$likelihood == "negbin") {
        expect_true(is.finite(r$preparation$alphacroot))
        expect_false(r$config$fixalpha)
        expect_false("alpha_fixed" %in% names(r$config))
      }
    }
    path <- tempfile(fileext = ".rds")
    saveRDS(r, path)
    expect_identical(readRDS(path), r)
    unlink(path)
  })
}

test_that("unit IDs use explicit IDs, complete names, then positions", {
  units <- list(list(y = 1, X = matrix(1)), list(y = 2, X = matrix(2)))
  data <- list(regdata = units)
  expect_identical(.hart_restart_data(data, "regdata")$unit_id, c("1", "2"))
  names(data$regdata) <- c("b", "a")
  expect_identical(.hart_restart_data(data, "regdata")$unit_id, c("b", "a"))
  data$unit_id <- c(8L, 3L)
  expect_identical(.hart_restart_data(data, "regdata")$unit_id, c("8", "3"))
  data$unit_id <- NULL
  for (nm in list(c("", "a"), c("a", "a"), c(NA_character_, "a"))) {
    names(data$regdata) <- nm
    expect_identical(.hart_restart_data(data, "regdata")$unit_id, c("1", "2"))
  }
})

for (likelihood in c("mnl", "linear", "negbin")) {
  test_that(paste("unit ID validation and absent Z:", likelihood), {
    data <- .guard_data(likelihood)
    units <- if (likelihood == "mnl") "lgtdata" else "regdata"
    n <- length(data[[units]])
    for (ids in list(seq_len(n - 1L), c(NA, seq_len(n - 1L)), rep("x", n))) {
      data$unit_id <- ids
      expect_error(.guard_fit(likelihood, data, .guard_priors$none, 42L), "Data\\$unit_id")
    }
    data$unit_id <- NULL
    data$Z <- NULL
    names(data[[units]]) <- paste0("unit", rev(seq_len(n)))
    fit <- .guard_fit(likelihood, data, .guard_priors$none, 42L)
    expect_identical(fit$restart$data$unit_id, names(data[[units]]))
    expect_true("Z" %in% names(fit$restart$data))
    expect_null(fit$restart$data$Z)
    expect_false(fit$restart$config$drawdelta)
    expect_equal(fit$restart$config$nz, 0L)
    expect_equal(fit$restart$state$completed_sweeps, 30L)
    expect_equal(fit$restart$state$betas, unclass(fit$betadraw)[, , 30L])
    expect_named(fit$restart$state$backend, c("delta", "prob", "ind"))
  })
}

test_that("linear restart retains variance defaults and supplied scales", {
  data <- .guard_data("linear")
  data$regdata[[1]]$y[] <- 2
  data$regdata[[2]]$y <- data$regdata[[2]]$y[1]
  data$regdata[[2]]$X <- data$regdata[[2]]$X[1, , drop = FALSE]
  fit <- .guard_fit("linear", data, .guard_priors$none, 42L)
  expected <- sapply(data$regdata, function(u) {
    v <- var(u$y)
    if (is.na(v) || v <= 0) 1 else v
  })
  expect_equal(fit$restart$config$prior$ssq, expected)
  expect_equal(fit$restart$config$prior$nu.e, 3)
  prior <- .guard_priors$none
  prior$ssq <- rep(7, length(data$regdata))
  fit <- .guard_fit("linear", data, prior, 42L)
  expect_identical(fit$restart$config$prior$ssq, prior$ssq)
})

test_that("MNL restart retains final sign-dependent priors", {
  data <- .guard_data("mnl")
  nvar <- ncol(data$lgtdata[[1]]$X)
  prior <- .guard_priors$heteroskedastic
  prior$SignRes <- c(-1, rep(0, nvar - 1L))
  fit <- .guard_fit("mnl", data, prior, 42L)
  r <- fit$restart$config
  expect_equal(r$prior$nu, nvar + 15)
  expect_equal(r$prior$V, (nvar + 15) * diag(c(.2, rep(1, nvar - 1L))))
  expect_equal(r$prior$mubar, matrix(2 * abs(prior$SignRes), nrow = 1))
  expect_equal(r$prior$Amu, .1 * diag(nvar))
  expect_identical(r$prior$SignRes, prior$SignRes)
  expect_equal(r$vartree$v0, c(.1, rep(1, nvar - 1L)))
  internal <- unclass(fit$betadraw)[, , 30L]
  internal[, prior$SignRes != 0] <- log(abs(internal[, prior$SignRes != 0]))
  expect_equal(fit$restart$state$betas, internal, tolerance = 1e-12)
})

test_that("NB restart distinguishes fixed and sampled alpha", {
  data <- .guard_data("negbin")
  for (fixed in c(TRUE, FALSE)) {
    mcmc <- .guard_mcmc$negbin
    mcmc$alpha <- 2.5
    mcmc$fixalpha <- fixed
    before <- serialize(list(data, mcmc), NULL)
    capture.output(fit <- rhierNegbinRw(data, .guard_priors$none, mcmc,
                                         r_verbose = FALSE))
    expect_identical(fit$restart$config$fixalpha, fixed)
    if (fixed) expect_identical(fit$restart$config$alpha_fixed, 2.5)
    else expect_false("alpha_fixed" %in% names(fit$restart$config))
    expect_equal(fit$restart$state$alpha, as.numeric(fit$alphadraw[30L]))
    expect_equal(fit$restart$config$prior$a, .5)
    expect_equal(fit$restart$config$prior$b, .1)
    expect_equal(fit$restart$config$prior$a_mix, 5)
    expect_identical(serialize(list(data, mcmc), NULL), before)
  }
})

for (case in .guard_cases()) {
  if (case$likelihood != "linear" && case$model != "heteroskedastic") next
  test_that(paste("terminal state ignores thinning and tree history:", case$name), {
    data <- .guard_data(case$likelihood)
    prior <- .guard_priors[[case$model]]
    if (case$model != "none") prior$store_trees <- FALSE
    if (case$likelihood == "mnl") prior$SignRes <- c(-1, 0)
    sampler <- switch(case$likelihood, mnl = rhierMnlRwMixture,
                      linear = rhierLinearMixture, negbin = rhierNegbinRw)
    set.seed(case$seed)
    capture.output(fit <- sampler(data, prior, list(R = 11L, keep = 3L, nprint = 0L),
                                   r_verbose = FALSE))
    expect_equal(fit$restart$state$completed_sweeps, 11L)
    last <- unclass(fit$betadraw)[, , 3L]
    if (case$likelihood == "mnl") last[, 1] <- log(abs(last[, 1]))
    expect_false(isTRUE(all.equal(fit$restart$state$betas, last)))
    set.seed(case$seed)
    capture.output(full <- sampler(data, prior, list(R = 11L, keep = 1L, nprint = 0L),
                                    r_verbose = FALSE))
    expect_equal(fit$restart$state, full$restart$state)
    if (case$model != "none") {
      expect_null(fit$bart_models)
      mean <- fit$restart$state$backend$mean
      if (case$model == "normalized") {
        expect_true(all(vapply(mean, function(x) nzchar(x$trees), logical(1))))
      } else {
        expect_true(nzchar(mean$trees))
      }
      if (case$model == "heteroskedastic") {
        cov <- fit$restart$state$backend$covariance
        expect_null(fit$var_models)
        expect_null(fit$phi_models)
        expect_true(all(vapply(cov$variance, function(x) nzchar(x$trees), logical(1))))
        expect_true(nzchar(cov$phi[[2]][[1]]$trees))
      }
    }
  })
}

test_that("terminal DART flags follow each ensemble's burn", {
  data <- .guard_data("linear")
  for (model in c("normalized", "heteroskedastic")) {
    prior <- .guard_priors[[model]]
    prior$bart$sparse <- TRUE
    prior$bart$burn <- 2L
    if (model == "heteroskedastic") {
      prior$vartree$sparse <- TRUE
      prior$vartree$burn <- 4L
      prior$phitree$sparse <- TRUE
      prior$phitree$burn <- 6L
    }
    for (R in 3:8) {
      set.seed(412)
      capture.output(fit <- rhierLinearMixture(data, prior,
        list(R = R, keep = 1L, nprint = 0L), r_verbose = FALSE))
      b <- fit$restart$state$backend
      if (model == "normalized") {
        expect_identical(vapply(b$mean, function(x) x$dartOn, logical(1)),
                         rep(R > 3L, length(b$mean)))
      } else {
        expect_identical(b$mean$dartOn, R > 3L)
        expect_identical(vapply(b$covariance$variance, function(x) x$dartOn, logical(1)),
                         rep(R > 5L, length(b$covariance$variance)))
        expect_identical(b$covariance$phi[[2]][[1]]$dartOn, R > 7L)
      }
    }
  }
})
