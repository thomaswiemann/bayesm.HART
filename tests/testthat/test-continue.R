# Unchanged-data continuation, terminal state and local draw storage.

library(testthat)

continue_sampler <- function(likelihood) {
  switch(likelihood, mnl = rhierMnlRwMixture,
         linear = rhierLinearMixture, negbin = rhierNegbinRw)
}#CONTINUE_SAMPLER

continue_trees <- function(stream, kept = NULL) {
  lines <- strsplit(stream, "\n", fixed = TRUE)[[1]]
  header <- scan(text = lines[1], quiet = TRUE)
  pos <- 2L
  out <- vector("list", header[1])
  for (i in seq_len(header[1])) {
    out[[i]] <- vector("list", header[2])
    for (j in seq_len(header[2])) {
      n <- as.integer(lines[pos])
      nodes <- do.call(rbind, lapply(lines[pos + seq_len(n)], function(x)
        scan(text = x, quiet = TRUE)))
      out[[i]][[j]] <- list(topology = nodes[, 1:3, drop = FALSE],
                            leaves = nodes[, -(1:3), drop = FALSE])
      pos <- pos + n + 1L
    }
  }
  if (is.null(kept)) out else out[kept]
}#CONTINUE_TREES

continue_tree_models <- function(x, kept = NULL) {
  if (!is.list(x)) return(x)
  for (name in names(x)) {
    if (name == "trees") x[[name]] <- continue_trees(x[[name]], kept)
    else x[[name]] <- continue_tree_models(x[[name]], kept)
  }
  if (is.null(names(x))) x <- lapply(x, continue_tree_models, kept = kept)
  x
}#CONTINUE_TREE_MODELS

continue_draws <- function(fit, kept) {
  out <- list()
  for (name in c("betadraw", "varcount", "varprob", "var_varcount", "var_varprob")) {
    x <- fit[[name]]
    if (is.null(x)) next
    index <- lapply(dim(x), seq_len)
    index[[length(index)]] <- kept
    out[[name]] <- do.call(`[`, c(list(unclass(x)), index, list(drop = FALSE)))
  }
  for (name in c("loglike", "Deltadraw", "taudraw", "alphadraw", "mu_draw",
                 "mean_scale_draw")) {
    x <- fit[[name]]
    if (!is.null(x)) out[[name]] <- if (is.null(dim(x))) x[kept]
      else unclass(x)[kept, , drop = FALSE]
  }
  if (!is.null(fit$nmix)) {
    out$nmix <- list(probdraw = fit$nmix$probdraw[kept, , drop = FALSE],
                     compdraw = fit$nmix$compdraw[kept], zdraw = fit$nmix$zdraw)
  }
  for (name in c("bart_models", "var_models", "phi_models")) {
    if (!is.null(fit[[name]]))
      out[[name]] <- continue_tree_models(fit[[name]], kept)
  }
  cache <- attr(fit, "hart_cache")
  for (name in names(cache)) {
    x <- cache[[name]]
    if (grepl("draws", name)) {
      index <- lapply(dim(x), seq_len)
      index[[length(index)]] <- kept
      x <- do.call(`[`, c(list(x), index, list(drop = FALSE)))
    }
    out[[paste0("cache.", name)]] <- x
  }
  out
}#CONTINUE_DRAWS

continue_difference <- function(x, y, path = "") {
  if (!identical(dim(x), dim(y)) || length(x) != length(y))
    stop(paste(path, "dimension mismatch"))
  if (is.list(x) && is.list(y)) {
    if (!identical(names(x), names(y))) stop(paste(path, "field mismatch"))
    if (!length(x)) return(0)
    return(max(vapply(seq_along(x), function(i)
      continue_difference(x[[i]], y[[i]],
        paste0(path, "/", if (is.null(names(x))) i else names(x)[i])),
      numeric(1))))
  }
  if (is.numeric(x) && is.numeric(y)) {
    if (grepl("topology", path, fixed = TRUE) && !identical(x, y))
      stop(paste(path, "tree topology mismatch"))
    if (!length(x)) return(0)
    return(max(abs(x - y)))
  }
  if (!identical(x, y)) stop(paste(path, "value mismatch"))
  0
}#CONTINUE_DIFFERENCE

continue_split <- function(case, k = 8L, keep = 1L, prior = NULL,
                            no_z = FALSE, extra = list(), label = case$name, data = NULL) {
  if (is.null(data)) data <- .guard_data(case$likelihood)
  if (no_z) data$Z <- NULL
  if (is.null(prior)) prior <- .guard_priors[[case$model]]
  sampler <- continue_sampler(case$likelihood)
  controls <- c(list(R = 2L * k, keep = 1L, nprint = 0L), extra)
  set.seed(case$seed)
  capture.output(full <- sampler(data, prior, controls, r_verbose = FALSE))
  controls$R <- k
  controls$keep <- keep
  set.seed(case$seed)
  capture.output(a <- sampler(data, prior, controls, r_verbose = FALSE))
  before <- serialize(a, NULL)
  capture.output(b <- sampler(Mcmc = list(R = k), previous = a,
                              r_verbose = FALSE))
  expect_identical(serialize(a, NULL), before)
  expect_identical(b$restart$data, a$restart$data)
  expect_identical(b$restart$config, a$restart$config)
  expect_identical(class(b), class(a))
  expect_equal(b$restart$state$completed_sweeps, 2L * k)
  got <- continue_draws(b, seq_len(k %/% keep))
  want <- continue_draws(full, k + seq_len(k %/% keep) * keep)
  got$state <- continue_tree_models(b$restart$state)
  want$state <- continue_tree_models(full$restart$state)
  expect_identical(names(got), names(want))
  delta <- vapply(names(got), function(name)
    continue_difference(got[[name]], want[[name]], name), numeric(1))
  worst <- names(delta)[which.max(delta)]
  cat(sprintf("SPLIT\t%s\t%s\t%.17g\n", label, worst, max(delta)))
  expect_lte(max(delta), 1e-8, label = paste(label, worst))
  invisible(b)
}#CONTINUE_SPLIT

for (case in .guard_cases()) {
  if (case$likelihood != "linear" && case$model != "heteroskedastic") next
  test_that(paste("split run:", case$name), {
    continue_split(case)
  })
  if (case$likelihood == "linear") test_that(
    paste("unkept terminal sweep:", case$name), {
    prior <- .guard_priors[[case$model]]
    if (case$model != "none") prior$store_trees <- FALSE
    continue_split(case, k = 7L, keep = 3L, prior = prior,
                   label = paste0(case$name, "-thinned-no-history"))
  })
  if (case$likelihood == "linear" && case$model == "none") {
    test_that(paste("no Z split run:", case$name), {
      continue_split(case, no_z = TRUE, label = paste0(case$name, "-noZ"))
    })
  } else if (case$likelihood == "linear" &&
             case$model %in% c("normalized", "heteroskedastic")) {
    test_that(paste("DART split boundaries:", case$name), {
      prior <- .guard_priors[[case$model]]
      prior$bart$sparse <- TRUE
      prior$bart$burn <- 2L
      if (case$model == "heteroskedastic") {
        prior$vartree$sparse <- TRUE
        prior$vartree$burn <- 3L
        prior$phitree$sparse <- TRUE
        prior$phitree$burn <- 4L
      }
      for (k in 2:6) continue_split(case, k = k, prior = prior,
        label = paste0(case$name, "-dart-k", k))
    })
  }
  if (case$likelihood == "mnl") {
    test_that(paste("signed MNL split run:", case$model), {
      prior <- .guard_priors[[case$model]]
      prior$SignRes <- c(-1, 0)
      continue_split(case, prior = prior, label = paste0(case$name, "-signs"))
    })
  }
  if (case$likelihood == "negbin") {
    test_that(paste("fixed alpha split run:", case$model), {
      continue_split(case, extra = list(alpha = 2.5),
                     label = paste0(case$name, "-fixed"))
    })
  }
}

for (likelihood in c("linear", "negbin")) {
  test_that(paste("repeated, serialized continuation and prediction:", likelihood), {
    data <- .guard_data(likelihood)
    sampler <- continue_sampler(likelihood)
    set.seed(201)
    capture.output(a <- sampler(data, .guard_priors$heteroskedastic,
      list(R = 8L, keep = 2L, nprint = 0L), r_verbose = FALSE))
    path <- tempfile(fileext = ".rds")
    saveRDS(a, path)
    disk <- readRDS(path)
    unlink(path)
    before <- serialize(disk, NULL)
    capture.output(b <- sampler(Data = NULL, Mcmc = list(R = 6L),
      previous = disk, r_verbose = FALSE))
    expect_identical(serialize(disk, NULL), before)
    before_b <- serialize(b, NULL)
    capture.output(c <- sampler(Mcmc = list(R = 5L, keep = 1L),
      previous = b, r_verbose = FALSE))
    expect_identical(serialize(b, NULL), before_b)
    expect_equal(b$restart$state$completed_sweeps, 14L)
    expect_equal(c$restart$state$completed_sweeps, 19L)
    expect_equal(dim(b$betadraw)[3], 3L)
    expect_equal(dim(c$betadraw)[3], 5L)
    expect_identical(c$restart$data, a$restart$data)
    expect_identical(c$restart$config, a$restart$config)
    expect_identical(class(c), class(a))
    field <- if (likelihood == "linear") "taudraw" else "alphadraw"
    expect_equal(attr(b[[field]], "mcpar"), c(1, 6, 2))
    expect_equal(attr(c[[field]], "mcpar"), c(1, 5, 1))
    prediction <- predict(c, newdata = list(Z = data$Z), type = "DeltaZ+mu",
                           burn = 0L, r_verbose = FALSE)
    expect_true(all(is.finite(prediction)))
    expect_equal(dim(prediction)[3], 5L)
  })

  if (likelihood == "linear") test_that(
    paste("continuation argument errors:", likelihood), {
    sampler <- continue_sampler(likelihood)
    data <- .guard_data(likelihood)
    capture.output(a <- sampler(data, .guard_priors$none,
      list(R = 4L, keep = 1L, nprint = 0L), r_verbose = FALSE))
    expect_error(sampler(previous = a), "Requires R")
    expect_error(sampler(previous = a, Mcmc = list(keep = 1L)), "Requires R")
    old <- a
    old$restart <- NULL
    expect_error(sampler(previous = old, Mcmc = list(R = 4)), "missing restart")
    old <- a
    old$restart$state <- NULL
    expect_error(sampler(previous = old, Mcmc = list(R = 4)), "missing restart state")
    old <- a
    old$restart$schema <- 2L
    expect_error(sampler(previous = old, Mcmc = list(R = 4)), "schema 1")
    old$restart$schema <- 1L
    old$restart$likelihood <- "wrong"
    expect_error(sampler(previous = old, Mcmc = list(R = 4)), "wrong likelihood")
    expect_error(sampler(Data = data, previous = a, Mcmc = list(R = 4)), "Data")
    expect_error(sampler(Prior = list(ncomp = 1), previous = a,
                         Mcmc = list(R = 4)), "Prior")
  })
}

for (likelihood in c("mnl", "linear", "negbin")) {
  test_that(paste("continuation rejects proposal controls:", likelihood), {
    sampler <- continue_sampler(likelihood)
    data <- .guard_data(likelihood)
    capture.output(a <- sampler(data, .guard_priors$none,
      list(R = 4L, keep = 1L, nprint = 0L), r_verbose = FALSE))
    fields <- switch(likelihood, mnl = "s", linear = "w",
                      negbin = c("s_beta", "alpha"))
    for (field in fields) {
      controls <- list(R = 4L)
      controls[[field]] <- 1
      expect_error(sampler(previous = a, Mcmc = controls),
                   "With previous, Mcmc accepts only R, keep and nprint")
    }
  })
}

test_that("NB continuation skips pooled subsampling on unchanged data", {
  case <- .guard_cases()[[3]]
  data <- .guard_data("negbin")
  data$regdata <- lapply(data$regdata, function(u) {
    index <- rep(seq_along(u$y), 12L)
    list(y = u$y[index], X = u$X[index, , drop = FALSE])
  })
  expect_gt(sum(vapply(data$regdata, function(u) length(u$y), integer(1))), 1000L)
  continue_split(case, k = 4L, data = data, label = "negbin-none-subsample")
})
