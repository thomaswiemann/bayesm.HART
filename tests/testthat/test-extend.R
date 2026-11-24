# Appended observations, stable unit order and population checkpoints.

library(testthat)

extend_fit <- function(likelihood, ...) {
  sampler <- switch(likelihood, mnl = rhierMnlRwMixture,
                    linear = rhierLinearMixture, negbin = rhierNegbinRw)
  capture.output(fit <- sampler(..., r_verbose = FALSE))
  fit
}#EXTEND_FIT

extend_batch <- function(saved, ids = c("new-a", "2", "new-b", "new-c")) {
  units <- if (saved$likelihood == "mnl") "lgtdata" else "regdata"
  data <- saved$data
  batch <- list(unit_id = ids)
  batch[[units]] <- lapply(seq_along(ids), function(i) {
    u <- data[[units]][[i]]
    rows <- if (units == "lgtdata") seq_len(data$p) else 1L
    list(y = if (units == "lgtdata") 1 else u$y[1],
         X = u$X[rows, , drop = FALSE])
  })
  if (!is.null(data$Z)) {
    batch$Z <- data$Z[seq_along(ids), , drop = FALSE]
    for (i in seq_along(ids)) {
      j <- match(ids[i], data$unit_id)
      if (!is.na(j)) batch$Z[i, ] <- data$Z[j, ]
      else if (i == 3L) batch$Z[i, ] <- colMeans(data$Z[1:2, , drop = FALSE])
      else if (i > 1L) batch$Z[i, ] <- batch$Z[i, ] + 5 * i
    }
  }
  batch
}#EXTEND_BATCH

extend_cutpoints <- function(x) {
  if (!is.list(x)) return(NULL)
  if (!is.null(x$cutpoints)) return(x$cutpoints)
  lapply(x, extend_cutpoints)
}#EXTEND_CUTPOINTS

for (likelihood in names(.test_mcmc)) {
  for (model in names(.test_priors)) {
    if (likelihood != "linear" && model != "heteroskedastic") next
    local({
      family <- likelihood
      regime <- model
      test_that(paste("mixed extension preserves inputs and cutpoints", family, regime), {
        data <- .test_data(family)
        units <- if (family == "mnl") "lgtdata" else "regdata"
        prior <- .test_priors[[regime]]
        if (family == "mnl") prior$SignRes <- c(-1, 0)
        a <- extend_fit(family, Data = data, Prior = prior,
                        Mcmc = list(R = 6, keep = 1, nprint = 0))
        batch <- extend_batch(a$restart)
        before <- unserialize(serialize(a, NULL))
        merged <- .hart_merge_data(a$restart, batch)
        expect_identical(a, before)
        expect_identical(merged$unit_id,
                         c(a$restart$data$unit_id, "new-a", "new-b", "new-c"))
        expect_equal(merged[[units]][[2]]$y,
                     c(data[[units]][[2]]$y, batch[[units]][[2]]$y))
        expect_equal(merged[[units]][[2]]$X,
                     rbind(data[[units]][[2]]$X, batch[[units]][[2]]$X))
        expect_identical(merged[[units]][[1]], a$restart$data[[units]][[1]])
        expect_equal(merged$Z, rbind(data$Z, batch$Z[c(1, 3, 4), , drop = FALSE]))
        for (j in seq_len(3))
          expect_identical(merged[[units]][[12 + j]], batch[[units]][[c(1, 3, 4)[j]]])
        before <- unserialize(serialize(a, NULL))
        b <- extend_fit(family, Data = batch, previous = a, Mcmc = list(R = 4))
        expect_identical(a, before)
        expect_identical(b$restart$data, merged)
        expect_equal(dim(b$betadraw)[1], 15)
        expect_equal(nrow(b$restart$state$betas), 15)
        expect_equal(b$restart$state$completed_sweeps, 10)
        expect_identical(extend_cutpoints(b$restart$state$backend),
                         extend_cutpoints(a$restart$state$backend))
        if (family == "linear") {
          expect_equal(b$restart$config$prior$ssq[1:12], a$restart$config$prior$ssq)
          expect_equal(b$restart$config$prior$ssq[13:15], rep(1, 3))
          expect_length(b$restart$state$tau, 15)
        }
        if (regime == "none") expect_length(b$restart$state$backend$ind, 15)
        if (regime != "none") {
          for (type in c("DeltaZ", "SigmaZ")) {
            cache <- predict(b, newdata = list(Z = merged$Z), type = type,
                             burn = 1, r_verbose = FALSE)
            trees <- predict(b, newdata = list(Z = merged$Z), type = type,
                             burn = 1, force_tree_eval = TRUE, r_verbose = FALSE)
            expect_equal(cache, trees, tolerance = 1e-6, ignore_attr = TRUE)
          }
        }
        batch2 <- extend_batch(b$restart, c("new-b", "new-d"))
        cfit <- extend_fit(family, Data = batch2, previous = b, Mcmc = list(R = 2))
        expect_identical(cfit$restart$data, .hart_merge_data(b$restart, batch2))
        expect_equal(cfit$restart$state$completed_sweeps, 12)
        expect_equal(nrow(cfit$restart$state$betas), 16)
      })
    })
  }
}

test_that("old-only and new-only batches retain exact order and observations", {
  for (family in names(.test_mcmc)) {
    a <- extend_fit(family, Data = .test_data(family), Prior = list(ncomp = 1),
                    Mcmc = list(R = 4, keep = 1, nprint = 0))
    for (ids in list(c("4", "1"), c("new-z", "new-y"))) {
      batch <- extend_batch(a$restart, ids)
      merged <- .hart_merge_data(a$restart, batch)
      b <- extend_fit(family, Data = batch, previous = a, Mcmc = list(R = 2))
      expect_identical(b$restart$data, merged)
      units <- if (family == "mnl") "lgtdata" else "regdata"
      for (j in seq_along(a$restart$data$unit_id)) {
        i <- match(a$restart$data$unit_id[j], ids)
        expected <- a$restart$data[[units]][[j]]
        if (!is.na(i)) {
          expected$y <- c(expected$y, batch[[units]][[i]]$y)
          expected$X <- rbind(expected$X, batch[[units]][[i]]$X)
        }
        expect_identical(merged[[units]][[j]], expected)
      }
      expect_identical(merged$unit_id, c(a$restart$data$unit_id,
                                        ids[!ids %in% a$restart$data$unit_id]))
    }
  }
})

test_that("linear ssq retains old values and defaults for no-Z extension", {
  data <- .test_data("linear")
  data$Z <- NULL
  a <- extend_fit("linear", Data = data, Prior = list(ncomp = 1),
                  Mcmc = list(R = 4, keep = 1, nprint = 0))
  batch <- extend_batch(a$restart, c("2", "new"))
  ssq <- c(a$restart$config$prior$ssq, 1)
  b <- extend_fit("linear", Data = batch, previous = a,
                  Mcmc = list(R = 2))
  expect_equal(b$restart$config$prior$ssq, ssq)
  expect_null(b$restart$data$Z)
  expect_equal(dim(b$betadraw)[1], 13)
  expect_equal(nrow(b$restart$state$betas), 13)
  expect_length(b$restart$state$tau, 13)
  batch$Z <- matrix(0, 2, 2)
  expect_error(.hart_merge_data(a$restart, batch), "without Z")
})

test_that("NB fixed and sampled alpha extend above the subsampling threshold", {
  for (fixed in c(FALSE, TRUE)) {
    a <- extend_fit("negbin", Data = .test_data("negbin"),
                    Prior = list(ncomp = 1),
                    Mcmc = list(R = 4, keep = 1, nprint = 0,
                                fixalpha = fixed, alpha = 2))
    batch <- extend_batch(a$restart, "2")
    u <- a$restart$data$regdata[[2]]
    batch$regdata[[1]] <- list(y = rep(u$y, 130),
                               X = u$X[rep(seq_len(nrow(u$X)), 130), , drop = FALSE])
    b <- extend_fit("negbin", Data = batch, previous = a, Mcmc = list(R = 2))
    expect_gt(sum(lengths(lapply(b$restart$data$regdata, `[[`, "y"))), 1000)
    expect_identical(b$restart$config$fixalpha, fixed)
    if (fixed) expect_equal(b$restart$state$alpha, 2)
    else expect_true(is.finite(b$restart$state$alpha))
  }
})

test_that("batch validation rejects incompatible IDs, coordinates and priors", {
  for (family in names(.test_mcmc)) {
    a <- extend_fit(family, Data = .test_data(family), Prior = list(ncomp = 1),
                    Mcmc = list(R = 4, keep = 1, nprint = 0))
    units <- if (family == "mnl") "lgtdata" else "regdata"
    batch <- extend_batch(a$restart, c("2", "new"))
    invalid <- list()
    invalid[[1]] <- batch; invalid[[1]]$unit_id <- NULL
    invalid[[2]] <- batch; invalid[[2]]$unit_id <- c("2", "2")
    invalid[[3]] <- batch; invalid[[3]]$unit_id[1] <- NA
    invalid[[4]] <- batch; invalid[[4]][[units]][[1]]$X <- matrix(0, 1, 3)
    invalid[[5]] <- batch; invalid[[5]]$Z[1, 1] <- batch$Z[1, 1] + 1
    invalid[[6]] <- batch; invalid[[6]]$Z <- NULL
    invalid[[7]] <- batch; invalid[[7]]$unit_id <- "2"
    if (family == "mnl") {
      invalid[[8]] <- batch; invalid[[8]]$p <- 4
    }
    for (bad in invalid) {
      rng <- .Random.seed
      expect_error(extend_fit(family, Data = bad, previous = a,
                              Mcmc = list(R = 2)), "Data")
      expect_identical(.Random.seed, rng)
    }
    expect_error(extend_fit(family, Data = batch, previous = a,
                            Prior = list(ncomp = 1), Mcmc = list(R = 2)),
                 "Prior is inherited with previous")
    named <- a$restart
    for (i in seq_along(named$data[[units]]))
      colnames(named$data[[units]][[i]]$X) <- c("a", "b")
    colnames(batch[[units]][[1]]$X) <- c("b", "a")
    expect_error(.hart_merge_data(named, batch), "column names")
  }
})

test_that("state extension preserves old live values and initializes new rows", {
  a <- extend_fit("linear", Data = .test_data("linear"), Prior = list(ncomp = 1),
                  Mcmc = list(R = 4, keep = 1, nprint = 0))
  beta <- matrix(7, 14, 2)
  state <- .hart_extend_state(a$restart, beta, rep(3, 14))
  expect_equal(state$betas[1:12, ], a$restart$state$betas)
  expect_equal(state$betas[13:14, ], beta[13:14, ])
  expect_equal(state$tau[1:12], as.vector(a$restart$state$tau))
  expect_equal(state$tau[13:14], c(3, 3))
  expect_equal(tail(state$backend$ind, 2), c(1, 1))
  expect_identical(state$backend$delta, a$restart$state$backend$delta)
})


test_that("new linear units use their variance and matrix responses append by row", {
  a <- extend_fit("linear", Data = .test_data("linear"), Prior = list(ncomp = 1),
                  Mcmc = list(R = 4, keep = 1, nprint = 0))
  batch <- extend_batch(a$restart, "new")
  batch$regdata[[1]] <- a$restart$data$regdata[[1]]
  batch$regdata[[1]]$y <- seq_along(batch$regdata[[1]]$y)
  b <- extend_fit("linear", Data = batch, previous = a, Mcmc = list(R = 2))
  expect_equal(b$restart$config$prior$ssq[13], var(batch$regdata[[1]]$y))
  saved <- a$restart
  saved$data$regdata[[1]]$y <- matrix(saved$data$regdata[[1]]$y, ncol = 1)
  batch$unit_id <- "1"
  merged <- .hart_merge_data(saved, batch)
  expect_equal(merged$regdata[[1]]$y,
               rbind(saved$data$regdata[[1]]$y, matrix(batch$regdata[[1]]$y, ncol = 1)))
  colnames(saved$data$Z) <- c("z1", "z2")
  colnames(batch$Z) <- c("z2", "z1")
  expect_error(.hart_merge_data(saved, batch), "Z column names")
})

# Normalized data log-likelihood at each kept draw.
extend_loglike <- function(fit, family) {
  data <- fit$restart$data
  units <- if (family == "mnl") data$lgtdata else data$regdata
  vapply(seq_len(dim(fit$betadraw)[3]), function(r) {
    sum(vapply(seq_along(units), function(i) {
      u <- units[[i]]
      eta <- as.vector(u$X %*% fit$betadraw[i, , r])
      if (family == "linear")
        return(sum(dnorm(u$y, eta, sqrt(fit$taudraw[r, i]), log = TRUE)))
      if (family == "negbin")
        return(sum(dnbinom(u$y, size = fit$alphadraw[r], mu = exp(eta), log = TRUE)))
      utility <- matrix(eta, ncol = data$p, byrow = TRUE)
      maximum <- apply(utility, 1, max)
      sum(utility[cbind(seq_along(u$y), u$y)] - maximum -
            log(rowSums(exp(utility - maximum))))
    }, numeric(1)))
  }, numeric(1))
}#EXTEND_LOGLIKE

for (family in names(.test_mcmc)) {
  local({
    likelihood <- family
    test_that(paste("kept likelihood accounts for every appended observation", likelihood), {
      a <- extend_fit(likelihood, Data = .test_data(likelihood),
                      Prior = .test_priors$heteroskedastic,
                      Mcmc = list(R = 6, keep = 1, nprint = 0))
      b <- extend_fit(likelihood, Data = extend_batch(a$restart),
                      previous = a, Mcmc = list(R = 8, keep = 1))
      difference <- max(abs(extend_loglike(b, likelihood) - as.vector(b$loglike)))
      cat(sprintf("\nACCOUNTING %s max_abs=%.17g\n", likelihood, difference))
      expect_lt(difference, 1e-8)
    })
  })
}
