# Warm and cold posterior means with independent batch-means errors.

library(testthat)

target_fit <- function(likelihood, ...) {
  sampler <- switch(likelihood, mnl = rhierMnlRwMixture,
                    linear = rhierLinearMixture, negbin = rhierNegbinRw)
  capture.output(fit <- sampler(..., r_verbose = FALSE))
  fit
}#TARGET_FIT

target_cutpoints <- function(x) {
  if (!is.list(x)) return(NULL)
  if (!is.null(x$cutpoints)) return(x$cutpoints)
  lapply(x, target_cutpoints)
}#TARGET_CUTPOINTS

# Fifty batches of 100 draws after 3000 settling sweeps.
target_summary <- function(fit, Z) {
  retained <- 3001:8000
  delta <- predict(fit, newdata = list(Z = Z), type = "DeltaZ+mu",
                   burn = 3000, r_verbose = FALSE)
  delta <- t(matrix(delta, nrow = prod(dim(delta)[1:2])))
  beta <- t(apply(fit$betadraw[, , retained, drop = FALSE], c(2, 3), mean))
  values <- cbind(delta, beta)
  colnames(values) <- c(paste0("population", seq_len(ncol(delta))),
                        paste0("beta", seq_len(ncol(beta))))
  if (!is.null(fit$taudraw))
    values <- cbind(values, tau = rowMeans(fit$taudraw[retained, , drop = FALSE]))
  if (!is.null(fit$alphadraw))
    values <- cbind(values, alpha = fit$alphadraw[retained])
  batches <- rowsum(values, rep(seq_len(50), each = 100)) / 100
  list(mean = colMeans(values), se = apply(batches, 2, sd) / sqrt(50))
}#TARGET_SUMMARY

# Initial ranges contain the new rows and both designs have mean-zero Z.
target_data <- function(family) {
  data <- switch(family,
    mnl = sim_hier_mnl(nlgt = 36, nT = 12, p = 3, nz = 2, nXa = 2,
                       nXd = 0, const = FALSE, seed = 811),
    linear = sim_hier_linear(nreg = 36, nobs = 12, nvar = 1, nz = 2,
                             const = TRUE, het_observed = "linear", seed = 812),
    negbin = sim_hier_negbin(nreg = 36, nobs = 12, nvar = 1, nz = 2,
                             const = TRUE, het_observed = "linear", seed = 813))
  units <- if (family == "mnl") "lgtdata" else "regdata"
  Z <- scale(data$Z[1:32, , drop = FALSE], scale = FALSE)
  attributes(Z) <- attributes(Z)[c("dim", "dimnames")]
  initial <- list(unit_id = as.character(1:32), Z = Z)
  initial[[units]] <- lapply(data[[units]][1:32], function(u) {
    list(y = u$y[1:8], X = u$X[seq_len(8 * if (family == "mnl") 3 else 1), , drop = FALSE])
  })
  batch <- list(unit_id = as.character(c(2, 1, 33:36)),
                 Z = rbind(Z[c(2, 1), , drop = FALSE], matrix(0, 4, 2)))
  batch[[units]] <- lapply(c(2, 1, 33:36), function(i) {
    u <- data[[units]][[i]]
    keep <- if (i <= 32) 9:12 else seq_along(u$y)
    rows <- if (family == "mnl") as.vector(sapply(keep, function(j) (3*j-2):(3*j))) else keep
    list(y = u$y[keep], X = u$X[rows, , drop = FALSE])
  })
  if (family == "mnl") initial$p <- 3
  list(initial = initial, batch = batch)
}#TARGET_DATA

cases <- list(c("mnl", "none"), c("mnl", "heteroskedastic"),
              c("linear", "none"), c("linear", "heteroskedastic"),
              c("negbin", "none"), c("negbin", "heteroskedastic"),
              c("linear", "normalized"), c("linear", "constant"))
for (case in seq_along(cases)) {
  local({
    index <- case
    family <- cases[[case]][1]
    model <- cases[[case]][2]
    test_that(paste("warm and cold posterior target", family, model), {
      skip_if_not(identical(Sys.getenv("BAYESM_HART_MC"), "true"),
                  "set BAYESM_HART_MC=true for the target tests")
      data <- target_data(family)
      prior <- .test_priors[[model]]
      if (family == "linear") prior$ssq <- rep(1, 32)
      set.seed(21000 + index)
      a <- target_fit(family, Data = data$initial, Prior = prior,
                      Mcmc = list(R = 500, keep = 1, nprint = 0))
      combined <- .hart_merge_data(a$restart, data$batch)
      expect_true(all(apply(data$batch$Z, 2, min) >= apply(data$initial$Z, 2, min)))
      expect_true(all(apply(data$batch$Z, 2, max) <= apply(data$initial$Z, 2, max)))
      set.seed(22000 + index)
      args <- list(Data = data$batch, previous = a, Mcmc = list(R = 8000))
      warm <- do.call(target_fit, c(list(likelihood = family), args))
      if (family == "linear") prior$ssq <- warm$restart$config$prior$ssq
      set.seed(23000 + index)
      cold <- target_fit(family, Data = combined, Prior = prior,
                         Mcmc = list(R = 8000, keep = 1, nprint = 0))
      same <- identical(target_cutpoints(warm$restart$state$backend),
                         target_cutpoints(cold$restart$state$backend))
      expect_true(same)
      if (!same) stop("Warm and cold cutpoint grids differ")
      w <- target_summary(warm, combined$Z[c(1, 10, 20), , drop = FALSE])
      c <- target_summary(cold, combined$Z[c(1, 10, 20), , drop = FALSE])
      combined_se <- sqrt(w$se^2 + c$se^2)
      ratio <- abs(w$mean - c$mean) / combined_se
      for (j in seq_along(ratio))
        cat(sprintf("\nTARGET %s %s %s warm=%.9g cold=%.9g warm_se=%.9g cold_se=%.9g ratio=%.9g\n",
                    family, model, names(ratio)[j], w$mean[j], c$mean[j],
                    w$se[j], c$se[j], ratio[j]))
      expect_true(all(is.finite(ratio)))
      expect_lt(max(ratio), 4)
    })
  })
}
