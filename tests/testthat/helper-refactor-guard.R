# Shared definitions for the refactor guard (test-refactor-guard.R) and the
# fixture generator (tools/refactor-guard-fixtures.R).  Each case pins one
# likelihood and one heterogeneity model on a fixed data set and seed; the
# recorded quantities are every draw the sampler writes plus the predictions
# derived from them.

# Resolves from the test directory and from the package root (the generator).
.guard_dir <- function() {
  root <- file.path("tests", "testthat", "fixtures", "refactor-guard")
  if (dir.exists(dirname(dirname(root)))) root
  else file.path("fixtures", "refactor-guard")
}#.GUARD_DIR

# Mean trees are few and shallow: the chain is short, and birth/death still
# fires often enough to pin the tree strings.
.guard_priors <- list(
  none       = list(ncomp = 1L),
  normalized = list(ncomp = 1L, hart_model = "normalized",
                    bart = list(num_trees = 5L)),
  constant   = list(ncomp = 1L, hart_model = "direct",
                    hart_cov = "constant", bart = list(num_trees = 5L)),
  heteroskedastic = list(ncomp = 1L, hart_model = "direct",
                         hart_cov = "heteroskedastic",
                         bart = list(num_trees = 5L),
                         vartree = list(num_trees = 3L),
                         phitree = list(num_trees = 3L)))

.guard_mcmc <- list(mnl    = list(R = 30L, keep = 1L, nprint = 0L),
                    linear = list(R = 30L, keep = 1L, nprint = 0L),
                    negbin = list(R = 30L, keep = 1L, nprint = 0L))

# One data set per likelihood, so a case differs only by the model.
.guard_data <- function(likelihood) {
  switch(likelihood,
    mnl = sim_hier_mnl(nlgt = 12L, nT = 5L, p = 3L, nz = 2L, nXa = 2L,
                       nXd = 0L, const = FALSE, beta_func_type = "linear",
                       seed = 11L),
    linear = sim_hier_linear(nreg = 12L, nobs = 8L, nvar = 1L, nz = 2L,
                             const = TRUE, het_observed = "linear",
                             seed = 12L),
    negbin = sim_hier_negbin(nreg = 12L, nobs = 8L, nvar = 1L, nz = 2L,
                             const = TRUE, het_observed = "linear",
                             seed = 13L))
}#.GUARD_DATA

.guard_fit <- function(likelihood, data, Prior, seed) {
  sampler <- switch(likelihood, mnl = rhierMnlRwMixture,
                    linear = rhierLinearMixture, negbin = rhierNegbinRw)
  set.seed(seed)
  utils::capture.output(
    fit <- sampler(Data = data, Prior = Prior,
                   Mcmc = .guard_mcmc[[likelihood]], r_verbose = FALSE))
  fit
}#.GUARD_FIT

# Prior and posterior predictive draws.  The posterior runs over all
# estimation units, the prior over the first `npred`.  Both consume RNG, so
# each call gets its own seed.
.guard_predictive <- function(fit, likelihood, data, npred = 4L) {
  out  <- list()
  seed_of <- function(k) 7000L + 10L * match(likelihood, names(.guard_mcmc)) + k
  units <- if (likelihood == "mnl") data$lgtdata else data$regdata
  Xall  <- lapply(units, function(u) list(X = u$X))
  Zp    <- data$Z[seq_len(npred), , drop = FALSE]
  Xp    <- lapply(seq_len(npred), function(i) units[[i]]$X)
  if (likelihood == "mnl") {
    type <- "choice_probs"
    post <- list(nlgtdata = Xall, p = data$p)
    pri  <- list(Z = Zp, X = Xp, p = data$p)
  } else {
    type <- "response"
    post <- list(regdata = Xall)
    pri  <- list(Z = Zp, X = Xp)
  }
  set.seed(seed_of(1L))
  out$posterior <- tryCatch(
    predict(fit, newdata = post, mode = "posterior", type = type, burn = 10L,
            r_verbose = FALSE),
    error = function(e) paste("error:", conditionMessage(e)))
  set.seed(seed_of(2L))
  out$prior <- tryCatch(
    predict(fit, newdata = pri, mode = "prior", type = type, burn = 10L,
            nsim = 3L, r_verbose = FALSE),
    error = function(e) paste("error:", conditionMessage(e)))
  # Z off the training rows: the unique-Z cache misses and the ensembles are
  # evaluated instead.
  pri$Z <- Zp + 0.013
  set.seed(seed_of(3L))
  out$prior_oos <- tryCatch(
    predict(fit, newdata = pri, mode = "prior", type = type, burn = 10L,
            nsim = 3L, r_verbose = FALSE),
    error = function(e) paste("error:", conditionMessage(e)))
  out
}#.GUARD_PREDICTIVE

# Every draw the fit carries, plus predictions at the training Z.  Absent
# fields are dropped: the guard compares the recorded ones and tolerates
# fields a later schema adds.
.guard_record <- function(fit, Z) {
  out <- list(betadraw = fit$betadraw, loglike = fit$loglike,
              Deltadraw = fit$Deltadraw, mu_draw = fit$mu_draw,
              mean_scale_draw = fit$mean_scale_draw,
              varcount = fit$varcount, var_varcount = fit$var_varcount,
              taudraw = fit$taudraw, alphadraw = fit$alphadraw,
              acceptrbeta = fit$acceptrbeta)
  if (!is.null(fit$nmix$compdraw)) {
    out$nmix_mu <- t(vapply(fit$nmix$compdraw, function(x) x[[1]]$mu,
                            numeric(dim(fit$betadraw)[2])))
    out$nmix_rooti <- vapply(fit$nmix$compdraw, function(x) x[[1]]$rooti,
                             matrix(0, dim(fit$betadraw)[2],
                                    dim(fit$betadraw)[2]))
  }
  # One vector-leaf ensemble (direct) or one scalar ensemble per coordinate.
  out$trees <- if (!is.null(fit$bart_models$treedraws))
    fit$bart_models$treedraws$trees
  else if (!is.null(fit$bart_models))
    vapply(fit$bart_models, function(b) b$treedraws$trees, character(1))
  out$var_trees <- if (!is.null(fit$var_models))
    vapply(fit$var_models, function(b) b$treedraws$trees, character(1))
  nd <- list(Z = Z)
  out$DeltaZ_mu <- predict(fit, newdata = nd, type = "DeltaZ+mu", burn = 0L,
                           r_verbose = FALSE)
  out$DeltaZ_trees <- tryCatch(
    predict(fit, newdata = nd, type = "DeltaZ", burn = 0L,
            force_tree_eval = TRUE, r_verbose = FALSE),
    error = function(e) NULL)
  out$SigmaZ <- tryCatch(
    predict(fit, newdata = nd, type = "SigmaZ", burn = 0L, r_verbose = FALSE),
    error = function(e) NULL)
  Filter(Negate(is.null), out)
}#.GUARD_RECORD

# (likelihood, model) -> one case; the seed encodes the pair.
.guard_cases <- function() {
  grid <- expand.grid(likelihood = c("mnl", "linear", "negbin"),
                      model = names(.guard_priors),
                      stringsAsFactors = FALSE)
  lapply(seq_len(nrow(grid)), function(i) {
    lik <- grid$likelihood[i]; mod <- grid$model[i]
    list(name = paste(lik, mod, sep = "-"), likelihood = lik, model = mod,
         seed = 1000L * match(lik, c("mnl", "linear", "negbin")) +
                match(mod, names(.guard_priors)))
  })
}#.GUARD_CASES

.guard_run <- function(case) {
  data <- .guard_data(case$likelihood)
  fit  <- .guard_fit(case$likelihood, data, .guard_priors[[case$model]],
                     case$seed)
  c(.guard_record(fit, data$Z),
    .guard_predictive(fit, case$likelihood, data))
}#.GUARD_RUN
