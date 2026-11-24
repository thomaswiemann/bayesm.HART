# Small fixed cases shared by the continuation tests: one data set per
# likelihood and one prior list per heterogeneity model.

# Mean trees are few and shallow: the chain is short.
.test_priors <- list(
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

.test_mcmc <- list(mnl    = list(R = 30L, keep = 1L, nprint = 0L),
                   linear = list(R = 30L, keep = 1L, nprint = 0L),
                   negbin = list(R = 30L, keep = 1L, nprint = 0L))

# One data set per likelihood, so a case differs only by the model.
.test_data <- function(likelihood) {
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
}#.TEST_DATA

.test_fit <- function(likelihood, data, Prior, seed) {
  sampler <- switch(likelihood, mnl = rhierMnlRwMixture,
                    linear = rhierLinearMixture, negbin = rhierNegbinRw)
  set.seed(seed)
  utils::capture.output(
    fit <- sampler(Data = data, Prior = Prior,
                   Mcmc = .test_mcmc[[likelihood]], r_verbose = FALSE))
  fit
}#.TEST_FIT

# (likelihood, model) -> one case; the seed encodes the pair.
.test_cases <- function() {
  grid <- expand.grid(likelihood = c("mnl", "linear", "negbin"),
                      model = names(.test_priors),
                      stringsAsFactors = FALSE)
  lapply(seq_len(nrow(grid)), function(i) {
    lik <- grid$likelihood[i]; mod <- grid$model[i]
    list(name = paste(lik, mod, sep = "-"), likelihood = lik, model = mod,
         seed = 1000L * match(lik, c("mnl", "linear", "negbin")) +
                match(mod, names(.test_priors)))
  })
}#.TEST_CASES
