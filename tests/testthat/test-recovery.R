# Recovery of Delta(Z) and Sigma(Z) on small designs with D = 2 coefficients.
# Gated: run with BAYESM_HART_MC=true; the file takes about three minutes.
# Under a constant true covariance the fitted log variances still wobble
# across units (sd 0.4-0.7 on these designs), hence the loose sd_logv bounds.

library(testthat)

skip_unless_mc <- function()
  skip_if_not(identical(Sys.getenv("BAYESM_HART_MC"), "true"),
              "set BAYESM_HART_MC=true for the recovery tests")

# v_1(Z) = s1 exp(b1 Z_1), v_2 = s2, correlation rho_max tanh(2 Z_2).
sigma_fn <- function(s1, s2, b1 = 0, rho_max = 0) function(Zi) {
  v <- c(s1 * exp(b1 * Zi[1]), s2)
  r <- rho_max * tanh(2 * Zi[2]) * sqrt(v[1] * v[2])
  matrix(c(v[1], r, r, v[2]), 2, 2)
}#SIGMA_FN

direct_prior <- function(...)
  list(ncomp = 1L, hart_model = "direct", bart = list(num_trees = 50L),
       vartree = list(num_trees = 20L), phitree = list(num_trees = 20L), ...)

mcmc <- list(R = 1500L, keep = 1L, nprint = 0L)
burn <- 500L

# Posterior means of Delta(Z) + mu and Sigma(Z) at the training Z.
posterior_summary <- function(fit, Z) {
  dz <- predict(fit, newdata = list(Z = Z), type = "DeltaZ+mu", burn = burn,
                r_verbose = FALSE)
  out <- list(delta = apply(dz, c(1, 2), mean))
  if (inherits(fit, "bayesm.HART.HeterCov")) {
    sz <- predict(fit, newdata = list(Z = Z), type = "SigmaZ", burn = burn,
                  r_verbose = FALSE)
    out$sigma <- apply(sz, c(1, 2, 3), mean)
  }
  out
}#POSTERIOR_SUMMARY

# Correlations across units with the truth (NA where the truth is constant),
# the bias of the log variances, and the spread of the fitted log variances.
metrics <- function(ps, betabar, Sigma_true) {
  cor_or_na <- function(a, b) if (sd(b) > 1e-12) cor(a, b) else NA_real_
  m <- list(cor_mean  = sapply(1:2, function(k)
                          cor_or_na(ps$delta[, k], betabar[, k])),
            rmse_mean = sqrt(colMeans((ps$delta - betabar)^2)))
  if (is.null(ps$sigma)) return(m)
  logv_hat  <- log(cbind(ps$sigma[, 1, 1], ps$sigma[, 2, 2]))
  logv_true <- log(cbind(Sigma_true[, 1, 1], Sigma_true[, 2, 2]))
  rho_hat   <- ps$sigma[, 1, 2] / sqrt(ps$sigma[, 1, 1] * ps$sigma[, 2, 2])
  rho_true  <- Sigma_true[, 1, 2] /
               sqrt(Sigma_true[, 1, 1] * Sigma_true[, 2, 2])
  c(m, list(cor_logv  = sapply(1:2, function(k)
                        cor_or_na(logv_hat[, k], logv_true[, k])),
            bias_logv = colMeans(logv_hat - logv_true),
            sd_logv   = apply(logv_hat, 2, sd),
            cor_rho   = cor_or_na(rho_hat, rho_true),
            rmse_rho  = sqrt(mean((rho_hat - rho_true)^2))))
}#METRICS

report <- function(label, m, secs) {
  f <- function(x) paste(round(x, 2), collapse = ",")
  cat(sprintf("\n[%s] %.0fs cor_mean=%s rmse_mean=%s", label, secs,
              f(m$cor_mean), f(m$rmse_mean)))
  if (!is.null(m$cor_logv))
    cat(sprintf(" cor_logv=%s bias_logv=%s sd_logv=%s cor_rho=%s rmse_rho=%s",
                f(m$cor_logv), f(m$bias_logv), f(m$sd_logv), f(m$cor_rho),
                f(m$rmse_rho)))
  cat("\n")
}#REPORT

constant_sigma <- function(n, s) {
  out <- array(0, dim = c(n, 2, 2))
  for (i in seq_len(n)) out[i, , ] <- s
  out
}#CONSTANT_SIGMA

# Linear ----------------------------------------------------------------------

linear_data <- function(seed, het, target_var_betabar, Sigma_func = NULL)
  sim_hier_linear(nreg = 200, nobs = 10, nvar = 1, nz = 2, const = TRUE,
                  het_observed = het, target_var_betabar = target_var_betabar,
                  target_var_eps = 0.25, sigma_sq = 0.25,
                  Sigma_func = Sigma_func, seed = seed)

fit_linear <- function(sim, Prior) {
  t0 <- Sys.time()
  capture.output(fit <- rhierLinearMixture(
    Data = list(regdata = sim$regdata, Z = sim$Z), Prior = Prior,
    Mcmc = mcmc, r_verbose = FALSE))
  list(fit = fit, secs = as.numeric(Sys.time() - t0, units = "secs"))
}#FIT_LINEAR

test_that("linear: mean-only heterogeneity", {
  skip_unless_mc()
  sim <- linear_data(601, "step", 1)
  r <- fit_linear(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               constant_sigma(200, sim$true_values$Sigma_eps))
  report("linear mean-only", m, r$secs)
  expect_gt(m$cor_mean[1], 0.9)
  expect_lt(m$rmse_mean[1], 0.35)
  expect_lt(max(abs(m$bias_logv)), 0.6)
  expect_lt(max(m$sd_logv), 0.8)
  expect_lt(m$rmse_rho, 0.35)
})#TEST_THAT

test_that("linear: variance-only heterogeneity", {
  skip_unless_mc()
  sim <- linear_data(602, "linear", 0, sigma_fn(0.25, 0.25, b1 = 1.2))
  r <- fit_linear(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               sim$true_values$Sigma_true)
  report("linear variance-only", m, r$secs)
  expect_gt(m$cor_logv[1], 0.7)
  expect_lt(max(m$rmse_mean), 0.3)
  expect_lt(m$rmse_rho, 0.35)
})#TEST_THAT

test_that("linear: joint heterogeneity with cross-covariance", {
  skip_unless_mc()
  sim <- linear_data(603, "step", 1,
                     sigma_fn(0.25, 0.25, b1 = 1.2, rho_max = 0.7))
  r <- fit_linear(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               sim$true_values$Sigma_true)
  report("linear joint", m, r$secs)
  expect_gt(m$cor_mean[1], 0.85)
  expect_gt(m$cor_logv[1], 0.6)
  expect_gt(m$cor_rho, 0.5)
})#TEST_THAT

# Negbin ----------------------------------------------------------------------

negbin_data <- function(seed, het, target_var_betabar, Sigma_func = NULL)
  sim_hier_negbin(nreg = 150, nobs = 20, nvar = 1, nz = 2, const = TRUE,
                  het_observed = het, target_var_betabar = target_var_betabar,
                  target_var_eps = 0.05, alpha = 5, Sigma_func = Sigma_func,
                  seed = seed)

fit_negbin <- function(sim, Prior) {
  t0 <- Sys.time()
  capture.output(fit <- rhierNegbinRw(
    Data = list(regdata = sim$regdata, Z = sim$Z), Prior = Prior,
    Mcmc = mcmc, r_verbose = FALSE))
  list(fit = fit, secs = as.numeric(Sys.time() - t0, units = "secs"))
}#FIT_NEGBIN

test_that("negbin: mean-only heterogeneity", {
  skip_unless_mc()
  sim <- negbin_data(611, "step", 0.25)
  r <- fit_negbin(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               constant_sigma(150, sim$true_values$Sigma_eps))
  report("negbin mean-only", m, r$secs)
  expect_gt(m$cor_mean[1], 0.85)
  expect_lt(max(m$sd_logv), 0.9)
})#TEST_THAT

test_that("negbin: variance-only heterogeneity", {
  skip_unless_mc()
  sim <- negbin_data(612, "linear", 0, sigma_fn(0.05, 0.05, b1 = 1.2))
  r <- fit_negbin(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               sim$true_values$Sigma_true)
  report("negbin variance-only", m, r$secs)
  expect_gt(m$cor_logv[1], 0.6)
  expect_lt(max(m$rmse_mean), 0.3)
})#TEST_THAT

test_that("negbin: joint heterogeneity with cross-covariance", {
  skip_unless_mc()
  sim <- negbin_data(613, "step", 0.25,
                     sigma_fn(0.05, 0.05, b1 = 1.2, rho_max = 0.7))
  r <- fit_negbin(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               sim$true_values$Sigma_true)
  report("negbin joint", m, r$secs)
  expect_gt(m$cor_mean[1], 0.8)
  expect_gt(m$cor_logv[1], 0.5)
  expect_gt(m$cor_rho, 0.4)
})#TEST_THAT

# MNL -------------------------------------------------------------------------

step_mean <- function(Zi) c(ifelse(Zi[1] > 0, 1, -1), 0.5)

mnl_data <- function(seed, Sigma_func = NULL, step = TRUE)
  sim_hier_mnl(nlgt = 150, nT = 10, p = 3, nz = 2, nXa = 2, nXd = 0,
               const = FALSE, beta_func_type = "custom",
               beta_func_args = list(func = if (step) step_mean
                                            else function(Zi) c(0, 0.5)),
               sigma_inv_diag = 2, Sigma_func = Sigma_func, seed = seed)

fit_mnl <- function(sim, Prior) {
  t0 <- Sys.time()
  capture.output(fit <- rhierMnlRwMixture(Data = sim, Prior = Prior,
                                          Mcmc = mcmc, r_verbose = FALSE))
  list(fit = fit, secs = as.numeric(Sys.time() - t0, units = "secs"))
}#FIT_MNL

test_that("mnl: mean-only heterogeneity, direct and normalized", {
  skip_unless_mc()
  sim <- mnl_data(621)
  truth_sigma <- constant_sigma(150, diag(0.5, 2))
  r <- fit_mnl(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               truth_sigma)
  report("mnl mean-only direct", m, r$secs)
  expect_gt(m$cor_mean[1], 0.85)
  expect_lt(max(m$sd_logv), 0.9)
  rn <- fit_mnl(sim, list(ncomp = 1L, hart_model = "normalized",
                          bart = list(num_trees = 50L)))
  mn <- metrics(posterior_summary(rn$fit, sim$Z), sim$true_values$betabar_true,
                truth_sigma)
  report("mnl mean-only normalized", mn, rn$secs)
  expect_gt(mn$cor_mean[1], 0.85)
})#TEST_THAT

test_that("mnl: variance-only heterogeneity", {
  skip_unless_mc()
  sim <- mnl_data(622, sigma_fn(0.5, 0.5, b1 = 1.2), step = FALSE)
  r <- fit_mnl(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               sim$true_values$Sigma_true)
  report("mnl variance-only", m, r$secs)
  expect_gt(m$cor_logv[1], 0.5)
  expect_lt(max(m$rmse_mean), 0.4)
})#TEST_THAT

test_that("mnl: joint heterogeneity with cross-covariance", {
  skip_unless_mc()
  sim <- mnl_data(623, sigma_fn(0.5, 0.5, b1 = 1.2, rho_max = 0.7))
  r <- fit_mnl(sim, direct_prior())
  m <- metrics(posterior_summary(r$fit, sim$Z), sim$true_values$betabar_true,
               sim$true_values$Sigma_true)
  report("mnl joint", m, r$secs)
  expect_gt(m$cor_mean[1], 0.8)
  expect_gt(m$cor_logv[1], 0.4)
  expect_gt(m$cor_rho, 0.3)
})#TEST_THAT
