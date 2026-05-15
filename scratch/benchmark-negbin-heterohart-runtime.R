library(bayesm.HART)

cat("heterohart_benchmark_start\n")

set.seed(20260515)

sim <- bayesm.HART::sim_hier_negbin(
  nreg = 100,
  nobs = 10,
  nvar = 3,
  nz = 3,
  const = TRUE,
  het_observed = "linear",
  alpha = 5.0
)

Z <- sim$Z
if (is.null(Z)) {
  Z <- matrix(rnorm(length(sim$regdata) * 3), ncol = 3)
}
Z <- scale(Z, center = TRUE, scale = FALSE)

Data <- list(
  regdata = sim$regdata,
  Z = as.matrix(Z)
)

Prior <- list(
  ncomp = 1L,
  bart = list(num_trees = 15),
  vartree = list(num_trees = 10)
)

Mcmc <- list(
  R = 250L,
  keep = 1L,
  nprint = 0L,
  s_beta = 2.93 / sqrt(3),
  s_alpha = 2.93
)

t1 <- system.time(
  fit1 <- bayesm.HART::rhierNegbinRw(
    Data = Data,
    Prior = Prior,
    Mcmc = Mcmc,
    r_verbose = FALSE
  )
)

t2 <- system.time(
  fit2 <- bayesm.HART::rhierNegbinRw(
    Data = Data,
    Prior = Prior,
    Mcmc = Mcmc,
    r_verbose = FALSE
  )
)

cat(sprintf("heterohart_negbin_runtime_sec_run1=%.3f\n", unname(t1[["elapsed"]])))
cat(sprintf("heterohart_negbin_runtime_sec_run2=%.3f\n", unname(t2[["elapsed"]])))
cat(sprintf("heterohart_acceptrbeta_run1=%.3f\n", fit1$acceptrbeta))
cat(sprintf("heterohart_acceptralpha_run1=%.3f\n", fit1$acceptralpha))

cat("heterohart_benchmark_done\n")
