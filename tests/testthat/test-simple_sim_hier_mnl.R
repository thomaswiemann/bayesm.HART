set.seed(12345)
TEST_NLGT <- 5000
TEST_NT <- 5
TEST_NZ <- 10
TOLERANCE_VAR <- 0.15
TOLERANCE_MEAN <- 0.1

.extract_betas <- function(sim_data) {
  nlg <- length(sim_data$lgtdata)
  nc <- sim_data$true_values$dimensions$ncoef 
  expect_equal(nc, 5, info="ncoef calculated in simulation is not 5")
  
  beta_mat <- sim_data$true_values$beta_true
  betabar_mat <- sim_data$true_values$betabar_true
  
  expect_equal(dim(beta_mat), c(nlg, nc))
  expect_equal(dim(betabar_mat), c(nlg, nc))
  
  return(list(beta = beta_mat, betabar = betabar_mat))
}

# Variance and Mean Tests

test_that("Variance and mean targets are met for 'none'", {
  target_eps_var <- 0.25
  sim_data <- simple_sim_hier_mnl(nlgt = TEST_NLGT, nT = TEST_NT, nz = TEST_NZ,
                                  het_observed = "none", 
                                  target_var_betabar = 1.0,
                                  target_var_eps = target_eps_var, 
                                  seed = 1)
  
  betas <- .extract_betas(sim_data)
  eps_mat <- betas$beta - betas$betabar
  
  expect_equal(mean(abs(betas$betabar)), 0)
  expect_equal(var(betas$betabar[,1]), 0)
  
  vars_eps <- apply(eps_mat, 2, var)
  expect_equal(vars_eps, rep(target_eps_var, 5), tolerance = TOLERANCE_VAR * target_eps_var)
})

test_that("Variance and mean targets are met for 'linear'", {
  target_bb_var <- 2.0
  target_eps_var <- 0.5
  sim_data <- simple_sim_hier_mnl(nlgt = TEST_NLGT, nT = TEST_NT, nz = TEST_NZ,
                                  het_observed = "linear", 
                                  target_var_betabar = target_bb_var,
                                  target_var_eps = target_eps_var, 
                                  seed = 2)
  
  betas <- .extract_betas(sim_data)
  eps_mat <- betas$beta - betas$betabar
  
  var_bb1 <- var(betas$betabar[, 1])
  mean_bb1 <- mean(betas$betabar[, 1])
  expect_equal(var_bb1, target_bb_var, tolerance = TOLERANCE_VAR * target_bb_var)
  expect_equal(mean_bb1, 0, tolerance = TOLERANCE_MEAN)
  
  expect_equal(var(betas$betabar[, 2]), 0) 
  expect_equal(betas$betabar[1, 2:5], c(-1, 1, -1, 1))
  
  vars_eps <- apply(eps_mat, 2, var)
  expect_equal(vars_eps, rep(target_eps_var, 5), tolerance = 5 * TOLERANCE_VAR * target_eps_var)
})

test_that("Variance and mean targets are met for 'step'", {
  target_bb_var <- 0.8
  target_eps_var <- 0.1
  sim_data <- simple_sim_hier_mnl(nlgt = TEST_NLGT, nT = TEST_NT, nz = TEST_NZ,
                                  het_observed = "step", 
                                  target_var_betabar = target_bb_var,
                                  target_var_eps = target_eps_var, 
                                  seed = 3)
  
  betas <- .extract_betas(sim_data)
  eps_mat <- betas$beta - betas$betabar
  
  var_bb1 <- var(betas$betabar[, 1])
  mean_bb1 <- mean(betas$betabar[, 1])
  expect_equal(var_bb1, target_bb_var, tolerance = TOLERANCE_VAR * target_bb_var)
  expect_equal(mean_bb1, 0, tolerance = TOLERANCE_MEAN)
  
  expect_equal(var(betas$betabar[, 2]), 0) 
  expect_equal(betas$betabar[1, 2:5], c(-1, 1, -1, 1))

  vars_eps <- apply(eps_mat, 2, var)
  expect_equal(vars_eps, rep(target_eps_var, 5), tolerance = 5 * TOLERANCE_VAR * target_eps_var)
})

test_that("Variance and mean targets are met for 'friedman'", {
  target_bb_var <- 1.5
  target_eps_var <- 1.0
  sim_data <- simple_sim_hier_mnl(nlgt = TEST_NLGT, nT = TEST_NT, nz = TEST_NZ,
                                  het_observed = "friedman", 
                                  target_var_betabar = target_bb_var,
                                  target_var_eps = target_eps_var, 
                                  seed = 4)
  
  betas <- .extract_betas(sim_data)
  eps_mat <- betas$beta - betas$betabar
  
  var_bb1 <- var(betas$betabar[, 1])
  mean_bb1 <- mean(betas$betabar[, 1])
  expect_equal(var_bb1, target_bb_var, tolerance = TOLERANCE_VAR * target_bb_var)
  expect_equal(mean_bb1, 0, tolerance = TOLERANCE_MEAN)
  
  expect_equal(var(betas$betabar[, 2]), 0) 
  expect_equal(betas$betabar[1, 2:5], c(-1, 1, -1, 1))

  vars_eps <- apply(eps_mat, 2, var)
  expect_equal(vars_eps, rep(target_eps_var, 5), tolerance = 5 * TOLERANCE_VAR * target_eps_var)
})

# Edge Case Tests

test_that("target_var_betabar = 0 works", {
  sim_data <- simple_sim_hier_mnl(nlgt = 1000, nT = TEST_NT, nz = TEST_NZ,
                                  het_observed = "linear", 
                                  target_var_betabar = 0,
                                  target_var_eps = 0.5, 
                                  seed = 5)
  betas <- .extract_betas(sim_data)
  expect_equal(var(betas$betabar[, 1]), 0, tolerance = 1e-6)
  expect_equal(mean(betas$betabar[, 1]), 0, tolerance = 1e-6)
})

test_that("target_var_eps = 0 works", {
  sim_data <- simple_sim_hier_mnl(nlgt = 1000, nT = TEST_NT, nz = TEST_NZ,
                                  het_observed = "linear", 
                                  target_var_betabar = 1.0,
                                  target_var_eps = 0, 
                                  seed = 6)
  betas <- .extract_betas(sim_data)
  eps_mat <- betas$beta - betas$betabar
  vars_eps <- apply(eps_mat, 2, var)
  expect_true(all(vars_eps < 1e-6))
})

test_that("nz = 0 works", {
  sim_data <- simple_sim_hier_mnl(nlgt = 1000, nT = TEST_NT, nz = 0,
                                  het_observed = "linear",
                                  target_var_betabar = 1.0,
                                  target_var_eps = 0.5, 
                                  seed = 7)
  expect_null(sim_data$Z)
  betas <- .extract_betas(sim_data)
  expect_equal(mean(abs(betas$betabar)), 0)
  eps_mat <- betas$beta - betas$betabar
  vars_eps <- apply(eps_mat, 2, var)
  expect_equal(vars_eps, rep(0.5, 5), tolerance = 5 * TOLERANCE_VAR * 0.5)
}) 
