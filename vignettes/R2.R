# Simple R script demonstrating linear hierarchical MNL simulation
# Using simple_sim_hier_mnl function

# Load required library
library(bayesm.HART)

# Set parameters for simulation
set.seed(123)
nlgt <- 200      # Number of individuals
nT <- 10         # Number of observations per individual  
nz <- 5          # Number of Z variables
target_var_betabar <- 2.0  # Target variance for betabar (observed heterogeneity)
target_var_eps <- 0.5      # Target variance for eps (unobserved heterogeneity)

# Generate simulation data with linear hierarchical prior
cat("Generating simulation data with linear hierarchical prior...\n")
sim_data <- simple_sim_hier_mnl(
  nlgt = nlgt,
  nT = nT, 
  nz = nz,
  het_observed = "linear",
  target_var_betabar = target_var_betabar,
  target_var_eps = target_var_eps,
  seed = 123
)

# Display basic information about the simulation
cat("\n=== Simulation Results ===\n")
cat("Number of individuals:", length(sim_data$lgtdata), "\n")
cat("Number of Z variables:", ifelse(is.null(sim_data$Z), 0, ncol(sim_data$Z)), "\n")
cat("Number of coefficients:", sim_data$true_values$dimensions$ncoef, "\n")

# Extract beta components for analysis
beta_true <- sim_data$true_values$beta_true
betabar_true <- sim_data$true_values$betabar_true
eps_true <- beta_true - betabar_true

# Check variance targets were met
cat("\n=== Variance Check ===\n")
cat("Target variance for betabar (coef 1):", target_var_betabar, "\n")
cat("Actual variance for betabar (coef 1):", round(var(betabar_true[, 1]), 3), "\n")
cat("Target variance for eps (coef 1):", target_var_eps, "\n")
cat("Actual variance for eps (coef 1):", round(var(eps_true[, 1]), 3), "\n")

# Calculate R-squared for first coefficient
total_var_k1 <- var(beta_true[, 1])
observed_var_k1 <- var(betabar_true[, 1])
R2_k1 <- observed_var_k1 / total_var_k1
cat("R-squared for coefficient 1:", round(R2_k1, 3), "\n")

# Display structure of individual data
cat("\n=== Individual Data Structure ===\n")
cat("First individual data structure:\n")
str(sim_data$lgtdata[[1]])

# Simple plots if running interactively
if (interactive()) {
  # Plot observed vs unobserved components for first coefficient
  par(mfrow = c(2, 2))
  
  # Histogram of betabar (observed heterogeneity)
  hist(betabar_true[, 1], main = "Observed Heterogeneity (betabar)", 
       xlab = "betabar coefficient 1", col = "lightblue")
  
  # Histogram of eps (unobserved heterogeneity)  
  hist(eps_true[, 1], main = "Unobserved Heterogeneity (eps)",
       xlab = "eps coefficient 1", col = "lightcoral")
  
  # Scatter plot of Z vs betabar for first coefficient
  if (!is.null(sim_data$Z)) {
    plot(sim_data$Z[, 1], betabar_true[, 1], 
         main = "Z vs betabar (Linear Relationship)",
         xlab = "Z variable 1", ylab = "betabar coefficient 1",
         pch = 16, col = "darkblue")
  }
  
  # Scatter plot of betabar vs total beta
  plot(betabar_true[, 1], beta_true[, 1],
       main = "Observed vs Total Heterogeneity", 
       xlab = "betabar coefficient 1", ylab = "beta coefficient 1",
       pch = 16, col = "darkgreen")
  abline(0, 1, col = "red", lty = 2)
  
  par(mfrow = c(1, 1))
}

cat("\n=== Simulation Complete ===\n")

# === COMPARISON WITH DIFFERENT NUMBERS OF Z VARIABLES ===

# Test different numbers of Z variables
nz_values <- c(2, 5, 10, 15)
comparison_results <- list()

for (nz_test in nz_values) {
  cat(paste("\n--- Testing with nz =", nz_test, "---\n"))
  
  # Generate data with current nz value
  sim_data_nz <- simple_sim_hier_mnl(
    nlgt = 100,  # Smaller sample for faster testing
    nT = 8,
    nz = nz_test,
    het_observed = "linear",
    target_var_betabar = 1.5,
    target_var_eps = 0.3,
    seed = 123 + nz_test  # Different seed for each nz
  )
  
  # MCMC parameters (shorter chains for comparison)
  R <- 500
  burn <- 50
  keep <- 1
  Mcmc <- list(R = R, keep = keep, nprint = 0)
  
  # Fit HART model
  cat("Fitting HART model...\n")
  out_hart <- bayesm.HART::rhierMnlRwMixture(
    Data = sim_data_nz, 
    Mcmc = Mcmc, 
    Prior = list(
      ncomp = 1, 
      bart = list(num_trees = 20) # HART prior parameters
    ),
    r_verbose = FALSE
  )
  
  # Fit linear hierarchical model
  cat("Fitting linear hierarchical model...\n") 
  out_lin <- bayesm.HART::rhierMnlRwMixture(
    Data = sim_data_nz, 
    Mcmc = Mcmc,
    Prior = list(ncomp = 1),  # No bart = linear hierarchical
    r_verbose = FALSE
  )
  
  # Calculate effective sample sizes and other diagnostics
  # Remove burn-in period
  beta_draws_hart <- out_hart$betadraw[, , (burn+1):dim(out_hart$betadraw)[3]]
  beta_draws_lin <- out_lin$betadraw[, , (burn+1):dim(out_lin$betadraw)[3]]
  
  # Calculate mean squared error vs true values
  true_beta <- sim_data_nz$true_values$beta_true
  
  # MSE for HART
  beta_mean_hart <- apply(beta_draws_hart, c(1,2), mean)
  mse_hart <- mean((beta_mean_hart - true_beta)^2)
  
  # MSE for Linear
  beta_mean_lin <- apply(beta_draws_lin, c(1,2), mean)
  mse_lin <- mean((beta_mean_lin - true_beta)^2)
  
  # Store results
  comparison_results[[paste0("nz_", nz_test)]] <- list(
    nz = nz_test,
    mse_hart = mse_hart,
    mse_lin = mse_lin,
    loglike_hart = mean(out_hart$loglike[(burn+1):length(out_hart$loglike)]),
    loglike_lin = mean(out_lin$loglike[(burn+1):length(out_lin$loglike)]),
    true_R2 = var(sim_data_nz$true_values$betabar_true[,1]) / var(sim_data_nz$true_values$beta_true[,1])
  )
  
  cat("HART MSE:", round(mse_hart, 4), "\n")
  cat("Linear MSE:", round(mse_lin, 4), "\n")
  cat("HART Log-likelihood:", round(comparison_results[[paste0("nz_", nz_test)]]$loglike_hart, 2), "\n")
  cat("Linear Log-likelihood:", round(comparison_results[[paste0("nz_", nz_test)]]$loglike_lin, 2), "\n")
  cat("True R-squared:", round(comparison_results[[paste0("nz_", nz_test)]]$true_R2, 3), "\n")
}

# === SUMMARY OF RESULTS ===
cat("\n=== SUMMARY OF COMPARISON RESULTS ===\n")
cat("nz\tHART_MSE\tLin_MSE\tHART_LL\tLin_LL\tTrue_R2\n")
cat("--\t--------\t-------\t-------\t------\t-------\n")

for (result in comparison_results) {
  cat(sprintf("%d\t%.4f\t\t%.4f\t\t%.1f\t%.1f\t%.3f\n", 
              result$nz, result$mse_hart, result$mse_lin, 
              result$loglike_hart, result$loglike_lin, result$true_R2))
}

# Create a simple visualization if running interactively
if (interactive()) {
  nz_vals <- sapply(comparison_results, function(x) x$nz)
  mse_hart_vals <- sapply(comparison_results, function(x) x$mse_hart)
  mse_lin_vals <- sapply(comparison_results, function(x) x$mse_lin)
  
  # Plot MSE comparison
  plot(nz_vals, mse_hart_vals, type = "b", col = "blue", pch = 16,
       xlab = "Number of Z Variables", ylab = "Mean Squared Error",
       main = "HART vs Linear Hierarchical: MSE Comparison",
       ylim = range(c(mse_hart_vals, mse_lin_vals)))
  lines(nz_vals, mse_lin_vals, type = "b", col = "red", pch = 17)
  legend("topright", legend = c("HART", "Linear"), 
         col = c("blue", "red"), pch = c(16, 17), lty = 1)
}

cat("\n=== Analysis Complete ===\n")
