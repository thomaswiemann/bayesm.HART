#' Simulate Hierarchical MNL Data with Predefined Settings & Var Control
#'
#' A wrapper around sim_hier_mnl (R/sim_hier_mnl.R) providing simpler control
#' over common DGPs using a custom betabar function. Fixes p=3, nXa=1, nXd=1,
#' const=TRUE, ncomp=1. Allows specifying target variances for betabar_i (with
#' mean zero for the first coef) and eps_i.
#'
#' @param nlgt Number of individuals.
#' @param nT Number of observations per individual.
#' @param nz Number of Z variables.
#' @param het_observed Character. Functional form of observed heterogeneity
#'   (betabar_i = f(Z_i)). Options: "none", "linear", "step", "friedman".
#' @param target_var_betabar Numeric. Target variance for the **first coefficient**
#'   of the observed component betabar_i = f(Z_i). Must be non-negative.
#'   The first coefficient will be centered to have mean zero. The remaining
#'   coefficients (2 to ncoef) are fixed constants [-1, 1, -1, 1].
#'   Ignored if het_observed = "none".
#' @param target_var_eps Numeric. Target variance for each coefficient in the
#'   unobserved component eps_i. Must be non-negative. Controls `sigma_inv_diag`.
#' @param seed Integer. Optional random seed.
#'
#' @return Output list from sim_hier_mnl.
#' @importFrom stats runif var
#' @export
#' @examples
#' \dontrun{
#' # --- Simulate data with linear observed heterogeneity ---
#' sim_lin <- simple_sim_hier_mnl(nlgt = 100, nT = 10, nz = 5,
#'                                  het_observed = "linear",
#'                                  target_var_betabar = 2.0,
#'                                  target_var_eps = 0.5,
#'                                  seed = 123)
#'
#' # Check dimensions and structure
#' print(names(sim_lin))
#' print(sim_lin$true_values$dimensions)
#' print(str(sim_lin$lgtdata[[1]])) # Structure for one individual
#' if (!is.null(sim_lin$Z)) print(dim(sim_lin$Z))
#'
#' # Plot observed vs unobserved components for first coefficient
#' beta_1 <- sim_lin$true_values$beta_true[, 1]
#' betabar_1 <- sim_lin$true_values$betabar_true[, 1]
#' eps_1 <- beta_1 - betabar_1
#'
#' print(paste("Var(betabar_k1) =", round(var(betabar_1), 2)))
#' print(paste("Var(eps_k1) =", round(var(eps_1), 2)))
#'
#' # --- Simulate data with Friedman observed heterogeneity ---
#' sim_fried <- simple_sim_hier_mnl(nlgt = 100, nT = 10, nz = 5,
#'                                    het_observed = "friedman",
#'                                    target_var_betabar = 1.5,
#'                                    target_var_eps = 0.2,
#'                                    seed = 456)
#'
#' # Plot Z vs betabar for the first coefficient
#' if (!is.null(sim_fried$Z) && ncol(sim_fried$Z) > 0) {
#'   plot(sim_fried$Z[, 1], sim_fried$true_values$betabar_true[, 1],
#'        xlab = "Z1", ylab = "betabar_k1 (Friedman)")
#' }
#'
#' # --- Simulate data with no observed heterogeneity ---
#' sim_none <- simple_sim_hier_mnl(nlgt = 100, nT = 10, nz = 5,
#'                                   het_observed = "none",
#'                                   target_var_betabar = 0, # Should be ignored
#'                                   target_var_eps = 1.0,
#'                                   seed = 789)
#'
#' # Check if betabar is zero
#' print(all(sim_none$true_values$betabar_true == 0))
#' print(paste("Var(eps_k1) =", round(var(sim_none$true_values$beta_true[, 1]), 2)))
#' }
simple_sim_hier_mnl <- function(nlgt, nT, nz, 
                                het_observed = c("none", "linear", "step", "friedman"),
                                target_var_betabar = 1.0, 
                                target_var_eps = 0.5, 
                                seed = NULL) {

  # --- Input Validation ---
  het_observed <- match.arg(het_observed)
  if (!is.numeric(target_var_betabar) || target_var_betabar < 0 || length(target_var_betabar) != 1) {
    stop("target_var_betabar must be a single non-negative number.")
  }
  if (!is.numeric(target_var_eps) || target_var_eps < 0 || length(target_var_eps) != 1) {
    stop("target_var_eps must be a single non-negative number.")
  }
  if (het_observed == "friedman" && nz < 5) {
      stop("Friedman function requires nz >= 5.")
  }

  # --- Fixed Parameters ---
  p <- 3
  nXa <- 1
  nXd <- 1
  const <- TRUE
  ncomp <- 1
  ncoef <- const * (p - 1) + (p - 1) * nXd + nXa # Should be 5
  if(ncoef != 5) warning("Calculated ncoef is not 5, fixed values for coefs 2-5 might be wrong.")

  # --- Pre-calculated Baseline Values (from data-raw script) ---
  baseline_means <- list(
    linear = 0,
    step = 0.0006,
    friedman = 3.127588
  )
  baseline_vars <- list(
    linear = 27.109929,
    step = 1.00001,
    friedman = 13.942078
  )
  # -----------------------------------------------------------

  # --- Determine sigma_inv_diag for eps_i ---
  if (target_var_eps <= 0) {
    # Effectively no error variance, use very high precision
    final_sigma_inv_diag <- 1e10 
  } else {
    final_sigma_inv_diag <- 1 / target_var_eps
  }

  # --- Handle betabar generation ---
  if (het_observed == "none" || nz == 0) {
    # No observed heterogeneity or no Z vars to drive it.
    # Use a custom function that always returns zero.
    final_beta_func_type <- "custom"
    zero_generator <- function(Zi) { rep(0, ncoef) }
    final_beta_func_args <- list(func = zero_generator) 

  } else {
    # Use custom function for scaling and centering
    final_beta_func_type <- "custom"

    # Retrieve baseline values
    bl_mean <- baseline_means[[het_observed]]
    bl_var <- baseline_vars[[het_observed]]

    # Calculate scaling and centering factors for the first coefficient
    if (bl_var <= 1e-10) { # Handle cases where baseline variance is essentially zero
      scale_factor <- 0
    } else {
      scale_factor <- sqrt(max(0, target_var_betabar) / bl_var) # Ensure sqrt >= 0
    }
    mean_shift <- bl_mean * scale_factor
    
    # Define the custom generator function (needs access to parent env)
    # Includes duplicated logic/helpers for initial value generation
    
    # Duplicate Friedman helper locally to avoid exporting sim_hier_mnl internal
    ._simple_beta_Z_friedman <- function(Zi, 
                                         coefs = c(10, 20, 10, 5), 
                                         shift = 0, scale = 8) {
      term1 <- coefs[1] * sin(pi * Zi[1] * Zi[2])
      term2 <- coefs[2] * (Zi[3] - 0.5)^2 
      term3 <- coefs[3] * Zi[4]
      term4 <- coefs[4] * Zi[5]
      return((term1 + term2 + term3 + term4 + shift) / scale) 
    }#._SIMPLE_BETA_Z_FRIEDMAN
    
    custom_betabar_generator <- function(Zi) {
      betabar_i_final <- numeric(ncoef) # Initialize vector (ncoef=5)
      
      # Generate initial scalar value for the FIRST coefficient
      betabar_i1_init <- 0
      if (het_observed == "linear") {
        Delta_vec <- numeric(nz)
        nz_used <- min(nz, 5)
        if (nz_used > 0) {
          fixed_coefs <- c(3, -3, 2, -2, 1)
          Delta_vec[1:nz_used] <- fixed_coefs[1:nz_used]
        }
        betabar_i1_init <- sum(Zi * Delta_vec) # Vector product
      } else if (het_observed == "step") {
        beta1_val <- -1; beta2_val <- 1; Z_index <- 1; cutoff <- 0
        betabar_i1_init <- ifelse(Zi[Z_index] > cutoff, beta1_val, beta2_val)
      } else if (het_observed == "friedman") {
        # Use locally defined helper
        betabar_i1_init <- ._simple_beta_Z_friedman(Zi)
      } 
      
      # Scale and de-mean the FIRST coefficient
      betabar_i_final[1] <- (betabar_i1_init * scale_factor) - mean_shift
      
      # Set fixed values for remaining coefficients
      if (ncoef > 1) {
          betabar_i_final[2:ncoef] <- c(-1, 1, -1, 1)
      }
      
      return(betabar_i_final)
    }#CUSTOM_BETABAR_GENERATOR

    final_beta_func_args <- list(func = custom_betabar_generator)

  }#END if het_observed != "none"

  # --- Calculate and report implied R-squared for first coefficient ---
  total_var_k1 <- target_var_betabar + target_var_eps
  implied_R2_k1 <- if (total_var_k1 > 1e-10) {
      target_var_betabar / total_var_k1
  } else {
      NA # R-squared is undefined if total variance is zero
  }
  message(paste(" -> Implied R2 for coef 1 (Var(betabar)/Var(beta)):", 
                round(implied_R2_k1, 3)))

  # --- Call sim_hier_mnl --- 
  # Assumes sim_hier_mnl is available in the environment
  sim_data <- sim_hier_mnl(
      nlgt = nlgt, nT = nT, p = p, nz = nz, nXa = nXa, nXd = nXd,
      const = const,
      # Z generation defaults are used from sim_hier_mnl
      beta_func_type = final_beta_func_type, 
      beta_func_args = final_beta_func_args,
      ncomp = ncomp,
      # mixture_comps = NULL (default used)
      sigma_inv_diag = final_sigma_inv_diag,
      # X generation defaults are used from sim_hier_mnl
      seed = seed
  )

  return(sim_data)

}#SIMPLE_SIM_HIER_MNL 