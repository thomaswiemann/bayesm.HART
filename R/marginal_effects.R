#' Calculate Marginal Effects for Hierarchical Models
#'
#' Computes the posterior distribution of average marginal effects by
#' varying a target covariate over a grid.
#'
#' @param object A fitted model object.
#' @param ... Other arguments passed to methods.
#'
#' @details This is a generic function. The main implementation for this package
#'   is `marginal_effects.rhierMnlRwMixture`.
#'
#' @return The result depends on the method implementation.
#' @export
marginal_effects <- function(object, ...) {
  UseMethod("marginal_effects")
}#MARGINAL_EFFECTS

#' @rdname marginal_effects
#' @method marginal_effects rhierMnlRwMixture
#'
#' @param object A fitted `rhierMnlRwMixture` object.
#' @param z_values A numeric matrix where each row represents a specific grid
#'   point for evaluation. Columns correspond to the covariates in the model's
#'   original `Z` matrix. Non-NA values in a row fix the corresponding
#'   covariate to that value for the grid point. NA values indicate that the
#'   covariate should take its value from the corresponding row in the base `Z`
#'   matrix provided.
#' @param Z A numeric matrix representing the base population or context over
#'   which the effects are averaged. It must have the same number of columns
#'   as required by the model (`nz`). For each grid point (row) in `z_values`,
#'   a counterfactual matrix is constructed based on `Z`, and the average
#'   effect across the rows of this counterfactual matrix is calculated.
#' @param burn Integer, the number of initial MCMC draws to discard.
#' @param verbose Logical, whether to print progress messages.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list object of class `"marginal_effects"` containing:
#'   \item{z_values}{The user-provided matrix of grid points used for evaluation.}
#'   \item{avg_betabar_draws}{A list where each element `[[i]]` is an
#'     `ncoef x ndraws_use` matrix representing the posterior draws of the
#'     average `betabar` (DeltaZ + mu), evaluated at the `i`-th grid point
#'     (row of `z_values`) and averaged over the rows of the base `Z` matrix.}
#'   \item{call}{The matched call to the function.}
#'   \item{burn}{The number of burn-in draws discarded.}
#' @export
#' @examples
#' # --- Simulate Data (requires bayesm.HART package) ---
#' if (requireNamespace("bayesm.HART", quietly = TRUE)) {
#'   sim_data <- bayesm.HART::sim_hier_mnl(nlgt = 50, nT = 10, p = 3, nz = 3, 
#'                                     nXa = 1, nXd = 0, seed = 123)
#'   Data <- list(p = sim_data$p, lgtdata = sim_data$lgtdata, Z = sim_data$Z)
#'   ncoef <- sim_data$true_values$dimensions$ncoef
#
#'   # --- Fit Model (minimal run for example) ---
#'   # Note: Use much larger R and keep for real analysis!
#'   Prior <- list(ncomp = 1)
#'   Mcmc <- list(R = 100, keep = 1)
#'   # Use try() to avoid errors stopping the example build if model fails
#'   fit <- try(bayesm.HART::rhierMnlRwMixture(Data = Data, Prior = Prior, Mcmc = Mcmc, 
#'                                         r_verbose = FALSE), silent = TRUE)
#
#'   if (!inherits(fit, "try-error")) {
#'     # --- Define Grid for Marginal Effects (vary Z1) ---
#'     target_z_index <- 1
#'     grid_z1 <- seq(min(Data$Z[, target_z_index]), 
#'                    max(Data$Z[, target_z_index]), 
#'                    length.out = 5)
#'     z_grid <- matrix(NA, nrow = length(grid_z1), ncol = ncol(Data$Z))
#'     z_grid[, target_z_index] <- grid_z1
#
#'     # --- Calculate Marginal Effects ---
#'     mfx_result <- marginal_effects(fit, 
#'                                      z_values = z_grid, 
#'                                      Z = Data$Z, 
#'                                      burn = 20, # Discard first 20 draws
#'                                      verbose = FALSE)
#
#'     print(names(mfx_result))
#'     print(dim(mfx_result$avg_betabar_draws[[1]])) # Check dimensions
#
#'     # --- Summarize Effects (see example for summary.marginal_effects) ---
#'
#'     # --- Plot Effects (see example for plot.summary.marginal_effects) ---
#'   }
#' }
marginal_effects.rhierMnlRwMixture <- function(
  object,
  z_values,
  Z,
  burn = 0,
  verbose = TRUE,
  ...
) {
  # --- 1. Input Checks & Initial Setup ---
  
  # Validate object class
  if (!inherits(object, "rhierMnlRwMixture")) {
    stop("Input 'object' must be of class 'rhierMnlRwMixture'.")
  }
  # Validate Z
  if (missing(Z) || !is.matrix(Z) || !is.numeric(Z)) {
    stop("'Z' must be a numeric matrix.")
  }
  # Validate z_values
  if (missing(z_values) || !is.matrix(z_values) || !is.numeric(z_values)) {
    stop("'z_values' must be a numeric matrix.")
  }
  # Ensure Z and z_values have compatible columns (nz)
  nz <- ncol(Z)
  if (ncol(z_values) != nz) {
      stop(paste0("Number of columns in 'z_values' (", ncol(z_values),
                  ") must match number of columns in 'Z' (", nz, ")."))
  }
  nlgt <- nrow(Z)
  n_grid_points <- nrow(z_values)
  if (n_grid_points == 0) stop("'z_values' cannot have zero rows.")

  # Validate column-wise NA constraint on z_values
  for (j in 1:nz) {
      na_count <- sum(is.na(z_values[, j]))
      if (na_count > 0 && na_count < n_grid_points) {
          stop(paste0("Column ", j, " in 'z_values' contains a mix of NA and non-NA values. ",
                      "Each column must be either entirely NA or contain no NAs."))
      }
  }#FOR j

  # Identify columns specified with non-NA grid values
  non_na_col_indices <- which(colSums(is.na(z_values)) == 0)

  # Ensure at least one column defines grid values
  if (length(non_na_col_indices) == 0) {
      stop("'z_values' must contain at least one non-NA column to define grid points.")
  }

  # Extract dimensions and handle burn-in
  if (is.null(object$betadraw)) stop("Invalid object: betadraw missing.")
  ncoef <- dim(object$betadraw)[2]
  ndraws_total <- dim(object$betadraw)[3]
  
  if (!is.numeric(burn) || length(burn) != 1 || burn < 0 || burn >= ndraws_total) {
      stop(paste0("burn must be a single non-negative integer less than the total number of draws (", ndraws_total, ")."))
  }#IF
  ndraws_use <- ndraws_total - burn

  # --- 2. Initialize Storage ---
  avg_draws_list <- vector(mode = "list", length = n_grid_points)

  # --- 3. Iterate Through Grid Points (Rows of z_values) ---
  for (i in 1:n_grid_points) {
    # Create counterfactual Z matrix
    Z_counterfactual <- Z

    # Overwrite relevant columns with grid values using vector assignment
    grid_values_vec <- z_values[i, non_na_col_indices]
    # Assign values; R recycles the single row of grid values down the columns
    Z_counterfactual[, non_na_col_indices] <- grid_values_vec

    # Prepare list for predict()
    newdata_counterfactual <- list(Z = Z_counterfactual)
    
    # Get posterior draws of betabar (DeltaZ+mu)
    # Result dim: [nlgt x ncoef x ndraws_use]
    betabar_draws <- predict(object, newdata = newdata_counterfactual, 
                           type = "DeltaZ+mu", burn = burn, 
                           r_verbose = verbose)
    
    # Calculate the mean over individuals (dim 1) for each draw
    # Result dim: [ncoef x ndraws_use]
    avg_betabar_draws_matrix <- apply(betabar_draws, c(2, 3), mean, na.rm = TRUE)
    
    # Store the result
    avg_draws_list[[i]] <- avg_betabar_draws_matrix
    
    # Optional progress indicator
    if (verbose && n_grid_points > 1) {
      cat("  Marginal effects progress:", i, "/", n_grid_points, "grid points calculated.", fill = TRUE)
      utils::flush.console()
    }#IF
  }#FOR i
  
  # --- 4. Package Results ---
  output <- list(
    z_values = z_values,
    avg_betabar_draws = avg_draws_list,
    call = match.call(),
    burn = burn
  )
  
  # --- 5. Return Final List ---
  class(output) <- "marginal_effects"
  return(output)

}#MARGINAL_EFFECTS.RHIERMNLRWMIXTURE


# ==============================================================================
# Summary Method for Marginal Effects
# ==============================================================================

#' Summarize Marginal Effects Object
#'
#' Computes posterior means and quantiles from the draws stored in a
#' `marginal_effects` object.
#'
#' @param object An object of class `"marginal_effects"` created by
#'   `marginal_effects.rhierMnlRwMixture`.
#' @param probs Numeric vector of quantile probabilities (between 0 and 1)
#'   to compute.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list object of class `"summary.marginal_effects"` containing:
#'   \item{summary_df}{A data frame with columns `coefficient_index`,
#'     columns for each non-NA column in the input `z_values` (e.g.,
#'     `z_col_1`, `z_col_2`, ... or using original variable names if available),
#'     `mean`, and columns for each requested quantile (e.g., `q2.5`, `q50`,
#'     `q97.5`). Each row represents a specific coefficient at a specific grid
#'     point defined by a row in `z_values`.
#'   }
#'   \item{z_values}{The original `z_values` matrix provided to `marginal_effects`.}
#'   \item{probs}{The numeric vector of quantile probabilities used.}
#'   \item{call}{The matched call to the `summary` function.}
#' @export
#' @examples
#' # --- Assumes `mfx_result` exists from marginal_effects.rhierMnlRwMixture example ---
#' if (exists("mfx_result")) {
#'   mfx_summary <- summary(mfx_result, probs = c(0.05, 0.5, 0.95))
#'   print(mfx_summary$summary_df)
#
#'   # --- Plotting example follows in plot.summary.marginal_effects ---
#'
#' } else {
#'   message("Run the example for 'marginal_effects.rhierMnlRwMixture' first.")
#' }
summary.marginal_effects <- function(
  object,
  probs = c(0.025, 0.05, 0.5, 0.95, 0.975),
  ...
) {

  # --- 1. Input Validation ---
  if (!inherits(object, "marginal_effects")) {
    stop("Input 'object' must be of class 'marginal_effects'.")
  }
  required_elements <- c("z_values", "avg_betabar_draws")
  if (!all(required_elements %in% names(object))) {
    stop(paste("Invalid 'marginal_effects' object. Missing elements:",
               paste(setdiff(required_elements, names(object)), collapse=", ")))
  }
  if (!is.numeric(probs) || any(probs < 0) || any(probs > 1)) {
    stop("'probs' must be a numeric vector with values between 0 and 1.")
  }
  probs <- sort(unique(probs))

  # --- 2. Extract Info ---
  z_values <- object$z_values
  avg_betabar_draws <- object$avg_betabar_draws

  if (!is.matrix(z_values)) stop("'object$z_values' must be a matrix.")
  if (!is.list(avg_betabar_draws) || length(avg_betabar_draws) == 0) {
      stop("'object$avg_betabar_draws' must be a non-empty list.")
  }

  n_grid_points <- nrow(z_values)
  if (n_grid_points == 0) stop("'object$z_values' cannot have zero rows.")
  if(length(avg_betabar_draws) != n_grid_points) {
      stop("Length mismatch between z_values and avg_betabar_draws.")
  }

  # Use first element to get ncoef
  first_draw_matrix <- avg_betabar_draws[[1]]
  if (!is.matrix(first_draw_matrix)) {
      stop("Elements of 'avg_betabar_draws' must be matrices.")
  }
  ncoef <- nrow(first_draw_matrix)
  if (ncoef == 0) stop("Cannot calculate summary with ncoef = 0.")

  # Identify non-NA columns in z_values and get their names for output df
  nz <- ncol(z_values)
  non_na_col_indices <- which(colSums(is.na(z_values)) < n_grid_points)
  z_col_names <- paste0("z_val_", non_na_col_indices)

  # Prepare quantile names
  quantile_names <- paste0("q", gsub("\\.", "", probs * 100))

  # --- 3. Iterate and Build List of Rows ---
  results_list <- list()

  for (i in 1:n_grid_points) {
    draws_matrix_i <- avg_betabar_draws[[i]]
    if (!is.matrix(draws_matrix_i) || nrow(draws_matrix_i) != ncoef) {
        warning(paste("Inconsistent matrix dimensions at grid point index", i,". Skipping."))
        next
    }

    # Get the specific grid point values (non-NA columns only)
    z_point_values_list <- list()
    z_slice <- z_values[i, non_na_col_indices, drop = FALSE]
    z_point_values_list <- as.list(z_slice)
    names(z_point_values_list) <- z_col_names

    for (k in 1:ncoef) {
      # Extract draws for this coefficient and grid point
      coef_draws <- draws_matrix_i[k, ]

      # Calculate summaries
      mean_val <- mean(coef_draws, na.rm = TRUE)
      quantile_vals <- stats::quantile(coef_draws, probs = probs, na.rm = TRUE)
      names(quantile_vals) <- quantile_names

      # Combine into a named list for this row
      row_data <- c(
          list(coefficient_index = k),
          z_point_values_list, # Grid point values
          list(mean = mean_val),
          as.list(quantile_vals)
      )
      results_list[[length(results_list) + 1]] <- row_data

    }#FOR k
  }#FOR i

  # --- 4. Combine Results ---
  if (length(results_list) == 0) {
      warning("No results were calculated. Returning empty data frame.")

      # Define column names for the empty data frame
      all_col_names <- c("coefficient_index", z_col_names, "mean", quantile_names)
      # Create a list of zero-length vectors with the correct names and types
      empty_cols <- structure(vector("list", length(all_col_names)), names = all_col_names)
      for (col in all_col_names) {
          empty_cols[[col]] <- if(col == "coefficient_index") integer(0) else double(0)
      }
      summary_df <- data.frame(empty_cols, check.names = FALSE)

  } else {
      # Convert list of lists to data frame efficiently
      summary_df <- do.call(rbind, lapply(results_list, data.frame, check.names = FALSE))
      # Ensure column order is consistent
      final_col_order <- c("coefficient_index", z_col_names, "mean", quantile_names)
      # Subset columns that actually exist
      final_col_order <- final_col_order[final_col_order %in% names(summary_df)]
      summary_df <- summary_df[, final_col_order, drop = FALSE]
  }

  # --- 5. Package Output ---
  output <- list(
    summary_df = summary_df,
    z_values = z_values,
    probs = probs,
    call = match.call()
  )

  # --- 6. Assign Class & Return ---
  class(output) <- "summary.marginal_effects"
  return(output)

}#SUMMARY.MARGINAL_EFFECTS 

# ==============================================================================
# Plot Method for Summarized Marginal Effects
# ==============================================================================

#' Plot Summarized Marginal Effects
#'
#' Creates a ggplot visualization of the summarized marginal effects from one or
#' more `summary.marginal_effects` objects.
#'
#' @param x An object of class `"summary.marginal_effects"`.
#' @param ... Additional *named* objects of class `"summary.marginal_effects"`
#'   to be plotted alongside `x`. The names will be used in the legend.
#' @param coef_index Integer, the index of the coefficient (`beta`) to plot.
#'   Defaults to 1.
#' @param ci_level Numeric (0 < ci_level < 1), the credible interval level
#'   to display as a ribbon (e.g., 0.95 for a 95% CI). Requires the
#'   corresponding quantiles (e.g., q2.5 and q97.5 or q25 and q975) to be
#'   present in the `summary_df` of the input object(s). Defaults to 0.95.
#' @param plot_axis_var Character, the name of the column in `summary_df` to use
#'   on the x-axis (must be one of the `z_val_*` columns derived from the
#'   non-NA columns of the original `z_values` input to `marginal_effects`).
#'   If `NULL` (default), the function will use the first column found in
#'   `x$summary_df` that matches the pattern `^z_val_[0-9]+$`.
#' @param xlab Character, custom label for the x-axis. If `NULL` (default),
#'   the label will be the value of `plot_axis_var`.
#' @param ylab Character, custom label for the y-axis. Defaults to describing
#'   the mean effect of the selected coefficient.
#' @param title Character, custom plot title. If `NULL` (default), the title
#'   will indicate the coefficient and the variable plotted on the x-axis.
#' @param color_values Named character vector for custom colors (e.g.,
#'   `c("Model A" = "blue", "Model B" = "red")`). If `NULL`, default ggplot
#'   colors are used.
#' @param linetype_values Named character vector for custom linetypes (e.g.,
#'   `c("Model A" = "solid", "Model B" = "dashed")`). If `NULL`, default ggplot
#'   linetypes are used.
#' @param fill_values Named character vector for custom ribbon fills (e.g.,
#'   `c("Model A" = "blue", "Model B" = "red")`). If `NULL`, default ggplot
#'   fills derived from colors are used.
#' @param show_legend Logical, should the legend be displayed? Defaults to `TRUE`.
#' @param x_name Character, optional name to assign to the primary object `x`
#'   in the legend. If `NULL` (default), the function attempts to infer the name.
#'
#' @details
#' This function plots the specified coefficient's posterior mean and credible
#' interval against the variable specified by `plot_axis_var`.
#'
#' If the original `marginal_effects` call used a `z_values` matrix that varied
#' across multiple dimensions (resulting in multiple `z_val_*` columns in the
#' `summary_df`), this plot will show the relationship against the chosen
#' `plot_axis_var`, **overlaying** the results from all combinations of the
#' *other* varying dimensions. For instance, if `z_val_1` and `z_val_2` exist
#' and `plot_axis_var = "z_val_1"`, the plot will show lines for each distinct
#' value of `z_val_2` (implicitly, as they will be overlaid).
#'
#' To visualize the effect along one dimension conditional on specific values of
#' other dimensions, you should filter the `summary_df` within the
#' `summary.marginal_effects` object *before* passing it to this function.
#'
#' @return A `ggplot` object.
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line scale_color_manual
#'   scale_linetype_manual scale_fill_manual labs theme_classic theme element_blank
#'   waiver .data
#' @importFrom dplyr bind_rows filter select rename
#' @importFrom rlang list2 `!!!`
#' @export
#' @examples
#' # --- Full Example Sequence for Plotting ---
#' # Requires ggplot2 and bayesm.HART
#' if (requireNamespace("bayesm.HART", quietly = TRUE) && 
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'
#'   # 1. Simulate Data (using a step function for beta_i)
#'   # Define simulation parameters
#'   nlgt_sim <- 200; nT_sim <- 20; p_sim <- 3; nz_sim <- 2
#'   nXa_sim <- 1; nXd_sim <- 0; const_sim <- TRUE
#'   # Calculate expected ncoef based on parameters
#'   ncoef_sim <- const_sim*(p_sim - 1) + (p_sim - 1)*nXd_sim + nXa_sim
#'   # Define arguments for the step function (using defaults from sim_hier_mnl)
#'   step_args_ex <- list(
#'     cutoff = 0,               # Default cutoff
#'     beta_1 = rep(-2, ncoef_sim), # Value above cutoff (default in sim_hier_mnl)
#'     beta_2 = rep(2, ncoef_sim),  # Value below cutoff (default in sim_hier_mnl)
#'     Z_index = 1               # Step based on Z1
#'   )
#'   sim_data <- bayesm.HART::sim_hier_mnl(nlgt = nlgt_sim, nT = nT_sim, p = p_sim, nz = nz_sim,
#'                                     nXa = nXa_sim, nXd = nXd_sim, const = const_sim,
#'                                     seed = 123,
#'                                     beta_func_type = "step",
#'                                     beta_func_args = step_args_ex # Pass the full list
#'                                     )
#'   Data <- list(p = sim_data$p, lgtdata = sim_data$lgtdata, Z = sim_data$Z)
#'   # Use actual ncoef from simulation output for consistency
#'   ncoef <- sim_data$true_values$dimensions$ncoef 
#'
#'   # 2. Fit Model (minimal run for example)
#'   Prior <- list(ncomp = 1,
#'                 bart = list(num_trees = 10,
#'                             num_cut = 10))
#'   Mcmc <- list(R = 500, keep = 1, nprint = 0)
#'   fit <- try(bayesm.HART::rhierMnlRwMixture(Data = Data, Prior = Prior, 
#'                                         Mcmc = Mcmc, 
#'                                         r_verbose = FALSE), silent = TRUE)
#'
#'   if (!inherits(fit, "try-error")) {
#'     # 3. Define Grid (Vary Z1, which drives the step function)
#'     target_z_index <- 1
#'     grid_z1 <- sort(c(seq(min(Data$Z[, target_z_index]), 
#'                            max(Data$Z[, target_z_index]), 
#'                            length.out = 6), 0)) # Use more points for step
#'     z_grid <- matrix(NA, nrow = length(grid_z1), ncol = ncol(Data$Z))
#'     z_grid[, target_z_index] <- grid_z1
#'
#'     # 4. Calculate Marginal Effects
#'     mfx_result <- marginal_effects(fit, 
#'                                      z_values = z_grid, 
#'                                      Z =Data$Z, 
#'                                      burn = 200, 
#'                                      verbose = FALSE)
#'
#'     # 5. Summarize Marginal Effects
#'     mfx_summary <- summary(mfx_result, probs = c(0.025, 0.5, 0.975))
#'
#'     # 6. Plot the Summary (showing effect of Z1 on coef 1)
#'     # Ensure the axis variable name matches the column in summary_df
#'     plot_var_name <- paste0("z_val_", target_z_index) 
#'     try(plot(mfx_summary, coef_index = 1, plot_axis_var = plot_var_name), silent = TRUE)
#'
#'   } else {
#'     message("Model fitting failed in example, skipping plotting.")
#'   }
#' } else {
#'   message("Requires bayesm.HART and ggplot2 packages for examples.")
#' }
plot.summary.marginal_effects <- function(
    x,
    ...,
    coef_index = 1,
    ci_level = 0.95,
    plot_axis_var = NULL,
    xlab = NULL,
    ylab = NULL,
    title = NULL,
    color_values = NULL,
    linetype_values = NULL,
    fill_values = NULL,
    show_legend = TRUE,
    x_name = NULL
) {

  # --- 1. Gather Inputs ---
  dot_args <- rlang::list2(...)
  
  # Determine name for the 'x' object (first argument)
  if (!is.null(x_name) && is.character(x_name) && length(x_name) == 1) {
      obj_name_x <- x_name
  } else {
      obj_name_x <- deparse(substitute(x))
      # Try to improve default name if deparse gives something like ".1"
      if (grepl("^\\\\.+[0-9]+$", obj_name_x)) { 
          obj_name_x <- "Model 1" 
      }
  }

  # Construct the list of summary objects, named by model
  if (obj_name_x %in% names(dot_args)) {
    # If x was passed via dots with the same name, use only the list from dots
    summary_list <- dot_args
  } else {
    # Add x with its determined name, then add the rest from dots
    summary_list <- c(rlang::list2(!!obj_name_x := x), dot_args)
  }

  if (length(summary_list) == 0) {
    stop("No summary.marginal_effects objects provided.")
  }
  
  # Assign default names to any unnamed objects in the final list
  list_names <- names(summary_list)
  unnamed_idx <- which(list_names == "")
  if (any(unnamed_idx)) {
      warning("Unnamed objects found in '...'. Using default names like 'Model 2', 'Model 3', etc. for them.")
      default_names <- paste("Model", seq_along(summary_list))
      list_names[unnamed_idx] <- default_names[unnamed_idx]
      names(summary_list) <- list_names
  }
  # Ensure the first object definitely has a name if unnamed
  if (is.null(names(summary_list)) || names(summary_list)[1] == "") {
      warning("Could not determine a name for the first object. Using 'Model 1'.")
      names(summary_list)[1] <- "Model 1"
  }


  # --- 2. Input Validation & X-axis determination ---
  
  # Validate ci_level range
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("'ci_level' must be a single number between 0 and 1 (exclusive).")
  }
  
  # Determine required quantile column names based on ci_level
  lower_q <- (1 - ci_level) / 2
  upper_q <- 1 - lower_q
  # Generate names like q2.5 -> q25, q5 -> q50
  lower_q_name <- paste0("q", gsub("\\.", "", lower_q * 100))
  upper_q_name <- paste0("q", gsub("\\.", "", upper_q * 100))
  required_q_cols <- c(lower_q_name, upper_q_name)

  # --- Determine plot_axis_var if NULL ---
  if (is.null(plot_axis_var)) {
    # Find columns starting with 'z_val_' in the first object's summary_df
    z_val_cols <- grep("^z_val_[0-9]+$", names(x$summary_df), value = TRUE)
    if (length(z_val_cols) == 0) {
      stop("'plot_axis_var' was not specified and no 'z_val_*' columns found in summary_df to use as default.")
    }
    plot_axis_var <- z_val_cols[1]
    message("Using default x-axis variable: ", plot_axis_var)
  } else {
    # Basic check if user provided it
    if (!is.character(plot_axis_var) || length(plot_axis_var) != 1 || plot_axis_var == "") {
      stop("'plot_axis_var' must be a single non-empty character string.")
    }
  }

  # Check structure and consistency of each object
  for (i in seq_along(summary_list)) {
    obj_name <- names(summary_list)[i]
    obj <- summary_list[[i]]

    if (!inherits(obj, "summary.marginal_effects")) {
      stop("Object '", obj_name, "' is not of class 'summary.marginal_effects'.")
    }
    if (is.null(obj$summary_df) || !is.data.frame(obj$summary_df)) {
      stop("Object '", obj_name, "' is missing the 'summary_df' data frame.")
    }
    if (is.null(obj$z_values)) {
        stop("Object '", obj_name, "' is missing 'z_values'.")
    }

    # Check for required quantile columns needed for CI ribbon
    if (!all(required_q_cols %in% names(obj$summary_df))) {
      stop("Object '", obj_name, "'$summary_df is missing required quantile columns ('",
           paste(required_q_cols, collapse = "', '"), "') for ci_level = ", ci_level, ".")
    }
    
    # Check if the determined plot_axis_var exists in this object
    if (!(plot_axis_var %in% names(obj$summary_df))) {
      stop("The specified 'plot_axis_var' ('", plot_axis_var, 
           "') was not found in the summary_df of object '", obj_name, "'.")
    }
    
    # Check if plot_axis_var is numeric
    if (!is.numeric(obj$summary_df[[plot_axis_var]])){
      stop("The column '", plot_axis_var, "' in object '", obj_name, 
           "' must be numeric for plotting.")
    }
    
  }

  # Check coef_index validity based on first object
  # (Assume consistent coef indices across objects)
  first_obj_df <- summary_list[[1]]$summary_df
  max_coef <- max(first_obj_df$coefficient_index, na.rm = TRUE)
  if (!is.numeric(coef_index) || length(coef_index) != 1 ||
      coef_index < 1 || coef_index > max_coef) {
    stop("'coef_index' must be a single integer between 1 and ", max_coef, ".")
  }

  # --- 3. Combine Data ---
  plot_data <- dplyr::bind_rows(lapply(summary_list, `[[`, "summary_df"),
                                .id = "model_name")

  # Ensure model_name is a factor with levels matching input order for legend
  plot_data$model_name <- factor(plot_data$model_name, levels = names(summary_list))


  # --- 4. Filter Data ---
  plot_data <- dplyr::filter(plot_data, .data$coefficient_index == !!coef_index)

  if (nrow(plot_data) == 0) {
      stop("No data remaining after filtering for coef_index = ", coef_index)
  }

  # --- 5. Prepare Plotting ---
  # Set default labels if not provided
  if (is.null(xlab)) {
    xlab <- plot_axis_var # Default x-label to the variable being plotted
  }
  if (is.null(ylab)) {
    ylab <- paste("Avg. Effect (Coefficient", coef_index, ")")
  }
  if (is.null(title)) {
    title <- paste("Marginal Effect on Coefficient", coef_index,
                   "vs", plot_axis_var)
  }

  # --- 6. Create Plot (ggplot2) ---
  # Use !!rlang::sym() or .data[[ ]] to handle variable column name
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[plot_axis_var]]))

  # Add CI ribbon
  p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[lower_q_name]],
                   ymax = .data[[upper_q_name]],
                   fill = .data$model_name),
      alpha = 0.2)

  # Add mean line
  p <- p + ggplot2::geom_line(
      ggplot2::aes(y = .data$mean,
                   color = .data$model_name,
                   linetype = .data$model_name))

  # Apply manual scales if provided
  if (!is.null(color_values)) {
    p <- p + ggplot2::scale_color_manual(values = color_values)
  }
  if (!is.null(linetype_values)) {
    p <- p + ggplot2::scale_linetype_manual(values = linetype_values)
  }
  if (!is.null(fill_values)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_values)
  }

  # Set labels and theme
  p <- p + ggplot2::labs(title = title, x = xlab, y = ylab,
                         color = "Model", linetype = "Model", fill = "Model") +
           ggplot2::theme_classic()

  # Control legend visibility
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  # --- 7. Return ---
  return(p)

}#PLOT.SUMMARY.MARGINAL_EFFECTS 