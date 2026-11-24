#' Calculate Marginal Effects for Hierarchical Models
#'
#' Computes the posterior distribution of average marginal effects by
#' varying a target covariate over a grid.
#'
#' @param object A fitted hierarchical model object (e.g.,
#'   `rhierMnlRwMixture`, `rhierLinearMixture`, `rhierNegbinRw`).
#' @param z_values A numeric matrix of grid values for the unit-level
#'   covariates `Z`. Each row defines one counterfactual `Z*`. Columns that
#'   are entirely `NA` are held at their training values; columns with no
#'   `NA`s are swept to the supplied grid value.
#' @param Z A numeric matrix of unit-level covariates from the training
#'   sample (typically `Data$Z` from the original fit). Must have the same
#'   number of columns as `z_values`.
#' @param burn Non-negative integer. Number of initial MCMC draws to drop
#'   before averaging (default `0`).
#' @param r_verbose Logical. Print progress per grid point (default `TRUE`).
#' @param ... Other arguments passed to methods.
#'
#' @details This is a generic function. Method implementations are provided
#'   for `rhierMnlRwMixture`, `rhierLinearMixture`, and `rhierNegbinRw` (the
#'   direct-model subclasses inherit dispatch through their base classes).
#'
#' @return An object of class `"marginal_effects"`; see
#'   `summary.marginal_effects()` for downstream summarization.
#' @export
marginal_effects <- function(object, z_values, Z, burn = 0, r_verbose = TRUE,
                             ...) {
  UseMethod("marginal_effects")
}#MARGINAL_EFFECTS

# Grid loop shared by the marginal_effects methods.
.marginal_effects_impl <- function(object, z_values, Z, burn, r_verbose, ...) {
  if (missing(Z) || !is.matrix(Z) || !is.numeric(Z)) {
    stop("'Z' must be a numeric matrix.")
  }
  if (missing(z_values) || !is.matrix(z_values) || !is.numeric(z_values)) {
    stop("'z_values' must be a numeric matrix.")
  }
  nz <- ncol(Z)
  if (ncol(z_values) != nz) {
      stop(paste0("Number of columns in 'z_values' (", ncol(z_values),
                  ") must match number of columns in 'Z' (", nz, ")."))
  }
  n_grid_points <- nrow(z_values)
  if (n_grid_points == 0) stop("'z_values' cannot have zero rows.")

  # A column of z_values is all NA (Z left as observed) or has no NA.
  for (j in 1:nz) {
      na_count <- sum(is.na(z_values[, j]))
      if (na_count > 0 && na_count < n_grid_points) {
          stop(paste0("Column ", j, " in 'z_values' contains a mix of NA and non-NA values. ",
                      "Each column must be either entirely NA or contain no NAs."))
      }
  }#FOR j

  non_na_col_indices <- which(colSums(is.na(z_values)) == 0)
  if (length(non_na_col_indices) == 0) {
      stop("'z_values' must contain at least one non-NA column to define grid points.")
  }

  if (is.null(object$betadraw)) stop("Invalid object: betadraw missing.")
  ndraws_total <- dim(object$betadraw)[3]

  if (!is.numeric(burn) || length(burn) != 1 || burn < 0 || burn >= ndraws_total) {
      stop(paste0("burn must be a single non-negative integer less than the total number of draws (", ndraws_total, ")."))
  }#IF

  # mu is added when the fit carries it (nmix or mu_draw).
  has_mu <- !is.null(object$nmix) || !is.null(object$mu_draw)
  predict_type <- if (has_mu) "DeltaZ+mu" else "DeltaZ"

  avg_draws_list <- vector(mode = "list", length = n_grid_points)

  for (i in 1:n_grid_points) {
    # Z with the grid columns set to row i of z_values.
    Z_counterfactual <- Z
    grid_values_vec <- z_values[i, non_na_col_indices]
    Z_counterfactual[, non_na_col_indices] <- grid_values_vec

    # npred x ncoef x ndraws
    betabar_draws <- predict(object, newdata = list(Z = Z_counterfactual),
                           mode = "coefficients",
                           type = predict_type, burn = burn,
                           r_verbose = r_verbose)

    # Average over units: ncoef x ndraws.
    avg_betabar_draws_matrix <- apply(betabar_draws, c(2, 3), mean, na.rm = TRUE)

    avg_draws_list[[i]] <- avg_betabar_draws_matrix

    if (r_verbose && n_grid_points > 1) {
      cat("  Marginal effects progress:", i, "/", n_grid_points, "grid points calculated.", fill = TRUE)
      utils::flush.console()
    }#IF
  }#FOR i

  output <- list(
    z_values = z_values,
    avg_betabar_draws = avg_draws_list,
    call = match.call(call = sys.call(-1)),
    burn = burn
  )
  class(output) <- "marginal_effects"
  return(output)
}#.MARGINAL_EFFECTS_IMPL

# Methods

#' @rdname marginal_effects
#' @method marginal_effects rhierMnlRwMixture
#' @export
marginal_effects.rhierMnlRwMixture <- function(
  object,
  z_values,
  Z,
  burn = 0,
  r_verbose = TRUE,
  ...
) {
  .marginal_effects_impl(object, z_values, Z, burn, r_verbose, ...)
}#MARGINAL_EFFECTS.RHIERMNLRWMIXTURE

#' @rdname marginal_effects
#' @method marginal_effects rhierLinearMixture
#' @export
marginal_effects.rhierLinearMixture <- function(
  object,
  z_values,
  Z,
  burn = 0,
  r_verbose = TRUE,
  ...
) {
  .marginal_effects_impl(object, z_values, Z, burn, r_verbose, ...)
}#MARGINAL_EFFECTS.RHIERLINEARMIXTURE

#' @rdname marginal_effects
#' @method marginal_effects rhierNegbinRw
#' @export
marginal_effects.rhierNegbinRw <- function(
  object,
  z_values,
  Z,
  burn = 0,
  r_verbose = TRUE,
  ...
) {
  .marginal_effects_impl(object, z_values, Z, burn, r_verbose, ...)
}#MARGINAL_EFFECTS.RHIERNEGBINRW

# ==============================================================================
# Summary Method for Marginal Effects
# ==============================================================================

#' Summarize Marginal Effects Object
#'
#' Computes posterior means and quantiles from the draws stored in a
#' `marginal_effects` object.
#'
#' @param object An object of class `"marginal_effects"` created by
#'   `marginal_effects()`.
#' @param probs Numeric vector of quantile probabilities (between 0 and 1)
#'   to compute.
#' @param ... Ignored.
#'
#' @return A list object of class `"summary.marginal_effects"` containing:
#'   \item{summary_df}{A data frame with one row per coefficient and grid
#'     point: `coefficient_index`, one column `z_val_<j>` for each grid
#'     column `j` of `z_values`, `mean`, and one column per quantile named
#'     `q` followed by the percentage without the decimal point (`q25` for
#'     0.025, `q50`, `q975`).}
#'   \item{z_values}{The original `z_values` matrix provided to `marginal_effects`.}
#'   \item{probs}{The numeric vector of quantile probabilities used.}
#'   \item{call}{The matched call to the `summary` function.}
#' @seealso [plot.summary.marginal_effects()] for a worked example.
#' @export
summary.marginal_effects <- function(
  object,
  probs = c(0.025, 0.05, 0.5, 0.95, 0.975),
  ...
) {

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

  first_draw_matrix <- avg_betabar_draws[[1]]
  if (!is.matrix(first_draw_matrix)) {
      stop("Elements of 'avg_betabar_draws' must be matrices.")
  }
  ncoef <- nrow(first_draw_matrix)
  if (ncoef == 0) stop("Cannot calculate summary with ncoef = 0.")

  # Columns of z_values that define the grid.
  nz <- ncol(z_values)
  non_na_col_indices <- which(colSums(is.na(z_values)) < n_grid_points)
  z_col_names <- paste0("z_val_", non_na_col_indices)

  quantile_names <- paste0("q", gsub("\\.", "", probs * 100))

  results_list <- list()

  for (i in 1:n_grid_points) {
    draws_matrix_i <- avg_betabar_draws[[i]]
    if (!is.matrix(draws_matrix_i) || nrow(draws_matrix_i) != ncoef) {
        warning(paste("Inconsistent matrix dimensions at grid point index", i,". Skipping."))
        next
    }

    z_point_values_list <- list()
    z_slice <- z_values[i, non_na_col_indices, drop = FALSE]
    z_point_values_list <- as.list(z_slice)
    names(z_point_values_list) <- z_col_names

    for (k in 1:ncoef) {
      coef_draws <- draws_matrix_i[k, ]

      mean_val <- mean(coef_draws, na.rm = TRUE)
      quantile_vals <- stats::quantile(coef_draws, probs = probs, na.rm = TRUE)
      names(quantile_vals) <- quantile_names

      row_data <- c(
          list(coefficient_index = k),
          z_point_values_list,
          list(mean = mean_val),
          as.list(quantile_vals)
      )
      results_list[[length(results_list) + 1]] <- row_data

    }#FOR k
  }#FOR i

  if (length(results_list) == 0) {
      warning("No results were calculated. Returning empty data frame.")

      # Empty frame with the output columns.
      all_col_names <- c("coefficient_index", z_col_names, "mean", quantile_names)
      empty_cols <- structure(vector("list", length(all_col_names)), names = all_col_names)
      for (col in all_col_names) {
          empty_cols[[col]] <- if(col == "coefficient_index") integer(0) else double(0)
      }
      summary_df <- data.frame(empty_cols, check.names = FALSE)

  } else {
      summary_df <- do.call(rbind, lapply(results_list, data.frame, check.names = FALSE))
      final_col_order <- c("coefficient_index", z_col_names, "mean", quantile_names)
      final_col_order <- final_col_order[final_col_order %in% names(summary_df)]
      summary_df <- summary_df[, final_col_order, drop = FALSE]
  }

  output <- list(
    summary_df = summary_df,
    z_values = z_values,
    probs = probs,
    call = match.call()
  )

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
#' @param ... Additional named objects of class `"summary.marginal_effects"`
#'   to be plotted alongside `x`. The names will be used in the legend.
#' @param coef_index Integer, the index of the coefficient (`beta`) to plot.
#'   Defaults to 1.
#' @param ci_level Numeric (0 < ci_level < 1), the credible interval level
#'   to display as a ribbon. The matching quantile columns (`q25` and
#'   `q975` for 0.95) must be in `summary_df`. Defaults to 0.95.
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
#'   in the legend. Default: the deparsed argument.
#'
#' @details
#' This function plots the specified coefficient's posterior mean and credible
#' interval against the variable specified by `plot_axis_var`.
#'
#' If `z_values` varied more than one column, the other grid dimensions are
#' overlaid. To condition on them, filter `summary_df` before plotting.
#'
#' @return A `ggplot` object.
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line scale_color_manual
#'   scale_linetype_manual scale_fill_manual labs theme_classic theme element_blank
#'   waiver .data
#' @importFrom dplyr bind_rows filter select rename
#' @importFrom rlang list2 `!!!` `:=`
#' @importFrom stats predict
#' @export
#' @examples
#' \donttest{
#' # Simulate
#' sim <- sim_hier_mnl(nlgt = 200, nT = 20, p = 3, nz = 2, nXa = 1, nXd = 0,
#'                     beta_func_type = "step", seed = 123)
#' Data <- list(p = sim$p, lgtdata = sim$lgtdata, Z = sim$Z)
#'
#' # Fit
#' Prior <- list(ncomp = 1, hart_model = "normalized",
#'               bart = list(num_trees = 10, numcut = 10))
#' fit <- rhierMnlRwMixture(Data = Data, Prior = Prior,
#'                          Mcmc = list(R = 500, keep = 1, nprint = 0),
#'                          r_verbose = FALSE)
#'
#' # Grid over Z1
#' z_grid <- matrix(NA, nrow = 7, ncol = ncol(Data$Z))
#' z_grid[, 1] <- sort(c(seq(min(Data$Z[, 1]), max(Data$Z[, 1]),
#'                           length.out = 6), 0))
#'
#' # Marginal effects, summary, plot
#' mfx <- marginal_effects(fit, z_values = z_grid, Z = Data$Z, burn = 200,
#'                         r_verbose = FALSE)
#' mfx_summary <- summary(mfx, probs = c(0.025, 0.5, 0.975))
#' plot(mfx_summary, coef_index = 1, plot_axis_var = "z_val_1")
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

  dot_args <- rlang::list2(...)
  
  # Label for x: x_name, else the deparsed argument.
  if (!is.null(x_name) && is.character(x_name) && length(x_name) == 1) {
      obj_name_x <- x_name
  } else {
      obj_name_x <- deparse(substitute(x))
      if (grepl("^\\\\.+[0-9]+$", obj_name_x)) { 
          obj_name_x <- "Model 1" 
      }
  }

  if (obj_name_x %in% names(dot_args)) {
    summary_list <- dot_args
  } else {
    summary_list <- c(rlang::list2(!!obj_name_x := x), dot_args)
  }

  if (length(summary_list) == 0) {
    stop("No summary.marginal_effects objects provided.")
  }
  
  list_names <- names(summary_list)
  unnamed_idx <- which(list_names == "")
  if (any(unnamed_idx)) {
      warning("Unnamed objects found in '...'. Using default names like 'Model 2', 'Model 3', etc. for them.")
      default_names <- paste("Model", seq_along(summary_list))
      list_names[unnamed_idx] <- default_names[unnamed_idx]
      names(summary_list) <- list_names
  }
  if (is.null(names(summary_list)) || names(summary_list)[1] == "") {
      warning("Could not determine a name for the first object. Using 'Model 1'.")
      names(summary_list)[1] <- "Model 1"
  }

  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("'ci_level' must be a single number between 0 and 1 (exclusive).")
  }
  
  # Interval quantile columns, named as in summary().
  lower_q <- (1 - ci_level) / 2
  upper_q <- 1 - lower_q
  lower_q_name <- paste0("q", gsub("\\.", "", lower_q * 100))
  upper_q_name <- paste0("q", gsub("\\.", "", upper_q * 100))
  required_q_cols <- c(lower_q_name, upper_q_name)

  # Default axis: the first z_val_* column.
  if (is.null(plot_axis_var)) {
    z_val_cols <- grep("^z_val_[0-9]+$", names(x$summary_df), value = TRUE)
    if (length(z_val_cols) == 0) {
      stop("'plot_axis_var' was not specified and no 'z_val_*' columns found in summary_df to use as default.")
    }
    plot_axis_var <- z_val_cols[1]
    message("Using default x-axis variable: ", plot_axis_var)
  } else {
    if (!is.character(plot_axis_var) || length(plot_axis_var) != 1 || plot_axis_var == "") {
      stop("'plot_axis_var' must be a single non-empty character string.")
    }
  }

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

    if (!all(required_q_cols %in% names(obj$summary_df))) {
      stop("Object '", obj_name, "'$summary_df is missing required quantile columns ('",
           paste(required_q_cols, collapse = "', '"), "') for ci_level = ", ci_level, ".")
    }
    
    if (!(plot_axis_var %in% names(obj$summary_df))) {
      stop("The specified 'plot_axis_var' ('", plot_axis_var, 
           "') was not found in the summary_df of object '", obj_name, "'.")
    }
    
    if (!is.numeric(obj$summary_df[[plot_axis_var]])){
      stop("The column '", plot_axis_var, "' in object '", obj_name, 
           "' must be numeric for plotting.")
    }
    
  }

  # coef_index is checked against the first object.
  first_obj_df <- summary_list[[1]]$summary_df
  max_coef <- max(first_obj_df$coefficient_index, na.rm = TRUE)
  if (!is.numeric(coef_index) || length(coef_index) != 1 ||
      coef_index < 1 || coef_index > max_coef) {
    stop("'coef_index' must be a single integer between 1 and ", max_coef, ".")
  }

  plot_data <- dplyr::bind_rows(lapply(summary_list, `[[`, "summary_df"),
                                .id = "model_name")

  # Legend in input order.
  plot_data$model_name <- factor(plot_data$model_name, levels = names(summary_list))

  plot_data <- dplyr::filter(plot_data, .data$coefficient_index == !!coef_index)

  if (nrow(plot_data) == 0) {
      stop("No data remaining after filtering for coef_index = ", coef_index)
  }

  if (is.null(xlab)) {
    xlab <- plot_axis_var
  }
  if (is.null(ylab)) {
    ylab <- paste("Avg. Effect (Coefficient", coef_index, ")")
  }
  if (is.null(title)) {
    title <- paste("Marginal Effect on Coefficient", coef_index,
                   "vs", plot_axis_var)
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[plot_axis_var]]))

  p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[lower_q_name]],
                   ymax = .data[[upper_q_name]],
                   fill = .data$model_name),
      alpha = 0.2)

  p <- p + ggplot2::geom_line(
      ggplot2::aes(y = .data$mean,
                   color = .data$model_name,
                   linetype = .data$model_name))

  if (!is.null(color_values)) {
    p <- p + ggplot2::scale_color_manual(values = color_values)
  }
  if (!is.null(linetype_values)) {
    p <- p + ggplot2::scale_linetype_manual(values = linetype_values)
  }
  if (!is.null(fill_values)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_values)
  }

  p <- p + ggplot2::labs(title = title, x = xlab, y = ylab,
                         color = "Model", linetype = "Model", fill = "Model") +
           ggplot2::theme_classic()

  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)

}#PLOT.SUMMARY.MARGINAL_EFFECTS 