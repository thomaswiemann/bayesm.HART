# Summarize Marginal Effects Object

Computes posterior means and quantiles from the draws stored in a
`marginal_effects` object.

## Usage

``` r
# S3 method for class 'marginal_effects'
summary(object, probs = c(0.025, 0.05, 0.5, 0.95, 0.975), ...)
```

## Arguments

- object:

  An object of class `"marginal_effects"` created by
  `marginal_effects.rhierMnlRwMixture`.

- probs:

  Numeric vector of quantile probabilities (between 0 and 1) to compute.

- ...:

  Additional arguments (currently ignored).

## Value

A list object of class `"summary.marginal_effects"` containing:

- summary_df:

  A data frame with columns `coefficient_index`, columns for each non-NA
  column in the input `z_values` (e.g., `z_col_1`, `z_col_2`, ... or
  using original variable names if available), `mean`, and columns for
  each requested quantile (e.g., `q2.5`, `q50`, `q97.5`). Each row
  represents a specific coefficient at a specific grid point defined by
  a row in `z_values`.

- z_values:

  The original `z_values` matrix provided to `marginal_effects`.

- probs:

  The numeric vector of quantile probabilities used.

- call:

  The matched call to the `summary` function.

## Examples

``` r
# --- Assumes `mfx_result` exists from marginal_effects.rhierMnlRwMixture example ---
if (exists("mfx_result")) {
  mfx_summary <- summary(mfx_result, probs = c(0.05, 0.5, 0.95))
  print(mfx_summary$summary_df)
  # --- Plotting example follows in plot.summary.marginal_effects ---

} else {
  message("Run the example for 'marginal_effects.rhierMnlRwMixture' first.")
}
#> Run the example for 'marginal_effects.rhierMnlRwMixture' first.
```
