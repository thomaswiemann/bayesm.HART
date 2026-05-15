# Calculate Marginal Effects for Hierarchical Models

Computes the posterior distribution of average marginal effects by
varying a target covariate over a grid.

## Usage

``` r
marginal_effects(object, z_values, Z, burn = 0, verbose = TRUE, ...)

# S3 method for class 'rhierMnlRwMixture'
marginal_effects(object, z_values, Z, burn = 0, verbose = TRUE, ...)

# S3 method for class 'rhierLinearMixture'
marginal_effects(object, z_values, Z, burn = 0, verbose = TRUE, ...)

# S3 method for class 'rhierNegbinRw'
marginal_effects(object, z_values, Z, burn = 0, verbose = TRUE, ...)
```

## Arguments

- object:

  A fitted hierarchical model object (e.g., `rhierMnlRwMixture`,
  `rhierLinearMixture`, `rhierNegbinRw`).

- z_values:

  A numeric matrix of grid values for the unit-level covariates `Z`.
  Each row defines one counterfactual `Z*`. Columns that are entirely
  `NA` are held at their training values; columns with no `NA`s are
  swept to the supplied grid value.

- Z:

  A numeric matrix of unit-level covariates from the training sample
  (typically `Data$Z` from the original fit). Must have the same number
  of columns as `z_values`.

- burn:

  Non-negative integer. Number of initial MCMC draws to drop before
  averaging (default `0`).

- verbose:

  Logical. Print progress per grid point (default `TRUE`).

- ...:

  Other arguments passed to methods.

## Value

An object of class `"marginal_effects"`; see
[`summary.marginal_effects()`](https://thomaswiemann.com/bayesm.HART/reference/summary.marginal_effects.md)
for downstream summarization.

## Details

This is a generic function. Method implementations are provided for
`rhierMnlRwMixture`, `rhierLinearMixture`, and `rhierNegbinRw` (the
heter-cov subclasses inherit dispatch through their base classes).
