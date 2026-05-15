# Simulate Hierarchical Linear Model Data

Generates simulated data from a hierarchical linear model with optional
observed heterogeneity (linear or nonlinear via Z covariates).

## Usage

``` r
sim_hier_linear(
  nreg = 100,
  nobs = 10,
  nvar = 2,
  nz = 3,
  const = TRUE,
  het_observed = c("none", "linear", "step", "friedman"),
  target_var_betabar = 1,
  target_var_eps = 0.5,
  sigma_sq = 1,
  seed = NULL
)
```

## Arguments

- nreg:

  Number of cross-sectional units.

- nobs:

  Number of observations per unit.

- nvar:

  Number of X variables (excluding intercept if const=TRUE).

- nz:

  Number of Z variables. Set to 0 for no observed heterogeneity.

- const:

  Logical. Include an intercept in X? Default TRUE.

- het_observed:

  Character. Functional form of observed heterogeneity. Options: "none",
  "linear", "step", "friedman".

- target_var_betabar:

  Numeric. Target variance for the first coefficient of the observed
  component betabar_i = f(Z_i). Default 1.0.

- target_var_eps:

  Numeric. Target variance for each coefficient in the unobserved
  component eps_i. Default 0.5.

- sigma_sq:

  Numeric. Error variance for y_i = X_i \* beta_i + e_i. Default 1.0.

- seed:

  Integer. Optional random seed.

## Value

A list containing:

- regdata:

  List of length nreg. Each element is a list with y, X, XpX, Xpy.

- Z:

  nreg x nz matrix of covariates (NULL if nz=0).

- true_values:

  List with beta_true, betabar_true, eps_true, sigma_sq, Delta (if
  linear).
