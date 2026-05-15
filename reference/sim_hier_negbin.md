# Simulate Hierarchical Negative Binomial Data

Generates simulated data from a hierarchical negative binomial model
with optional observed heterogeneity via Z covariates.

## Usage

``` r
sim_hier_negbin(
  nreg = 100,
  nobs = 50,
  nvar = 2,
  nz = 3,
  const = TRUE,
  het_observed = c("none", "linear", "step"),
  target_var_betabar = 0.5,
  target_var_eps = 0.25,
  alpha = 5,
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
  "linear", "step".

- target_var_betabar:

  Numeric. Target variance for the first coefficient of the observed
  component. Default 0.5.

- target_var_eps:

  Numeric. Target variance for each coefficient in the unobserved
  component eps_i. Default 0.25.

- alpha:

  Numeric. Overdispersion parameter for negative binomial. Default 5.

- seed:

  Integer. Optional random seed.

## Value

A list containing:

- regdata:

  List of length nreg. Each element has y, X.

- hessdata:

  List of length nreg. Each element has hess (approx Hessian).

- Z:

  nreg x nz matrix of covariates (NULL if nz=0).

- true_values:

  List with beta_true, betabar_true, alpha, etc.
