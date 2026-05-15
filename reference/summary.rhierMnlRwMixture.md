# Summary Method for Hierarchical Models

Computes model diagnostics including acceptance rates, posterior mean
betas, BART variable importance, and Bayesian R².

## Usage

``` r
# S3 method for class 'rhierMnlRwMixture'
summary(object, Z = NULL, burn = 0, coefs = NULL, r_verbose = FALSE, ...)

# S3 method for class 'rhierLinearMixture'
summary(object, Z = NULL, burn = 0, coefs = NULL, r_verbose = FALSE, ...)

# S3 method for class 'rhierNegbinRw'
summary(object, Z = NULL, burn = 0, coefs = NULL, r_verbose = FALSE, ...)
```

## Arguments

- object:

  A fitted hierarchical model object.

- Z:

  Optional matrix of unit-level characteristics for R² computation. If
  `NULL` and cached DeltaZ draws exist, in-sample R² is still computed.

- burn:

  Number of initial draws to discard (thinned units).

- coefs:

  Integer vector of coefficient indices for R². Default: all.

- r_verbose:

  Print progress? Default FALSE.

- ...:

  Ignored.

## Value

An object of class `summary.rhierModel` with print method.
