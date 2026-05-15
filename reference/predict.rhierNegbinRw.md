# Predict Method for rhierNegbinRw Objects

Computes coefficient-level or predictive count-response draws.

## Usage

``` r
# S3 method for class 'rhierNegbinRw'
predict(
  object,
  newdata = NULL,
  type = "DeltaZ+mu",
  burn = 0,
  mode = "coefficients",
  nsim = 10,
  r_verbose = TRUE,
  force_tree_eval = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted `rhierNegbinRw` object.

- newdata:

  A list whose required fields depend on `mode`:

  - For `mode = "coefficients"`: requires `newdata$Z` (if the model was
    fit with Z).

  - For `mode = "posterior"`: requires `newdata$regdata`, a list with
    per-unit design matrices `X` matching the fitted units.

  - For `mode = "prior"`: requires `newdata$X`, a list of design
    matrices for prediction units, and `newdata$Z` when the model was
    fit with Z.

- type:

  Character; interpretation depends on `mode`:

  - For `mode = "coefficients"`: one of `"DeltaZ"`, `"DeltaZ+mu"`,
    `"SigmaZ"`.

  - For `mode %in% c("prior","posterior")`: must be `"response"`.

- burn:

  Integer, number of initial MCMC draws to discard.

- mode:

  Prediction mode: `"coefficients"` (default), `"posterior"`, or
  `"prior"`.

- nsim:

  Integer, number of Monte Carlo draws for prior predictive response
  mixing.

- r_verbose:

  Logical, print progress updates?

- ...:

  Additional arguments passed to `pwbart` for BART models.

## Value

Depends on `mode` / `type`:

- For `mode = "coefficients"`, `type %in% c("DeltaZ", "DeltaZ+mu")`: 3D
  array `[npred, ncoef, ndraws_out]` of predicted betabar values.

- For `mode = "coefficients"`, `type = "SigmaZ"`: 4D array
  `[npred, ncoef, ncoef, ndraws_out]` of covariance draws at each
  prediction unit.

- For `mode %in% c("prior","posterior")`, `type = "response"`: list of
  matrices with simulated count draws per unit (`nobs_i x ndraws_out`).
