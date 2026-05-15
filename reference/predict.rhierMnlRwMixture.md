# Predict Method for rhierMnlRwMixture Objects

Predict Method for rhierMnlRwMixture Objects

## Usage

``` r
# S3 method for class 'rhierMnlRwMixture'
predict(
  object,
  newdata = NULL,
  type = "DeltaZ+mu",
  burn = 0,
  nsim = 10,
  mode = "coefficients",
  r_verbose = TRUE,
  force_tree_eval = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted rhierMnlRwMixture object.

- newdata:

  Optional list containing data for prediction. Structure depends on
  `mode` and `type`:

  - For `mode = "coefficients"` with
    `type %in% c("DeltaZ", "DeltaZ+mu", "SigmaZ")`: requires
    `newdata$Z`, a matrix with `npred` rows for prediction units (if
    model was fit with Z).

  - For `mode = "posterior"` with `type = "choice_probs"`: requires
    `newdata$nlgtdata`, a list of length `nlgt` (original number of
    units). Each element `\\[[i]]` must contain `$X`, the design matrix
    `(T_i*p) x nvar` for unit `i`. Also requires `newdata$p`, the number
    of alternatives.

  - For `mode = "prior"` with `type = "choice_probs"`: requires
    `newdata$Z` (if model fit with Z, determining `npred`), `newdata$p`,
    and `newdata$X` (a list of length `npred`, each element `\\[[i]]`
    having design matrix `(T_i*p) x nvar`).

- type:

  Type of prediction within the selected `mode`:

  - `"DeltaZ"`: Expected part-worths of the representative consumer,
    \\\Delta(Z)\\.

  - `"DeltaZ+mu"`: Expected part-worths plus the mean of the unobserved
    heterogeneity component, \\\Delta(Z) + \mu_1\\. The package supports
    only `ncomp = 1`.

  - `"choice_probs"`: Predictive choice probabilities (use with
    `mode = "posterior"` or `mode = "prior"`).

  - `"SigmaZ"`: Draws of the heteroscedastic covariance matrix
    \\\Sigma(Z)\\. Available only for models fit with `Prior$vartree`
    (class marker `"bayesm.HART.HeterCov"`).

- burn:

  Integer, number of initial MCMC draws to discard.

- nsim:

  Integer, number of draws from the heterogeneity distribution per
  posterior draw for `mode = "prior"` and `type = "choice_probs"`.

- mode:

  Prediction mode:

  - `"coefficients"`: coefficient-level outputs (`DeltaZ`, `DeltaZ+mu`,
    `SigmaZ`).

  - `"posterior"`: posterior predictive output
    (`type = "choice_probs"`).

  - `"prior"`: prior predictive output (`type = "choice_probs"`).

- r_verbose:

  Logical, print progress updates?

- ...:

  Additional arguments passed to underlying prediction functions (e.g.,
  `mc.cores`, `verbose` for BART `DeltaZ` predictions via `pwbart`).

## Value

Depends on `type`:

- For `type %in% c("DeltaZ", "DeltaZ+mu")`: 3D array
  `[npred, nvar, ndraws_out]` of predicted expected part-worths.

- For `type = "SigmaZ"`: 4D array `[npred, nvar, nvar, ndraws_out]` of
  covariance draws at each prediction unit.

- For `mode = "posterior", type = "choice_probs"`: List of length
  `nlgt`. Each element `\\[[i]]` is a 3D array `[T_i, p, ndraws_out]` of
  posterior predictive choice probabilities for unit `i`.

- For `mode = "prior", type = "choice_probs"`: List of length `npred`.
  Each element `\\[[i]]` is a 3D array `[T_i, p, ndraws_out]` of prior
  predictive choice probabilities for prediction unit `i`.
