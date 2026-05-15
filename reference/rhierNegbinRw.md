# Bayesian Hierarchical Negative Binomial Model with HART Prior

`rhierNegbinRw` implements an MCMC algorithm for a Bayesian hierarchical
negative binomial regression model with a Hierarchical Additive
Regression Trees (HART) prior. HART is a hierarchical nonparametric
prior that allows for flexible modeling of the representative unit as a
function of potentially many observed characteristics. Per-unit
\\\beta_i\\ are drawn via Random Walk Metropolis-Hastings; the
dispersion \\\alpha\\ is drawn on the log scale. `Prior$bart` allows for
modeling the conditional mean of the normal prior.

## Usage

``` r
rhierNegbinRw(Data, Prior, Mcmc, r_verbose = TRUE)
```

## Arguments

- Data:

  A list containing:

  - `regdata`: A list of length `nreg`. Each element `regdata[[i]]` must
    be a list with:

    - `y`: `n_i x 1` vector of count responses.

    - `X`: `n_i x nvar` design matrix (including intercept column).

  - `Z` (optional): `nreg x nz` matrix of unit-level covariates. If
    omitted, `drawdelta` is set to FALSE and no second-stage covariates
    are used.

- Prior:

  A list containing prior parameters:

  - `ncomp` (required): Number of components. Must be 1.

  - `nu` (optional): Degrees of freedom for IW prior on the normal prior
    covariance (default: nvar + 3).

  - `V` (optional): `nvar x nvar` location matrix for IW prior (default:
    nu \* I).

  - `a` (optional): Shape parameter for gamma prior on alpha (default:
    0.5). Note: This differs from `rhierLinearMixture` where `a` is the
    Dirichlet prior.

  - `b` (optional): Rate parameter for gamma prior on alpha (default:
    0.1).

  - `mubar` (optional): `1 x nvar` prior mean vector for the normal
    prior (default: 0).

  - `Amu` (optional): `1 x 1` prior precision for the normal prior
    (default: 0.01).

  - `Ad` (optional): `nvar*nz x nvar*nz` prior precision for vec(Delta)
    (default: 0.01 \* I). Only used when `drawdelta = TRUE` and
    `useBART = FALSE`.

  - `deltabar` (optional): `nvar*nz` prior mean for vec(Delta) (default:
    0). Only used when `drawdelta = TRUE` and `useBART = FALSE`.

  - `bart` (optional): List of HART prior parameters. See
    `rhierLinearMixture` for details.

  - `vartree` (optional, *experimental extension beyond Wiemann 2025*):
    List of parameters enabling heteroscedastic covariance
    \\\Sigma(Z_i)\\ via product-of-trees variance models on the Modified
    Cholesky diagonal \\d_j(\cdot)\\. Requires `bart` to also be
    specified and `ncomp = 1`. When `nvar > 1`, the package
    automatically promotes to the full-Cholesky structure described
    under `phitree`. See "Heteroscedastic Covariance" in Details.

  - `phitree` (optional, *experimental extension beyond Wiemann 2025*):
    List of parameters enabling sum-of-trees regression models on the
    Modified Cholesky off-diagonals \\\phi\_{jk}(\cdot)\\. Requires
    `vartree`. Auto-enabled when `vartree` is supplied and `nvar > 1`.

- Mcmc:

  A list containing MCMC parameters:

  - `R`: Number of MCMC iterations (required).

  - `keep` (optional): Thinning parameter (default: 1).

  - `nprint` (optional): Print progress every `nprint` draws (default:
    100, 0 for none).

  - `s_beta` (optional): RW scaling for beta (default: 2.93/sqrt(nvar)).

  - `s_alpha` (optional): RW scaling for alpha (default: 2.93).

  - `w` (optional): Fractional likelihood weight used in candidate
    Hessian construction (default: 0.1).

  - `alpha` (optional): Fixed alpha value. If provided, alpha is not
    sampled (sets `fixalpha = TRUE`).

  - `fixalpha` (optional): Logical, whether to fix alpha (default:
    FALSE).

- r_verbose:

  Logical. Print startup messages? Default TRUE.

## Value

A list of class `"rhierNegbinRw"` containing:

- `betadraw`: `nreg x nvar x (R/keep)` array of unit-level beta draws.

- `alphadraw`: `(R/keep) x 1` vector of overdispersion draws.

- `loglike`: `(R/keep) x 1` vector of log-likelihoods.

- `nmix`: Legacy list containing prior draws (omitted under heter-cov).

- `acceptrbeta`, `acceptralpha`: Metropolis acceptance rates (percent).

- If `drawdelta` and non-BART: `Deltadraw`.

- If BART: `bart_models`, `varcount`, `varprob`.

- If heter-cov (additional class `"rhierNegbinRwHeterCov"`):
  `var_models`, `phi_models` (jagged or `NULL`), `mu_draw`,
  `var_varcount`, `var_varprob`.

## Details

### Model Specification

\\(y_i\|\lambda_i, \alpha) \sim NegBin(\lambda_i, \alpha)\\,
\\\ln(\lambda_i) = X_i \beta_i\\

The unit-level coefficients are modeled as: \\\beta_i = Z_i \Delta +
u_i\\  
\\u_i \sim N(0, \Sigma_i)\\  
Note: This documentation describes `ncomp = 1` behavior as higher-order
components are not supported.  

### HART Prior Details

If `Prior$bart` is a list, \\\Delta(Z_i)\\ is modeled via a sum-of-trees
(BART) prior instead of a linear hierarchical specification. See
`rhierMnlRwMixture` for complete details and hyperparameter definitions.

### Heteroscedastic Covariance \\\Sigma(Z_i)\\ (experimental extension)

When `Prior$vartree` is supplied, the homoscedastic prior covariance
\\\Sigma\\ is replaced with a unit-specific modified Cholesky
decomposition where the components are modeled as tree ensembles. See
`rhierMnlRwMixture` for complete details.

When this extension is active, the returned object additionally inherits
classes `"rhierNegbinRwHeterCov"` and `"bayesm.HART.HeterCov"`.
[`predict()`](https://rdrr.io/r/stats/predict.html) dispatches on these
classes to evaluate \\\Sigma(Z^\*)\\ at any new \\Z^\*\\.

## Note

Currently, only `ncomp = 1` is supported.

## See also

[`predict.rhierNegbinRw()`](https://thomaswiemann.com/bayesm.HART/reference/predict.rhierNegbinRw.md),
[`marginal_effects.rhierNegbinRw()`](https://thomaswiemann.com/bayesm.HART/reference/marginal_effects.md)

## Author

Sridhar Narayanan, Peter Rossi (original bayesm code), Thomas Wiemann
(HART modifications).

## Examples

``` r
# \donttest{
set.seed(20260513)
sim <- bayesm.HART::sim_hier_negbin(
  nreg = 30, nobs = 12, nvar = 1, nz = 2,
  const = TRUE, het_observed = "linear",
  target_var_betabar = 0.5, target_var_eps = 0.25,
  alpha = 5.0)
Prior <- list(ncomp = 1L,
              bart    = list(num_trees = 20),
              vartree = list(num_trees = 20))
Mcmc  <- list(R = 200L, keep = 1L, nprint = 0L)
fit   <- bayesm.HART::rhierNegbinRw(
  Data = list(regdata = sim$regdata, Z = sim$Z),
  Prior = Prior, Mcmc = Mcmc, r_verbose = FALSE)
str(fit, max.level = 1)
#> List of 13
#>  $ mu_draw     : 'bayesm.mat' num [1:200, 1:2] 0.2497 0.0736 0.2043 0.0518 0.0238 ...
#>   ..- attr(*, "mcpar")= num [1:3] 1 200 1
#>  $ varcount    : num [1:2, 1:2, 1:200] 7 8 3 7 7 7 6 9 7 10 ...
#>  $ varprob     : num [1:2, 1:2, 1:200] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
#>  $ var_varcount: num [1:2, 1:2, 1:200] 9 6 7 8 8 7 8 9 9 7 ...
#>  $ var_varprob : num [1:2, 1:2, 1:200] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
#>  $ bart_models :List of 2
#>  $ var_models  :List of 2
#>  $ phi_models  :List of 2
#>  $ betadraw    : 'bayesm.hcoef' num [1:30, 1:2, 1:200] 0.537 0.181 0.161 0.917 0.557 ...
#>  $ loglike     : num [1:200, 1] -576 -583 -584 -580 -580 ...
#>  $ acceptrbeta : num 36.6
#>  $ acceptralpha: num 74.5
#>  $ alphadraw   : 'bayesm.mat' num [1:200, 1] 0.87 0.851 0.851 0.954 0.954 ...
#>   ..- attr(*, "mcpar")= num [1:3] 1 200 1
#>  - attr(*, "class")= chr [1:3] "rhierNegbinRwHeterCov" "bayesm.HART.HeterCov" "rhierNegbinRw"
#>  - attr(*, "hart_cache")=List of 6
# }
```
