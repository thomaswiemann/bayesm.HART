# Bayesian Multinomial Logit Model with HART Prior

`rhierMnlRwMixture` implements an MCMC algorithm for a Bayesian
hierarchical multinomial logit model with a Hierarchical Additive
Regression Trees (HART) prior. HART is a hierarchical nonparametric
prior that allows for flexible modeling of the representative consumer
as a function of potentially many observed characteristics. `Prior$bart`
allows for modeling the conditional mean of the normal prior.

## Usage

``` r
rhierMnlRwMixture(Data, Prior, Mcmc, r_verbose = TRUE)
```

## Arguments

- Data:

  A list containing:

  - `p`: Number of choice alternatives (integer).

  - `lgtdata`: A list of length `nlgt`. Each element `lgtdata[[i]]` must
    be a list with:

    - `y`: `n_i x 1` vector of multinomial outcomes (1 to `p`).

    - `X`: `(n_i * p) x nvar` matrix of alternative-specific attributes.

  - `Z` (optional): `nlgt x nz` matrix of observed characteristics for
    each unit. Should NOT contain an intercept and should be centered.

- Prior:

  A list containing prior parameters:

  - `ncomp` (required): Number of components. Must be 1.

  - `deltabar` (optional): `nz * nvar x 1` prior mean for `vec(Delta)`
    (default: 0). Ignored if HART is used.

  - `Ad` (optional): Prior precision matrix for `vec(Delta)` (default:
    `0.01 * I`). Ignored if HART is used.

  - `mubar` (optional): `nvar x 1` prior mean vector for the normal
    prior (default: 0 if unrestricted, 2 if restricted).

  - `Amu` (optional): Prior precision for the normal prior (default:
    0.01 if unrestricted, 0.1 if restricted).

  - `nu` (optional): Degrees of freedom for IW prior on component
    `Sigma` (default: `nvar+3` if unrestricted, `nvar+15` if
    restricted).

  - `V` (optional): Location matrix for IW prior on component `Sigma`
    (default: `nu * I` or scaled based on restriction).

  - `SignRes` (optional): `nvar x 1` vector of sign restrictions. Must
    contain values of 0, -1, or 1. The value 0 means no restriction, -1
    ensures the coefficient is negative, and 1 ensures the coefficient
    is positive. For example, `SignRes = c(0,1,-1)` means the first
    coefficient is unconstrained, the second will be positive, and the
    third will be negative. Default: `rep(0, nvar)`.

  - `bart` (optional): List of parameters for the HART prior. If
    specified, this models the representative consumer \\\Delta(Z)\\ as
    a scaled sum-of-trees factor model. See Details.

  - `vartree` (optional, *experimental extension beyond Wiemann 2025*):
    List of parameters enabling heteroscedastic covariance
    \\\Sigma(Z_i)\\ via product-of-trees variance models on the Modified
    Cholesky diagonal \\d_j(\cdot)\\. Requires `bart` to also be
    specified and `ncomp = 1`. Compatible with `SignRes` (see
    "Heteroscedastic Covariance" in Details). When `nvar > 1`, the
    package automatically promotes to the full-Cholesky structure
    described under `phitree` (diagonal-only \\\Sigma(Z_i)\\ would force
    every conditional cross-correlation to zero, which is essentially
    never a believable posterior).

  - `phitree` (optional, *experimental extension beyond Wiemann 2025*):
    List of parameters enabling sum-of-trees regression models on the
    Modified Cholesky off-diagonals \\\phi\_{jk}(\cdot)\\. Requires
    `vartree`. Adds the full-Cholesky structure on top of the diagonal
    `vartree` model. Auto-enabled when `vartree` is supplied and
    `nvar > 1`.

- Mcmc:

  A list containing MCMC parameters:

  - `R`: Number of MCMC iterations (required).

  - `keep` (optional): Thinning parameter (default: 1).

  - `nprint` (optional): Print progress every `nprint` draws (default:
    100, 0 for none).

  - `s` (optional): RW Metropolis scaling parameter (default:
    `2.93 / sqrt(nvar)`).

  - `w` (optional): Fractional likelihood weighting parameter (default:
    0.1).

- r_verbose:

  Logical. Print startup messages? Default TRUE.

## Value

A list containing:

- `Deltadraw`: If `Z` provided and `bart=NULL`, `(R/keep) x (nz * nvar)`
  matrix of `vec(Delta)` draws.

- `betadraw`: `nlgt x nvar x (R/keep)` array of unit-level `beta_i`
  draws.

- `nmix`: Legacy list containing prior draws (omitted when `vartree` is
  supplied).

- `loglike`: `(R/keep) x 1` vector of log-likelihood values at kept
  draws.

- `SignRes`: `nvar x 1` vector of sign restrictions used.

- `acceptrbeta`: Metropolis acceptance rate (percent) for the unit-level
  MNL random-walk updates of `beta_i`.

- `bart_models`: If HART used, list of length `nvar` containing tree
  structures and related parameters for the mean trees
  \\\delta_j(\cdot)\\.

- `var_models` (only with `vartree`): list of length `nvar` of
  variance-tree ensembles for \\d_j(\cdot)\\ (product of trees with
  \\\chi^{-2}\\ leaves).

- `phi_models` (only with `vartree` + `phitree`): jagged list of length
  `nvar`. `phi_models[[1]]` is `NULL`; for `j > 1`, `phi_models[[j]]` is
  a list of length `j-1` containing the sum-of-trees ensemble for
  \\\phi\_{jk}(\cdot)\\.

- `mu_draw` (only with `vartree`): `(R/keep) x nvar` matrix of \\\mu\\
  draws (single-component path).

## Details

### Model Specification

\\y_i \sim MNL(X_i, \beta_i)\\ for unit \\i = 1, ..., nlgt\\. The
unit-level coefficients (part-worths) \\\beta_i\\ are modeled as:
\$\$\beta_i = \Delta(Z_i) + u_i\$\$ where \\\Delta(Z_i)\\ is the
*representative consumer* component, which depends on observed
characteristics \\Z_i\\, and \\u_i\\ is the unobserved heterogeneity
component.

The representative consumer component is specified as:

- If `Z` is provided and `Prior$bart` is `NULL`: \\\Delta(Z_i) = Z_i
  \Delta\\ where \\\Delta\\ is an `nz x nvar` matrix (linear
  hierarchical model).

- If `Z` is provided and `Prior$bart` is a list: \\\Delta(Z_i)\\ is
  modeled with a HART prior (scaled sum-of-trees factor model).

- If `Z` is `NULL`: \\\Delta(Z_i) = 0\\.

With `ncomp = 1` (currently required), the unobserved heterogeneity
component follows: \$\$u_i \sim N(\mu_1, \Sigma_1)\$\$

### Prior Specifications

- **Linear model**: \\\delta = vec(\Delta) \sim N(deltabar, A_d^{-1})\\

- **Component means**: \\\mu_1 \sim N(mubar, \Sigma_1 \otimes
  Amu^{-1})\\ (covariance scaled by \\\Sigma_1\\)

- **Component covariance**: \\\Sigma_1 \sim IW(\nu, V)\\

- **HART model**: A sum-of-trees prior is placed on each factor of the
  scaled sum-of-trees model (see HART details below).

### HART Prior Details

If `Prior$bart` is a list, it specifies a HART prior for the
representative consumer \\\Delta(Z)\\. This replaces the conventional
linear hierarchical specification. The HART prior models the
representative consumer using a scaled vector of `nvar` sum-of-trees
models.

**HART Parameters** (defaults used if not specified in `Prior$bart`):

- `num_trees`: Number of trees H in each sum-of-trees model (default:
  200).

- `power`, `base`: Parameters for the tree structure prior. The
  probability of a node at depth `q` splitting is
  \\\alpha(1+q)^{-\beta}\\, where `base`=\\\alpha\\ and
  `power`=\\\beta\\. Defaults are `base=0.95`, `power=2`, which strongly
  favors shallow trees.

- `tau`: Parameter controlling the prior variance of terminal leaf
  coefficients. The default is \\\tau = 1/\sqrt{H}\\ where
  \\\lambda\_{dhg} \sim N(0, \tau^2)\\ for terminal leaf coefficients.

- `numcut`: Number of grid points for proposing splitting rules for
  continuous variables (default: 100).

- `sparse`: If `TRUE`, use the Dirichlet HART prior to induce sparsity
  in variable selection (default: `FALSE`).

**Dirichlet HART** (`sparse = TRUE`): The Dirichlet HART model augments
the HART prior to induce sparsity in variable selection, following
Linero (2018). Instead of uniform probability for selecting splitting
variables, the selection probabilities \\\tau = (\tau^{(1)}, \ldots,
\tau^{(K)})\\ are given a sparse Dirichlet prior: \\(\tau^{(1)}, \ldots,
\tau^{(K)}) \sim Dirichlet(\theta/K, \ldots, \theta/K)\\, where K is the
number of characteristics. The concentration parameter \\\theta\\ is
given a hierarchical prior: \\\theta/(\theta+\rho) \sim Beta(a,b)\\.

- `a`, `b`: Shape parameters for the Beta hyperprior. The default
  (`a=0.5, b=1`) induces sparsity where few variables have high
  selection probabilities.

- `rho`: Parameter influencing sparsity. Default is the number of
  characteristics K. Reducing rho below K encourages greater sparsity.

- `theta`: When used, sets Dirichlet concentration parameter without
  additional hyper-prior (default: 0.0).

- `burn`: Number of internal burn-in iterations for the Dirichlet HART
  sampler before variable selection is allowed (default: 100).

### Heteroscedastic Covariance \\\Sigma(Z_i)\\ (experimental extension)

Wiemann (2025) treats \\\Sigma\\ as global with an inverse-Wishart
prior, used as a fixed factor-model loading \\\Sigma^{1/2}\\. When
`Prior$vartree` is supplied, the package replaces this with
unit-specific \$\$\Sigma(Z_i)^{-1} = L(Z_i)^\top D(Z_i)^{-1} L(Z_i)\$\$
where \\L(Z_i)\\ is unit lower-triangular with \\L\_{jk} =
-\phi\_{jk}(Z_i)\\ for \\k \< j\\ and \\D(Z_i) = \mathrm{diag}(d_1(Z_i),
\ldots, d_D(Z_i))\\, \\d_j \> 0\\. Each \\d_j(\cdot)\\ is modeled as a
product of trees with \\\chi^{-2}\\ leaves (Pratola et al., 2020); each
\\\phi\_{jk}(\cdot)\\ (when `Prior$phitree` is supplied) is modeled as a
sum of trees with \\N(0, \tau^2)\\ leaves. Only the single-component
path (`ncomp = 1`) is supported in this mode.

**`Prior$vartree` parameters** (defaults shown):

- `num_trees` (40): Number of trees per dimension \\m'\\.

- `nu` (10), `lambda` (auto-calibrated from \\\mathrm{var}(\theta_i)\\):
  Baseline parameters of the \\\chi^{-2}\\ prior. Pratola per-tree
  calibration is applied internally.

- `power` (2), `base` (0.95), `numcut` (100): Tree-structure prior
  parameters.

- DART hyperparameters (`sparse`, `a`, `b`, `rho`, `theta`, `omega`,
  `aug`, `burn`): same meaning as `Prior$bart`.

**`Prior$phitree` parameters** (defaults shown):

- `num_trees` (40): Number of trees per \\(j,k)\\ pair \\m''\\.

- `tau` (\\1/\sqrt{m''}\\): Prior standard deviation of leaf
  \\\lambda\_{jkh}\\.

- `power` (2), `base` (0.95), `numcut` (100): Tree-structure prior
  parameters.

- `nmin` (2), `ess_min` (5): Leaf-admissibility floors. The default
  `ess_min = 5` requires the effective sample size \\\sum\_{i \in \ell}
  (\theta_i^{(k)} - \mu^{(k)})^2 / d_j(Z_i)\\ of each candidate leaf to
  exceed 5; `nmin = 2` is the corresponding raw-count floor.

- DART hyperparameters: same meaning as `Prior$bart`.

When this extension is active, the returned object additionally inherits
class `"rhierMnlRwMixtureHeterCov"` and contains slots `var_models`,
`phi_models` (jagged, lower-triangular; `NULL` if `nvar == 1`), and
`mu_draw` (single-component posterior mean draws).
[`predict()`](https://rdrr.io/r/stats/predict.html) dispatches on this
class to evaluate \\\Sigma(Z^\*)\\ at new \\Z^\*\\.

### Sign Restrictions

If `SignRes[k]` is non-zero, the k-th coefficient \\\beta\_{ik}\\ is
modeled as \$\$\beta\_{ik} = SignRes\[k\] \cdot
\exp(\beta^\*\_{ik}).\$\$ The `betadraw` output contains the draws for
\\\beta\_{ik}\\ (with the restriction applied). The `nmix` output
contains draws for the *unrestricted* prior covariance and mean.

**Note:** Care should be taken when selecting priors on any sign
restricted coefficients.

## Note

Currently, only `ncomp = 1` is supported.

## References

Chipman, Hugh A., Edward I. George, and Robert E. McCulloch (2010).
"BART: Bayesian Additive Regression Trees." Annals of Applied Statistics
4.1.

Linero, Antonio R. (2018). "Bayesian regression trees for
high-dimensional prediction and variable selection." Journal of the
American Statistical Association 113.522, pp. 626-636.

Pratola, M. T., Chipman, H. A., George, E. I., and McCulloch, R. E.
(2020). "Heteroscedastic BART via Multiplicative Regression Trees."
Journal of Computational and Graphical Statistics 29.2, pp. 405-417.

Rossi, Peter E., Greg M. Allenby, and Robert McCulloch (2009). Bayesian
Statistics and Marketing. Reprint. Wiley Series in Probability and
Statistics. Chichester: Wiley.

Rossi, Peter (2023). bayesm: Bayesian Inference for
Marketing/Micro-Econometrics. Comprehensive R Archive Network.

Wiemann, Thomas (2025). "Personalization with HART." Working paper.

## See also

[`predict.rhierMnlRwMixture()`](https://thomaswiemann.com/bayesm.HART/reference/predict.rhierMnlRwMixture.md)

## Author

Peter Rossi (original bayesm code), Thomas Wiemann (HART modifications).
