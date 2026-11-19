# bayesm.HART 0.1.0

* Adds `previous` to the three hierarchical fitters for continuation on retained data or batches of observations and units, with terminal state stored in `restart`. Each result contains only the new run's draws.
* Adds `Prior$hart_model` to select the preference model: `"none"` (the bayesm paths), `"normalized"` (HART prior on the representative unit), or `"direct"` (HART mean trees in coefficient units with covariate-dependent covariance).
* Adds the direct HART model for `rhierMnlRwMixture()`, `rhierLinearMixture()`, and `rhierNegbinRw()`: one vector-leaf mean ensemble with learned scales (`Prior$bart$a_tau`, `tau0`), product-of-trees variances (`Prior$vartree$tau`, `v0`), and sum-of-trees Cholesky off-diagonals (`Prior$phitree$tau`). Outputs add `mean_scale_draw`; `predict()` supports `type = "SigmaZ"`.
* Adds `Prior$hart_cov` for the direct model: `"heteroskedastic"` (default) keeps the covariance ensembles, `"constant"` draws one `Sigma ~ IW(nu, V)` for all units under the same vector-leaf mean ensemble, so the mean and covariance models can be assessed separately. Both emit `mu_draw`; `"constant"` also carries bayesm's `nmix`.
* `predict(type = "SigmaZ")` works for every model: under `hart_cov = "heteroskedastic"` the draws vary with `Z`, otherwise the single `Sigma` draw is repeated on each row of `newdata$Z`.
* Under a HART model the prior on `mu` is the independent `N(mubar, Amu^{-1})`; `Prior$Amu` accepts a scalar or an `nvar x nvar` matrix.
* `tau` is the ensemble scale in every tree list: the prior standard deviation of each `f*_d(Z)` (`bart`, leaf sd `tau / sqrt(num_trees)`), of each `phi_jk(Z)` (`phitree`), and of each `log d_j(Z)` (`vartree`, leaf shape from `num_trees * trigamma(a_v) = tau^2`).
* Adds `nmin`, the minimum number of units per child leaf, to every tree list.
* Adds `Sigma_func` to `sim_hier_mnl()`, `sim_hier_linear()`, and `sim_hier_negbin()` for covariate-dependent heterogeneity.
* Fixes the birth-proposal reverse probability inherited from BART, which under-accepted births from single-split trees.
* Fixes the normalized-model covariance update, which targeted neither the independent nor the conditional prior on `mu`.
* Compiles with Armadillo's internal OpenMP disabled (`-DARMA_DONT_USE_OPENMP`); BART's prediction functions keep their own OpenMP. With 320 or more units, the elementwise `exp()` in the mixture label draw started a thread team every iteration, which made the bayesm paths an order of magnitude slower when several fits run concurrently.
* Prior predictive draws (`predict(mode = "prior")`) now come from one code path for every model, `mu + Sigma(Z)^{1/2} z`. The draws are distributionally unchanged but differ for a given seed except under `hart_cov = "heteroskedastic"` with the multinomial logit.
* Tree prediction goes through one templated C++ entry for scalar and vector leaves; vector-leaf prediction now uses OpenMP across draws.
* Removes unused internals: BART's standalone sampler and MNL tree predictor, the `heterbartfuns` kernels (the one live score moved into `vcbartfuns.h`), and the BART debug logging.
* Breaking changes:
    - A tree list in `Prior` without `Prior$hart_model` is an error.
    - Direct-model outputs: `bart_models` is one vector-leaf `treedraws` with `leaf_dim`, and `varcount`/`varprob` are `nz x (R/keep)`.
    - `Prior$vartree$nu`, `Prior$vartree$lambda`, and `Prior$phitree$ess_min` are removed; `Prior$nu` and `Prior$V` are not used under `hart_model = "direct"`.
    - Fits saved by earlier versions are not supported by `predict()`.
