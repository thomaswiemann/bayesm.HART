# bayesm.HART 0.1.0

* `Prior$hart_model` selects the heterogeneity model: `"none"` (bayesm), `"normalized"`, or `"direct"`. Tree lists in `Prior` require it.
* Direct HART: vector-leaf mean trees with learned scales and covariate-dependent covariance (`Prior$vartree`, `Prior$phitree`); `Prior$hart_cov = "constant"` gives one shared covariance.
* `previous` continues a fit, optionally appending observations and units.
* `predict(type = "SigmaZ")` for every model.
* `marginal_effects()` takes `r_verbose`, as the samplers do.
* `tau` is the ensemble scale in every tree list; `nmin` sets the minimum units per leaf.
* `Sigma_func` in the simulators for covariate-dependent heterogeneity.
* Fixes the birth-proposal reverse probability inherited from BART.
* Fixes the covariance update of the normalized model.
* Fits saved by earlier versions are not supported by `predict()`.
