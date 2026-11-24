# bayesm.HART 0.1.0

* `Prior$hart_model`: `"none"` (bayesm), `"normalized"`, or `"direct"`. Required with tree lists.
* Direct HART: vector-leaf mean trees and covariance trees (`vartree`, `phitree`); `hart_cov = "constant"` for one shared covariance.
* `previous`: continue a fit, optionally with new observations and units.
* `predict(type = "SigmaZ")` for all models.
* `r_verbose` in `marginal_effects()`.
* Tree lists: `tau` (ensemble scale), `nmin` (minimum units per leaf).
* `Sigma_func` in the simulators.
* Fixed: birth-proposal reverse probability inherited from BART.
* Fixed: covariance update of the normalized model.
* Fixed: `summary()` R-squared under a constant covariance used the wrong residual variance.
* Breaking: `predict()` does not support fits from earlier versions.
