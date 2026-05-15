# ==============================================================================
# Heteroscedastic-Covariance Prior Parsing & Validation
# ==============================================================================
# Shared R-side helpers for the heter-cov extension (modified-Cholesky tree
# ensembles for Sigma(Z_i)).  Used by rhierMnlRwMixture and (in Steps 2-3)
# rhierNegbinRw + rhierLinearMixture.  All functions are package-internal.
#
# All defaults documented in the model-specific roxygen blocks (see
# `R/rhierMnlRwMixture_rcpp.R`'s "Heteroscedastic Covariance" details).

# ------------------------------------------------------------------------------
# Flag parsing.  Returns a 3-element list: useHeterCov, useFullCov, auto_full.
#
# Diagonal-only Sigma(Z) is honest only when nvar == 1 (where there is one
# variance to model and no off-diagonals exist).  For nvar > 1, auto-promote
# to full Cholesky even when the user did not supply Prior$phitree -- diagonal
# Sigma(Z) for multivariate beta would force every conditional cross-correlation
# to hard zero at every Z, which is essentially never a believable posterior.
# ------------------------------------------------------------------------------
.parse_heter_cov_flags <- function(Prior, nvar) {
  useHeterCov <- !is.null(Prior$vartree)
  user_phi    <- useHeterCov && !is.null(Prior$phitree)
  useFullCov  <- useHeterCov && (user_phi || nvar > 1)
  auto_full   <- useFullCov && !user_phi
  list(useHeterCov = useHeterCov,
       useFullCov  = useFullCov,
       auto_full   = auto_full)
}

# ------------------------------------------------------------------------------
# Validation.  Called as soon as flags are known (model-specific code may add
# further checks afterwards).  Pass `SignRes = NULL` for models without sign
# restrictions (Negbin, Linear).
# ------------------------------------------------------------------------------
.validate_heter_cov <- function(useHeterCov, useBART, ncomp, drawdelta,
                                SignRes = NULL) {
  if (!useHeterCov) return(invisible(NULL))
  # Order matters: report the most diagnostic missing input first.  Wrappers
  # auto-disable useBART when drawdelta is FALSE (no Z was passed), so the
  # "no Z" check has to fire BEFORE the "no bart" check or users without Z
  # would be told they need bart even when bart was actually supplied.
  if (!drawdelta)
    pandterm("Prior$vartree requires Z (drawdelta = TRUE).")
  if (!useBART)
    pandterm("Prior$vartree requires Prior$bart (heter-cov needs mean trees).")
  if (ncomp != 1)
    pandterm("Prior$vartree requires ncomp == 1 (no mixture under heteroscedastic Sigma).")
  # Sign restrictions are compatible with heter-cov: Sigma(Z) is the
  # heteroscedastic covariance of the unconstrained scale beta* = log|beta|/SignRes,
  # exactly as in the homoscedastic HART path.  The post-loop transform
  # beta = SignRes * exp(beta*) is unconditional in the C++ kernel.
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# var_params: chi^{-2} prior + DART hyperparameters for the diagonal d_j(.) trees.
#
# `oldbetas` is a baseline matrix of per-unit coefficient estimates (MLE for
# MNL/Negbin, OLS for Linear) used solely to auto-calibrate `lambda` when the
# user did not supply one.  `nz` sets the default `rho` (Dirichlet HART
# concentration spread across all characteristics).
# ------------------------------------------------------------------------------
.parse_var_params <- function(prior_vartree, oldbetas, nz) {
  lambda_default <- mean(apply(oldbetas, 2, var))
  if (!is.finite(lambda_default) || lambda_default <= 0) lambda_default <- 1.0
  num_var_trees <- ifelse(is.null(prior_vartree$num_trees), 40L,
                          as.integer(prior_vartree$num_trees))
  list(
    num_trees = num_var_trees,
    power     = ifelse(is.null(prior_vartree$power),  2.0,           prior_vartree$power),
    base      = ifelse(is.null(prior_vartree$base),   0.95,          prior_vartree$base),
    nu        = ifelse(is.null(prior_vartree$nu),     10.0,          prior_vartree$nu),
    lambda    = ifelse(is.null(prior_vartree$lambda), lambda_default,prior_vartree$lambda),
    numcut    = ifelse(is.null(prior_vartree$numcut), 100L,          prior_vartree$numcut),
    sparse    = ifelse(is.null(prior_vartree$sparse), FALSE,         prior_vartree$sparse),
    theta     = ifelse(is.null(prior_vartree$theta),  0.0,           prior_vartree$theta),
    omega     = ifelse(is.null(prior_vartree$omega),  1.0,           prior_vartree$omega),
    a         = ifelse(is.null(prior_vartree$a),      0.5,           prior_vartree$a),
    b         = ifelse(is.null(prior_vartree$b),      1.0,           prior_vartree$b),
    rho       = ifelse(is.null(prior_vartree$rho),    nz,            prior_vartree$rho),
    aug       = ifelse(is.null(prior_vartree$aug),    FALSE,         prior_vartree$aug),
    burn      = ifelse(is.null(prior_vartree$burn),   100L,          prior_vartree$burn)
  )
}

# ------------------------------------------------------------------------------
# phi_params: N(0, tau^2)-leaf prior + DART hyperparameters for the off-diagonal
# phi_{jk}(.) trees.  Returns NULL `tau` if user supplied; otherwise auto-sets to
# 1/sqrt(num_trees) per the standard BART convention.
# ------------------------------------------------------------------------------
.parse_phi_params <- function(prior_phitree, nz, nu, nvar) {
  num_phi_trees <- ifelse(is.null(prior_phitree$num_trees), 40L,
                          as.integer(prior_phitree$num_trees))
  
  # Exact tau calibration matching Inverse-Wishart(nu)
  if (nvar > 1) {
    i_vec <- 2:nvar
    V_phi_bar <- (2 / (nvar * (nvar - 1))) * sum((i_vec - 1) / (nu - nvar + i_vec - 2))
    tau_default <- sqrt(V_phi_bar / num_phi_trees)
  } else {
    tau_default <- 1.0 / sqrt(num_phi_trees)
  }

  list(
    num_trees = num_phi_trees,
    power     = ifelse(is.null(prior_phitree$power),   2.0,                   prior_phitree$power),
    base      = ifelse(is.null(prior_phitree$base),    0.95,                  prior_phitree$base),
    tau       = ifelse(is.null(prior_phitree$tau),     tau_default,           prior_phitree$tau),
    numcut    = ifelse(is.null(prior_phitree$numcut),  100L,                  prior_phitree$numcut),
    sparse    = ifelse(is.null(prior_phitree$sparse),  FALSE,                 prior_phitree$sparse),
    theta     = ifelse(is.null(prior_phitree$theta),   0.0,                   prior_phitree$theta),
    omega     = ifelse(is.null(prior_phitree$omega),   1.0,                   prior_phitree$omega),
    a         = ifelse(is.null(prior_phitree$a),       0.5,                   prior_phitree$a),
    b         = ifelse(is.null(prior_phitree$b),       1.0,                   prior_phitree$b),
    rho       = ifelse(is.null(prior_phitree$rho),     nz,                    prior_phitree$rho),
    aug       = ifelse(is.null(prior_phitree$aug),     FALSE,                 prior_phitree$aug),
    nmin      = ifelse(is.null(prior_phitree$nmin),    2L,                    prior_phitree$nmin),
    ess_min   = ifelse(is.null(prior_phitree$ess_min), 5.0,                   prior_phitree$ess_min)
  )
}

# ------------------------------------------------------------------------------
# Verbose summary printer.  `lambda_user_supplied` drives the "(auto / user)"
# annotation on the lambda line; `auto_full` triggers a one-line note when the
# full Cholesky was synthesized rather than user-requested.
# ------------------------------------------------------------------------------
.print_heter_cov_summary <- function(useFullCov, var_params, phi_params,
                                     lambda_user_supplied,
                                     auto_full = FALSE) {
  cat("\nHeteroscedastic Sigma(Z) enabled (",
      if (useFullCov) "full Cholesky" else "diagonal", ")\n", sep = "")
  cat(sprintf("  vartree:   num_trees=%d  nu=%g  lambda=%g (%s)\n",
              var_params$num_trees, var_params$nu, var_params$lambda,
              ifelse(lambda_user_supplied, "user", "auto-calibrated")))
  if (useFullCov) {
    cat(sprintf("  phitree:   num_trees=%d  tau=%g  nmin=%d  ess_min=%g%s\n",
                phi_params$num_trees, phi_params$tau,
                phi_params$nmin, phi_params$ess_min,
                if (auto_full) "  (auto-enabled: nvar>1)" else ""))
  }
  invisible(NULL)
}
