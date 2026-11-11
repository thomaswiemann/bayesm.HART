# ==============================================================================
# HART Prior Configuration
# ==============================================================================
# One validated configuration for the HART preference models, shared by
# rhierMnlRwMixture, rhierLinearMixture, and rhierNegbinRw.  All functions are
# package-internal.
#
#   Prior$hart_model = "none"        bayesm paths: constant or linear Delta(Z)
#                    = "normalized"  theta_i = mu + S {f*(Z_i) + eps_i}
#                    = "direct"      theta_i = mu + f(Z_i) + S_i eps_i
#
# The selector is explicit.  Supplying a tree list (`bart`, `vartree`,
# `phitree`) without a HART value of `hart_model` is an error, and every tree
# list is optional under a HART model because all hyperparameters have
# defaults.
#
# Status (discussions/2026-11-10-hart-package-update-plan.md): "normalized"
# runs on NormalizedHartBackend (stage 1); "direct" still runs on the pre-plan
# HeterCovBackend until stage 4.  Variance- and phi-tree parsing still needs
# the warm-start coefficients (data-calibrated lambda), so the wrappers call
# .parse_var_params() / .parse_phi_params() after the warm start and store the
# result back into cfg$vartree / cfg$phitree.
#
# Under a HART model the prior on mu is the independent N(mubar, Amu^{-1});
# .hart_Amu() builds the nvar x nvar precision from Prior$Amu.

.hart_models <- c("none", "normalized", "direct")

# Configuration ----------------------------------------------------------------

# Build the validated HART configuration.
#
# @param Prior The user's Prior list.
# @param nvar Number of unit-level coefficients.
# @param nz Number of characteristics (0 when Data$Z is absent).
# @param drawdelta TRUE iff Data$Z was supplied.
# @param coef_names Column names of X, recorded as the coefficient ordering.
# @return A list with model, nvar, nz, coef_names, store_trees, full_cov, and
#   the tree lists bart, vartree, phitree (NULL when unused).
.hart_config <- function(Prior, nvar, nz, drawdelta, coef_names = NULL) {
  model <- Prior$hart_model %||% "none"
  if (!is.character(model) || length(model) != 1L || !(model %in% .hart_models))
    pandterm(paste0("Prior$hart_model must be one of: ",
                    paste0("'", .hart_models, "'", collapse = ", "), "."))
  .validate_hart(Prior, model, nvar, drawdelta)

  cfg <- list(model       = model,
              nvar        = as.integer(nvar),
              nz          = as.integer(nz),
              coef_names  = coef_names,
              store_trees = isTRUE(Prior$store_trees %||% TRUE),
              full_cov    = FALSE,
              bart        = NULL,
              vartree     = NULL,
              phitree     = NULL)
  if (model == "none") return(cfg)

  cfg$bart <- .parse_bart_params(Prior$bart %||% list(), nz)
  cfg$bart$store_trees <- cfg$store_trees
  if (model == "direct") {
    # Off-diagonal ensembles exist iff there is more than one coefficient.
    cfg$full_cov <- nvar > 1L
    cfg$vartree  <- Prior$vartree %||% list()   # parsed after the warm start
    cfg$phitree  <- if (cfg$full_cov) Prior$phitree %||% list() else NULL
  }
  cfg
}#.HART_CONFIG

# Reject inconsistent selector / tree-list combinations.  Most diagnostic
# check first: a tree list without a HART selector, then the data requirement,
# then model-specific lists.
.validate_hart <- function(Prior, model, nvar, drawdelta) {
  has <- c(bart    = !is.null(Prior$bart),
           vartree = !is.null(Prior$vartree),
           phitree = !is.null(Prior$phitree))
  if (model == "none") {
    if (has["bart"])
      pandterm("Prior$bart requires Prior$hart_model = 'normalized' or 'direct'.")
    if (has["vartree"] || has["phitree"])
      pandterm("Prior$vartree and Prior$phitree require Prior$hart_model = 'direct'.")
    return(invisible(NULL))
  }
  if (!drawdelta)
    pandterm(paste0("Prior$hart_model = '", model, "' requires Data$Z."))
  if (model == "normalized" && (has["vartree"] || has["phitree"]))
    pandterm("Prior$vartree and Prior$phitree require Prior$hart_model = 'direct'.")
  if (model == "direct" && has["phitree"] && nvar == 1L)
    pandterm("Prior$phitree requires nvar > 1 (no off-diagonal terms when nvar == 1).")
  invisible(NULL)
}#.VALIDATE_HART

# Prior precision of mu under a HART model: independent mu ~ N(mubar, Amu^{-1}).
# Accepts NULL (default * I), a positive scalar (scalar * I), or an SPD
# nvar x nvar matrix.
.hart_Amu <- function(Amu, nvar, default) {
  msg <- "Prior$Amu must be a positive scalar or an nvar x nvar SPD matrix."
  if (is.null(Amu)) return(default * diag(nvar))
  if (length(Amu) == 1L) {
    if (!is.numeric(Amu) || !is.finite(Amu) || Amu <= 0) pandterm(msg)
    return(as.numeric(Amu) * diag(nvar))
  }
  Amu <- unname(as.matrix(Amu))
  if (!is.numeric(Amu) || !all(dim(Amu) == nvar)) pandterm(msg)
  if (!isSymmetric(Amu, tol = 1e-8) ||
      any(eigen(Amu, symmetric = TRUE, only.values = TRUE)$values <= 0))
    pandterm("Prior$Amu must be symmetric positive definite.")
  Amu
}#.HART_AMU

# Hyperparameters --------------------------------------------------------------

# Mean-tree hyperparameters with defaults.  `nz` sets the default DART `rho`.
.parse_bart_params <- function(bart_list, nz) {
  num_trees <- bart_list$num_trees %||% 200
  list(
    num_trees = num_trees,
    power     = bart_list$power  %||% 2.0,
    base      = bart_list$base   %||% 0.95,
    tau       = bart_list$tau    %||% 1.0,   # ensemble scale; leaf sd tau/sqrt(H)
    numcut    = bart_list$numcut %||% 100,
    nmin      = bart_list$nmin   %||% 5L,
    sparse    = bart_list$sparse %||% FALSE,
    theta     = bart_list$theta  %||% 0.0,
    omega     = bart_list$omega  %||% 1.0,
    a         = bart_list$a      %||% 0.5,
    b         = bart_list$b      %||% 1.0,
    rho       = bart_list$rho    %||% nz,
    aug       = bart_list$aug    %||% FALSE,
    burn      = bart_list$burn   %||% 100
  )
}#.PARSE_BART_PARAMS

# Variance-tree hyperparameters: chi^{-2} leaf prior plus DART.  `oldbetas`
# (per-unit warm-start coefficients) calibrates `lambda` when the user did not
# supply one; stage 4 replaces this with the appendix's fixed (a_v, v0).
.parse_var_params <- function(prior_vartree, oldbetas, nz) {
  lambda_default <- mean(apply(oldbetas, 2, var))
  if (!is.finite(lambda_default) || lambda_default <= 0) lambda_default <- 1.0
  list(
    num_trees = as.integer(prior_vartree$num_trees %||% 40L),
    power     = prior_vartree$power  %||% 2.0,
    base      = prior_vartree$base   %||% 0.95,
    nu        = prior_vartree$nu     %||% 10.0,
    lambda    = prior_vartree$lambda %||% lambda_default,
    numcut    = prior_vartree$numcut %||% 100L,
    nmin      = prior_vartree$nmin   %||% 5L,
    sparse    = prior_vartree$sparse %||% FALSE,
    theta     = prior_vartree$theta  %||% 0.0,
    omega     = prior_vartree$omega  %||% 1.0,
    a         = prior_vartree$a      %||% 0.5,
    b         = prior_vartree$b      %||% 1.0,
    rho       = prior_vartree$rho    %||% nz,
    aug       = prior_vartree$aug    %||% FALSE,
    burn      = prior_vartree$burn   %||% 100L
  )
}#.PARSE_VAR_PARAMS

# Off-diagonal (phi) tree hyperparameters: N(0, tau^2) leaves plus DART.  The
# default leaf sd matches the implied inverse-Wishart(nu) scale of the
# off-diagonal Cholesky entries; stage 4 replaces this with the appendix's
# ensemble-scale convention and drops `ess_min`.
.parse_phi_params <- function(prior_phitree, nz, nu, nvar) {
  num_trees <- as.integer(prior_phitree$num_trees %||% 40L)
  if (nvar > 1) {
    i_vec <- 2:nvar
    V_phi_bar <- (2 / (nvar * (nvar - 1))) *
      sum((i_vec - 1) / (nu - nvar + i_vec - 2))
    tau_default <- sqrt(V_phi_bar / num_trees)
  } else {
    tau_default <- 1.0 / sqrt(num_trees)
  }
  list(
    num_trees = num_trees,
    power     = prior_phitree$power   %||% 2.0,
    base      = prior_phitree$base    %||% 0.95,
    tau       = prior_phitree$tau     %||% tau_default,
    numcut    = prior_phitree$numcut  %||% 100L,
    sparse    = prior_phitree$sparse  %||% FALSE,
    theta     = prior_phitree$theta   %||% 0.0,
    omega     = prior_phitree$omega   %||% 1.0,
    a         = prior_phitree$a       %||% 0.5,
    b         = prior_phitree$b       %||% 1.0,
    rho       = prior_phitree$rho     %||% nz,
    aug       = prior_phitree$aug     %||% FALSE,
    nmin      = prior_phitree$nmin    %||% 2L,
    ess_min   = prior_phitree$ess_min %||% 5.0
  )
}#.PARSE_PHI_PARAMS

# Printing ---------------------------------------------------------------------

# Verbose summary of the mean-tree prior, printed with the wrappers' setup
# information when a HART model is selected.
.print_hart_prior <- function(cfg) {
  b <- cfg$bart
  cat(sprintf("\nHART prior (hart_model = '%s'):\n", cfg$model))
  cat(sprintf("num_trees = %d\n", as.integer(b$num_trees)))
  cat(sprintf("power = %g\n", b$power))
  cat(sprintf("base = %g\n", b$base))
  cat(sprintf("tau = %g\n", b$tau))
  cat(sprintf("numcut = %d\n", as.integer(b$numcut)))
  cat(sprintf("nmin = %d\n", as.integer(b$nmin)))
  cat(sprintf("sparse = %d\n", b$sparse))
  if (isTRUE(b$sparse)) {
    cat(sprintf("theta = %g\n", b$theta))
    cat(sprintf("omega = %g\n", b$omega))
    cat(sprintf("a = %g\n", b$a))
    cat(sprintf("b = %g\n", b$b))
    cat(sprintf("rho = %g\n", b$rho))
    cat(sprintf("aug = %g\n", b$aug))
    cat(sprintf("burn = %d\n", as.integer(b$burn)))
  }
  invisible(NULL)
}#.PRINT_HART_PRIOR

# Verbose summary of the covariance-tree priors (direct model), printed once
# the variance- and phi-tree lists are parsed.
.print_hart_covariance <- function(cfg, lambda_user_supplied) {
  cat("\nHeteroscedastic Sigma(Z) enabled (",
      if (cfg$full_cov) "full Cholesky" else "diagonal", ")\n", sep = "")
  v <- cfg$vartree
  cat(sprintf("  vartree:   num_trees=%d  nu=%g  lambda=%g (%s)\n",
              v$num_trees, v$nu, v$lambda,
              if (lambda_user_supplied) "user" else "auto-calibrated"))
  if (cfg$full_cov) {
    p <- cfg$phitree
    cat(sprintf("  phitree:   num_trees=%d  tau=%g  nmin=%d  ess_min=%g\n",
                p$num_trees, p$tau, p$nmin, p$ess_min))
  }
  invisible(NULL)
}#.PRINT_HART_COVARIANCE
