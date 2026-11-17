# ==============================================================================
# HART Prior Configuration
# ==============================================================================
# Internal configuration of the HART preference models, shared by
# rhierMnlRwMixture, rhierLinearMixture, and rhierNegbinRw.
#
#   Prior$hart_model = "none"        bayesm paths: constant or linear Delta(Z)
#                    = "normalized"  theta_i = mu + S {f*(Z_i) + eps_i}
#                    = "direct"      theta_i = mu + f(Z_i) + S_i eps_i
#
# The tree lists bart, vartree, phitree are optional and require a HART value
# of hart_model.  Under a HART model mu ~ N(mubar, Amu^{-1}), independent of
# Sigma.  `tau` is an ensemble scale (leaf sd tau / sqrt(num_trees)).  Under
# the direct model the mean scales are tau_d^2 ~ IG(a_tau, (a_tau - 1) tau0_d^2)
# and the variance leaves IG(a_v, (a_v - 1) v0_j^{1/H_v}) with H_v trigamma(a_v)
# = tau^2, so E[v_j(Z)] = v0_j and sd log v_j(Z) = tau whatever H_v is.
# User lists are read with `[[` (exact names: tau vs tau0, a vs a_v).

.hart_models <- c("none", "normalized", "direct")
.hart_covs   <- c("heteroskedastic", "constant")

# Configuration ----------------------------------------------------------------

# Build the validated HART configuration.
#
# @param Prior The user's Prior list.
# @param nvar Number of unit-level coefficients.
# @param nz Number of characteristics (0 when Data$Z is absent).
# @param drawdelta TRUE iff Data$Z was supplied.
# @param coef_names Column names of X, recorded as the coefficient ordering.
# @param v0_default Prior mean of v_j(Z) per coefficient (direct model).
# @return A list with model, nvar, nz, coef_names, store_trees, full_cov, and
#   the tree lists bart, vartree, phitree (NULL when unused).
.hart_config <- function(Prior, nvar, nz, drawdelta, coef_names = NULL,
                         v0_default = 1) {
  model <- Prior[["hart_model"]] %||% "none"
  if (!is.character(model) || length(model) != 1L || !(model %in% .hart_models))
    pandterm(paste0("Prior$hart_model must be one of: ",
                    paste0("'", .hart_models, "'", collapse = ", "), "."))
  cov <- Prior[["hart_cov"]] %||% "heteroskedastic"
  if (!is.character(cov) || length(cov) != 1L || !(cov %in% .hart_covs))
    pandterm(paste0("Prior$hart_cov must be one of: ",
                    paste0("'", .hart_covs, "'", collapse = ", "), "."))
  .validate_hart(Prior, model, cov, nvar, drawdelta)

  cfg <- list(model       = model,
              cov         = if (model == "direct") cov else NULL,
              nvar        = as.integer(nvar),
              nz          = as.integer(nz),
              coef_names  = coef_names,
              store_trees = isTRUE(Prior[["store_trees"]] %||% TRUE),
              full_cov    = FALSE,
              bart        = NULL,
              vartree     = NULL,
              phitree     = NULL)
  if (model == "none") return(cfg)

  cfg$bart <- .parse_bart_params(Prior[["bart"]] %||% list(), nz, model, nvar)
  cfg$bart$store_trees <- cfg$store_trees
  if (model == "direct" && cov == "heteroskedastic") {
    # Off-diagonal ensembles exist iff there is more than one coefficient.
    cfg$full_cov <- nvar > 1L
    cfg$vartree  <- .parse_var_params(Prior[["vartree"]] %||% list(), nz,
                                      nvar, v0_default)
    cfg$phitree  <- if (cfg$full_cov)
      .parse_phi_params(Prior[["phitree"]] %||% list(), nz) else NULL
  }
  cfg
}#.HART_CONFIG

# Reject inconsistent selector / tree-list combinations.
.validate_hart <- function(Prior, model, cov, nvar, drawdelta) {
  bart <- Prior[["bart"]]; vartree <- Prior[["vartree"]]
  phitree <- Prior[["phitree"]]
  msg_cov <- paste("Prior$vartree and Prior$phitree require",
                   "Prior$hart_model = 'direct'.")
  if (!is.null(Prior[["hart_cov"]]) && model != "direct")
    pandterm("Prior$hart_cov requires Prior$hart_model = 'direct'.")
  if (model == "none") {
    if (!is.null(bart))
      pandterm(paste("Prior$bart requires Prior$hart_model = 'normalized'",
                     "or 'direct'."))
    if (!is.null(vartree) || !is.null(phitree)) pandterm(msg_cov)
    return(invisible(NULL))
  }
  if (!drawdelta)
    pandterm(paste0("Prior$hart_model = '", model, "' requires Data$Z."))
  if (model == "normalized") {
    if (!is.null(vartree) || !is.null(phitree)) pandterm(msg_cov)
    if (!is.null(bart[["a_tau"]]) || !is.null(bart[["tau0"]]))
      pandterm(paste("Prior$bart$a_tau and Prior$bart$tau0 require",
                     "Prior$hart_model = 'direct'; use Prior$bart$tau."))
  }
  if (model == "direct") {
    if (cov == "constant" && (!is.null(vartree) || !is.null(phitree)))
      pandterm(paste("Prior$vartree and Prior$phitree are not used under",
                     "Prior$hart_cov = 'constant'."))
    if (!is.null(phitree) && nvar == 1L)
      pandterm("Prior$phitree requires nvar > 1 (no off-diagonal terms).")
    if (!is.null(bart[["tau"]]))
      pandterm(paste("Prior$bart$tau is not used under Prior$hart_model =",
                     "'direct'; set Prior$bart$a_tau and Prior$bart$tau0."))
    if (!is.null(vartree[["nu"]]) || !is.null(vartree[["lambda"]]) ||
        !is.null(vartree[["a_v"]]))
      pandterm(paste("Prior$vartree takes tau (prior sd of log d_j(Z)) and",
                     "v0 (prior mean of d_j(Z)); nu, lambda, a_v are not",
                     "used."))
    if (!is.null(phitree[["ess_min"]]))
      pandterm(paste("Prior$phitree$ess_min is not used; leaf admissibility",
                     "is Prior$phitree$nmin."))
    if (cov == "heteroskedastic" &&
        (!is.null(Prior[["nu"]]) || !is.null(Prior[["V"]])))
      message(paste("Prior$nu and Prior$V are not used under hart_model =",
                    "'direct' with hart_cov = 'heteroskedastic'."))
  }
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

# Tree-structure and DART fields shared by every tree list.
.parse_tree_params <- function(x, nz, num_trees, nmin) {
  list(
    num_trees = as.integer(x[["num_trees"]] %||% num_trees),
    power     = x[["power"]]  %||% 2.0,
    base      = x[["base"]]   %||% 0.95,
    numcut    = as.integer(x[["numcut"]] %||% 100L),
    nmin      = as.integer(x[["nmin"]]   %||% nmin),
    sparse    = isTRUE(x[["sparse"]] %||% FALSE),
    theta     = x[["theta"]]  %||% 0.0,
    omega     = x[["omega"]]  %||% 1.0,
    a         = x[["a"]]      %||% 0.5,
    b         = x[["b"]]      %||% 1.0,
    rho       = x[["rho"]]    %||% nz,
    aug       = isTRUE(x[["aug"]] %||% FALSE),
    burn      = as.integer(x[["burn"]] %||% 100L)
  )
}#.PARSE_TREE_PARAMS

# Mean-tree hyperparameters: `tau` under the normalized model, the scale
# hyperprior (a_tau, tau0) under the direct model.
.parse_bart_params <- function(bart_list, nz, model, nvar) {
  out <- .parse_tree_params(bart_list, nz, num_trees = 200L, nmin = 5L)
  if (model == "direct") {
    out$a_tau <- .check_shape(bart_list[["a_tau"]] %||% 3, "Prior$bart$a_tau")
    out$tau0  <- .check_scale(bart_list[["tau0"]] %||% 1, nvar,
                              "Prior$bart$tau0")
  } else {
    out$tau <- bart_list[["tau"]] %||% 1.0
  }
  out
}#.PARSE_BART_PARAMS

# Variance-tree hyperparameters: v_j(Z) = prod_h s_jh, s ~ IG(a_v, b_vj), with
# `tau` the prior sd of log v_j(Z) (H_v trigamma(a_v) = tau^2) and `v0` its
# prior mean.
.parse_var_params <- function(prior_vartree, nz, nvar, v0_default = 1) {
  out <- .parse_tree_params(prior_vartree, nz, num_trees = 40L, nmin = 5L)
  out$tau <- .check_scale(prior_vartree[["tau"]] %||% 1, 1L,
                          "Prior$vartree$tau")
  out$a_v <- .ig_shape(out$tau, out$num_trees)
  out$v0  <- .check_scale(prior_vartree[["v0"]] %||% v0_default, nvar,
                          "Prior$vartree$v0")
  out
}#.PARSE_VAR_PARAMS

# Off-diagonal tree hyperparameters: leaves N(0, tau^2 / H_phi).
.parse_phi_params <- function(prior_phitree, nz) {
  out <- .parse_tree_params(prior_phitree, nz, num_trees = 40L, nmin = 5L)
  out$tau <- prior_phitree[["tau"]] %||% 1.0
  out
}#.PARSE_PHI_PARAMS

# Inverse-gamma leaf shape giving a product of H leaves a prior sd of `tau` on
# the log scale: H trigamma(a) = tau^2.
.ig_shape <- function(tau, H) {
  stats::uniroot(function(a) trigamma(a) - tau^2 / H, c(1 + 1e-6, 1e8),
                 tol = 1e-10)$root
}#.IG_SHAPE

# Inverse-gamma shape; > 1 so the prior mean exists.
.check_shape <- function(a, name) {
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a <= 1)
    pandterm(paste0(name, " must be a scalar greater than 1."))
  as.numeric(a)
}#.CHECK_SHAPE

# Positive scale of length 1 (recycled) or nvar.
.check_scale <- function(x, nvar, name) {
  if (!is.numeric(x) || !(length(x) %in% c(1L, nvar)) ||
      any(!is.finite(x)) || any(x <= 0))
    pandterm(paste0(name, " must be positive, of length 1 or nvar."))
  rep_len(as.numeric(x), nvar)
}#.CHECK_SCALE

# Printing ---------------------------------------------------------------------

# Print the mean-tree prior (verbose wrapper output).
.print_hart_prior <- function(cfg) {
  b <- cfg$bart
  cat(sprintf("\nHART prior (hart_model = '%s'%s):\n", cfg$model,
              if (cfg$model == "direct")
                sprintf(", hart_cov = '%s'", cfg$cov) else ""))
  cat(sprintf("num_trees = %d\n", b$num_trees))
  cat(sprintf("power = %g\n", b$power))
  cat(sprintf("base = %g\n", b$base))
  if (cfg$model == "direct") {
    cat(sprintf("a_tau = %g\n", b$a_tau))
    cat("tau0 =", format(b$tau0), "\n")
  } else {
    cat(sprintf("tau = %g\n", b$tau))
  }
  cat(sprintf("numcut = %d\n", b$numcut))
  cat(sprintf("nmin = %d\n", b$nmin))
  cat(sprintf("sparse = %d\n", b$sparse))
  if (b$sparse) {
    cat(sprintf("theta = %g\n", b$theta))
    cat(sprintf("omega = %g\n", b$omega))
    cat(sprintf("a = %g\n", b$a))
    cat(sprintf("b = %g\n", b$b))
    cat(sprintf("rho = %g\n", b$rho))
    cat(sprintf("aug = %d\n", b$aug))
    cat(sprintf("burn = %d\n", b$burn))
  }
  invisible(NULL)
}#.PRINT_HART_PRIOR

# Print the covariance-tree priors (direct model, verbose wrapper output).
.print_hart_covariance <- function(cfg) {
  v <- cfg$vartree
  cat("\nCovariance trees (",
      if (cfg$full_cov) "full Cholesky" else "diagonal", "):\n", sep = "")
  cat(sprintf("  vartree: num_trees = %d, tau = %g (a_v = %.3g), v0 = %s\n",
              v$num_trees, v$tau, v$a_v, paste(format(v$v0), collapse = " ")))
  if (cfg$full_cov) {
    p <- cfg$phitree
    cat(sprintf("  phitree: num_trees = %d, tau = %g, nmin = %d\n",
                p$num_trees, p$tau, p$nmin))
  }
  invisible(NULL)
}#.PRINT_HART_COVARIANCE
