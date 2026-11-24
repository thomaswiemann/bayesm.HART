# ==============================================================================
# Restart metadata
# ==============================================================================

# Raw observations and unit IDs in fitting order.
.hart_restart_data <- function(Data, units) {
  unit_data <- Data[[units]]
  ids <- Data[["unit_id"]]
  if (!is.null(ids)) {
    if (length(ids) != length(unit_data) || anyNA(ids))
      pandterm("Data$unit_id must have one non-missing ID per unit")
    ids <- as.character(ids)
    if (anyDuplicated(ids))
      pandterm("Data$unit_id must contain unique non-missing IDs")
  } else {
    ids <- names(unit_data)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      ids <- as.character(seq_along(unit_data))
  }
  out <- list(unit_id = ids)
  out[[units]] <- lapply(unit_data, function(u) list(y = u$y, X = u$X))
  out["Z"] <- list(Data[["Z"]])
  if (units == "lgtdata") out$p <- Data$p
  out
}#.HART_RESTART_DATA

# Resolved fitting recipe without terminal sampler state.
.hart_restart <- function(likelihood, data, cfg, ncomp, drawdelta, prior,
                          mcmc, preparation = list()) {
  config <- cfg[c("model", "cov", "nvar", "nz", "coef_names")]
  config$ncomp <- ncomp
  config$drawdelta <- drawdelta
  config$prior <- prior
  for (name in c("bart", "vartree", "phitree")) {
    params <- cfg[[name]]
    if (!is.null(params)) {
      params$z_index <- NULL
      params$n_unique <- NULL
      config[[name]] <- params
    }
  }
  list(schema = 1L, likelihood = likelihood, data = data, config = config,
       mcmc = mcmc, preparation = preparation)
}#.HART_RESTART

# A compatible terminal checkpoint for the called likelihood.
.hart_previous <- function(previous, likelihood) {
  if (!is.list(previous) || is.null(previous$restart))
    pandterm("previous is missing restart; refit with restart support")
  saved <- previous$restart
  if (!is.list(saved) || !identical(saved$schema, 1L))
    pandterm("previous requires restart schema 1")
  if (!identical(saved$likelihood, likelihood))
    pandterm("previous has the wrong likelihood")
  for (name in c("state", "data", "config", "mcmc", "preparation")) {
    if (is.null(saved[[name]]))
      pandterm(paste("previous is missing restart", name))
  }
  saved
}#.HART_PREVIOUS

# Omitted run controls inherit; R counts additional sweeps.
.hart_continue_controls <- function(Mcmc, saved) {
  if (is.null(Mcmc$R)) pandterm("Requires R in Mcmc for continuation")
  if (any(!names(Mcmc) %in% c("R", "keep", "nprint")))
    pandterm("With previous, Mcmc accepts only R, keep and nprint")
  out <- saved$mcmc
  for (name in names(Mcmc)) {
    if (!is.null(Mcmc[[name]])) out[[name]] <- Mcmc[[name]]
  }
  for (name in c("R", "keep", "nprint")) {
    x <- out[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != floor(x) || x < if (name == "nprint") 0 else 1)
      pandterm(paste("Mcmc", name, "must be a valid integer"))
  }
  if (out$R < out$keep) pandterm("Mcmc$R must retain at least one draw")
  if (saved$likelihood == "negbin") {
    out$fixalpha <- saved$config$fixalpha
    if (out$fixalpha) out$alpha <- saved$config$alpha_fixed
  }
  out
}#.HART_CONTINUE_CONTROLS

# Resolved tree settings and output flags without public prior parsing.
.hart_previous_config <- function(saved) {
  cfg <- saved$config
  cfg$Amu <- cfg$prior$Amu
  cfg$store_trees <- cfg$bart$store_trees %||% TRUE
  cfg$full_cov <- !is.null(cfg$phitree)
  cfg
}#.HART_PREVIOUS_CONFIG

# Batch observations append in saved unit order, with fixed Z coordinates.
.hart_merge_data <- function(saved, batch) {
  units <- if (saved$likelihood == "mnl") "lgtdata" else "regdata"
  old <- saved$data
  entries <- batch[[units]]
  ids <- batch$unit_id
  if (is.null(ids) || !is.list(entries) || !length(entries) ||
      length(ids) != length(entries) || anyNA(ids))
    pandterm("Data$unit_id must have one non-missing ID per batch entry")
  ids <- as.character(ids)
  if (anyDuplicated(ids)) pandterm("Data$unit_id must be unique")
  if (units == "lgtdata" && !is.null(batch$p) &&
      !isTRUE(all.equal(batch$p, old$p, tolerance = 0)))
    pandterm("Data$p must equal the saved p")
  positions <- match(ids, old$unit_id)
  fresh <- which(is.na(positions))
  for (i in seq_along(entries)) {
    u <- entries[[i]]
    X <- u$X
    reference <- old[[units]][[if (is.na(positions[i])) 1L else positions[i]]]$X
    if (!is.matrix(X) || ncol(X) != ncol(reference))
      pandterm("Data batch X must have the saved number of columns")
    if (!is.null(colnames(X)) && !is.null(colnames(reference)) &&
        !identical(colnames(X), colnames(reference)))
      pandterm("Data batch X column names differ from the saved names")
    rows <- length(u$y) * if (units == "lgtdata") old$p else 1L
    if (!length(u$y) || nrow(X) != rows)
      pandterm("Data batch y and X have incompatible dimensions")
  }
  if (is.null(old$Z)) {
    if (!is.null(batch$Z)) pandterm("Data Z supplied to a fit without Z")
  } else {
    Z <- batch$Z
    if (!is.matrix(Z) || nrow(Z) != length(ids) ||
        ncol(Z) != ncol(old$Z))
      pandterm("Data batch Z must have one row per entry and saved columns")
    if (!is.null(colnames(Z)) && !is.null(colnames(old$Z)) &&
        !identical(colnames(Z), colnames(old$Z)))
      pandterm("Data batch Z column names differ from the saved names")
    for (i in which(!is.na(positions))) {
      if (!isTRUE(all.equal(Z[i, ], old$Z[positions[i], ],
                           tolerance = 0, check.attributes = FALSE)))
        pandterm("Data batch Z differs from the stored row for an old unit")
    }
    if (length(fresh)) old$Z <- rbind(old$Z, Z[fresh, , drop = FALSE])
  }
  for (i in seq_along(entries)) {
    u <- list(y = entries[[i]]$y, X = entries[[i]]$X)
    j <- positions[i]
    if (is.na(j)) {
      old[[units]][[length(old[[units]]) + 1L]] <- u
    } else {
      unit <- old[[units]][[j]]
      u$y <- if (is.matrix(unit$y))
        rbind(unit$y, matrix(u$y, ncol = 1L)) else c(unit$y, u$y)
      u$X <- rbind(unit$X, u$X)
      old[[units]][[j]] <- u
    }
  }
  old$unit_id <- c(old$unit_id, ids[fresh])
  old
}#.HART_MERGE_DATA

# Old live values overwrite initializers; population state persists.
.hart_extend_state <- function(saved, betas, tau = NULL) {
  state <- saved$state
  old <- seq_along(saved$data$unit_id)
  betas[old, ] <- state$betas
  state$betas <- betas
  if (!is.null(tau)) {
    tau[old] <- state$tau
    state$tau <- tau
  }
  if (saved$config$model == "none")
    state$backend$ind <- c(state$backend$ind,
                          rep(1L, nrow(betas) - length(old)))
  state
}#.HART_EXTEND_STATE
