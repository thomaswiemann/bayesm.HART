#' Simulate Hierarchical Multinomial Logit Data
#'
#' Generates simulated data suitable for testing hierarchical multinomial logit
#' models, particularly those involving individual-specific covariates (Z)
#' influencing coefficients (beta_i). Supports various functional forms for the
#' Z-beta relationship and mixture models for residual heterogeneity.
#'
#' @param nlgt Integer. Number of individuals or cross-sectional units.
#' @param nT Integer. Number of choice observations per individual.
#' @param p Integer. Number of choice alternatives (including outside option if any).
#' @param nz Integer. Number of demographic/individual-specific variables in Z.
#'   If `nz = 0`, no Z matrix is generated, `betabar_true` is set to zero, 
#'   and `beta_func_type`/`beta_func_args` are ignored.
#' @param nXa Integer. Number of alternative-specific variables in X.
#' @param nXd Integer. Number of choice-invariant variables in X (e.g., price).
#' @param const Logical. Include p-1 intercepts in the model?
#' @param z_dist_func Function. A function to generate the Z matrix.
#'   Must accept arguments `n` (nlgt) and `d` (nz) and return an n x d matrix.
#'   Default: `function(n, d) matrix(runif(n*d, -1, 1), n, d)`.
#' @param standardize_Z Logical. Standardize the generated Z matrix (mean 0, sd 1)?
#' @param beta_func_type Character. Specifies the functional form mapping Z to
#'   the systematic component of beta (`betabar_i`). Ignored if `nz = 0`.
#'   Options:
#'   * `"linear"`: Linear function `betabar_i = Z_i %*% Delta`. Requires `Delta` in `beta_func_args`.
#'   * `"step"`: Step function based on one Z variable. Requires `cutoff`, `beta_1`, `beta_2`, `Z_index` in `beta_func_args`.
#'   * `"friedman"`: Friedman benchmark function (modified) based on first 5 Z variables. Requires `coef_index` in `beta_func_args` to specify which coefficient it applies to (others are zero).
#'   * `"custom"`: A user-defined function provided in `beta_func_args$func`.
#' @param beta_func_args List. Arguments needed for the chosen `beta_func_type`.
#'   Ignored if `nz = 0`.
#'   * For `"linear"`: `list(Delta = matrix(runif(ncoef * nz), nrow=nz))`. Delta is `nz x ncoef`.
#'   * For `"step"`: `list(cutoff = 0, beta_1 = rep(-1, ncoef), beta_2 = rep(1, ncoef), Z_index = 1)`. `beta_1`/`beta_2` are vectors of length `ncoef`. `Z_index` is the column of Z to use.
#'   * For `"friedman"`: `list(coef_index = 1)`. `nz` must be >= 5. The function is applied to `betabar_i[coef_index]`, others are 0.
#'   * For `"custom"`: `list(func = function(Zi) { ... })`. The function must take a vector `Zi` (a row of Z) and return a vector `betabar_i` of length `ncoef`.
#' @param ncomp Integer. Number of components in the normal mixture for residual heterogeneity (`eps_i`).
#' @param mixture_comps List. Optional pre-specified mixture components. A list of
#'   length `ncomp`, where each element is `list(mu = ..., rooti = ...)`. `mu` is
#'   the mean vector (length `ncoef`), `rooti` is the upper Cholesky factor of the
#'   inverse covariance matrix (`ncoef x ncoef`). If `NULL`, components are generated
#'   based on `sigma_inv_diag`.
#' @param sigma_inv_diag Numeric. Diagonal value for the inverse covariance matrix
#'   (precision) of mixture components if `mixture_comps` is `NULL`. Assumes identity
#'   covariance scaled by this.
#' @param Xa_dist_func Function. Function to generate alternative-specific variables `Xa`.
#'   Takes `n` (nT), `p`, `na` (nXa) and returns a matrix (usually `n x (p*na)` or similar structure expected by `createX`).
#'   Default: `function(n, p, na) matrix(runif(n*p*na, -1, 1), ncol=p*na)`.
#' @param Xd_dist_func Function. Function to generate choice-invariant variables `Xd`.
#'   Takes `n` (nT), `nd` (nXd) and returns an `n x nd` matrix.
#'   Default: `function(n, nd) matrix(rnorm(n*nd), ncol=nd)`.
#' @param seed Integer. Optional random seed for reproducibility.
#'
#' @return A list suitable for direct use as the `Data` argument in
#'   `rhierMnlRwMixture`, containing:
#'   * `p`: Number of alternatives.
#'   * `lgtdata`: List of length `nlgt`. Each element `i` is `list(y=y_i, X=X_i, beta=beta_i, betabar=betabar_i)`.
#'   * `Z`: The `nlgt x nz` matrix of individual-specific covariates (standardized if requested).
#'   Additionally, the list contains `true_values`:
#'   * `true_values$beta_true`: `nlgt x ncoef` matrix of true `beta_i`.
#'   * `true_values$betabar_true`: `nlgt x ncoef` matrix of true `betabar_i = f(Z_i)`.
#'   * `true_values$true_params`: List containing parameters used for generation (`beta_func_type`, `beta_func_args`, `mixture_comps`, `pvec`).
#'   * `true_values$dimensions`: List containing key dimensions used (`p`, `nlgt`, `nT`, `nz`, `ncoef`, etc.).
#'
#' @export
#' @importFrom stats rnorm runif rmultinom optim quantile
#' @importFrom utils flush.console
#' @importFrom bayesm rdirichlet rmixture
#' @examples
#' # Simple linear example
#' sim_data_linear <- sim_hier_mnl(nlgt = 50, nT = 5, p = 3, nz = 2, nXa = 1, nXd = 0,
#'                                beta_func_type = "linear", seed = 123)
#' str(sim_data_linear)
#'
#' # Step function example
#' sim_data_step <- sim_hier_mnl(nlgt = 50, nT = 5, p = 3, nz = 2, nXa = 1, nXd = 0,
#'                               beta_func_type = "step",
#'                               beta_func_args = list(Z_index = 1),
#'                               seed = 456)
#' str(sim_data_step$true_values$true_params) # Check args used
#' plot(sim_data_step$Z[,1], sim_data_step$true_values$betabar_true[,1]) # Visualize step
#'
sim_hier_mnl <- function(
  # Dimensions
  nlgt = 300,       # Number of individuals/groups
  nT = 10,          # Number of choices per individual
  p = 3,            # Number of choice alternatives (incl. outside)
  nz = 5,           # Number of demographic variables (in Z)
  nXa = 2,          # Number of alternative-specific X vars
  nXd = 1,          # Number of choice-invariant X vars
  const = TRUE,     # Include intercepts?

  # Z Generation
  z_dist_func = function(n, d) matrix(stats::runif(n*d, -1, 1), n, d), # Function to generate Z
  standardize_Z = TRUE, # Standardize Z matrix?

  # Beta Generation (betabar_i = f(Z_i))
  beta_func_type = "linear", # Type: "linear", "step", "friedman", "custom"
  beta_func_args = list(),   # List of arguments for the chosen beta function

  # Beta Generation (eps_i ~ Mixture)
  ncomp = 1,             # Number of mixture components for eps_i
  mixture_comps = NULL,  # Optional: Pre-specify list of components 
  sigma_inv_diag = 1,    # Diagonal value for Sigma_inv if mixture_comps is NULL

  # X Generation
  Xa_dist_func = function(n, p, na) matrix(stats::runif(n*p*na, -1, 1), ncol=p*na), # Func for Xa
  Xd_dist_func = function(n, nd) matrix(stats::rnorm(n*nd), ncol=nd), # Func for Xd

  # Other
  seed = NULL     # Random seed
) {
  # ============================================================================
  # 1. Input Checks & Setup
  # ============================================================================
  if (!is.null(seed)) set.seed(seed)

  # Calculate number of coefficients
  ncoef <- const * (p - 1) + (p - 1) * nXd + nXa
  if (ncoef == 0) stop("Model must have at least one coefficient.")

  # Validate beta_func_type
  valid_beta_funcs <- c("linear", "step", "friedman", "custom")
  if (nz > 0 && !(beta_func_type %in% valid_beta_funcs)) {
    stop(paste("Invalid beta_func_type. Choose from:",
               paste(valid_beta_funcs, collapse=", ")))
  }

  if (nz > 0 && beta_func_type == "friedman" && nz < 5) {
      stop("Friedman function requires nz >= 5.")
  }
  if (nz > 0 && beta_func_type == "custom" && !is.function(beta_func_args$func)) {
      stop("Custom beta function requires 'func' in beta_func_args.")
  }

  # ============================================================================
  # 2. Generate Z
  # ============================================================================
  Z <- NULL # Initialize Z as NULL
  if (nz > 0) {
      Z <- z_dist_func(nlgt, nz)
      if (!is.matrix(Z) || nrow(Z) != nlgt || ncol(Z) != nz) {
          stop("z_dist_func did not return a matrix with correct dimensions.")
      }

      # Standardize Z if requested
      if (standardize_Z) {
          # Use base::scale which should dispatch to stats::scale.default
          Z <- base::scale(Z)
      }
      # Handle potential NaN if sd is 0 for a column
      Z[is.nan(Z)] <- 0
  }

  # ============================================================================
  # 3. Generate betabar_true = f(Z_i)
  # ============================================================================
  betabar_true <- matrix(0, nrow = nlgt, ncol = ncoef)
  final_beta_func_args <- list() # Default empty list

  if (nz > 0) {
      # Only calculate betabar from Z if nz > 0
      args <- beta_func_args # Local copy for modification

      # --- Set default arguments based on type ---
      if (beta_func_type == "linear") {
        if (is.null(args$Delta)) {
          # Default Delta: nz x ncoef matrix of ~U(0,1)
          args$Delta <- matrix(stats::runif(nz * ncoef), nrow = nz, ncol = ncoef)
        } else if (!is.matrix(args$Delta) || nrow(args$Delta) != nz ||
                   ncol(args$Delta) != ncoef) {
            stop("Provided Delta matrix has incorrect dimensions.")
        }
        # Vectorized calculation
        betabar_true <- Z %*% args$Delta

      } else { # Non-vectorized cases: Loop through individuals
        if (beta_func_type == "step") {
            if (is.null(args$cutoff)) args$cutoff <- 0
            if (is.null(args$beta_1)) args$beta_1 <- rep(-1, ncoef)
            if (is.null(args$beta_2)) args$beta_2 <- rep(1, ncoef)
            if (is.null(args$Z_index)) args$Z_index <- 1
            if (args$Z_index > nz) stop("Z_index exceeds number of Z variables.")
            if (length(args$beta_1) != ncoef || length(args$beta_2) != ncoef) {
                stop("beta_1 and beta_2 must have length ncoef.")
            }
        } else if (beta_func_type == "friedman") {
            if (is.null(args$coef_index)) args$coef_index <- 1
            if (args$coef_index > ncoef) stop("coef_index exceeds ncoef.")
            # Set defaults for other Friedman parameters if not provided
            if (is.null(args$coefs)) args$coefs <- c(10, 20, 10, 5)
            if (is.null(args$shift)) args$shift <- 0
            if (is.null(args$scale)) args$scale <- 8
        }

        # Apply function row-by-row
        for (i in 1:nlgt) {
            Zi <- Z[i, ]
            if (beta_func_type == "step") {
                betabar_true[i, ] <- .beta_Z_step(Zi, args$Z_index, args$cutoff,
                                                  args$beta_1, args$beta_2)
            } else if (beta_func_type == "friedman") {
                # Applies Friedman value to coef_index, others are 0
                # Now uses args which are guaranteed to have defaults
                friedman_val <- .beta_Z_friedman(Zi, args$coefs,
                                                 args$shift,
                                                 args$scale)
                betabar_true[i, args$coef_index] <- friedman_val
            } else if (beta_func_type == "custom") {
                res <- args$func(Zi)
                if (length(res) != ncoef) {
                    stop("Custom function did not return vector of length ncoef.")
                }
                betabar_true[i, ] <- res
            }
        }
      }

      # Store the actual arguments used (including defaults)
      final_beta_func_args <- args
  }

  # ============================================================================
  # 4. Define Mixture Components for eps_i = beta_i - betabar_i
  # ============================================================================
  if (is.null(mixture_comps)) {
    # Generate default components if not provided
    mixture_comps <- vector("list", ncomp)
    base_mu <- rep(0, ncoef)
    base_rooti <- diag(sqrt(sigma_inv_diag), ncoef)

    if (ncomp == 1) {
      mixture_comps[[1]] <- list(mu = base_mu, rooti = base_rooti)
    } else {
      # Simple example: spread means along first dimension if possible
      # More sophisticated generation could be added here.
      spread <- if(ncoef > 0) seq(-1, 1, length.out = ncomp) else rep(0, ncomp)
      for (c in 1:ncomp) {
        comp_mu <- base_mu
        if(ncoef > 0) comp_mu[1] <- comp_mu[1] + spread[c]
        mixture_comps[[c]] <- list(mu = comp_mu, rooti = base_rooti)
      }
    }
  } else {
    # Validate provided mixture_comps structure
    if (!is.list(mixture_comps) || length(mixture_comps) != ncomp) {
      stop("Provided mixture_comps must be a list of length ncomp.")
    }
    for(c in 1:ncomp) {
      comp <- mixture_comps[[c]]
      if (!is.list(comp) || !all(c("mu", "rooti") %in% names(comp))) {
        stop(paste("Component", c, "missing 'mu' or 'rooti'."))
      }
      if (length(comp$mu) != ncoef) {
        stop(paste("Component", c, "mu has wrong length."))
      }
      if (!is.matrix(comp$rooti) || !all(dim(comp$rooti) == c(ncoef, ncoef))) {
        stop(paste("Component", c, "rooti has wrong dimensions."))
      }
    }
  }

  # Mixture weights (equal by default, could be argument later)
  # Using bayesm::rdirichlet requires bayesm in Imports
  pvec <- bayesm::rdirichlet(rep(2, ncomp))

  # ============================================================================
  # 5. Generate beta_true = betabar_true + eps_i
  # ============================================================================
  beta_true <- matrix(0, nrow = nlgt, ncol = ncoef)
  # Using bayesm::rmixture requires bayesm in Imports
  eps_draws <- bayesm::rmixture(nlgt, pvec, mixture_comps)
  # Ensure eps_draws$x is treated as a matrix even if nlgt=1 or ncoef=1
  eps_matrix <- matrix(eps_draws$x, nrow = nlgt, ncol = ncoef)
  beta_true <- betabar_true + eps_matrix

  # ============================================================================
  # 6. Generate Choice Data (lgtdata)
  # ============================================================================
  simlgtdata <- vector(mode = "list", length = nlgt)
  ni <- rep(nT, nlgt) # Number of observations for each unit i

  for (i in 1:nlgt) {
    ## ---- 6a. Simulate X covariates for unit i ----
    Xa <- Xa_dist_func(ni[i], p, nXa)
    Xd <- if(nXd > 0) Xd_dist_func(ni[i], nXd) else NULL

    ## ---- 6b. Create design matrix X ----
    # Assumes createX function is available (now specifically bayesm::createX)
    X <- bayesm::createX(p = p, na = nXa, nd = nXd, Xa = Xa, Xd = Xd, INT = const,
                         base = 1) # Assuming base alternative is 1

    ## ---- 6c. Simulate choices y_i ----
    # Assumes simmnlwX function is available (e.g., from this package)
    # Need to handle case where nT = 1 correctly within simmnlwX if necessary
    out_choices <- .simmnlwX(ni[i], X, beta_true[i, ])

    ## ---- 6d. Store results for unit i ----
    simlgtdata[[i]] <- list(y = out_choices$y,
                            X = X,
                            beta = beta_true[i, ],
                            betabar = betabar_true[i, ])
  }

  # ============================================================================
  # 7. Prepare Output List
  # ============================================================================
  output <- list(
    # Data structure matching rhierMnlRwMixture input
    p = p,
    lgtdata = simlgtdata,
    Z = Z,

    # Additional details about the simulation
    true_values = list(
      beta_true = beta_true,
      betabar_true = betabar_true,
      true_params = list(
        beta_func_type = beta_func_type,
        beta_func_args = final_beta_func_args, # Store args used (empty if nz=0)
        mixture_comps = mixture_comps,
        pvec = as.vector(pvec) # Store as vector
      ),
      dimensions = list(
        nlgt = nlgt, nT = nT, p = p, nz = nz, nXa = nXa, nXd = nXd,
        ncoef = ncoef, const = const, ncomp = ncomp
      )
    )
  )

  return(output)
}

# --- Helper Functions (to be potentially moved or refined) ---

# Placeholder for simmnlwX if not sourced/imported elsewhere
# simmnlwX <- function(n,X,beta) { ... }

# Placeholder for createX if not sourced/imported elsewhere
# createX <- function(p, na, nd, Xa, Xd, INT = TRUE, base = 1) { ... }

# --- Specific Beta Generation Functions ---

# These will be called within sim_hier_mnl based on beta_func_type

.beta_Z_const <- function(Zi, Delta) {
  # Zi: vector for one individual
  # Delta: nz x ncoef matrix
  as.vector(Zi %*% Delta)
}

.beta_Z_step <- function(Zi, Z_index, cutoff, beta_1, beta_2) {
  # Zi: vector for one individual
  # Z_index: column index of Z to use for the step
  # cutoff: threshold value
  # beta_1: coefficient vector if Zi[Z_index] > cutoff
  # beta_2: coefficient vector if Zi[Z_index] <= cutoff
  if (Zi[Z_index] > cutoff) {
    return(beta_1)
  } else {
    return(beta_2)
  }
}

# Note: Assumes nz >= 5
.beta_Z_friedman <- function(Zi, coefs = c(10, 20, 10, 5), shift = 0, scale = 8) {
  # Zi: vector for one individual (uses first 5 elements)
  # Returns a scalar value based on Friedman benchmark function
  term1 <- coefs[1] * sin(pi * Zi[1] * Zi[2])
  term2 <- coefs[2] * (Zi[3] - 0.5)^2 # Centered
  term3 <- coefs[3] * Zi[4]
  term4 <- coefs[4] * Zi[5]
  return((term1 + term2 + term3 + term4 + shift) / scale) # Scaled
}

# ==============================================================================
# Internal Helper: Simulate from MNL model conditional on X matrix
# ==============================================================================
.simmnlwX <- function(n, X, beta) {
  # n: number of observations (e.g., nT)
  # X: design matrix (stacked alternatives), (n*p) x k
  # beta: coefficient vector, k x 1
  k <- length(beta)
  if (nrow(X) %% n != 0) {
      stop("Internal error: nrow(X) not divisible by n in .simmnlwX")
  }
  p <- nrow(X) / n # Number of alternatives implied by X
  Xbeta <- X %*% beta

  # Reshape Xbeta into n x p matrix
  Xbeta <- matrix(Xbeta, byrow = TRUE, ncol = p)

  # Calculate probabilities (using safe softmax)
  max_utility <- apply(Xbeta, 1, max) # Find max utility for each observation
  Xbeta_scaled <- Xbeta - max_utility # Subtract max for numerical stability
  Prob <- exp(Xbeta_scaled)
  iota <- c(rep(1, p))
  denom <- Prob %*% iota
  Prob <- Prob / as.vector(denom)

  y <- vector("double", n)
  ind <- 1:p
  for (i in 1:n) {
    # Check for invalid probabilities (should not happen with safe softmax)
    if(any(is.na(Prob[i,])) || any(Prob[i,] < 0)) {
       warning(paste("Invalid probabilities generated for observation", i,
                   "in .simmnlwX. Check inputs."))
       # Handle potential error - e.g., assign a default choice or stop
       # For now, assign first choice
       y[i] <- 1
       next
    }
    # Ensure probabilities sum to 1 (within tolerance)
    prob_sum <- sum(Prob[i,])
    if (abs(prob_sum - 1.0) > 1e-6) {
        # Attempt to re-normalize if slightly off
        Prob[i,] <- Prob[i,] / prob_sum
    }

    yvec <- stats::rmultinom(1, 1, Prob[i, ])
    y[i] <- ind %*% yvec
  }#FOR i
  # Return only y, as X and beta are already known outside
  return(list(y = y))
} 