# Exact check of the tree-structure kernels vecbd and vcbd on a tiny grid.
#
# With one characteristic and `numcut` cutpoints every tree structure can be
# enumerated.  The structure posterior is pi(R) prod_leaves M_g with pi(R) the
# Chipman-George-McCulloch prior as implemented in bartfuns: grow with
# alpha (1 + depth)^-beta when a split is possible, variable uniform over the
# splittable variables, cutpoint uniform over the node's remaining range, and
# support truncated to trees whose leaves all hold at least nmin units.  A
# long single-tree chain must visit structures with these frequencies.

library(testthat)

# All structures over cutpoint indices [lo, hi] (0-based, inclusive).
enum_trees <- function(lo, hi) {
  out <- list(list(split = NULL))
  if (hi < lo) return(out)
  for (c in lo:hi)
    for (L in enum_trees(lo, c - 1))
      for (R in enum_trees(c + 1, hi))
        out[[length(out) + 1]] <- list(split = c, left = L, right = R)
  out
}#ENUM_TREES

# Interior nodes (nid, cut), log prior, and leaf unit sets of one structure.
walk <- function(tr, x, cuts, lo, hi, depth, units, alpha, beta, nid = 1) {
  pg <- if (hi >= lo) alpha * (1 + depth)^(-beta) else 0
  if (is.null(tr$split))
    return(list(nodes = NULL, logprior = log(1 - pg), leaves = list(units)))
  c  <- tr$split
  lu <- units[x[units] <  cuts[c + 1]]
  ru <- units[x[units] >= cuts[c + 1]]
  L <- walk(tr$left,  x, cuts, lo, c - 1, depth + 1, lu, alpha, beta,
            2 * nid)
  R <- walk(tr$right, x, cuts, c + 1, hi, depth + 1, ru, alpha, beta,
            2 * nid + 1)
  list(nodes    = rbind(c(nid, c), L$nodes, R$nodes),
       logprior = log(pg) - log(hi - lo + 1) + L$logprior + R$logprior,
       leaves   = c(L$leaves, R$leaves))
}#WALK

key_of <- function(nodes) {
  if (is.null(nodes)) return("")
  nodes <- nodes[order(nodes[, 1]), , drop = FALSE]
  paste0(nodes[, 1], ":0:", nodes[, 2], ";", collapse = "")
}#KEY_OF

# Normalised structure posterior; `leaf_score(u)` is a leaf's log marginal for
# unit indices u (empty leaves score zero).
exact_posterior <- function(x, cuts, alpha, beta, nmin, leaf_score) {
  trees <- enum_trees(0, length(cuts) - 1)
  rows <- lapply(trees, function(tr) {
    w <- walk(tr, x, cuts, 0, length(cuts) - 1, 0, seq_along(x), alpha, beta)
    sizes <- vapply(w$leaves, length, integer(1))
    loglik <- sum(vapply(w$leaves, function(u)
      if (length(u) == 0) 0 else leaf_score(u), numeric(1)))
    data.frame(key = key_of(w$nodes), logpost = w$logprior + loglik,
               admissible = all(sizes >= nmin), stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df <- df[df$admissible, ]
  df$prob <- exp(df$logpost - max(df$logpost))
  df$prob <- df$prob / sum(df$prob)
  df
}#EXACT_POSTERIOR

leaf_stats <- function(Y, Qs) {
  D <- ncol(Y)
  P <- matrix(0, D, D); b <- numeric(D)
  for (i in seq_len(nrow(Y))) { P <- P + Qs[[i]]; b <- b + Qs[[i]] %*% Y[i, ] }
  list(P = P, b = as.vector(b))
}#LEAF_STATS

# Six units on a line, three equally spaced cutpoints (as makexinfo builds).
grid <- function() {
  x <- 1:6; numcut <- 3L
  list(x = x, numcut = numcut,
       cuts = min(x) + seq_len(numcut) * (max(x) - min(x)) / (numcut + 1))
}#GRID

# Vector-leaf problem: D = 2, step signal, per-unit precisions scaled by qscale
# (1e-6 gives a flat likelihood).
vec_problem <- function(seed = 7L, qscale = 1) {
  g <- grid(); x <- g$x; D <- 2
  set.seed(seed)
  Y <- cbind(ifelse(x > 3.5, 0.6, -0.6), ifelse(x > 2.3, 0.3, -0.3)) +
       matrix(rnorm(6 * D, sd = 0.5), 6, D)
  Qs <- lapply(seq_along(x), function(i) {
    A <- matrix(rnorm(D * D), D, D); qscale * (crossprod(A) / D + 0.6 * diag(D))
  })
  Q <- array(unlist(Qs), c(D, D, length(x))); C <- 0.5 * diag(D)
  c(g, list(
    leaf_score = function(u) {
      st <- leaf_stats(Y[u, , drop = FALSE], Qs[u])
      veclh_R(st$P, st$b, length(u), C)
    },
    chain = function(nmin, alpha, beta, niter)
      vecbart_chain_R(matrix(x, 1), Y, Q, m = 1L, C = C, alpha = alpha,
                      beta = beta, nmin = nmin, numcut = g$numcut,
                      niter = niter)))
}#VEC_PROBLEM

# Varying-coefficient problem: step phi, regressors a, variances v scaled by
# vscale (1e6 gives a flat likelihood).
vc_problem <- function(seed = 9L, vscale = 1) {
  g <- grid(); x <- g$x
  set.seed(seed)
  a <- rnorm(6); v <- vscale * runif(6, 0.5, 1.5); tau <- 0.6
  y <- a * ifelse(x > 3.5, 0.7, -0.7) + rnorm(6, sd = sqrt(v))
  c(g, list(
    leaf_score = function(u)
      vclh_R(sum(a[u]^2 / v[u]), sum(a[u] * y[u] / v[u]), tau),
    chain = function(nmin, alpha, beta, niter)
      vcbart_chain_R(matrix(x, 1), y, a, v, m = 1L, tau = tau, alpha = alpha,
                     beta = beta, nmin = nmin, numcut = g$numcut,
                     niter = niter)))
}#VC_PROBLEM

check_kernel <- function(pb, nmin, alpha = 0.95, beta = 1, niter = 150000L) {
  exact <- exact_posterior(pb$x, pb$cuts, alpha, beta, nmin, pb$leaf_score)
  out <- pb$chain(nmin, alpha, beta, niter)
  expect_equal(out$treedraws$cutpoints[[1]], pb$cuts, tolerance = 1e-12)
  keys <- out$keys[-seq_len(1000)]
  # Support: the chain never leaves the admissible set.
  expect_true(all(keys %in% exact$key))
  freq <- as.numeric(table(factor(keys, levels = exact$key))) / length(keys)
  expect_lt(max(abs(freq - exact$prob)), 0.02)
  # Every structure with appreciable mass is visited.
  expect_true(all(freq[exact$prob > 0.02] > 0))
  invisible(list(exact = exact, freq = freq))
}#CHECK_KERNEL

test_that("vecbd reproduces the tree prior under a flat likelihood", {
  res <- check_kernel(vec_problem(qscale = 1e-6), nmin = 1L)
  expect_equal(nrow(res$exact), 15L)
})#TEST_THAT

test_that("vecbd reproduces the exact structure posterior (nmin = 1)", {
  res <- check_kernel(vec_problem(), nmin = 1L)
  expect_equal(nrow(res$exact), 15L)
})#TEST_THAT

test_that("vecbd respects the nmin truncation (nmin = 2)", {
  res <- check_kernel(vec_problem(seed = 8L), nmin = 2L)
  expect_lt(nrow(res$exact), 15L)
})#TEST_THAT

test_that("vcbd reproduces the tree prior under a flat likelihood", {
  res <- check_kernel(vc_problem(vscale = 1e6), nmin = 1L)
  expect_equal(nrow(res$exact), 15L)
})#TEST_THAT

test_that("vcbd reproduces the exact structure posterior (nmin = 1)", {
  res <- check_kernel(vc_problem(), nmin = 1L)
  expect_equal(nrow(res$exact), 15L)
})#TEST_THAT

test_that("vcbd respects the nmin truncation (nmin = 2)", {
  res <- check_kernel(vc_problem(seed = 10L), nmin = 2L)
  expect_lt(nrow(res$exact), 15L)
})#TEST_THAT
