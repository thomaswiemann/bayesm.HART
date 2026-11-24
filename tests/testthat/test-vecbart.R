# Tests for the vector-leaf tree machinery (src/vectree.cpp,
# src/vecbartfuns.cpp, src/vecbd.cpp, src/vecbart.cpp) and its predictor.

library(testthat)

rand_spd <- function(D, scale = 1) {
  A <- matrix(rnorm(D * D), D, D)
  crossprod(A) / D + scale * diag(D)
}#RAND_SPD

# Marginal of the stacked response under y_i | m ~ N(m, Q_i^{-1}), m ~ N(0, C):
# N(0, J C J' + blockdiag(Q_i^{-1})) with J = 1_n (x) I_D.
log_marginal_stacked <- function(Y, Qs, C) {
  n <- nrow(Y); D <- ncol(Y)
  J <- kronecker(rep(1, n), diag(D))
  V <- J %*% C %*% t(J)
  for (i in seq_len(n)) {
    idx <- (i - 1) * D + seq_len(D)
    V[idx, idx] <- V[idx, idx] + solve(Qs[[i]])
  }
  yv <- as.vector(t(Y))
  -0.5 * n * D * log(2 * pi) -
    0.5 * as.numeric(determinant(V, logarithm = TRUE)$modulus) -
    0.5 * drop(t(yv) %*% solve(V, yv))
}#LOG_MARGINAL_STACKED

leaf_stats <- function(Y, Qs) {
  D <- ncol(Y)
  P <- matrix(0, D, D); b <- numeric(D)
  for (i in seq_len(nrow(Y))) { P <- P + Qs[[i]]; b <- b + Qs[[i]] %*% Y[i, ] }
  list(P = P, b = as.vector(b))
}#LEAF_STATS

test_that("veclh matches the stacked-Gaussian marginal up to unit constants", {
  set.seed(1)
  D <- 2; n <- 4
  Qs <- replicate(n, rand_spd(D, 0.5), simplify = FALSE)
  Y  <- matrix(rnorm(n * D), n, D)
  C  <- rand_spd(D, 0.3)
  st <- leaf_stats(Y, Qs)
  # Constants dropped by veclh: per-unit normalising terms and y_i' Q_i y_i.
  unit_const <- sum(vapply(seq_len(n), function(i)
    -0.5 * D * log(2 * pi) +
      0.5 * as.numeric(determinant(Qs[[i]], logarithm = TRUE)$modulus) -
      0.5 * drop(t(Y[i, ]) %*% Qs[[i]] %*% Y[i, ]), numeric(1)))
  expect_equal(veclh_R(st$P, st$b, n, C),
               log_marginal_stacked(Y, Qs, C) - unit_const, tolerance = 1e-10)
  # An empty leaf has marginal one.
  expect_equal(veclh_R(matrix(0, D, D), rep(0, D), 0L, C), 0)
})#TEST_THAT

test_that("veclh with D = 1 is the scalar BART leaf marginal plus log sigma", {
  n <- 7; sy <- 2.3; sigma <- 0.8; tau <- 0.4
  lh <- -0.5 * log(n * tau^2 + sigma^2) +
        tau^2 * sy^2 / (2 * sigma^2 * (n * tau^2 + sigma^2))
  expect_equal(veclh_R(matrix(n / sigma^2), sy / sigma^2, n, matrix(tau^2)),
               lh + log(sigma), tolerance = 1e-12)
})#TEST_THAT

test_that("vecdrawnodemu draws from the leaf conditional", {
  set.seed(2)
  D <- 3
  P <- rand_spd(D, 2); b <- rnorm(D); C <- rand_spd(D, 0.5)
  Pt <- solve(C) + P
  m <- solve(Pt, b); V <- solve(Pt)
  draws <- vecdrawnodemu_R(P, b, C, 20000L)
  se <- sqrt(diag(V) / 20000)
  expect_true(all(abs(colMeans(draws) - m) < 4 * se))
  expect_equal(cov(draws), V, tolerance = 0.1)
})#TEST_THAT

chain_problem <- function(seed = 3L, n = 40L, D = 2L, p = 2L) {
  set.seed(seed)
  Zt <- matrix(runif(p * n, -1, 1), p, n)
  f  <- cbind(ifelse(Zt[1, ] > 0, 1, -1), 0.5 * Zt[2, ])
  Y  <- f + matrix(rnorm(n * D, sd = 0.3), n, D)
  Q  <- array(0, c(D, D, n))
  for (i in seq_len(n)) Q[, , i] <- rand_spd(D, 5)
  list(Zt = Zt, Y = Y, Q = Q, D = D, C = 0.2 * diag(D))
}#CHAIN_PROBLEM

test_that("vecbart keeps split counts and leaf counts consistent", {
  pb <- chain_problem()
  out <- with(pb, vecbart_chain_R(Zt, Y, Q, m = 3L, C = C, alpha = 0.95,
                                  beta = 2, nmin = 5L, numcut = 20L,
                                  niter = 60L))
  final_leaves <- out$nleaves[length(out$nleaves)]
  # Every split adds one leaf, so interior nodes = leaves - trees.
  expect_equal(sum(out$nv), final_leaves - 3L)
  expect_true(all(out$nleaves >= 3L))
  expect_true(all(is.finite(out$leaf_sumsq)) && all(out$leaf_sumsq >= 0))
  expect_equal(dim(out$fits), c(40L, 2L, 60L))
  expect_true(all(is.finite(out$fits)))
})#TEST_THAT

test_that("serialized vector trees reproduce the in-sample fit", {
  pb <- chain_problem(seed = 4L)
  out <- with(pb, vecbart_chain_R(Zt, Y, Q, m = 4L, C = C, alpha = 0.95,
                                  beta = 1, nmin = 3L, numcut = 15L,
                                  niter = 40L))
  expect_equal(out$treedraws$leaf_dim, 2L)
  pred <- .pw_vec_trees(t(pb$Zt), out$treedraws, pb$D) # npred x nvar x ndraws
  expect_equal(dim(pred), c(40L, 2L, 40L))
  expect_equal(pred, out$fits, tolerance = 1e-7)
  # Transposed input and new rows evaluate without error and keep shape.
  Znew <- matrix(runif(2 * 5, -1, 1), 2, 5)
  pred_new <- .pw_vec_trees(Znew, out$treedraws, pb$D, transposed = TRUE)
  expect_equal(dim(pred_new), c(5L, 2L, 40L))
  expect_true(all(is.finite(pred_new)))
  expect_error(.pw_vec_trees(Znew[1, , drop = FALSE], out$treedraws, pb$D,
                             transposed = TRUE),
               regexp = "number of columns")
})#TEST_THAT

test_that("DART bookkeeping carries over to vector leaves", {
  pb <- chain_problem(seed = 5L)
  out <- with(pb, vecbart_chain_R(Zt, Y, Q, m = 3L, C = C, alpha = 0.95,
                                  beta = 2, nmin = 5L, numcut = 20L,
                                  niter = 40L, sparse = TRUE, burn = 2L))
  expect_equal(sum(out$pv), 1, tolerance = 1e-10)
  expect_true(all(out$pv > 0))
  expect_equal(sum(out$nv), out$nleaves[length(out$nleaves)] - 3L)
})#TEST_THAT
