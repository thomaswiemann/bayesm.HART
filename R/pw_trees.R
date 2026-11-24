# Evaluate a serialized ensemble at Z.  `product = TRUE` is the variance
# ensembles, whose draws multiply; otherwise the trees are summed.  Scalar
# leaves return an (R/keep) x npred matrix of draws, or the posterior mean
# over draws when `dodraws = FALSE`; vector leaves return npred x nvar x
# (R/keep).
.pw_trees <- function(x.test, treedraws, product = FALSE, mu = 0,
                      mc.cores = 1L, transposed = FALSE, dodraws = TRUE,
                      verbose = FALSE) {
  res <- .pw_eval(x.test, treedraws, 1L, product, mc.cores, transposed,
                  verbose)
  res <- t(res[, 1, ])                                   # ndraws x npred
  if (dodraws) res + mu else apply(res, 2, mean) + mu
}#.PW_TREES

.pw_vec_trees <- function(x.test, treedraws, leaf_dim, transposed = FALSE,
                          mc.cores = 1L, verbose = FALSE) {
  .pw_eval(x.test, treedraws, leaf_dim, FALSE, mc.cores, transposed, verbose)
}#.PW_VEC_TREES

.pw_eval <- function(x.test, treedraws, leaf_dim, product, mc.cores,
                     transposed, verbose) {
  if (!transposed) x.test <- t(bartModelMatrix(x.test))
  p <- length(treedraws$cutpoints)
  if (p != nrow(x.test))
    stop(paste0("The number of columns in x.test must be equal to ", p))
  cpwtrees(treedraws, x.test, as.integer(leaf_dim), as.integer(mc.cores),
           product, verbose)
}#.PW_EVAL
