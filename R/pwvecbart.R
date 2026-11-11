## Prediction for vector-leaf ensembles (vecbart).  Mirrors pwbart() but each
## tree returns one vector leaf, so the result is npred x nvar x ndraws.

pwvecbart <- function(
   x.test,                       # x matrix to predict at
   treedraws,                    # $treedraws from a vecbart entry
   transposed = FALSE
) {
   if (!transposed) x.test <- t(bartModelMatrix(x.test))
   p <- length(treedraws$cutpoints)
   if (p != nrow(x.test))
      stop(paste0("The number of columns in x.test must be equal to ", p))
   cpwvecbart(treedraws, x.test)
}
