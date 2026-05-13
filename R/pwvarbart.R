## Product-of-trees prediction for varbart ensembles.
## Mirrors the API of pwbart() but returns prod_l h_{jl}(Z*) rather than
## the BART sum-of-trees sum_l g_{jl}(Z*).

pwvarbart <- function(
   x.test,                       # x matrix to predict at
   treedraws,                    # $treedraws from a varbart entry
   mc.cores  = 1L,
   transposed = FALSE,
   dodraws   = TRUE,
   verbose   = FALSE
) {
   if (!transposed) x.test <- t(bartModelMatrix(x.test))
   p <- length(treedraws$cutpoints)
   if (p != nrow(x.test))
      stop(paste0("The number of columns in x.test must be equal to ", p))

   res <- .Call("cpwvarbart",
                treedraws,
                x.test,
                mc.cores,
                verbose)
   if (dodraws) return(res$yhat.test)
   else         return(apply(res$yhat.test, 2, mean))
}
