#' bayesm.HART: Hierarchical Logit Models with Bayesian Additive Regression Trees
#'
#' @description The `bayesm.HART` package provides a hierarchical Bayesian 
#'   machine learning approach for personalization. At its core, the package 
#'   implements a Metropolis-within-Gibbs sampler for a hierarchical logit 
#'   model with a Hierarchical Additive Regression Trees (HART) prior. HART is a
#'   nonparametric prior that models the "representative consumer" as a flexible 
#'   sum-of-trees function of observed characteristics. This extends 
#'   conventional hierarchical models that are often limited to linear functions 
#'   of few characteristics.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib bayesm.HART, .registration = TRUE
## usethis namespace: end
NULL
