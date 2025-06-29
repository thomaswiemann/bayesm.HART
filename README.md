
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesm.HART

`bayesm.HART` implements a Metropolis-within-Gibbs sampler for the
Hierarchical Additive Regression Trees (HART) logit model proposed in
Wiemann (2025).

HART generalizes conventional hierarchical models by defining the representative consumer as a flexible function of potentially many characteristics. Drawing on the Bayesian Additive Regression Trees (BART) of Chipman, George, and McCulloch (2010), HART specifies this function as a sum-of-trees factor model. HART's combination of a flexible nonparametric prior within the hierarchical model provides a coherent framework for (Bayes-) optimal managerial decisions that adapt to the firm's familiarity with the consumer: first, HART flexibly leverages observed characteristics for granular predictions about new consumers; second, as a consumer's choices accumulate, their individual-level preferences adaptively deviate from this representative unit.

See the corresponding working paper [Personalization with HART](https://thomaswiemann.com/assets/pdfs/wiemann_jmp.pdf) for further discussion and details.

## Installation

Install the latest development version from GitHub (requires
[devtools](https://github.com/r-lib/devtools) package):

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("thomaswiemann/bayesm.HART", dependencies = TRUE)
```

## Example

The following example applies the HART logit model to the canonical
conjoint dataset of Allenby and Ginter (1995) on credit card design,
which is included in the `bayesm` package. The first code block loads
the data and formats it into the list structure required by
`rhierMnlRwMixture`. It then fits the HART logit model using a
Metropolis-within-Gibbs sampler. We specify a sum-of-trees prior with 50
trees for the representative consumer and run the MCMC algorithm for
2,000 iterations, keeping every second draw.

``` r
library(bayesm.HART)
library(bayesm)

# Load and prepare data from the 'bank' dataset
data(bank)
choiceAtt <- bank$choiceAtt
hh <- levels(factor(choiceAtt$id))
nhh <- length(hh)
lgtdata <- vector("list", length = nhh)
for (i in 1:nhh) {
    y = 2 - choiceAtt[choiceAtt[,1]==hh[i], 2]
    nobs = length(y)
    X_temp = as.matrix(choiceAtt[choiceAtt[,1]==hh[i], c(3:16)])
    X = matrix(0, nrow = nrow(X_temp) * 2, ncol = ncol(X_temp))
    X[seq(1, nrow(X), by = 2), ] = X_temp
    lgtdata[[i]] = list(y=y, X=X)
}
Z <- as.matrix(bank$demo[, -1])
Z <- t(t(Z) - colMeans(Z))
Data <- list(lgtdata = lgtdata, Z = Z, p = 2)

# Set MCMC parameters and suppress sampler output
Mcmc <- list(R = 2000, keep = 2, nprint = 500)

# Fit the HART logit model
out <- bayesm.HART::rhierMnlRwMixture(
    Data = Data, Mcmc = Mcmc, 
    Prior = list(ncomp = 1, bart = list(num_trees = 50))
)
#> Table of Y values pooled over all units
#> ypooled
#>    1    2 
#> 6473 8326 
#>  
#> Starting MCMC Inference for Hierarchical Logit:
#>    Normal Mixture with 1 components for first stage prior
#>    2  alternatives;  14  variables in X
#>    for  946  cross-sectional units
#>  
#> Prior Parms: 
#> nu = 17
#> V 
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]   17    0    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    0   17    0    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    0   17    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    0    0   17    0    0    0    0    0     0     0     0     0
#>  [5,]    0    0    0    0   17    0    0    0    0     0     0     0     0
#>  [6,]    0    0    0    0    0   17    0    0    0     0     0     0     0
#>  [7,]    0    0    0    0    0    0   17    0    0     0     0     0     0
#>  [8,]    0    0    0    0    0    0    0   17    0     0     0     0     0
#>  [9,]    0    0    0    0    0    0    0    0   17     0     0     0     0
#> [10,]    0    0    0    0    0    0    0    0    0    17     0     0     0
#> [11,]    0    0    0    0    0    0    0    0    0     0    17     0     0
#> [12,]    0    0    0    0    0    0    0    0    0     0     0    17     0
#> [13,]    0    0    0    0    0    0    0    0    0     0     0     0    17
#> [14,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>       [,14]
#>  [1,]     0
#>  [2,]     0
#>  [3,]     0
#>  [4,]     0
#>  [5,]     0
#>  [6,]     0
#>  [7,]     0
#>  [8,]     0
#>  [9,]     0
#> [10,]     0
#> [11,]     0
#> [12,]     0
#> [13,]     0
#> [14,]    17
#> mubar 
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0     0
#> Amu 
#>      [,1]
#> [1,] 0.01
#> a 
#> [1] 5
#> 
#> BART prior parameters:
#> num_trees = 50
#> power = 2
#> base = 0.95
#> tau = 0.141421
#> numcut = 100
#> sparse = 0
#>  
#> MCMC Parms: 
#> s= 0.636  w=  0.1  R=  2000  keep=  2  nprint=  500
#> 
#> initializing Metropolis candidate densities for  946  units ...
#>   completed unit # 50
#>   completed unit # 100
#>   completed unit # 150
#>   completed unit # 200
#>   completed unit # 250
#>   completed unit # 300
#>   completed unit # 350
#>   completed unit # 400
#>   completed unit # 450
#>   completed unit # 500
#>   completed unit # 550
#>   completed unit # 600
#>   completed unit # 650
#>   completed unit # 700
#>   completed unit # 750
#>   completed unit # 800
#>   completed unit # 850
#>   completed unit # 900
#>  MCMC Iteration (est time to end - min) 
#>  500 (0.9)
#>  1000 (0.6)
#>  1500 (0.3)
#>  2000 (0.0)
#>  Total Time Elapsed: 1.22
```

With the MCMC draws from the fitted model, we can compute posterior
estimates for any quantity of interest. A key object is the
representative consumer, which represents the expected
part-worths for a consumer with characteristics. The following
code computes the posterior mean and standard deviation of these
expected part-worths for three credit card attributes.

``` r
DeltaZ_hat <- predict(out, newdata = Data, type = "DeltaZ+mu", 
                      burn = 250, r_verbose=F)

posterior_mean <- apply(DeltaZ_hat, 2, mean)
posterior_sd <- apply(DeltaZ_hat, 2, sd)

# Indices for the desired coefficients:
# 2: Interest Low Fixed
# 8: Annual Fee Low
# 10: Bank Out-of-State
selected_indices <- c(2, 8, 10)

# Create a matrix with means and sds as rows
results_matrix <- rbind(
  `Posterior Mean` = posterior_mean[selected_indices],
  `Posterior SD`   = posterior_sd[selected_indices]
)

# Convert to a data frame and set column names
results_df <- as.data.frame(results_matrix)
colnames(results_df) <- c("Interest Low Fixed", "Annual Fee Low", "Bank Out-of-State")

# Print the data frame to the console
print(results_df, digits = 3)
#>                Interest Low Fixed Annual Fee Low Bank Out-of-State
#> Posterior Mean              5.154           4.43             -3.44
#> Posterior SD                0.959           1.01              1.03
```

## Learn More about `bayesm.HART`

Check out our articles to learn more:

- `vignette("bayesm.HART")` provides a more detailed introduction
- `vignette("marginal-effects")` discusses how to compute and plot
  marginal effects.

## Acknowledgements

`bayesm.HART` originated as a fork of the `bayesm` and `BART` packages.
Its current implementation heavily leverages the codebase and
foundational work from both packages. I gratefully acknowledge the
contributions of their respective authors:

- [**`bayesm`**](https://cran.r-project.org/web/packages/bayesm/index.html):
  Peter Rossi
- [**`BART`**](https://cran.r-project.org/web/packages/BART/index.html):
  Robert McCulloch, Rodney Sparapani, Robert Gramacy, Matthew Pratola,
  Charles Spanbauer, Martyn Plummer, Nicky Best, Kate Cowles, Karen
  Vines

## References

Allenby, Greg M. and James L. Ginter (1995). “Using Extremes to Design
Products and Segment Markets.” *Journal of Marketing Research* 32.4,
pp. 392–403.

Chipman, Hugh A., Edward I. George, and Robert E. McCulloch (2010).
“BART: Bayesian Additive Regression Trees.” Annals of Applied Statistics
4.1.

Rossi, Peter E., Greg M. Allenby, and Robert McCulloch (2009). Bayesian
Statistics and Marketing. Reprint. Wiley Series in Probability and
Statistics. Chichester: Wiley.

Rossi, Peter (2023). bayesm: Bayesian Inference for
Marketing/Micro-Econometrics. Comprehensive R Archive Network.

Sparapani, Rodney, Charles Spanbauer, and Robert McCulloch (2021).
“Nonparametric Machine Learning and Efficient Computation with Bayesian
Additive Regression Trees: The BART R Package.” Journal of Statistical
Software 97, pp. 1–66.

Wiemann, Thomas (2025). “Personalization with HART.” Working paper.
