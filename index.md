# bayesm.HART

`bayesm.HART` implements MCMC routines for hierarchical models with
Hierarchical Additive Regression Trees (HART) priors, as developed in
Wiemann (2025).

In the hierarchical logit specification, HART models the representative
consumer as a flexible function of observed characteristics. This
generalizes conventional hierarchical specifications that use a linear
projection on a small set of characteristics. The same framework
supports prediction for new consumers and posterior updating for
consumers with accumulating choice histories.

See the corresponding working paper [Personalization with
HART](https://thomaswiemann.com/assets/pdfs/jmp_wiemann.pdf) for further
discussion and details.

`bayesm.HART` builds on `bayesm` and `BART`. The syntax for estimating
the HART logit model mirrors `bayesm` syntax. In many existing `bayesm`
workflows, the main change is the `Prior` specification.

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

The following example applies the HART logit model to the Allenby and
Ginter (1995) bank conjoint dataset. The first code block loads the
data, formats the list structure required by `rhierMnlRwMixture`, and
sets MCMC hyperparameters. This setup syntax is the same in `bayesm` and
`bayesm.HART`.

``` r

# Load bayesm.HART for the sampler; load bayesm for the bank data
library(bayesm.HART)
library(bayesm)

# Load and prepare data from the 'bank' dataset (same as bayesm)
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
Z <- as.matrix(bank$demo[, -1]) # omit id
Z <- t(t(Z) - colMeans(Z)) # de-mean covariates as required by bayesm

# Final data object (same as bayesm)
Data <- list(lgtdata = lgtdata, Z = Z, p = 2)

# MCMC hyperparameters (same as bayesm)
Mcmc <- list(R = 2000, keep = 2, nprint = 500)
```

To fit the HART logit model, call `rhierMnlRwMixture`. Relative to a
linear hierarchical specification, the main change is adding the `bart`
entry in `Prior`. For illustration, the code below uses 50 trees per
factor and leaves other hyperparameters at defaults. See
[`?rhierMnlRwMixture`](https://thomaswiemann.com/bayesm.HART/reference/rhierMnlRwMixture.md)
for details.

``` r

# Fit the HART logit model
if (!model_cache_ok) {
  out <- bayesm.HART::rhierMnlRwMixture(
      Data = Data, Mcmc = Mcmc, 
      Prior = list(
        ncomp = 1, 
        bart = list(num_trees = 50) # new HART prior parameters
        ),
      r_verbose = F # suppress R print output (optional)
  )
  saveRDS(list(out = out, generated_at = Sys.time()), model_cache_file)
}
```

With posterior draws from the fitted model, we can summarize any
estimand. A central object is the *representative respondent*, i.e.,
expected part-worths conditional on respondent characteristics. The
following code computes posterior means and standard deviations for
three credit card attributes.

``` r

DeltaZ_hat <- predict(out, newdata = list(Z = Z), type = "DeltaZ+mu", 
                      burn = 250, r_verbose = FALSE)

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
#> Posterior Mean              5.047          4.152            -3.378
#> Posterior SD                0.902          0.901             0.989
```

## Learn More about `bayesm.HART`

See the `bayesm.HART` vignettes:

- [`vignette("bayesm-HART")`](https://thomaswiemann.com/bayesm.HART/articles/bayesm-HART.html):
  Get started (hierarchical logit with HART)
- [`vignette("bayesm-HART-linear")`](https://thomaswiemann.com/bayesm.HART/articles/bayesm-HART-linear.html):
  Hierarchical linear model
- [`vignette("bayesm-HART-negbin")`](https://thomaswiemann.com/bayesm.HART/articles/bayesm-HART-negbin.html):
  Hierarchical negative binomial model
- [`vignette("bayesm-HART-heteroskedastic-hart-bank")`](https://thomaswiemann.com/bayesm.HART/articles/bayesm-HART-heteroskedastic-hart-bank.html):
  Bank conjoint – heteroskedastic HART

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
Products and Segment Markets.” Journal of Marketing Research 32.4,
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

Wiemann, Thomas (2025). “[Personalization with
HART](https://thomaswiemann.com/assets/pdfs/jmp_wiemann.pdf).” Working
paper.
