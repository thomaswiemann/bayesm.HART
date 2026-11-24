# summary()$r2 under a constant covariance: the residual variance is the
# diagonal of the Sigma that predict(type = "SigmaZ") returns.

test_that("homoskedastic R2 uses the diagonal of Sigma", {
  data <- .test_data("mnl")
  for (model in c("none", "normalized", "constant")) {
    fit <- .test_fit("mnl", data, .test_priors[[model]], seed = 21L)
    burn <- 10L
    r2 <- summary(fit, Z = data$Z, burn = burn)$r2
    Delta <- predict(fit, newdata = list(Z = data$Z), type = "DeltaZ",
                     burn = burn, r_verbose = FALSE)
    Sigma <- predict(fit, newdata = list(Z = data$Z[1, , drop = FALSE]),
                     type = "SigmaZ", burn = burn, r_verbose = FALSE)
    for (k in seq_len(dim(Delta)[2])) {
      explained <- apply(Delta[, k, ], 2, var)
      expect_equal(r2$R2_draws[, k], explained / (explained + Sigma[1, k, k, ]),
                   ignore_attr = TRUE, info = model)
    }
  }
})
