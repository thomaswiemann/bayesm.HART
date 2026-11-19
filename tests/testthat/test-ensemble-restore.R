# Terminal ensemble state, original cutpoints and reconstructed fit buffers.

library(testthat)

ensemble_problem <- function(kind) {
  z <- seq(-1, 1, length.out = 60)
  Zt <- rbind(z, sin(seq_along(z)))
  y <- matrix(ifelse(z < 0, -2, 2) + .1 * cos(seq_along(z)), ncol = 1)
  if (kind == "varbart") y <- exp(y)
  if (kind == "vecbart") y <- cbind(y, -y / 2 + Zt[2, ])
  list(Zt = Zt, y = y)
}#ENSEMBLE_PROBLEM

for (kind in c("bart", "varbart", "vcbart", "vecbart")) {
  for (active in c(FALSE, TRUE)) {
    test_that(paste("ensemble round trip:", kind, "DART", active), {
      dat <- ensemble_problem(kind)
      set.seed(412)
      trained <- ensemble_state_R(kind, dat$Zt, dat$y, 4L, 80L, active)
      saved <- trained$snapshot
      expect_named(saved, c("cutpoints", "trees", "nv", "pv", "lpv", "theta", "dartOn"))
      expect_gt(trained$nleaves, 4)
      expect_identical(saved$dartOn, active)
      expect_type(saved$nv, "double")
      if (active) {
        expect_false(isTRUE(all.equal(saved$pv, c(.5, .5))))
        expect_false(identical(saved$theta, 1))
        expect_length(saved$lpv, 2L)
      } else if (kind != "vecbart") {
        expect_length(saved$lpv, 0L)
      }
      seed <- .Random.seed
      restored <- ensemble_state_R(kind, dat$Zt, dat$y, 4L, 0L, !active, saved)
      expect_identical(.Random.seed, seed)
      expect_identical(restored$snapshot, saved)
      expect_equal(restored$fit, trained$fit, tolerance = 1e-12)

      Zmore <- cbind(dat$Zt, c(-3, .4), c(4, -.7), c(.3, .2))
      ymore <- rbind(dat$y, matrix(0, 3, ncol(dat$y)))
      path <- tempfile(fileext = ".rds")
      saveRDS(saved, path)
      disk <- readRDS(path)
      unlink(path)
      expect_identical(disk, saved)
      extended <- ensemble_state_R(kind, Zmore, ymore, 4L, 0L, !active, disk)
      expect_identical(.Random.seed, seed)
      expect_identical(extended$snapshot, saved)
      expect_equal(extended$fit[seq_len(ncol(dat$Zt)), , drop = FALSE],
                   trained$fit, tolerance = 1e-12)
      independent <- .pw_eval(Zmore, saved, ncol(dat$y), kind == "varbart",
                               1L, TRUE, FALSE)
      expect_equal(extended$fit, matrix(independent, ncol = ncol(dat$y)),
                   tolerance = 1e-12)
      if (kind == "varbart") {
        summed <- .pw_eval(Zmore, saved, 1L, FALSE, 1L, TRUE, FALSE)
        expect_true(all(extended$fit > 0))
        expect_gt(diff(range(extended$fit)), .1)
        expect_gt(max(abs(extended$fit - matrix(summed, ncol = 1))), .1)
      }
      if (kind == "vecbart") {
        expect_true(all(trained$leaf_sumsq > 0))
        expect_identical(restored$leaf_sumsq, trained$leaf_sumsq)
      }
    })
  }
}

for (kind in c("bart", "vecbart")) {
  test_that(paste("ensemble dimension and stream errors:", kind), {
    dat <- ensemble_problem(kind)
    set.seed(412)
    saved <- ensemble_state_R(kind, dat$Zt, dat$y, 4L, 10L, TRUE)$snapshot
    expect_error(ensemble_state_R(kind, dat$Zt, dat$y, 5L, 0L, FALSE, saved),
                 "tree count mismatch")
    expect_error(ensemble_state_R(kind, rbind(dat$Zt, 0), dat$y, 4L, 0L, FALSE, saved),
                 "cutpoint dimension mismatch")
    bad <- saved
    bad$cutpoints <- bad$cutpoints[-1]
    expect_error(ensemble_state_R(kind, dat$Zt, dat$y, 4L, 0L, FALSE, bad),
                 "cutpoint dimension mismatch")
    bad <- saved
    bad$trees <- "broken"
    expect_error(ensemble_state_R(kind, dat$Zt, dat$y, 4L, 0L, FALSE, bad),
                 "failed tree stream read")
    header <- strsplit(saved$trees, "\n", fixed = TRUE)[[1]][1]
    for (body in c("", "0\n", "3\n1 0")) {
      bad$trees <- paste0(header, "\n", body)
      expect_error(ensemble_state_R(kind, dat$Zt, dat$y, 4L, 0L, FALSE, bad),
                   "failed tree stream read")
    }
    if (kind == "vecbart") {
      expect_error(ensemble_state_R(kind, dat$Zt, cbind(dat$y, 0), 4L, 0L, FALSE, saved),
                   "leaf dimension mismatch")
    }
  })
}

test_that("DART log probabilities survive probability underflow", {
  for (kind in c("bart", "varbart", "vcbart", "vecbart")) {
    dat <- ensemble_problem(kind)
    set.seed(412)
    saved <- ensemble_state_R(kind, dat$Zt, dat$y, 4L, 10L, TRUE)$snapshot
    saved$pv[1] <- 0
    saved$lpv[1] <- -1000
    seed <- .Random.seed
    restored <- ensemble_state_R(kind, dat$Zt, dat$y, 4L, 0L, FALSE, saved)
    expect_identical(restored$snapshot, saved)
    expect_identical(.Random.seed, seed)
  }
})
