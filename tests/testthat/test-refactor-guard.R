# Pins every draw and prediction of the four heterogeneity models under the
# three likelihoods against fixtures, so a refactor that is meant to leave
# behaviour untouched is checked rather than assumed.  Gated: the fixtures
# are only valid while the RNG sequence is unchanged, and regenerating them
# is a deliberate act (tools/refactor-guard-fixtures.R).
#
# Fields the fixture does not carry are ignored, so a schema addition passes;
# a field it carries must match to 1e-10.

library(testthat)

skip_unless_guard <- function()
  skip_if_not(identical(Sys.getenv("BAYESM_HART_REFACTOR_GUARD"), "true"),
              "set BAYESM_HART_REFACTOR_GUARD=true for the refactor guard")

for (case in .guard_cases()) local({
  cs <- case
  test_that(paste("draws and predictions are unchanged:", cs$name), {
    skip_unless_guard()
    file <- file.path(.guard_dir(), paste0(cs$name, ".rds"))
    skip_if_not(file.exists(file), paste("missing fixture", basename(file)))
    want <- readRDS(file)
    got  <- .guard_run(cs)
    expect_true(all(names(want) %in% names(got)))
    for (nm in names(want))
      expect_equal(got[[nm]], want[[nm]], tolerance = 1e-10, info = nm)
  })#TEST_THAT
})
