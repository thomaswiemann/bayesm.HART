# Regenerates the refactor-guard fixtures from the installed package.
# Run from the package root:  Rscript tools/refactor-guard-fixtures.R
#
# The fixtures pin the RNG sequence, so regenerate only when a change is
# meant to alter it (the unified birth/death step), never to silence a
# failing guard.

library(bayesm.HART)
source(file.path("tests", "testthat", "helper-refactor-guard.R"))

dir.create(.guard_dir(), recursive = TRUE, showWarnings = FALSE)
for (case in .guard_cases()) {
  rec <- .guard_run(case)
  saveRDS(rec, file.path(.guard_dir(), paste0(case$name, ".rds")),
          compress = "xz", version = 2)
  cat(sprintf("%-28s %2d fields  %6.1f KB\n", case$name, length(rec),
              file.size(file.path(.guard_dir(),
                                  paste0(case$name, ".rds"))) / 1024))
}
