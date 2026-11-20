# Precompile vignettes (too slow for CRAN).
#
# Usage:
#   Rscript vignettes/precompile.R              # compile all
#   Rscript vignettes/precompile.R hart negbin  # compile only matching vignettes
#
# An argument is a case-insensitive regular expression matched against the
# paths below; every match is built.
#
# Each document is knitted in its own R process.  The convergence study forks
# one worker per core, and a shared session would hand each fork a copy of
# every earlier document's fitted draws.

all_vignettes <- c(
  "vignettes/bayesm-HART.Rmd.txt",
  "vignettes/bayesm-HART-linear.Rmd.txt",
  "vignettes/bayesm-HART-negbin.Rmd.txt",
  "vignettes/bayesm-HART-heteroskedastic-hart-bank.Rmd.txt",
  "vignettes/marginal-effects.Rmd.txt",
  "vignettes/simulation-convergence.Rmd.txt",
  "README.Rmd"
)

build_one <- function(src) {
  library(knitr)
  library(bayesm.HART)
  if (src == "README.Rmd") {
    message(">>> Rendering: ", src, " -> README.md")
    rmarkdown::render(src)
  } else {
    out <- sub("\\.txt$", "", src)
    message(">>> Knitting: ", src, " -> ", out)
    # Knit inside vignettes/ so the figure paths resolve there.
    owd <- setwd(dirname(src))
    on.exit(setwd(owd))
    knit(basename(src), basename(out))
  }
}#BUILD_ONE

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2L && args[1] == "--one") {
  build_one(args[2])
  quit(save = "no")
}

to_build <- all_vignettes
if (length(args) > 0) {
  to_build <- unique(unlist(lapply(args, function(a)
    grep(a, all_vignettes, value = TRUE, ignore.case = TRUE))))
  if (length(to_build) == 0)
    stop("No vignettes matched: ", paste(args, collapse = ", "),
         "\nAvailable: ", paste(all_vignettes, collapse = ", "))
}

rscript <- file.path(R.home("bin"), "Rscript")
for (src in to_build) {
  status <- system2(rscript, c("vignettes/precompile.R", "--one", src))
  if (status != 0)
    stop("Failed to build ", src, " (exit ", status, ").")
}

message("Done. Built ", length(to_build), " file(s).")
