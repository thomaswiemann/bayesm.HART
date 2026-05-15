# Precompile vignettes (too slow for CRAN).
#
# Usage:
#   Rscript vignettes/precompile.R              # compile all
#   Rscript vignettes/precompile.R hart negbin  # compile only matching vignettes
#
# Partial matching: any .Rmd.txt whose basename contains the argument is built.

library(knitr)

all_vignettes <- c(
  "vignettes/bayesm-HART.Rmd.txt",
  "vignettes/bayesm-HART-linear.Rmd.txt",
  "vignettes/bayesm-HART-negbin.Rmd.txt",
  "vignettes/bayesm-HART-heteroskedastic-hart-bank.Rmd.txt",
  "vignettes/marginal-effects.Rmd.txt",
  "README.Rmd"
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  matched <- unique(unlist(lapply(args, function(a) {
    grep(a, all_vignettes, value = TRUE, ignore.case = TRUE)
  })))
  if (length(matched) == 0) {
    stop("No vignettes matched: ", paste(args, collapse = ", "),
         "\nAvailable: ", paste(all_vignettes, collapse = ", "))
  }
  to_build <- matched
} else {
  to_build <- all_vignettes
}

library(bayesm.HART)

for (src in to_build) {
  if (src == "README.Rmd") {
    message(">>> Rendering: ", src, " -> README.md")
    rmarkdown::render(src)
  } else {
    out <- sub("\\.txt$", "", src)
    message(">>> Knitting: ", src, " -> ", out)
    
    # knit sets the working directory to the directory of the input file by default, 
    # so we can just pass the path.
    knit(src, out)
  }
}

message("Done. Built ", length(to_build), " file(s).")
