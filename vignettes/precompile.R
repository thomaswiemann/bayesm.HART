library(knitr)

setwd("vignettes")
# takes too long to run on cran...
knit("bayesm-HART.Rmd.txt",
     "bayesm-HART.Rmd")
