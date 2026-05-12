library(knitr)

setwd("vignettes")
# takes too long to run on cran...
knit("bayesm-HART.Rmd.txt",
     "bayesm-HART.Rmd")
knit("bayesm-HART-linear.Rmd.txt",
     "bayesm-HART-linear.Rmd")
knit("bayesm-HART-negbin.Rmd.txt",
     "bayesm-HART-negbin.Rmd")
