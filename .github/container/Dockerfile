# syntax=docker/dockerfile:1
# R 4.5 + bayesm.HART at the commit being built + the simulation stack.
# Built by .github/workflows/container.yml; the build context holds the
# package source under bayesm.HART/.
#
# Layers are ordered least-volatile first, so the CRAN compiles survive in the
# build cache across package commits. Labels carry the commit and go last.

FROM rocker/r-ver:4.5.3

# Apptainer bind-mounts $HOME, and R derives R_LIBS_USER from it, so a stray
# package in the cluster home would shadow the ones installed below.
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# --- System libraries (ggplot2's graphics stack) --------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
    libpng-dev libjpeg-dev libtiff-dev libfontconfig1-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev libcairo2-dev \
    && rm -rf /var/lib/apt/lists/*

# --- Compile toolchain ----------------------------------------------------
RUN install2.r --error --skipinstalled -n 4 \
    Rcpp RcppArmadillo \
    && rm -rf /tmp/downloaded_packages

# --- bayesm.HART Imports, from DESCRIPTION --------------------------------
RUN install2.r --error --skipinstalled -n 4 \
    rlang nnet bayesm dplyr ggplot2 \
    && rm -rf /tmp/downloaded_packages

# --- Simulation stack (edit here) -----------------------------------------
#   posterior  rhat and effective sample sizes
#   digest     run fingerprints
RUN install2.r --error --skipinstalled -n 4 \
    posterior digest \
    && rm -rf /tmp/downloaded_packages

# --- bayesm.HART ----------------------------------------------------------
COPY bayesm.HART /tmp/bayesm.HART
RUN R CMD INSTALL /tmp/bayesm.HART && rm -rf /tmp/bayesm.HART \
    && Rscript -e 'stopifnot(requireNamespace("bayesm.HART", quietly = TRUE)); \
                   cat("bayesm.HART", as.character(packageVersion("bayesm.HART")), "\n")'

ARG HART_SHA=
LABEL org.opencontainers.image.description="R 4.5 + bayesm.HART + simulation stack"
LABEL bayesm.hart.revision="${HART_SHA}"
