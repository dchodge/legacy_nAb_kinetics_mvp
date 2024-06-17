if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, purrr, ggplot2, here, posterior, tidybayes,
    rstan, deSolve, patchwork, Rcpp, devtools, coda, bayesplot,
    extraDistr, hmer, loo, ggh4x)


# install cmdstanr
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

#source(here::here("R/utils.R"))

inv.logit <- function(x) {
    exp(x)/(1+exp(x)) 
} 