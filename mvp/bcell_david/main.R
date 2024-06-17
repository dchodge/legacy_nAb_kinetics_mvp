if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, purrr, ggplot2, here, posterior, tidybayes,
    cmdstanr, rstan, deSolve, patchwork, Rcpp, devtools, ptmc, coda, bayesplot,
    extraDistr, hmer, loo, ggh4x)

#source(here::here("R/utils.R"))

inv.logit <- function(x) {
    exp(x)/(1+exp(x)) 
} 