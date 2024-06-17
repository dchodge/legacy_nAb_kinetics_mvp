# run main.R first to load the required packages

source(here::here("mvp", "bcell_david", "post", "post_utils.R"))

stanfit_name <- "nih_vac_wu_s_i"
stanfit_s <- readRDS(file = here::here("mvp", "bcell_david", "fit", paste0(stanfit_name, ".RDS")) )

plot_posteriors(stanfit_s, "nih_wu_s")
plot_marg_posteriors(stanfit_s, "nih_wu_s")
