# run main.R first to load the required packages


# Run HMC sampler in stan
# Change iter_warmup and iter_sampling to 1000 for a real run, just short here for testing
run_sim <- function(name, model, data_list, folder_name) {
    fit_ode_model <- model$sample(
        data = data_list,    # named list of data
        refresh = 1,
        iter_warmup = 10,          # number of warmup iterations per chain
        iter_sampling = 10,            # total number of iterations per chain
        chains = 4,            # number of cores (could use one per chain)
        parallel_chains = 4,
        adapt_delta = 0.8,
        max_treedepth = 10,
        init = 0
    )
    fit_ode_model$save_object(file = here::here("mvp", "bcell_david", "fit", paste0(name, ".RDS")))
}

#Load and compile the stan model
ode_model <- cmdstan_model(here::here("mvp", "bcell_david", "fit", "cauchy_rk45_waneLLPC_ind.stan"), compile = TRUE)
ode_model_m3_ind <- ode_model$compile(stanc_options = list("O1", "canonicalize"))

# Load the stan data
data_list_nih_s <- readRDS(here::here("mvp", "bcell_david", "clean", "nih_wu_s", "list_stan_data.RDS"))

# Run the HMC sampler
run_sim("nih_vac_wu_s_i", ode_model_m3_ind, data_list_nih_s)