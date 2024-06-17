#' @title Clean utils
#' @description Utility functions for cleaning the data
#' @author David Hodgson
#' 

#' Add groupings to the meta data
#' @description Adds groupings to the meta data
#' @param meta_pre Dataframe of the meta data
#' @return Dataframe with groupings
make_meta_groups <- function(meta_pre) {
    meta_pre %>% mutate(time_d1 = case_when(
        (d1_d2int < 28) ~"<28 days",
        (d1_d2int >= 28) ~"28+ days"
    )) %>%     
    mutate(age_cat = cut(age, breaks = c(0, 30, 40, 50, 60, 100), 
        labels = c("<30 years", "30–39 years", "40–49 years", "50–59 years", "60+ years") 
    ))
}

#' @param df_wide Dataframe in wide format
#' @return Dataframe in long format
#' @description Converts a dataframe from wide to long format
#' 
covert_data_long <- function(df_wide) {
    complete_time <- df_wide$VisitN %>% unique %>% sort
    df_wide %>%
        group_by(PID) %>% 
        complete(VisitN = complete_time) %>% 
        mutate(IC50hs = log10(IC50hs)) %>% 
        mutate(time_post_exp = case_when(
            VisitN == 0~0,VisitN != 0 ~  time_post_exp
        ))  %>% replace_na(list(Bmem = -1, pbBmem = -1, IC50hs = -1) )
}

#' Convert data to stan format 
#' @description Takes the data in long format and converts it to a list for use in stan
#' @param df_model_full Dataframe in long format
#' @param study Character string of the study (s or rbd)
#' @param antigen_str Character string of the antigen
#' @param T Number of time points
create_datalist_stan <- function(df_model_full, study, antigen_str, T) {

    init_av <- df_model_full %>% filter(VisitN == 0) %>% ungroup %>%
            filter(!is.na(SampleID)) %>% unique %>%
            summarise(Bmem = mean(Bmem, na.rm = TRUE), pbBmem = mean(pbBmem, na.rm = TRUE) ) %>% as.matrix

    init_av <- c(init_av, 0)

    init_list <- list()
    y_list <- list()
    t_vals <- list()

    PIDs <- df_model_full$PID %>% unique
    for(PID_i in PIDs) {
        PID_i_info <- df_model_full %>% filter(PID == PID_i) %>% ungroup %>% 
            select(time_post_exp, Bmem, pbBmem, IC50hs) %>%
            filter(!is.na(time_post_exp)) %>% replace_na( 
            list(Bmem = -1, pbBmem = -1, IC50hs = -1) 
            ) %>% as.matrix 
        t_vals[[PID_i]] <- as.numeric(PID_i_info[-1, 1])
        init_list_t <- c(PID_i_info[1, 2:4])
        if(init_list_t[3] == -1) {
            init_list_t[3] <- 0
        }
        init_list[[PID_i]] <- init_list_t
        y_list[[PID_i]] <-  matrix(PID_i_info[-1, 2:4], c(nrow(PID_i_info[-1, 2:4]), 3))
        T = dim(y_list[[PID_i]])[2]
        if (T == 1) {
            y_list[[PID_i]] <- y_list[[PID_i]] %>% t
        }
    }

    df_meta <- df_model_full %>% select(PID, vax_type, time_d1, age_cat) %>% filter(!is.na(vax_type)) %>% unique

    saveRDS(df_meta, here::here("mvp", "bcell_david", "clean", paste0(study, "_", antigen_str), paste0("df_meta.RDS")))

    meta_data_comp <- compose_data(df_meta)
    x_cov <- matrix(c(meta_data_comp$vax_type - 1,
    meta_data_comp$time_d1 - 1,
    meta_data_comp$age_cat - 1), ncol = 3)

    data_model_stan <- list(
        N = meta_data_comp$n,
        N_id = meta_data_comp$PID,
        N_vac = meta_data_comp$vax_type,
        x_cov = x_cov,
        t0 = 0,
        tf = 220,
        y0_init = init_list,
        y_obs = y_list,
        init_av = init_av,
        y_times = t_vals
    )

    N <- data_model_stan$N
    N_id <- data_model_stan$N_id
    x_vac <- data_model_stan$x_cov[, 1] + 1
    x_time <- data_model_stan$x_cov[, 2]  + 1
    x_age <- data_model_stan$x_cov[, 3]  + 1

    init_av <- data_model_stan$init_av

    prop_vac <- x_vac %>% table(.) %>% `/`(., sum(.)) %>% as.numeric
    prop_time <- x_time %>% table(.) %>% `/`(., sum(.)) %>% as.numeric 
    prop_age <- x_age %>% table(.) %>% `/`(., sum(.)) %>% as.numeric 
    
    y0_list <- data_model_stan$y0_init %>% unlist %>% matrix(c(N, 3), byrow = TRUE) %>% cbind(rep(0, N)) %>% cbind(rep(0, N)) %>% cbind(rep(0, N)) 
    for(i in 1:nrow(y0_list)) {
        if(y0_list[i, 1] == -1) {
            y0_list[i, ] <- c(init_av, 0, 0, 0)
        }
    }
    # Get priors overall
    y0_list_vac <- cbind(y0_list, x_vac) %>% as.data.frame %>%
        summarise(mean(V1), mean(V2), mean(V3), mean(V4), mean(V5), mean(V6), .by = x_vac) %>% arrange(x_vac) %>% .[2:7] %>% as.matrix

    # Get priors overall
    y0_list_time <- cbind(y0_list, x_time) %>% as.data.frame %>%
        summarise(mean(V1), mean(V2), mean(V3), mean(V4), mean(V5), mean(V6), .by = x_time) %>% arrange(x_time) %>% .[2:7] %>% as.matrix

    # Get priors overall
    y0_list_age <- cbind(y0_list, x_age) %>% as.data.frame %>%
        summarise(mean(V1), mean(V2), mean(V3), mean(V4), mean(V5), mean(V6), .by = x_age) %>% arrange(x_age) %>% .[2:7] %>% as.matrix

    y_times <- data_model_stan$y_times
    for(i in 1:length(data_model_stan$y_times)) {
        n <- length(y_times[[i]])
        y_times[[i]] <- c(y_times[[i]], rep(-1, T - n))
    }

    y_times <- y_times %>% unlist %>% matrix(c(N, T), byrow = TRUE)
    t_max <- data_model_stan$y_times %>% map(~max(.x)) %>% unlist %>% as.numeric
    t_length <- data_model_stan$y_times %>% map(~length(.x)) %>% unlist %>% as.numeric

    y_obs <- data_model_stan$y_obs
    N_likelihood <- 0
    for(i in seq_along(1:length(y_obs))) {
        N_likelihood <- N_likelihood + ((y_obs[[i]] > -1) %>% sum)
        tvals <- y_obs[[i]] %>% dim %>% .[1]
        tvals_j <- y_obs[[i]] %>% dim %>% .[2]
        if (tvals == T) {
            #y_obs[[i]][y_obs[[i]] == 0] <- -1
        } else {
            null_add <- replicate(T - tvals, c(-1, -1, -1), simplify = FALSE)
            y_obs[[i]] <- do.call(rbind, c(list(y_obs[[i]]), null_add))
            #y_obs[[i]][y_obs[[i]] == 0] <- -1
        }
    }

    data_list_stan_model <- list(

        N = N,
        N_id = N_id,

        N_vac = 2,
        x_vac = x_vac,

        N_time = 2,
        x_time = x_time,

        N_age = 5,
        x_age = x_age,

        prop_vac = prop_vac,
        prop_time = prop_time,
        prop_age = prop_age,

        T = T,
        y_times = y_times,
        t_max = t_max,
        t_length = t_length,

        y0_init = y0_list,
        y_obs = y_obs,

        y0_init_vac = as.matrix(y0_list_vac),
        y0_init_time = as.matrix(y0_list_time),
        y0_init_age = as.matrix(y0_list_age)

    )

    saveRDS(data_list_stan_model, here::here("mvp", "bcell_david", "clean", paste0(study, "_", antigen_str), "list_stan_data.RDS"))

}

#' Clean the NIH data
#' @description Takes the raw sVNT and B cell data and convets into one data frame, useful for later analysis. Also produces the data in stan format.
#' @param antigen_str Character string of the antigen (s or rbd)
#' @return Dataframe in long format
clean_nih_data <- function(antigen_str) {
    study <- "nih"

    df_svnt <- read.csv(here::here("mvp", "bcell_david", "clean", "dataraw", "paired B and sVNT_for_DHodgson.csv"))
    df_svnt_long <- read.csv(here::here("mvp", "bcell_david", "clean", "dataraw", "COVAX_longitudinal_analysis.csv"))
    df_bio_nih <- read.csv(here::here("mvp", "bcell_david", "clean", "dataraw",  "CoVax_probe_analysis.csv"))
    # Get meta
    df_bio_nih_meta <- df_bio_nih %>% select(pid = PID, sex = gender, age = age_screening, vax_type = Covax2_brand, d1_d2int) %>%
        unique %>%  mutate(vax_type = recode(vax_type, "AZ" = "ChAdOx1", "P" = "BNT162b2")) %>%
        make_meta_groups %>% rename(PID = pid)


    if (antigen_str == "wu_s") {
        df_cleaned_bcell_nih <- df_bio_nih %>%
            select(SampleID, PID, VisitN, 
                nInfections, d2bleedint, d1_d2int, WuSofBmem, pbOfWuS) %>%
            mutate(WuSpbofBmem = WuSofBmem * pbOfWuS / 100, WuSofBmem = WuSofBmem - WuSpbofBmem) %>% filter(is.na(nInfections)) %>% select(!c(d1_d2int)) %>% 
            mutate(time_post_exp = case_when(
                VisitN == 0~0,VisitN != 0 ~  d2bleedint
            )) %>% rename(Bmem = WuSofBmem, pb = pbOfWuS, pbBmem = WuSpbofBmem) %>% mutate(antigen = "Ancestral spike")
    } else if (antigen_str == "wu_rbd") {
        df_cleaned_bcell_nih <- df_bio_nih %>%
            select(SampleID, PID, VisitN, 
                nInfections, d2bleedint, d1_d2int, WuRBDofBmem, pbOfWuRBD) %>%
            mutate(WuRBDpbofBmem = WuRBDofBmem * pbOfWuRBD / 100, WuRBDofBmem = WuRBDofBmem - WuRBDpbofBmem) %>% filter(is.na(nInfections)) %>% select(!c(d1_d2int)) %>% 
            mutate(time_post_exp = case_when(
                VisitN == 0~0,VisitN != 0 ~  d2bleedint
            )) %>% rename(Bmem = WuRBDofBmem, pb = pbOfWuRBD, pbBmem = WuRBDpbofBmem) %>% mutate(antigen = "Ancestral RBD")
    }


    df_meta <- df_bio_nih_meta
    # Remove missing info
    PIDs <- df_cleaned_bcell_nih %>% pull(PID) %>% unique 
    PIDSlol <- df_svnt %>% pull(PID) %>% unique 
    PIDs <- intersect(PIDSlol, PIDs) # "WCH-801" "ALF-807" . are excluded as they have been previously infected

    df_svnt_trim <- 
        bind_rows(df_svnt %>% filter(PID %in% PIDs) %>%
            select(PID, Visit, IC50hs = sVNT_W_PV.IC50) %>% 
            mutate(VisitN = substr(Visit, 2, 3)),
            df_svnt %>% filter(PID %in% PIDs) %>% select(PID, IC50hs = sVNT_V220_2021.IC50 ) %>% 
                mutate(VisitN = "220" ) %>% filter(!is.na(IC50hs)) %>% unique #%>% left_join(timings_220)
        ) %>% select(!Visit) %>% mutate(VisitN = as.numeric(VisitN))

    joining_dates <- bind_rows(
        df_cleaned_bcell_nih %>% select(PID, VisitN, time_post_exp),
        df_svnt_long %>% select(Sample_ID, PID, VisitN = V, time_post_exp = DaysSinceD2) %>% filter(VisitN == 220)
    ) %>% select(!Sample_ID)


    df_bio_trim_join <- df_cleaned_bcell_nih %>% full_join(df_svnt_trim) %>% select(!time_post_exp) %>% left_join(joining_dates) %>% select(!c(nInfections, d2bleedint))
    df_wide_model_nih <- df_bio_trim_join %>% left_join(df_meta) %>% arrange(PID, time_post_exp) 

    # Make long
    df_model_nih <- df_wide_model_nih %>% covert_data_long
    dir.create(here::here("mvp", "bcell_david", "clean", paste0(study, "_", antigen_str)), recursive = TRUE, showWarnings = FALSE)
    saveRDS(df_model_nih, here::here("mvp", "bcell_david", "clean", paste0(study, "_", antigen_str), paste0("df_model_data.RDS")))

    # Convert to stan
    create_datalist_stan(df_model_nih, study, antigen_str, 3)

}