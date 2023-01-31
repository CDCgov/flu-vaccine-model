## One-way sensitivity analysis

library(flumodels)
library(flumodelsutil)
library(tidyverse)

nsim <- 10

## Load equations and fixed input parameters -----------------------------------

hd_frac <- 0.75

source("0_simulation_functions.R")
source("1_call_simulation.R")
source("2_get_inputs.R")
source("5a_sensitivity_functions.R")


## Get parameter ranges --------------------------------------------------------

# Run function
parlist_low <- get_pars(pars_old = c(hd_frac = 0.75, hd_frac_miss = 0.1,  hd_frac_gain = 0.1,  ndelay = 21, hdVE_mult = 1.225),
                        pars_min = c(hd_frac = 0.6,  hd_frac_miss = 0.0,  hd_frac_gain = 0.0,  ndelay = 0,  hdVE_mult = 1.075),
                        pars_max = c(hd_frac = 0.8,  hd_frac_miss = 0.15, hd_frac_gain = 0.15, ndelay = 42, hdVE_mult = 1.525),
                        nsim = nsim
                        )

parlist_high <- get_pars(pars_old = c(hd_frac = 0.75, hd_frac_miss = 0.1,  hd_frac_gain = 0.1,  ndelay = 21, hdVE_mult = 1.45),
                         pars_min = c(hd_frac = 0.6,  hd_frac_miss = 0.0,  hd_frac_gain = 0.0,  ndelay = 0,  hdVE_mult = 1.15),
                         pars_max = c(hd_frac = 0.8,  hd_frac_miss = 0.15, hd_frac_gain = 0.15, ndelay = 42, hdVE_mult = 2.05),
                         nsim = nsim
                        )


## Late/mild season -------------------------------------------

late_list <- list()

j <- 1

n_pars <- length(parlist_low)

for (p in 1:n_pars) {
  parmat <- parlist_low[[p]]
  parmat <- check_pars(parmat)
  
  for (i in 1:nsim) {
    
    parsample <- parmat[i,]
    
    late_list[[j]] <- get_comparison(season = "Low",
                                       ndelay = as.numeric(parsample["ndelay"]),
                                       pars = pars, diffpars = diffpars,
                                       diff_fracs = diff_fracs,
                                       chr = chr, hfr = hfr,
                                       fracSympt_by_age = fracSympt_by_age,
                                       hd_frac = as.numeric(parsample["hd_frac"]),
                                       hd_frac_miss = as.numeric(parsample["hd_frac_miss"]),
                                       hd_frac_gain = as.numeric(parsample["hd_frac_gain"]),
                                       hdVE_mult = as.numeric(parsample["hdVE_mult"]),
                                       rel_VE_by_age = rel_VE_by_age, rel_hdVE_by_age = rel_hdVE_by_age
    ) %>%
      mutate(hd_frac = parsample["hd_frac"],
             hd_frac_miss = parsample["hd_frac_miss"],
             hd_frac_gain = parsample["hd_frac_gain"],
             hdVE_mult = parsample["hdVE_mult"],
             ndelay = parsample["ndelay"],
             nsim = j)
    
    j <- j + 1
  }
}

late_all <- bind_rows(late_list)


## High severity --------------------

early_list <- list()

k <- 1

for (p in 1:n_pars) {
  parmat <- parlist_high[[p]]
  parmat <- check_pars(parmat)
  
  for (i in 1:nsim) {
    
    parsample <- parmat[i,]
    
    early_list[[k]] <- get_comparison(season = "High",
                                        ndelay = as.numeric(parsample["ndelay"]),
                                        pars = pars, diffpars = diffpars,
                                        diff_fracs = diff_fracs,
                                        chr = chr, hfr = hfr,
                                        fracSympt_by_age = fracSympt_by_age,
                                        hd_frac = as.numeric(parsample["hd_frac"]),
                                        hd_frac_miss = as.numeric(parsample["hd_frac_miss"]),
                                        hd_frac_gain = as.numeric(parsample["hd_frac_gain"]),
                                        hdVE_mult = as.numeric(parsample["hdVE_mult"]),
                                        rel_VE_by_age = rel_VE_by_age, rel_hdVE_by_age = rel_hdVE_by_age
    ) %>%
      mutate(hd_frac = parsample["hd_frac"],
             hd_frac_miss = parsample["hd_frac_miss"],
             hd_frac_gain = parsample["hd_frac_gain"],
             hdVE_mult = parsample["hdVE_mult"],
             ndelay = parsample["ndelay"],
             nsim = k)
    
    k <- k + 1
  }
}

early_all <- bind_rows(early_list)


## Save output -----------------------------------------------------------------

str_pars <- colnames(parlist_low[[1]])

lhs_labels <- c(hd_frac = "Baseline EV uptake", 
                hd_frac_miss = "Reduction in overall coverage",
                hd_frac_gain = "Increase in EV uptake",
                ndelay = "Delay in additional EVs",
                hdVE_mult = "Relative VE of EVs")

late_pars <- parlist_low %>% map_df(as_tibble) %>% 
  mutate(nsim = 1:nrow(.),
         which_var = rep(str_pars, each = !!(nsim)),
         which_var_labels = lhs_labels[which_var]) %>%
  mutate(season = "Low")

early_pars <- parlist_high %>% map_df(as_tibble) %>% 
  mutate(nsim = 1:nrow(.),
         which_var = rep(str_pars, each = !!(nsim)),
         which_var_labels = lhs_labels[which_var] ) %>%
  mutate(season = "High")

osa <- list(late = late_all, early = early_all, pars = rbind(late_pars,early_pars) )

#save(osa, file = "results/osa.RData")


