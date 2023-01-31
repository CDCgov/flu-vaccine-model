## Two-way sensitivity analyses: cost-scenario parameters

## Load packages ---------------------------------------------------------------

library(flumodels)
library(flumodelsutil)
library(tidyverse)


## Load equations and fixed input parameters -----------------------------------

hd_frac <- 0.75

source("0_simulation_functions.R")
source("1_call_simulation.R")
source("2_get_inputs.R")
source("5a_sensitivity_functions.R")


## Get parameter vectors -------------------------------------------------------

# Values for main text
 pars_old <- c(hd_frac = 0.75, hd_frac_gain = 0.10, hdVE_mult_low = 1.225, hdVE_mult_high = 1.45)

# Values for supplement text
# pars_old <- c(hd_frac = 0.75, hd_frac_gain = 0.15, hdVE_mult_low = 1.225, hdVE_mult_high = 1.45)

nsim <- 11

pars_new <- expand.grid(hd_frac_miss = seq(0, 0.1, length.out = nsim),
                        ndelay       = seq(0, 21 , length.out = nsim))


pars_new <- cbind(pars_new, t(pars_old))
pars_new <- check_pars(pars_new)


## Low severity season ---------------------------------------------------------

low_list <- list()

for (p in 1:nrow(pars_new)) {
  
  parsample <- pars_new[p,]
  
  low_list[[p]] <- get_comparison(season = "Low",
                                   ndelay = as.numeric(parsample["ndelay"]),
                                   pars = pars, diffpars = diffpars,
                                   diff_fracs = diff_fracs,
                                   chr = chr, hfr = hfr,
                                   fracSympt_by_age = fracSympt_by_age,
                                   hd_frac = as.numeric(parsample["hd_frac"]),
                                   hd_frac_miss = as.numeric(parsample["hd_frac_miss"]),
                                   hd_frac_gain = as.numeric(parsample["hd_frac_gain"]),
                                   hdVE_mult = as.numeric(parsample["hdVE_mult_low"]),
                                   rel_VE_by_age = rel_VE_by_age, rel_hdVE_by_age = rel_hdVE_by_age
  ) %>%
    mutate(hd_frac = as.numeric(parsample["hd_frac"]),
           hd_frac_miss = as.numeric(parsample["hd_frac_miss"]),
           hd_frac_gain = as.numeric(parsample["hd_frac_gain"]),
           hdVE_mult = as.numeric(parsample["hdVE_mult_low"]),
           ndelay = as.numeric(parsample["ndelay"]),
           nsim = p)
}

low_all <- bind_rows(low_list)


## High severity season --------------------------------------------------------

high_list <- list()

for (p in 1:nrow(pars_new)) {
  
  parsample <- pars_new[p,]
  
  high_list[[p]] <- get_comparison(season = "High",
                                    ndelay = as.numeric(parsample["ndelay"]),
                                    pars = pars, diffpars = diffpars,
                                    diff_fracs = diff_fracs,
                                    chr = chr, hfr = hfr,
                                    fracSympt_by_age = fracSympt_by_age,
                                    hd_frac = as.numeric(parsample["hd_frac"]),
                                    hd_frac_miss = as.numeric(parsample["hd_frac_miss"]),
                                    hd_frac_gain = as.numeric(parsample["hd_frac_gain"]),
                                    hdVE_mult = as.numeric(parsample["hdVE_mult_high"]),
                                    rel_VE_by_age = rel_VE_by_age, rel_hdVE_by_age = rel_hdVE_by_age
  ) %>%
    mutate(hd_frac = as.numeric(parsample["hd_frac"]),
           hd_frac_miss = as.numeric(parsample["hd_frac_miss"]),
           hd_frac_gain = as.numeric(parsample["hd_frac_gain"]),
           hdVE_mult = as.numeric(parsample["hdVE_mult_high"]),
           ndelay = as.numeric(parsample["ndelay"]),
           nsim = p)
}

high_all <- bind_rows(high_list)


all <- rbind(low_all, high_all)


## Calculate change in burden --------------------------------------------------

net_cases <- all %>% 
  filter(time == max(time)) %>%
  select(season = Scenario, scenario, time, val = Discharged6, hd_frac_miss, ndelay) %>% 
  mutate(val = (val/chr[6])) %>%
  group_by(`hd_frac_miss`, `ndelay`, season) %>%
  mutate(val = round((val[1] - val) * popsize) )  %>%
  summarize(val = sum(val)) %>%
  ungroup() 


net_hosps <- all %>% 
  filter(time == max(time)) %>%
  select(season = Scenario, scenario, time, val = Discharged6, hd_frac_miss, ndelay) %>% 
  group_by(hd_frac_miss, ndelay, season) %>%
  mutate(val = round((val[1] - val) * popsize) )  %>%
  summarize(val = sum(val)) %>%
  ungroup()


net_deaths <- all %>% 
  filter(time == max(time)) %>%
  select(season = Scenario, scenario, time, val = Died6, hd_frac_miss, ndelay) %>% 
  group_by(hd_frac_miss, ndelay, season) %>%
  mutate(val = round((val[1] - val) * popsize) ) %>%
  summarize(val = sum(val)) %>%
  ungroup() 

burden_tsa <- list(cases_tsa = net_cases, hosps_tsa = net_hosps, deaths_tsa = net_deaths)

save(file = paste0("results/tsa_main.RData"), burden_tsa)

#save(file = paste0("results/tsa_supplement.RData"), burden_tsa)

