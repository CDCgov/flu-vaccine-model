library(flumodels)
library(flumodelsutil)
library(tidyverse)


## Load equations and input parameters -----------

hd_frac <- 0.75

source("0_simulation_functions.R")
source("1_call_simulation.R")
source("2_get_inputs.R")


## Load extra data -------------------------------------------------------------

# cumulative ILI for Fig 1B comparison ('ili_cum')
load("data/ILInet.RData")


## Baseline model --------------------------------------------------------------

low0 <- simulate_model(season = "Low", 
                       pars = pars, diffpars = diffpars,
                       chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                       vaccineMultiplier = vaccineMultiplier, 
                       hd_frac = hd_frac, hd_frac_miss = 0,
                       daily_coverage = daily_coverage, 
                       hd_daily_coverage = hd_daily_coverage,
                       rel_VE_by_age = rel_VE_by_age, 
                       rel_hdVE_by_age = rel_hdVE_by_age)

high0 <- simulate_model(season = "High", 
                        pars = pars, diffpars = diffpars,
                        chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                        vaccineMultiplier = vaccineMultiplier, 
                        hd_frac = hd_frac, hd_frac_miss = 0,
                        daily_coverage = daily_coverage, 
                        hd_daily_coverage = hd_daily_coverage,
                        rel_VE_by_age = rel_VE_by_age, 
                        rel_hdVE_by_age = rel_hdVE_by_age)

all0 <- rbind(low0, high0)
