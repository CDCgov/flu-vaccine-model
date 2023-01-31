## For: ACIP analyses using LHS parameter sampling

## Run as array job on high-performance computing cluster
datestamp <- "220706"

library(flumodels)
library(flumodelsutil)
library(tidyverse)
library(lhs)


## Load equations and input parameters -----------

hd_frac <- 0.75

source("0_simulation_functions.R")
source("1_call_simulation.R")
source("2_get_inputs.R")


## Function to process LHS grid for low/high severity seasons ------------------

get_grid <- function(parmat, pars_min, pars_max, pargrid, nsim, nfixed) {
  
  parmat <- t( apply(parmat, 1, function(x) qunif(x, as.numeric(pars_min), as.numeric(pars_max)) ) )
  
  parmat <- parmat %>% data.frame() %>% 
    mutate(count = rep(nfixed * nfixed, nsim)) %>% 
    uncount(count) %>%
    mutate(ndelay       = rep(pargrid$delay, nsim),
           hd_frac_miss = rep(pargrid$miss,  nsim))
  
  totalEV <- parmat[,"hd_frac"] + parmat[,"hd_frac_gain"] + parmat[,"hd_frac_miss"] 
  
  gain_min <- 1 - 1e-9 - parmat[which(totalEV >= 1), "hd_frac"] - parmat[which(totalEV >= 1), "hd_frac_miss"]
  
  # if hd_frac_gain corrected value < 0, 
  #    hd_frac + hd_frac_miss > 1 and so hd_frac needs to be reduced by that overshoot amount
  # i.e. -ve gain values are corrected for in new hd_frac values
  hd_frac_min <- ifelse(gain_min < 0, 
                        parmat[which(totalEV >= 1), "hd_frac"] + gain_min, 
                        parmat[which(totalEV >= 1), "hd_frac"])
  
  parmat[which(totalEV >= 1), "hd_frac_gain"] <- pmax(0, gain_min)  
  parmat[which(totalEV >= 1), "hd_frac"]      <- hd_frac_min
  
  return(parmat)
}


## Get LHS sampling grid -------------------------------------------------------

# Get input arguments from .sh file
args <- as.numeric(commandArgs(trailingOnly = TRUE))

# set.seed so that the same parmat matrix is created each time
set.seed(1233456)

nsim <- args[1]
nfixed <- 3

delay     <- seq(0, 6*7, length.out = nfixed)
frac_miss <- seq(0, 0.2, length.out = nfixed)

pargrid <- expand.grid(delay = delay, miss = frac_miss)

# rVE of 5% in low  severity => 7.5% increase; rVE of 35% => 52.5% increase
# rVE of 5% in high severity => 15 % increase; rVE of 35% => 105 % increase

pars_min_low <- c(hd_frac = 0.6, hd_frac_gain = 0.0, hdVE_mult = 1.075)
pars_max_low <- c(hd_frac = 0.8, hd_frac_gain = 0.2, hdVE_mult = 1.525)

pars_min_high <- c(hd_frac = 0.6, hd_frac_gain = 0.0, hdVE_mult = 1.15)
pars_max_high <- c(hd_frac = 0.8, hd_frac_gain = 0.2, hdVE_mult = 2.05)

npars <- length(pars_min_low)

parmat <- randomLHS(nsim, npars)
colnames(parmat) <- names(pars_min_low)

# Different grids for different seasons
parmat_low  <- get_grid(parmat, pars_min_low,  pars_max_low,  pargrid, nsim, nfixed)
parmat_high <- get_grid(parmat, pars_min_high, pars_max_high, pargrid, nsim, nfixed)

parmat <- parmat_low[, c("hd_frac", "hd_frac_gain", "ndelay", "hd_frac_miss")]
parmat$hdVE_mult_low  <- parmat_low[, "hdVE_mult"]
parmat$hdVE_mult_high <- parmat_high[, "hdVE_mult"]


## Break up grid based on array ID ---------------------------------------------

tot_array <- args[2]
array_id  <- args[3]

select_rows <- rep((1:nsim %% tot_array) == (array_id - 1), each = nfixed * nfixed)

parmat  <- parmat[select_rows,]


## Simulate samples ------------------------------------------------------------

late <- early <- list()

for (i in 1:nrow(parmat)) {
  parsample  <- parmat[i,]
  
  print(parsample)
  
  ## Late/mild season
  late[[i]] <- compare_scenario(season = "Low",
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
    mutate(hd_frac = parsample["hd_frac"],
           hd_frac_miss = parsample["hd_frac_miss"],
           hd_frac_gain = parsample["hd_frac_gain"],
           hdVE_mult = parsample["hdVE_mult_low"],
           ndelay = parsample["ndelay"],
           nsim = i) %>%
    filter(time == max(time))
  
  
  ## Early/severe season
  early[[i]] <- compare_scenario(season = "High",
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
    mutate(hd_frac = parsample["hd_frac"],
           hd_frac_miss = parsample["hd_frac_miss"],
           hd_frac_gain = parsample["hd_frac_gain"],
           hdVE_mult = parsample["hdVE_mult_high"],
           ndelay = parsample["ndelay"],
           nsim = i) %>%
    filter(time == max(time))
}

late  <- bind_rows(late)
early <- bind_rows(early)

all <- list(late     = late, 
            early    = early, 
            pars     = parmat, 
            array_id = array_id, 
            sim_ids  = rep( c(1:nsim)[(1:nsim %% tot_array) == (array_id - 1)], each = nfixed) )

save(all, file = paste0("results/LHS_array", array_id, "of", tot_array, "_", datestamp, ".RData"))
