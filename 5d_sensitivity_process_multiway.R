library(tidyverse)
library(broom)
library(patchwork)
library(sensitivity)

datestring <- "_220708"

narrays <- 10


## Load parameter samples -------------------------------------------------

load(paste0("results/LHS_prcc_parmat", datestring, ".RData"))

parmat_sim <- data.frame(parmat) %>% mutate(nsim = 1:nrow(.)) 


## Load simulation results -----------------------------------------------------

late_list <- early_list <- parlist <- list()

for (i in 1:narrays) {
  load(paste0("results/LHS_prcc_array", i, "of", narrays, datestring, ".RData"))
  
  list2env(all, .GlobalEnv)
  
  late_list[[i]] <- late %>% mutate(nsim = sim_ids[nsim])
  early_list[[i]] <- early %>% mutate(nsim = sim_ids[nsim])
}

late <- bind_rows(late_list) %>% mutate(Season = "Low")
early <- bind_rows(early_list) %>% mutate(Season = "High")

both_seasons <- rbind(late, early)


## Calculate summary metric ----------------------------------------------------

net_change <- both_seasons %>% 
  select(Scenario, Season, nsim, time, "Discharged6") %>% 
  mutate(val = Discharged6) %>%
  select( - time, - Discharged6) %>%
  group_by(nsim, Season) %>%
  mutate(val = val[1] - val)  %>%
  summarize(val = sum(val)) %>%
  ungroup() %>% 
  left_join(parmat_sim) %>%
  mutate(hdVE_mult = ifelse(Season == "Low", hdVE_mult_low, hdVE_mult_high)) %>%
  select(-hdVE_mult_low, -hdVE_mult_high)


## PRCC ------------------------------------------------------------------------

get_prcc <- function(output, season, nboot = 100, rank = TRUE){
  
  tmp <- output %>% filter(Season == season)
  
  pcc(X = tmp[, c(colnames(parmat)[1:4], "hdVE_mult")], 
      y = tmp$val, 
      nboot = nboot, rank = rank)$PRCC %>%
    mutate(metric = "Change",
           param = rownames(.),
           Season = season)
}

prc_high <- get_prcc(net_change, "High")
prc_low <- get_prcc(net_change, "Low")

