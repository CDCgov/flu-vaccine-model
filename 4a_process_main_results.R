## Set parameters used in simulations ------------------------------------------

nsim <- 1000
narrays <- 10

nfixed <- 3

# Note: this code procesess results when:
# cost scenario parameters are discrete (best, intermediate, worst)
# hd_frac and relVE_hd are continuous
# run FOR all combinations of cost-scenarios (i.e. 3*3 = 9)

if (exists("distribution")) {
  distribution <- distribution
} else {
  distribution <- TRUE
}

namestring <- ""
datestring <- "220706"

if (distribution) {
  namestring <- "_distgrid"
  datestring <- "220711"
}


## Load simulation results -----------------------------------------------------

late_list <- early_list <- parlist <- list()

for (i in 1:narrays) {
  load(paste0("results/LHS", namestring, "_array", i, "of", narrays, "_", datestring, ".RData"))
  
  list2env(all, .GlobalEnv)
  
  sim_ids <- rep(sim_ids, each = nfixed)
  
  late_list[[i]] <- late %>% mutate(nsim = sim_ids[nsim])
  early_list[[i]] <- early %>% mutate(nsim = sim_ids[nsim])
}

late <- bind_rows(late_list) %>% mutate(Season = "Low")
early <- bind_rows(early_list) %>% mutate(Season = "High")

both_seasons <- rbind(late, early) %>% rowwise() %>% 
  mutate(Costs = paste0(ndelay/7, "wk delay\n", hd_frac_miss * 100, "% decrease")) %>%
  ungroup()


## Process different burden metrics (over 65s)  --------------------------------

cases <- both_seasons %>% 
  select(Scenario, Costs, Season, nsim, time, "Discharged6") %>% 
  mutate(val = (Discharged6/chr[6])) %>%
  select( - time, - Discharged6) %>%
  group_by(nsim, Season, Costs) %>%
  mutate(val = round((val[1] - val) * popsize) )  %>%
  ungroup() %>%
  group_by(Scenario, Season, Costs) %>%
  summarize(mean  = mean(val), 
            lower = quantile(val, prob = c(0.025)), 
            upper = quantile(val, prob = c(0.975)),
            min   = min(val),
            max   = max(val) ) %>%
  ungroup()


hosps <- both_seasons %>% 
  select(Scenario, Costs, Season, nsim, time, "Discharged6") %>% 
  mutate(val = (Discharged6)) %>%
  group_by(nsim, Season, Costs) %>%
  mutate(val = round((val[1] - val) * popsize) )  %>%
  ungroup() %>%
  group_by(Scenario, Season, Costs) %>%
  summarize(mean  = mean(val), 
            lower = quantile(val, prob = c(0.025)), 
            upper = quantile(val, prob = c(0.975)),
            min   = min(val),
            max   = max(val) ) %>%
  ungroup()


deaths <- both_seasons %>% 
  select(Scenario, Costs, Season, nsim, time, "Died6") %>% 
  mutate(val = (Died6)) %>%
  group_by(nsim, Season, Costs) %>%
  mutate(val = round((val[1] - val) * popsize) )  %>%
  ungroup() %>%
  group_by(Scenario, Season, Costs) %>%
  summarize(mean  = mean(val), 
            lower = quantile(val, prob = c(0.025)), 
            upper = quantile(val, prob = c(0.975)),
            min   = min(val),
            max   = max(val)) %>%
  ungroup()