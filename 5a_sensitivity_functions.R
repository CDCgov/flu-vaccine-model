## Functions to run sensitivity analyses

## Functions -------------------------------------------------------------------

## Get parameter ranges for simulation ----
get_pars <- function(pars_old, nsim = 10, pars_min, pars_max) {
  
  n_pars <- length(pars_old)
  str_pars <- names(pars_old)
  
  parlist <- lapply(str_pars, 
                    function(x) {
                      
                      tmp <- cbind( seq(pars_min[x], pars_max[x], 
                                        length.out = nsim),
                                    
                                    matrix(pars_old[str_pars != x], 
                                           ncol = n_pars - 1, nrow = nsim, 
                                           byrow = T)
                      )
                      colnames(tmp) <- c(x, str_pars[str_pars != x])
                      tmp
                    }
  )
  return(parlist)
}


## Check vaccine coverage proportions dont exceed 1 ----
check_pars <- function(parmat) {
  
  totalEV <- parmat[,"hd_frac"] + parmat[,"hd_frac_gain"] + parmat[,"hd_frac_miss"] 
  
  gain_min <- 1 - 1e-9 - parmat[which(totalEV >= 1), "hd_frac"] - parmat[which(totalEV >= 1), "hd_frac_miss"]
  
  # if hd_frac_gain corrected value < 0, 
  #    hd_frac + hd_frac_miss > 1, so hd_frac needs to be reduced by the overshoot amount
  # i.e. -ve gain values are corrected for in new hd_frac values
  hd_frac_min <- ifelse(gain_min < 0, 
                        parmat[which(totalEV >= 1), "hd_frac"] + gain_min, 
                        parmat[which(totalEV >= 1), "hd_frac"])
  
  parmat[which(totalEV >= 1), "hd_frac_gain"] <- pmax(0, gain_min)  
  parmat[which(totalEV >= 1), "hd_frac"]      <- hd_frac_min
  
  return(parmat)
}


## Compare baseline to intervention for each set of parameters ----

get_comparison <- function(season, ndelay, ndelay_baseline = 0,
                           pars, diffpars, diff_fracs,
                           chr, hfr, fracSympt_by_age,
                           R0_amp = 0.35, base_day = 315,
                           hd_frac, hd_frac_miss, hd_frac_gain,
                           hdVE_mult,
                           rel_VE_by_age, rel_hdVE_by_age) {
  
  pop_coverage <- diff_fracs %>% 
    distinct(Day, Age, age_grp, age_frac, mean) %>% 
    mutate(meanpop = mean * age_frac * c(1, 1, 1, 1, 1, 1 - hd_frac) /100,
           hd_meanpop = mean * age_frac * c(0, 0, 0, 0, 0, hd_frac) /100) %>% 
    group_by(Day) %>% 
    summarize(meanpop = sum(meanpop),
              hd_meanpop = sum(hd_meanpop)) %>% ungroup()
  
  vaccineMultiplier <- diff_fracs %>% 
    distinct(Day, mean, Relative_mean, Age, age_frac) %>%
    mutate(daily_mean = match(Day, unique(diff_fracs$Day)),
           daily_mean = mean/monthdays[daily_mean]) %>%
    left_join(pop_coverage) %>%
    split(., f = .$Day) 
  
  
  daily_coverage <- rep(pop_coverage$meanpop/monthdays, times = monthdays) * popsize
  
  hd_daily_coverage <- rep(pop_coverage$hd_meanpop/monthdays, times = monthdays) * popsize
  hd_daily_coverage <- c(rep(0, ndelay_baseline), hd_daily_coverage)
  
  # make data frame with all parameters that vary during high/low severity seasons
  diffpars <- R0pars %>% mutate(VE = c(VEp_mild, VEp_severe),
                                hdVE_mult = hdVE_mult)
  
  ## Baseline model -------------------
  baseline <- simulate_model(season = season, 
                             pars = pars, diffpars = diffpars,
                             chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                             vaccineMultiplier = vaccineMultiplier, 
                             hd_frac = hd_frac, hd_frac_miss = 0,
                             daily_coverage = daily_coverage, 
                             hd_daily_coverage = hd_daily_coverage,
                             rel_VE_by_age = rel_VE_by_age, 
                             rel_hdVE_by_age = rel_hdVE_by_age) %>% 
    mutate(scenario = "Baseline")
  
  ## Set up alternative scenario -------------------
  
  hd_frac_new <- hd_frac + hd_frac_gain
  
  ## Change in HD vaccinees
  
  hd_vaccineMultiplier <- diff_fracs %>%
    distinct(Day, mean, Relative_mean, Age, age_frac) %>%
    mutate(daily_mean = match(Day, unique(diff_fracs$Day)),
           daily_mean = mean/monthdays[daily_mean],
           Day = Day + ndelay) %>%
    left_join(pop_coverage) %>%
    split(., f = .$Day)
  
  
  pop_coverage_extra <- diff_fracs %>% 
    distinct(Day, Age, age_grp, age_frac, mean) %>% 
    mutate(hd_meanpop = mean * age_frac * c(0, 0, 0, 0, 0, hd_frac_gain) /100) %>% 
    group_by(Day) %>% 
    summarize(hd_meanpop = sum(hd_meanpop)) %>% ungroup()
  
  hd_daily_coverage_extra <- rep(pop_coverage_extra$hd_meanpop/monthdays, times = monthdays) * popsize
  
  hd_daily_coverage_alt <- c(hd_daily_coverage, rep(0, ndelay + ndelay_baseline)) + 
    c(rep(0, ndelay + ndelay_baseline), hd_daily_coverage_extra)
  
  
  ## Change in standard vaccinees
  
  pop_coverage_alt <- diff_fracs %>% 
    distinct(Day, Age, age_grp, age_frac, mean) %>% 
    mutate(meanpop = mean * age_frac * c(1, 1, 1, 1, 1, 1 - hd_frac_new - hd_frac_miss) /100,
           hd_meanpop = mean * age_frac * c(0, 0, 0, 0, 0, hd_frac_new) /100) %>% 
    group_by(Day) %>% 
    summarize(meanpop = sum(meanpop),
              hd_meanpop = sum(hd_meanpop)) %>% ungroup()
  
  daily_coverage_alt <- rep(pop_coverage_alt$meanpop/monthdays, times = monthdays) * popsize
  
  
  ## Run alternative scenario  -------------------
  alternative <- simulate_model(season = season, #ndelay = ndelay,
                                pars = pars, diffpars = diffpars,
                                chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                                vaccineMultiplier = vaccineMultiplier, 
                                #hd_vaccineMultiplier = hd_vaccineMultiplier,
                                hd_frac = hd_frac_new, 
                                hd_frac_miss = hd_frac_miss,
                                daily_coverage = daily_coverage_alt, 
                                hd_daily_coverage = hd_daily_coverage_alt,
                                rel_VE_by_age = rel_VE_by_age, 
                                rel_hdVE_by_age = rel_hdVE_by_age) %>% 
    mutate(scenario = "Alternative")
  
  joined <- rbind(baseline, alternative)
}