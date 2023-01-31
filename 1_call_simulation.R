## Functions to simulate the model (using all the 
## helper & equation functions defined in '0_simulation_functions.R')

## Function (1) to simulate one model run ------------------------------------------

simulate_model <- function(season, ndelay = 0,
                           pars, diffpars,
                           chr, hfr, fracSympt_by_age,
                           R0_amp = 0.35, base_day = 315,
                           vaccineMultiplier, hd_vaccineMultiplier = NULL,
                           hd_frac, hd_frac_miss = 0,
                           daily_coverage, hd_daily_coverage,
                           rel_VE_by_age, rel_hdVE_by_age) {

  if(is.null(hd_vaccineMultiplier)) {
    hd_vaccineMultiplier <- vaccineMultiplier
  }
  
  hdVE_mult <- diffpars$hdVE_mult[diffpars$scenario == season]
  
  VEp <- plyr::round_any(rel_VE_by_age * diffpars$VE[diffpars$scenario == season], 0.05)

  SEAIRV_2vaccines(
    # transmission parameters
    R0 = diffpars$base[diffpars$scenario == season], 
    Seasonality = TRUE,
    offset = diffpars$peak[diffpars$scenario == season],
    R0_amp = R0_amp,
    dayOfR0 = base_day,
    
    # natural history parameters
    latentPeriod = pars$value[pars$name == "incubation"], 
    infectiousPeriod = pars$value[pars$name == "postsympt"], 
    fractionLatentThatIsInfectious = 
      pars$value[pars$name == "presympt"] / 
      pars$value[pars$name == "incubation"],
    relativeInfectivityAsymptomatic = 
      pars$value[pars$name == "asympt_trans"],
    fractionSymptomatic = fracSympt_by_age,
    
    # initial conditions parameters
    seedInfections = 500,
    seedStartDay = 0,
    priorImmunity = pars$value[pars$name == "prop_immune"],
    simulationLength = 365,
    population = popsize,
    tolerance = 1e-10,
    
    # demography
    populationFractions = age_fracs,
    contactMatrix = contacts,
    
    # ignore antiviral treatment parameters
    fractionDiagnosedAndPrescribedOutpatient = 0,
    
    # vaccination parameters
    vaccineAvailabilityByDay = daily_coverage, 
    hd_vaccineAvailabilityByDay = hd_daily_coverage, 
    vaccineAdministrationRatePerDay = max(daily_coverage) * 10, 
    vaccineUptakeMultipliers = lapply(vaccineMultiplier, 
                                      function(x)  {
                                        x$daily_mean * c(1, 1, 1, 1, 1, 1 - hd_frac - hd_frac_miss)
                                      }),
    vaccineUptakeMultiplierStartDay = c(0, cumsum(monthdays[-length(monthdays)])) + c(0, rep(1, 10)),
    
    hd_vaccineUptakeMultipliers = lapply(hd_vaccineMultiplier, 
                                         function(x) {
                                           x$daily_mean * c(0, 0, 0, 0, 0, hd_frac)
                                         }),
    hd_vaccineUptakeMultiplierStartDay = c(0, cumsum(monthdays[-length(monthdays)])) + c(0, rep(1 + ndelay, 10)),
    
    VEs = VEp * pars$value[pars$name == "VE_inf"],
    VEi = VEp * pars$value[pars$name == "VE_trans"],
    VEp = VEp,
    
    hd_VEs = VEp * hdVE_mult * pars$value[pars$name == "VE_inf"],
    hd_VEi = VEp * hdVE_mult * pars$value[pars$name == "VE_trans"],
    hd_VEp = VEp * hdVE_mult,
    vaccineEfficacyDelay = pars$value[pars$name == "VE_delay"],
    
    # clinical severity parameters
    hospDuration = c(4, 2, 2, 2, 7, 14),
    caseHospitalizationRatio = chr,
    DelayOnsetToHosp = 3,
    hospFatalityRatio = hfr
  )$rawOutput %>% 
    data.frame() %>%
    mutate(Scenario = season)
}




## Function (2) to compare one trade-off scenario w/ the baseline model ------------

compare_scenario <- function(season, ndelay, 
                             pars, diffpars, diff_fracs,
                             chr, hfr, fracSympt_by_age,
                             R0_amp = 0.35, base_day = 315,
                             hd_frac, hd_frac_miss, hd_frac_gain,
                             hdVE_mult,
                             rel_VE_by_age, rel_hdVE_by_age) {

  diffpars$hdVE_mult <- hdVE_mult
  
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
    mutate(Scenario = "Baseline")
  
  
  ## Set up alternative scenario -------------------
  
  hd_frac_new <- hd_frac + hd_frac_gain
  
  ## Change in HD vaccinees
  
  pop_coverage_extra <- diff_fracs %>% 
    distinct(Day, Age, age_grp, age_frac, mean) %>% 
    mutate(hd_meanpop = mean * age_frac * c(0, 0, 0, 0, 0, hd_frac_gain) /100) %>% 
    group_by(Day) %>% 
    summarize(hd_meanpop = sum(hd_meanpop)) %>% ungroup()
  
  hd_daily_coverage_extra <- rep(pop_coverage_extra$hd_meanpop/monthdays, times = monthdays) * popsize
  
  hd_daily_coverage_alt <- c(hd_daily_coverage, rep(0, ndelay)) + c(rep(0, ndelay), hd_daily_coverage_extra)
  
  
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
  alternative <- simulate_model(season = season, 
                                pars = pars, diffpars = diffpars,
                                chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                                vaccineMultiplier = vaccineMultiplier, 
                                hd_frac = hd_frac_new, 
                                hd_frac_miss = hd_frac_miss,
                                daily_coverage = daily_coverage_alt, 
                                hd_daily_coverage = hd_daily_coverage_alt,
                                rel_VE_by_age = rel_VE_by_age, 
                                rel_hdVE_by_age = rel_hdVE_by_age) %>%
    mutate(Scenario = "Alternative")
  
  joined <- rbind(baseline, alternative)
}

