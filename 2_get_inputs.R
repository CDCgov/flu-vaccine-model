## Setting up the baseline model: input parameters


getPopulationFractions <- function(ages,
                                   year = 2020,
                                   population = flumodels_data$population.US) {
  
  # Sinead edit: focus on one month (summing over multiple will inflate #s)
  population <- filter(population, AGE!=999, YEAR == year, MONTH == 1) %>% 
    select("AGE","TOT_POP")
  
  # Original code
  if(nrow(population) == 0)
    stop("No population for given year")
  
  if(is.unsorted(ages))
    stop("Ages must be increasing order")
  
  if(sum(ages - round(ages)) != 0  | min(ages) < 1)
    stop("Ages must be positive integers")
  
  maxAge <- max(population[,"AGE"])
  
  if (max(ages) >= maxAge)
    stop("Final age bracket includes no people. The final element in ages is the starting point of the final age bin. There is no need to include an endpoint representing the max age in the population.")
  
  ages <- append(ages, maxAge + 1)
  
  population.fractions <- vapply(ages, function(x) {sum(population[population$AGE <= x, "TOT_POP"])},
                                 FUN.VALUE = numeric(1))
  population.fractions <- c(population.fractions[1],
                            population.fractions[2:length(population.fractions)] - population.fractions[1:(length(population.fractions)-1)] ) /
    sum(population[,"TOT_POP"])
  names(population.fractions) <- c(paste0("0-", ages[1]),
                                   if(length(ages) > 1) {
                                     vapply(2:length(ages), function(x) { paste0(ages[x-1]+1, ifelse(ages[x] >= maxAge, "+", paste0("-", ages[x])))}, FUN.VALUE = character(1))
                                   }
  )
  
  return(population.fractions)
}



## Natural history parameters --------------------------------------------------

pars <- data.frame(
  value = c(2, 1, 2, 0.5, 1, 0.15, 14, 0, 0, 0, 0),
  name = c("incubation", "presympt", "postsympt", "frac_sympt", 
                "asympt_trans", "prop_immune", "VE_delay",
                "VE_inf", "VE_trans", "hdVE_inf", "hdVE_trans")
)


## Demography: age groups and contacts -----------------------------------------

age_grps <- c("0-4", "5-12", "13-17", "18-49", "50-64", "65+")

load("data/contacts.RData")

list2env(contacts_list, .GlobalEnv)

popsize <- census %>% filter(AGE == 999, YEAR == 2020) %>% select(TOT_POP)
popsize <- mean(popsize$TOT_POP)

age_fracs <- getPopulationFractions(ages = c(4, 12, 17, 49, 64), 
                                    year = "2020", population = census)

contacts <- makeContactMatrix(ages = c(4, 12, 17, 49, 64),
                              originalContactMatrix = polymod,
                              originalContactMatrixAges = data.frame(age_grps_poly),
                              originalPopulationFractions = age_fracs_poly)


## Vaccine coverage ------------------------------------------------------------

# Load data from FluVaxView
load("data/fluVaxView.RData")

list2env(coverage_list, .GlobalEnv)

# For converting monthly coverage to daily coverage
monthdays <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31)
monthnames <- c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "June")


diff_fracs <- diff_ages %>% 
  mutate(age_grp = match(Age, unique(diff_ages$Age)),
         age_frac = age_fracs[age_grp]) %>%
  arrange(Day, Season)

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

# cumulative coverage (all ages) for Fig 1B comparison 
vax_cum <- coverage_list$all_ages

## Case Hospitalization Ratio and Hospitalization Fatality Ratio ---------------

load("data/hosp_rates.RData")

list2env(burden, .GlobalEnv)

chr <- burden_age$CHR[c(1, 2, 2, 3, 4, 5)]
hfr <- burden_age$HFR[c(1, 2, 2, 3, 4, 5)]

# frac symptoms by age (Cohen et al 2021, Lancet GH Table 3): 
fracSympt_by_age <- c(80, 50, 60, 40, 50, 55)/100


## Transmission seasonality ----------------------------------------------------

R0_amp <- 0.35                
base_day <- 315      

R0pars <- data.frame(base = c(1.15, 1.25), peak = c(160, 120), scenario = c("Low", "High"))



## Vaccine effectiveness -------------------------------------------------------

load("data/VEnetwork.RData")

VEp_severe <- 0.35
VEp_mild <- 0.5

VEtable <- VE_list$VEtable %>% filter(!(Season %in% c("2010-11", "2009-10", "2008-09")))

VE_by_age <- VE_list$VEage$VE[c(1, 1, 2:5)]
VE_by_age[2] <- mean(c(VE_by_age[1], VE_by_age[3]))

VE <- mean(VEtable$VE)
rel_VE_by_age <- rel_hdVE_by_age <- VE_by_age/ VE

rel_VE_early <- c(0.45, 0.4, 0.35, 0.3, 0.35, 0.25)/VEp_severe
rel_VE_late <- c(0.6, 0.55, 0.5, 0.45, 0.5, 0.4)/VEp_mild


## Relative VE -----------------------------------------------------------------

# rVE of 15% => absolute HD VE of 49.0% in low  severity => 22.5% increase
# rVE of 15% => absolute HD VE of 36.3% in high severity => 45% increase

hdVE_mult <- c(1.225, 1.45)

# make data frame with all parameters that vary during high/low severity seasons
diffpars <- R0pars %>% mutate(VE = c(VEp_mild, VEp_severe),
                              hdVE_mult = hdVE_mult)


