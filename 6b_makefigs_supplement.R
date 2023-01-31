## Load packages and prep code -------------------------------------------------

rm(list=ls())

library(patchwork)
library(viridis)
library(scico)

source("3a_baseline_simulation.R")


## Plot settings -------------------------------------------------------

errorwidth <- 0.25
errorsize  <- 0.5

linesize  <- 1.3
linesize4 <- 1

my_lapaz    <- scico(7, palette = 'lapaz')[c(1, 4)]
my_lapaz1   <- scico(7, palette = 'lapaz')
mylapaz_rep <- scico(6, palette = 'lapaz')

textsize <- 8
mytheme8 <- theme_bw() + theme(strip.text.x = element_text(size = textsize),
                                axis.text  = element_text(size = textsize - 1),
                                axis.title   = element_text(size = textsize),
                                title        = element_text(size = textsize),
                                legend.title = element_text(size = textsize),
                                legend.text  = element_text(size = textsize),
                                strip.placement  = "outside", 
                                strip.background = element_blank())

textsize <- 10
mytheme10 <- theme_bw() + theme(strip.text.x = element_blank(),
                               axis.text.y  = element_text(size = textsize - 1),
                               axis.text.x  = element_text(size = textsize - 1),
                               axis.title   = element_text(size = textsize),
                               title        = element_text(size = textsize),
                               legend.title = element_text(size = textsize),
                               legend.text  = element_text(size = textsize),
                               plot.title   = element_text(size = textsize))


textsize <- 10
mythemeS8 <- theme_bw() + theme(strip.text.x = element_text(size = textsize),
                                axis.text.x  = element_blank(), 
                                axis.ticks.x = element_blank(),
                                axis.title   = element_text(size = textsize),
                                title        = element_text(size = textsize - 1),
                                legend.title = element_text(size = textsize),
                                legend.text  = element_text(size = textsize),
                                strip.placement  = "outside", 
                                strip.background = element_blank())

textsize <- 12
mythemeS9 <- theme_bw() + theme(strip.text.x = element_text(size = textsize - 1),
                                axis.text.y  = element_text(size = textsize),
                                axis.title   = element_text(size = textsize),
                                title        = element_text(size = textsize),
                                legend.title = element_text(size = textsize),
                                legend.text  = element_text(size = textsize),
                                strip.placement  = "outside", 
                                strip.background = element_blank(),
                                axis.text.x  = element_blank(), 
                                axis.ticks.x = element_blank(),
                                plot.title   = element_text(size = textsize)) 


## Fig S2: Contact rates -------------------------------------------------------

Cairo::CairoPDF(file = "figures/FigS2.pdf", width = 5, height = 5)

gplots::heatmap.2(contacts, Rowv   = NA, Colv  = NA, 
                  labRow = age_grps, labCol = age_grps, 
                  cexRow = 1.2, cexCol = 1.2,
                  margins = c(6, 6), 
                  density.info = "none",
                  trace = "none",
                  colsep = 1:nrow(contacts), # Add vertical grid lines
                  rowsep = 1:nrow(contacts), # Add horizontal grid lines
                  sepcolor = "black",
                  col = viridis::viridis_pal(),
                  breaks = c(0:10),
                  
                  # legend arguments
                  key.title = "", 
                  key.xlab = "Scaled contacts",
                  keysize = 2.7,
                  key.par = list(cex = 1))
dev.off()


## Fig S3: Vaccine coverage ----------------------------------------------------

FigS3 <- ggplot(diff_ages, aes(x = Day, y = mean, colour = Age)) +
  geom_line(size = 0.8) + mytheme8 +
  scale_colour_manual("Age group", labels = age_grps, values = my_lapaz1) +
  labs(y = "Average monthly coverage \n(% within age group)") +
  scale_x_continuous("", labels = monthnames, breaks = c(0, cumsum(monthdays)))

ggsave(file = "figures/FigS3.pdf", FigS3, width = 5, height = 2)


## Fig S4: R0 seasonality ------------------------------------------------------------------

days <- 1:365
freq <- 2 * pi / 365

seas <- list()

for (b in 1:nrow(R0pars)) {
  seas[[b]] <- data.frame(day = days,
                          R0 = R0pars$base[b] - (R0_amp/2) * ( cos(freq * (base_day - R0pars$peak[b])) - 
                                                                 cos(freq * (days - R0pars$peak[b])) ),
                          scenario = R0pars$scenario[b]
  )
  
  print( paste("R0 =", round(mean(seas[[b]]$R0), 1), "/ Re =",  round(mean(seas[[b]]$R0) * 0.85, 1)) )
}

seas <- bind_rows(seas)

# Get proportion immune to plot Re
Re <- all0 %>% 
  mutate(R = R1 + R2 + R3 + R4 + R5 + R6 + Rv1 + Rv2 + Rv3 + Rv4 + Rv5 + Rv6 + Rvb6) %>% 
  select(day = time, R, scenario = Scenario) %>%
  left_join(seas) %>%
  filter(day < (365 - 61))

FigS4 <- ggplot(Re, aes(x = day, y = R0 * (1 - R), colour = scenario)) + 
  geom_line(size = 0.8) + 
  scale_colour_manual("Season\nseverity", values = my_lapaz) +
  scale_x_continuous("", labels = monthnames, breaks = c(0, cumsum(monthdays))) +
  mytheme8 +
  labs(y = "Effective reproduction number")

ggsave(file = "figures/FigS4.pdf", FigS4, width = 5, height = 2)


## Fig S5: Model calibration ---------------------------------------------------

## S5a: Timing of peak incidence ----------------------------------

sympt0 <- all0 %>% 
  select(time, Scenario, contains("Discharged")) %>% 
  gather(age, val, contains("Discharged")) %>%
  mutate(age_grp = str_extract(age, pattern = "[0-9]+"),
         val = val/chr[as.numeric(age_grp)]
  ) %>%
  group_by(age_grp, Scenario) %>%
  mutate(val = val * 100000,
         I_inc = diff(c(0, val), lag = 1)) %>%
  ungroup() %>%
  group_by(Scenario, time) %>%
  summarize(I_inc = sum(I_inc)) %>% ungroup()

inc_peak <- 
  ggplot(sympt0, aes(x = time, y = I_inc, colour = Scenario)) +
  geom_rect(data = NULL, aes(xmin = 150, xmax = 270, ymin = 0, ymax = Inf, fill = "past peaks"), 
            alpha = 0.03,  colour = NA) +
  geom_line(size = linesize) + 
  scale_colour_manual("Season severity", values = my_lapaz) +
  scale_fill_manual("Shaded regions", values =  "grey85") +
  labs(y = "Symptomatic cases\nper 100,000") +
  scale_x_continuous("", labels = monthnames, breaks = c(0, cumsum(monthdays))) +
  mytheme10 +
  theme(axis.text.x = element_text(angle = 90)) + 
  guides(fill = guide_legend(override.aes = list(alpha = 0.8)))


## S5b: Vacc timing vs cumulative cases ---------------------------------------------

# vaccination
vacc2 <- all0 %>% 
  mutate(V1 = V1 + Vb1, V2 = V2 + Vb2, V3 = V3 + Vb3,
         V4 = V4 + Vb4, V5 = V5 + Vb5, V6 = V6 + Vb6) %>% 
  select(time, Scenario, V1:V6) %>%
  mutate(vaccinated = rowSums(.[-c(1:2)], na.rm = TRUE),
         percent = vaccinated/max(vaccinated) * 100,
         type = "model") %>%
  select(Scenario, type, time, percent)

# cumulative ILI
sympt2 <- sympt0 %>% 
  group_by(Scenario) %>%
  mutate(percent = cumsum(I_inc)/sum(I_inc) * 100,
         type = "model") %>%
  ungroup()

# Sims start on July 1st = epiweek 26
# ILI seasons start on epiweek 39

ili_compare <- ili_cum %>% 
  mutate(time = week_rel * 7  + 7 * (39 - 26),
         Scenario = ifelse(season == "2011-12", "Low", "High"),
         type = "data") %>%
  select(Scenario, type, time, percent) %>%
  full_join(sympt2)

vacc_timing <- 
  ggplot() + 
  geom_ribbon(data = vacc2, aes(x = time, ymax = percent, ymin = 0, fill = "vaccination"), alpha = 0.35) +
  geom_line(data = sympt2, aes(x = time, y = percent, colour = Scenario), 
            size = linesize, linetype = "solid") +
  mytheme10 + labs(y = "Percentage of total (%)")  +
  scale_x_continuous("", labels = monthnames, breaks = c(0, cumsum(monthdays))) +
  scale_colour_manual( values = my_lapaz) +
  scale_fill_manual(NULL, values =  my_lapaz1[6]) +
  theme(axis.text.x = element_text(angle = 90), 
        legend.title = element_blank(),
        legend.spacing.y = unit(-0.65, 'cm')) + 
  guides(colour = "none")



## S5c: Age distribution ----------------------------------

# CDC burden estimates
obs <- burden_season_age %>% 
  filter(Season %in% c("2012-2013", "2011-2012")) %>%
  mutate(Scenario = ifelse(Season == "2012-2013", "High\n(12/13)", "Low\n(11/12)")) 

obs_ili <- obs %>% select(age_group, percent = ili, Scenario)
obs_hosps <- obs %>% select(age_group, percent = hosps, Scenario)
obs_deaths <- obs %>% select(age_group, percent = deaths, Scenario)

# Model predictions
ili_ages <- all0 %>% 
  group_by(Scenario) %>%
  filter(time == max(time)) %>%
  ungroup() %>%
  select(Scenario, time, contains("Discharged")) %>% 
  gather(age, val, contains("Discharged")) %>%
  mutate(age_ind = as.numeric(str_extract(age, "[0-9]+")),
         val = val/chr[age_ind],
         age_group = age_grps[age_ind]) %>%
  group_by(Scenario) %>%
  mutate(percent = val/sum(val) * 100) %>%
  ungroup() %>% 
  mutate(age_group = ifelse(substr(age_group, start = 1, stop = 2) == "5-", 
                            gsub(age_group, pattern = "5", replacement = "05"),
                            age_group),
         age_group = ifelse(age_group %in% c("05-12", "05-17", "13-17"), 
                            "05-17", age_group)) %>% 
  select(age_group, percent, Scenario) %>%
  mutate(Scenario = ifelse(Scenario == "High", "High\n(model)", "Low\n(model)")) %>%
  rbind(., obs_ili) %>%
  ggplot(aes(x = Scenario, y = percent, fill = age_group)) + 
  geom_bar(stat = "identity") + mytheme10 +
  scale_fill_manual(NULL, values = mylapaz_rep) +
  labs(title = "Symptomatic cases", 
       y = "Percentage of total (%)", 
       x = "Season severity")

hosps_ages <- all0 %>% 
  group_by(Scenario) %>%
  filter(time == max(time)) %>%
  ungroup() %>%
  select(Scenario, time, contains("Discharged")) %>% 
  gather(age, val, contains("Discharged")) %>%
  mutate(age_ind = as.numeric(str_extract(age, "[0-9]+")),
         age_group = age_grps[age_ind]) %>%
  group_by(Scenario) %>%
  mutate(percent = val/sum(val) * 100) %>%
  ungroup() %>% 
  mutate(age_group = ifelse(substr(age_group, start = 1, stop = 2) == "5-", 
                            gsub(age_group, pattern = "5", replacement = "05"),
                            age_group),
         age_group = ifelse(age_group %in% c("05-12", "05-17", "13-17"), 
                            "05-17", age_group)) %>%
  select(age_group, percent, Scenario) %>%
  mutate(Scenario = ifelse(Scenario == "High", "High\n(model)", "Low\n(model)")) %>%
  rbind(., obs_hosps) %>%
  ggplot(aes(x = Scenario, y = percent, fill = age_group)) + 
  geom_bar(stat = "identity") + mytheme10 +
  scale_fill_manual(NULL, values = mylapaz_rep) +
  labs(title = "Hospitalizations", 
       y = "Percentage of total (%)",
       x = "Season severity")

deaths_ages <- all0 %>% 
  group_by(Scenario) %>%
  filter(time == max(time)) %>%
  ungroup() %>%
  select(Scenario, time, contains("Died")) %>% 
  gather(age, val, contains("Died")) %>%
  mutate(age_ind = as.numeric(str_extract(age, "[0-9]+")),
         age_group = age_grps[age_ind]) %>%
  group_by(Scenario) %>%
  mutate(percent = val/sum(val) * 100) %>%
  ungroup() %>% 
  mutate(age_group = ifelse(substr(age_group, start = 1, stop = 2) == "5-", 
                            gsub(age_group, pattern = "5", replacement = "05"),
                            age_group),
         age_group = ifelse(age_group %in% c("05-12", "05-17", "13-17"), 
                            "05-17", age_group)) %>%
  select(age_group, percent, Scenario) %>%
  mutate(Scenario = ifelse(Scenario == "High", "High\n(model)", "Low\n(model)")) %>%
  rbind(., obs_deaths) %>%
  ggplot(aes(x = Scenario, y = percent, fill = age_group)) + 
  geom_bar(stat = "identity") + mytheme10 +
  scale_fill_manual(NULL, values = mylapaz_rep) +
  labs(title = "Deaths", 
       y = "Percentage of total (%)", 
       x = "Season severity")

age_dist <- ili_ages + hosps_ages + deaths_ages + plot_layout(guide = "collect", nrow = 1)


## S5: Main plot ---------------------------------------------------------------

FigS5 <- (inc_peak + vacc_timing + plot_layout(widths = c(1, 1), guides = "collect") ) / 
  age_dist + plot_annotation(tag_levels = c("A")) + 
  plot_layout(heights = c(0.8, 1)) &
  theme(plot.tag = element_text(size = 12))

ggsave("figures/FigS5.pdf", FigS5, width = 11, height = 6)


## Fig S6 & 7 simulations: hd_frac, hd_mult ranges & baseline ------------------

pargrid <- expand.grid( hd_frac   = c(0.6, 0.75, 0.8),
                        hdVE_mult = c(1, 2, 3) ) 

pargrid$hdVE_low  <- c(1.075, 1.225, 1.525)[pargrid$hdVE_mult]
pargrid$hdVE_high <- c(1.150, 1.450, 2.050)[pargrid$hdVE_mult]

low_list <- high_list <- list()

for (k in 1:nrow(pargrid) ) {
  diffpars$hdVE_mult <- c(pargrid$hdVE_low[k], pargrid$hdVE_high[k])
    
  low_list[[k]] <- simulate_model(season = "Low", 
                                  pars = pars, diffpars = diffpars,
                                  chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                                  vaccineMultiplier = vaccineMultiplier, 
                                  hd_frac = pargrid$hd_frac[k], 
                                  hd_frac_miss = 0,
                                  daily_coverage = daily_coverage, 
                                  hd_daily_coverage = hd_daily_coverage,
                                  rel_VE_by_age = rel_VE_by_age, 
                                  rel_hdVE_by_age = rel_hdVE_by_age) %>%
    mutate(sim = k, 
           hd_frac = pargrid$hd_frac[k], 
           hdVE_mult = pargrid$hdVE_mult[k])
  
  high_list[[k]] <- simulate_model(season = "High", 
                                   pars = pars, diffpars = diffpars,
                                   chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                                   vaccineMultiplier = vaccineMultiplier, 
                                   hd_frac = pargrid$hd_frac[k], 
                                   hd_frac_miss = 0,
                                   daily_coverage = daily_coverage, 
                                   hd_daily_coverage = hd_daily_coverage,
                                   rel_VE_by_age = rel_VE_by_age, 
                                   rel_hdVE_by_age = rel_hdVE_by_age) %>%
    mutate(sim = k, 
           hd_frac = pargrid$hd_frac[k], 
           hdVE_mult = pargrid$hdVE_mult[k])
}

high_range <- bind_rows(high_list)
low_range <- bind_rows(low_list)

all_range <- rbind(low_range, high_range)


## S6A: Incidence timing plot --------------------------------------------------

rVE_vals <- c(5, 15, 35)

simlabels <- paste0("Prop. HDAV = ", pargrid$hd_frac, ", rVE = ", rVE_vals[pargrid$hdVE_mult], "%")

sympt_range <- all_range %>% 
  select(time, Scenario, hdVE_mult, hd_frac, sim, contains("Discharged")) %>% 
  gather(age, val, contains("Discharged")) %>%
  mutate(age_grp = str_extract(age, pattern = "[0-9]+"),
         val = val/chr[as.numeric(age_grp)]
  ) %>%
  group_by(age_grp, Scenario, hdVE_mult, hd_frac, sim) %>%
  mutate(val = val * 100000,
         I_inc = diff(c(0, val), lag = 1)) %>%
  ungroup() %>%
  group_by(Scenario, time, hdVE_mult, hd_frac, sim) %>%
  summarize(I_inc = sum(I_inc)) %>% ungroup() 

inc_range <- 
  ggplot(sympt_range, aes(x = time, y = I_inc, colour = as.factor(sim))) +
  geom_rect(data = NULL, aes(xmin = 150, xmax = 270, ymin = 0, ymax = Inf), 
            alpha = 0.2, fill = "grey90", colour = NA) +
  geom_line(size = linesize4) + 
  facet_wrap(~ Scenario) +  
  scale_colour_viridis(NULL, discrete = TRUE, direction = -1, labels = simlabels) +
  labs(y = "Symptomatic cases\nper 100,000") +
  scale_x_continuous("", labels = monthnames, breaks = c(0, cumsum(monthdays))) +
  mytheme8 + 
  guides(colour = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "top")


## S6B: Vaccine timing plot ----------------------------------------------------

# vaccination
vacc_range2 <- all_range %>% 
  mutate(V1 = V1 + Vb1, V2 = V2 + Vb2, V3 = V3 + Vb3,
         V4 = V4 + Vb4, V5 = V5 + Vb5, V6 = V6 + Vb6) %>% 
  select(time, Scenario, hdVE_mult, hd_frac, sim, V1:V6) %>%
  mutate(vaccinated = rowSums(.[-c(1:5)], na.rm = TRUE),
         percent = vaccinated/max(vaccinated) * 100,
         type = "model") %>%
  select(Scenario, type, hdVE_mult, hd_frac, sim, time, percent)

# cumulative ILI
sympt_range2 <- sympt_range %>% 
  group_by(Scenario, hdVE_mult, hd_frac, sim) %>%
  mutate(percent = cumsum(I_inc)/sum(I_inc) * 100,
         type = "model") %>%
  ungroup()


sim_labeller <- function(k) {
  paste0("Prop. HDAV = "       , pargrid$hd_frac  [as.numeric(k)] ,
         ", rVE = ", rVE_vals[pargrid$hdVE_mult[as.numeric(k)]], "%")
}

vacc_range <- 
  ggplot() + 
  geom_ribbon(data = vacc_range2, aes(x = time, ymax = percent, ymin = 0), fill = "grey90") +
  geom_line(data = sympt_range2, aes(x = time, y = percent, colour = Scenario), 
            size = linesize4, linetype = "solid") +
  mytheme8 + labs(y = "Percentage of total (%)") +
  facet_wrap(~ sim, labeller = labeller(sim = sim_labeller)) +
  scale_x_continuous("", labels = monthnames, breaks = c(0, cumsum(monthdays))) +
  scale_colour_manual("Season\nseverity", values = my_lapaz) +
  theme(axis.text.x = element_text(angle = 90))


## S6: Main plot ---------------------------------------------------------------

FigS6 <- inc_range / vacc_range + 
  plot_layout(heights = c(0.4, 1)) +
  plot_annotation(tag_levels = c("A")) + 
  theme(plot.tag = element_text(size = 12))

ggsave(file = "figures/FigS6.pdf", FigS6, width = 7, height = 7)


## Table S3/Fig S7: calculate burden w/ baseline uncertainty -----------------------------

# Cases
cases <- sympt_range %>% 
  group_by(Scenario, hdVE_mult, hd_frac, sim) %>% 
  summarize(Total_ILI = sum(I_inc) / 100000 * popsize) %>%
  ungroup() %>%
  select(Scenario, hdVE_mult, hd_frac, val = Total_ILI)%>%
  mutate(burden = "Symptomatic cases")

# Hospitalizations
hosps <- all_range %>% 
  group_by(Scenario, hdVE_mult, hd_frac) %>%
  filter(time == max(time)) %>% 
  ungroup() %>%
  select(Scenario, hdVE_mult, hd_frac, contains("Discharged")) %>%
  mutate(val = rowSums(.[,-c(1:3)], na.rm = TRUE) * popsize)  %>%
  select(Scenario, hdVE_mult, hd_frac, val)%>%
  mutate(burden = "Hospitalizations")

# Deaths
deaths <- all_range %>% 
  group_by(Scenario, hdVE_mult, hd_frac) %>%
  filter(time == max(time)) %>% 
  ungroup() %>%
  select(Scenario, hdVE_mult, hd_frac, contains("Died")) %>%
  mutate(val = rowSums(.[,-c(1:3)], na.rm = TRUE) * popsize) %>%
  select(Scenario, hdVE_mult, hd_frac, val) %>%
  mutate(burden = "Deaths")

burden_range <- rbind(cases, hosps, deaths)


## Fig S7: Main plot -----------------------------------------------------------

# Plot
plot_func <- function(df, name = name) {
  ggplot(data = df, aes(x = as.factor(hdVE_mult), y = as.factor(hd_frac), fill = val)) +
    geom_tile() + 
    scale_fill_scico(name, labels = scales::comma, palette = "lapaz") +
    mytheme8 + 
    labs(x = "Relative VE of HDAVs", y = "Proportion receiving HDAV") +
    scale_x_discrete(labels = c("5%", "15%", "35%"))
}

nested_high <- burden_range %>% 
  filter(Scenario == "High") %>%
  group_by(burden) %>% 
  nest() %>% 
  mutate(plots = map2(data, burden, plot_func)) 

high_range <- gridExtra::grid.arrange(grobs = nested_high$plots)

nested_low <- burden_range %>% 
  filter(Scenario == "Low") %>%
  group_by(burden) %>% 
  nest() %>% 
  mutate(plots = map2(data, burden, plot_func)) 

low_range <- gridExtra::grid.arrange(grobs = nested_low$plots)


FigS7 <- ggpubr::ggarrange(high_range, low_range, nrow = 1, labels = c("A", "B"), font.label = list(size = 12, face = "plain"))

ggsave(file = "figures/FigS7.pdf", FigS7, width = 7, height = 6.5)


## Fig S8: 3x3 symptomatic cases and deaths averted ---------------

# FALSE for uniform rVE
distribution <- FALSE

source("4a_process_main_results.R")

p_cases <- cases %>% filter(Scenario != "Baseline") %>%
  ggplot(aes(x = Scenario, y = mean, fill = Season)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  mythemeS8 + facet_wrap(~Costs, ncol = 3) +
  labs(title = "Symptomatic cases", y = "Number averted", x = "") + 
  scale_fill_manual("Season\nseverity", values = my_lapaz) +
  scale_y_continuous(labels = scales::comma)

p_deaths <- deaths %>% filter(Scenario != "Baseline") %>%
  ggplot(aes(x = Scenario, y = mean, fill = Season)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  mythemeS8 + facet_wrap(~Costs, ncol = 3) +
  labs(title = "Deaths", y = "Number averted", x = "") + 
  scale_fill_manual("Season\nseverity", values = my_lapaz) +
  scale_y_continuous(labels = scales::comma)


FigS8 <- p_cases + p_deaths + 
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A")


ggsave("figures/FigS8.pdf", FigS8, width = 9.5, height = 6)


## Fig S9: Indirect effects (10%) ----------------------------------------------------

# FALSE for uniform rVE
distribution <- FALSE

# Main case: no indirect effects
source("4a_process_main_results.R") #-- used for FigS8, don't need to run again 

cases0  <- cases  %>% mutate(burden = "cases",            indirect = "none")
hosps0  <- hosps  %>% mutate(burden = "hospitalizations", indirect = "none")
deaths0 <- deaths %>% mutate(burden = "deaths",           indirect = "none")

burden0 <- both_seasons %>% 
  select(Scenario, Costs, Season, nsim, time, "Died6") %>% 
  mutate(val = (Died6)) %>%
  group_by(nsim, Season, Costs) %>%
  mutate(val = (val[1] - val)/val[1] * 100 )  %>%
  ungroup() %>%
  group_by(Scenario, Season, Costs) %>%
  summarize(mean  = mean(val), 
            lower = quantile(val, prob = c(0.025)), 
            upper = quantile(val, prob = c(0.975))) %>%
  ungroup() %>% mutate(indirect = "none")


# Indirect protection against onward transmission
indirect <- "trans"
datestring <- "_220707"

source("4b_process_indirect_results.R")

cases_trans  <- cases  %>% mutate(burden = "cases",            indirect = "transmission")
hosps_trans  <- hosps  %>% mutate(burden = "hospitalizations", indirect = "transmission")
deaths_trans <- deaths %>% mutate(burden = "deaths",           indirect = "transmission")

burden_trans <- both_seasons %>% 
  select(Scenario, Costs, Season, nsim, time, "Died6") %>% 
  mutate(val = (Died6)) %>%
  group_by(nsim, Season, Costs) %>%
  mutate(val = (val[1] - val)/val[1] * 100 )  %>%
  ungroup() %>%
  group_by(Scenario, Season, Costs) %>%
  summarize(mean  = mean(val), 
            lower = quantile(val, prob = c(0.025)), 
            upper = quantile(val, prob = c(0.975))) %>%
  ungroup() %>% mutate(indirect = "transmission")


# Indirect protection against infection
indirect <- "inf"
source("4b_process_indirect_results.R")

cases_inf  <- cases  %>% mutate(burden = "cases",            indirect = "infection")
hosps_inf  <- hosps  %>% mutate(burden = "hospitalizations", indirect = "infection")
deaths_inf <- deaths %>% mutate(burden = "deaths",           indirect = "infection")

burden_inf <- both_seasons %>% 
  select(Scenario, Costs, Season, nsim, time, "Died6") %>% 
  mutate(val = (Died6)) %>%
  group_by(nsim, Season, Costs) %>%
  mutate(val = (val[1] - val)/val[1] * 100 )  %>%
  ungroup() %>%
  group_by(Scenario, Season, Costs) %>%
  summarize(mean  = mean(val), 
            lower = quantile(val, prob = c(0.025)), 
            upper = quantile(val, prob = c(0.975))) %>%
  ungroup() %>% mutate(indirect = "infection")


# Indirect protection against onward transmission and infection
indirect <- "both"
source("4b_process_indirect_results.R")

cases_both  <- cases  %>% mutate(burden = "cases",            indirect = "both")
hosps_both  <- hosps  %>% mutate(burden = "hospitalizations", indirect = "both")
deaths_both <- deaths %>% mutate(burden = "deaths",           indirect = "both")

burden_both <- both_seasons %>% 
  select(Scenario, Costs, Season, nsim, time, "Died6") %>% 
  mutate(val = (Died6)) %>%
  group_by(nsim, Season, Costs) %>%
  mutate(val = (val[1] - val)/val[1] * 100 )  %>%
  ungroup() %>%
  group_by(Scenario, Season, Costs) %>%
  summarize(mean  = mean(val), 
            lower = quantile(val, prob = c(0.025)), 
            upper = quantile(val, prob = c(0.975))) %>%
  ungroup() %>% mutate(indirect = "both")


## Absolute numbers  ----------------------------------------------------------

all10 <- rbind(cases0,      hosps0,      deaths0,
               cases_trans, hosps_trans, deaths_trans,
               cases_inf,   hosps_inf,   deaths_inf,
               cases_both,  hosps_both,  deaths_both)


high10 <- all10 %>% filter(Scenario != "Baseline", 
                           Season == "High",
                           burden == "hospitalizations") %>%
  mutate(indirect = factor(indirect, 
                           levels = c("none", "transmission", "infection", "both"))) %>%
  ggplot(aes(x = indirect, y = mean, fill = indirect)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  mythemeS9 + facet_wrap(~Costs, ncol = 3) +
  labs(title = "High severity", y = "Number averted", x = "") + 
  scale_fill_manual("Indirect\nprotection", values = mylapaz_rep) +
  scale_y_continuous(labels = scales::comma)

low10 <- all10 %>% filter(Scenario != "Baseline", 
                          Season == "Low",
                          burden == "hospitalizations") %>%
  mutate(indirect = factor(indirect, 
                           levels = c("none", "transmission", "infection", "both"))) %>%
  ggplot(aes(x = indirect, y = mean, fill = indirect)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  mythemeS9 + facet_wrap(~Costs, ncol = 3) +
  labs(title = "Low severity", y = "Number averted", x = "") + 
  scale_fill_manual("Indirect\nprotection", values = mylapaz_rep) +
  scale_y_continuous(labels = scales::comma)


## Percentage & figure ----------------------------------------------------------

burden10 <- rbind(burden0, burden_trans, burden_inf, burden_both)


high10p <- burden10 %>% filter(Scenario != "Baseline", 
                               Season == "High") %>%
  mutate(indirect = factor(indirect, 
                           levels = c("none", "transmission", "infection", "both"))) %>%
  ggplot(aes(x = indirect, y = mean, fill = indirect)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  mythemeS9 + facet_wrap(~Costs, ncol = 3) +
  labs(title = "High severity", y = "Percentage averted", x = "") + 
  scale_fill_manual("Indirect\nprotection", values = mylapaz_rep) +
  scale_y_continuous(labels = scales::comma)


low10p <- burden10 %>% filter(Scenario != "Baseline", 
                              Season == "Low") %>%
  mutate(indirect = factor(indirect, 
                           levels = c("none", "transmission", "infection", "both"))) %>%
  ggplot(aes(x = indirect, y = mean, fill = indirect)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  mythemeS9 + facet_wrap(~Costs, ncol = 3) +
  labs(title = "Low severity", y = "Percentage averted", x = "") + 
  scale_fill_manual("Indirect\nprotection", values = mylapaz_rep) +
  scale_y_continuous(labels = scales::comma)


FigS9 <- (high10  + low10) /
  (high10p + low10p) + 
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = "A")

ggsave(FigS9, width = 10, height = 11, file = "figures/FigS9.pdf")



## Fig S10: PRCC -----------------------------------

source("5d_sensitivity_process_multiway.R")

lhs_labels <- c(hd_frac = "Baseline HDAV uptake", 
                hd_frac_miss = "Reduction in overall coverage",
                hd_frac_gain = "Increase in HDAV uptake",
                ndelay = "Delay in additional HDAVs",
                hdVE_mult = "Relative VE of HDAVs")

FigS10 <- rbind(prc_high, prc_low) %>%
  ggplot(aes(x = param, y = original, colour = Season)) +
  geom_point(position = position_dodge(.5), size = 2) +
  geom_errorbar(aes(ymin = `min. c.i.`, ymax = `max. c.i.`), 
                width = 0.35, size = 0.35, 
                position = position_dodge(.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_color_manual("Season\nseverity", values = my_lapaz) +
  scale_y_continuous("Partial rank\ncorrelation coefficient", limits = c(-1, 1)) +
  scale_x_discrete("", labels = str_wrap(lhs_labels[sort(names(lhs_labels))], width = 12)) +
  mytheme10

ggsave(FigS10, width = 7, height = 3, file = "figures/FigS10.pdf")


## Fig S11: TSA different parameters --------------------------------------------

load("results/tsa_main.RData")
burden_tsa1 <- burden_tsa

load("results/tsa_supplement.RData")
burden_tsa2 <- burden_tsa


zlims <- c( max = max(burden_tsa1$hosps_tsa$val, burden_tsa2$hosps_tsa$val),
            min = min(burden_tsa1$hosps_tsa$val, burden_tsa2$hosps_tsa$val))

tsa_plot1 <- ggplot(burden_tsa1$hosps_tsa, 
                    aes(x = hd_frac_miss, y = ndelay/7, fill = val)) +
  geom_tile() + 
  geom_contour(aes(z = val), breaks = c(0), colour = "black") + 
  facet_wrap(~ season) + 
  scale_fill_scico("Hospitalizations\naverted", 
                   limits = c(zlims["min"], zlims["max"]),
                   breaks = c(),
                   palette = "lapaz") +
  mytheme10 + 
  scale_x_continuous(labels = scales::label_percent()) +
  theme( strip.text.x     = element_text(size = 9),
         strip.background = element_blank() ) +
  labs(x = "Decrease in overall coverage", 
       y = "Delay in additional\nHDAV uptake (weeks)",
       title = "Benefits = 10% increase; baseline uptake = 75%; relative VE = 15%") 


tsa_plot2 <- ggplot(burden_tsa2$hosps_tsa, 
                    aes(x = hd_frac_miss, y = ndelay/7, fill = val)) +
  geom_tile() + 
  geom_contour(aes(z = val), breaks = 0, colour = "black") + 
  facet_wrap(~ season) + 
  scale_fill_scico("Hospitalizations\naverted", 
                   limits = c(zlims["min"], zlims["max"]),
                   palette = "lapaz") +
  mytheme10 + 
  scale_x_continuous(labels = scales::label_percent()) +
  theme( strip.text.x     = element_text(size = 9),
         strip.background = element_blank() ) +
  labs(x = "Decrease in overall coverage", 
       y = "Delay in additional\nHDAV uptake (weeks)",
       title = "Benefits = 15% increase; baseline uptake = 75%; relative VE = 15%") 


FigS11 <- tsa_plot1 / tsa_plot2 + 
  plot_layout(guides = "collect") & plot_annotation(tag_levels = "A")

ggsave(FigS11, width = 7, height = 6, file = "figures/FigS11.pdf")


## Fig S12 Compare distributions -----

## Uniform distribution ------------

distribution <- FALSE

source("4a_process_main_results.R")

cases1  <- cases  %>% mutate(distribution = "uniform")
hosps1  <- hosps  %>% mutate(distribution = "uniform")
deaths1 <- deaths %>% mutate(distribution = "uniform")


## Logit-normal distribution -------

distribution <- TRUE

source("4a_process_main_results.R")

cases2  <- cases  %>% mutate(distribution = "logit")
hosps2  <- hosps  %>% mutate(distribution = "logit")
deaths2 <- deaths %>% mutate(distribution = "logit")


## Main plot -------------

p_hosps2 <- rbind(hosps1, hosps2) %>% 
  filter(Scenario != "Baseline") %>%
  mutate(combined = paste0(Season, " (", distribution, ")")) %>%
  ggplot(aes(x = Scenario, y = mean, fill = combined)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  mytheme8 + theme(axis.text.x = element_blank(), 
                   axis.ticks.x = element_blank()) +
  facet_wrap(~Costs, ncol = 3) +
  labs(title = "Hospitalizations", y = "Number averted", x = "") + 
  scale_fill_manual("Season\nseverity", values = mylapaz_rep) +
  scale_y_continuous(labels = scales::comma)

ggsave("figures/FigS12.pdf", p_hosps2, width = 6, height = 6)
