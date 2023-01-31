## Check baseline output when indirect effects of vaccination are added

## Load packages and prep code -------------------------------------------------

rm(list=ls())

library(patchwork)
library(viridis)
library(scico)

source("3a_baseline_simulation.R")


## Plot settings ---------------------------------------------------------------

linesize4 <- 1

textsize <- 8
mytheme8 <- theme_bw() + theme(strip.text.x = element_text(size = textsize),
                               axis.text  = element_text(size = textsize - 1),
                               axis.title   = element_text(size = textsize),
                               title        = element_text(size = textsize),
                               legend.title = element_text(size = textsize),
                               legend.text  = element_text(size = textsize),
                               strip.placement  = "outside", 
                               strip.background = element_blank())


## Simulate changing levels of indirect effects --------------------------------

pargrid <- expand.grid( VEs_mult = c(0, 0.1, 0.5),
                        VEi_mult = c(0, 0.1, 0.5) ) 

low_list <- high_list <- list()

for (k in 1:nrow(pargrid) ) {
  
  pars$value[pars$name == "VE_inf"]   <- pargrid$VEs_mult[k]
  pars$value[pars$name == "VE_trans"] <- pargrid$VEi_mult[k]
  
  low_list[[k]] <- simulate_model(season = "Low", 
                                  pars = pars, diffpars = diffpars,
                                  chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                                  vaccineMultiplier = vaccineMultiplier, 
                                  hd_frac = 0.75, hd_frac_miss = 0,
                                  daily_coverage = daily_coverage, 
                                  hd_daily_coverage = hd_daily_coverage,
                                  rel_VE_by_age = rel_VE_by_age, 
                                  rel_hdVE_by_age = rel_hdVE_by_age) %>%
    mutate(sim = k, 
           VEs_mult = pargrid$VEs_mult[k], 
           VEi_mult = pargrid$VEi_mult[k])
  
  high_list[[k]] <- simulate_model(season = "High", 
                                   pars = pars, diffpars = diffpars,
                                   chr = chr, hfr = hfr, fracSympt_by_age = fracSympt_by_age,
                                   vaccineMultiplier = vaccineMultiplier, 
                                   hd_frac = 0.75, hd_frac_miss = 0,
                                   daily_coverage = daily_coverage, 
                                   hd_daily_coverage = hd_daily_coverage,
                                   rel_VE_by_age = rel_VE_by_age, 
                                   rel_hdVE_by_age = rel_hdVE_by_age) %>%
    mutate(sim = k, 
           VEs_mult = pargrid$VEs_mult[k], 
           VEi_mult = pargrid$VEi_mult[k])
}

high_range2 <- bind_rows(high_list)
low_range2 <- bind_rows(low_list)

all_range2 <- rbind(low_range2, high_range2)


## Create plot to check timing -------------------------------------------------

sympt_range2 <- all_range2 %>% 
  select(time, Scenario, VEs_mult, VEi_mult, sim, contains("Discharged")) %>% 
  gather(age, val, contains("Discharged")) %>%
  mutate(age_grp = str_extract(age, pattern = "[0-9]+"),
         val = val/chr[as.numeric(age_grp)]
  ) %>%
  group_by(age_grp, Scenario, VEs_mult, VEi_mult, sim) %>%
  mutate(val = val * 100000,
         I_inc = diff(c(0, val), lag = 1)) %>%
  ungroup() %>%
  group_by(Scenario, time, VEs_mult, VEi_mult, sim) %>%
  summarize(I_inc = sum(I_inc)) %>% ungroup() 

simlabels2 <- paste0("VEs % = ", pargrid$VEs_mult * 100, ", VEi % = ", pargrid$VEi * 100)

inc_range2 <- 
  ggplot(sympt_range2, aes(x = time, y = I_inc, colour = as.factor(sim))) +
  geom_rect(data = NULL, aes(xmin = 150, xmax = 270, ymin = 0, ymax = Inf), 
            alpha = 0.2, fill = "grey90", colour = NA) +
  geom_line(size = linesize4) + 
  facet_wrap(~ Scenario) +  
  scale_colour_viridis(NULL, discrete = TRUE, direction = -1, labels = simlabels2) +
  labs(y = "Symptomatic cases\nper 100,000") +
  scale_x_continuous("", labels = monthnames, breaks = c(0, cumsum(monthdays))) +
  mytheme8 + 
  guides(colour = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "top")


## Check total burden output ---------------------------------------------------

# Cases
cases <- sympt_range2 %>% 
  group_by(Scenario, VEs_mult, VEi_mult, sim) %>% 
  summarize(Total_ILI = sum(I_inc) / 100000 * popsize) %>%
  ungroup() %>%
  select(Scenario, VEs_mult, VEi_mult, val = Total_ILI)%>%
  mutate(burden = "Symptomatic cases")

# Hospitalizations
hosps <- all_range2 %>% 
  group_by(Scenario, VEs_mult, VEi_mult, sim) %>%
  filter(time == max(time)) %>% 
  ungroup() %>%
  select(Scenario, VEs_mult, VEi_mult, contains("Discharged")) %>%
  mutate(val = rowSums(.[,-c(1:3)], na.rm = TRUE) * popsize)  %>%
  select(Scenario, VEs_mult, VEi_mult, val)%>%
  mutate(burden = "Hospitalizations")

# Deaths
deaths <- all_range2 %>% 
  group_by(Scenario, VEs_mult, VEi_mult, sim) %>%
  filter(time == max(time)) %>% 
  ungroup() %>%
  select(Scenario, VEs_mult, VEi_mult, contains("Died")) %>%
  mutate(val = rowSums(.[,-c(1:3)], na.rm = TRUE) * popsize) %>%
  select(Scenario, VEs_mult, VEi_mult, val) %>%
  mutate(burden = "Deaths")

burden_range2 <- rbind(cases, hosps, deaths)


## Plot effects on burden output -----------------------------------------------

plot_func <- function(df, name = name) {
  ggplot(data = df, aes(x = as.factor(VEs_mult), y = as.factor(VEi_mult), fill = val)) +
    geom_tile() + 
    scale_fill_scico(name, labels = scales::comma, palette = "lapaz") +
    mytheme8 + 
    labs(x = "VE against infection", y = "VE against transmission") +
    scale_x_discrete(labels = c("0%", "10%", "50%"))
}

nested_high <- burden_range2 %>% 
  filter(Scenario == "High") %>%
  group_by(burden) %>% 
  nest() %>% 
  mutate(plots = map2(data, burden, plot_func)) 

high_range2 <- gridExtra::grid.arrange(grobs = nested_high$plots)

nested_low <- burden_range2 %>% 
  filter(Scenario == "Low") %>%
  group_by(burden) %>% 
  nest() %>% 
  mutate(plots = map2(data, burden, plot_func)) 

low_range2 <- gridExtra::grid.arrange(grobs = nested_low$plots)

fig_check <- ggpubr::ggarrange(high_range2, low_range2, 
                               nrow = 1, labels = c("A", "B"), 
                               font.label = list(size = 12, face = "plain"))

