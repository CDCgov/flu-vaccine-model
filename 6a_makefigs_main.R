## Load packages and prep code -------------------------------------------------

rm(list=ls())

library(tidyverse)
library(patchwork)
library(ggforce)
library(scico)

source("3a_baseline_simulation.R")


## Colour palettes -------------------------------------------------------

my_lapaz <- scico(7, palette = 'lapaz')[c(1, 4)]


## Plot settings -------------------------------------------------------

errorwidth <- 0.25
errorsize  <- 0.35

linesize <- 1.3

textsize <- 9
mytheme1 <- theme_bw() + theme(strip.text.x = element_text(size = textsize - 1),
                               strip.text.y = element_text(size = textsize - 1),
                               axis.text.y  = element_text(size = textsize - 1),
                               axis.title   = element_text(size = textsize),
                               title        = element_text(size = textsize),
                               legend.title = element_text(size = textsize),
                               legend.text  = element_text(size = textsize),
                               strip.placement  = "outside", 
                               strip.background = element_blank(),
                               axis.text.x  = element_blank(), 
                               axis.ticks.x = element_blank(),
                               plot.title   = element_text(size = textsize)) 

textsize <- 9
mytheme2 <- theme_bw() + theme(strip.text.x = element_blank(),
                               axis.text    = element_text(size = textsize - 1),
                               axis.title   = element_text(size = textsize),
                               title        = element_text(size = textsize),
                               legend.title = element_text(size = textsize),
                               legend.text  = element_text(size = textsize),
                               plot.title   = element_text(size = textsize))



## Fig 1: Main results ---------------------------------------------------------

# FALSE for uniform rVE
distribution <- FALSE

source("4a_process_main_results.R")

# absolute # averted
p_absolute <- hosps %>% filter(Scenario != "Baseline") %>%
  ggplot(aes(x = Scenario, y = mean, fill = Season)) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.85) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  mytheme1 + facet_wrap(~Costs, ncol = 3) +
  labs(y = "Number of hospitalizations averted", x = "") + 
  scale_fill_manual("Season\nseverity", values = my_lapaz) +
  scale_y_continuous(labels = scales::comma)

# percentage averted
burden <- both_seasons %>% 
  select(Scenario, Costs, Season, nsim, time, "Died6") %>% 
  mutate(val = (Died6)) %>%
  group_by(nsim, Season, Costs) %>%
  mutate(val = (val[1] - val)/val[1] * 100 )  %>%
  ungroup() %>%
  group_by(Scenario, Season, Costs) %>%
  summarize(mean  = mean(val), 
            lower = quantile(val, prob = c(0.025)), 
            upper = quantile(val, prob = c(0.975)),
            min   = min(val),
            max   = max(val)) %>%
  ungroup()

p_percentage <- burden %>% filter(Scenario != "Baseline") %>%
  ggplot(aes(x = Scenario, y = mean, fill = Season)) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.85) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = errorwidth, size = errorsize,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  mytheme1 + facet_wrap(~Costs, ncol = 3) +
  labs(y = "Percentage of hospitalizations averted", x = "") + 
  scale_fill_manual("Season\nseverity", values = my_lapaz) +
  scale_y_continuous()

fig1 <- p_absolute + p_percentage + 
          plot_layout(guides = "collect") &
          plot_annotation(tag_levels = "A")

ggsave("figures/Fig1.pdf", fig1, width = 7, height = 6)


## Fig 2: Sensitivity analyses -------------------------------------------------

## Fig 2A: One way -------------------

load("results/osa.RData")

osa_early <- osa$early
osa_late <- osa$late

osa_pars <- osa$pars %>% 
            mutate(which_var_labels = gsub(which_var_labels, pattern = "EV", replace = "HDAV"))

late_hosps <- osa_late %>% 
  filter(time == max(time)) %>%
  select(season = Scenario, nsim, time, "Discharged6") %>% 
  mutate(val = (Discharged6)) %>%
  group_by(nsim, season) %>%
  mutate(val = round((val[1] - val) * popsize) )  %>%
  summarize(val = sum(val)) %>%
  ungroup() %>% left_join(osa_pars) 

early_hosps <- osa_early %>% 
  filter(time == max(time)) %>%
  select(season = Scenario, nsim, time, "Discharged6") %>% 
  mutate(val = (Discharged6)) %>%
  group_by(nsim) %>%
  mutate(val = round((val[1] - val) * popsize) )  %>%
  summarize(val = sum(val)) %>%
  ungroup() %>% left_join(osa_pars) 

plotting_order <- early_hosps %>% group_by(which_var_labels) %>% 
  summarize(val1 = max(max(val) - min(val)) ) %>% arrange(val1) 

tornado_both <- 
  late_hosps %>% mutate(Season = "Low") %>%
  rbind(., {
    early_hosps %>% mutate(Season = "High")
  }) %>%
  ggplot(aes(y = val, x = which_var_labels, colour = Season)) + 
  geom_link2(size = 4, position = position_dodge(0.75), lineend = "round") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  mytheme2 +
  scale_colour_manual("Season\nseverity", values = my_lapaz) +
  labs(    y = "Hospitalizations averted", 
           x = "", 
           title = "One-way") +
  scale_x_discrete(labels = label_wrap_gen(width = 21), limits = plotting_order$which_var_labels ) +
  scale_y_continuous(labels = scales::comma) +
  coord_flip()


## Fig 2B: Two way -------------------

load("results/tsa_main.RData")

burden_tsa1 <- burden_tsa$hosps_tsa %>% 
  mutate(season = paste(season, "severity"))

tsa_plot1 <- ggplot(burden_tsa1 , 
                    aes(x = hd_frac_miss, y = ndelay/7, fill = val)) +
  geom_tile() + 
  geom_contour(aes(z = val), breaks = c(0), colour = "black") + 
  facet_wrap(~ season) + 
  scale_fill_scico("Hospitalizations\naverted", 
                   palette = "lapaz") + 
  mytheme2 + 
  scale_x_continuous(labels = scales::label_percent()) +
  theme( strip.text.x     = element_text(size = 9),
         strip.background = element_blank() ) +
  labs(x = "Decrease in overall coverage", 
       y = "Delay in additional\nHDAV uptake (weeks)",
       title = "Two-way") 


## Main fig -------------

fig2 <- tornado_both / 
  ( tsa_plot1  + theme(axis.title.y = element_text(margin = margin(r = -3.5, unit = "cm"))) ) +
  plot_layout(heights = c(1, 1)) & 
  plot_annotation(tag_levels = 'A') 

ggsave("figures/Fig2.pdf", fig2, width = 7, height = 5)
