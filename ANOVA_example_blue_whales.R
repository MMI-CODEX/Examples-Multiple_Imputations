# Multiple Imputations ANOVA 
#
# KC Bierlich 
# Center of Drone Excellence (CODEX), Marine Mammal Institute (MMI), Oregon State University
# 241025


# Incorporating uncertainty into an ANOVA model 

#
# Example: Examining differences in body area index (BAI) between different
# populations of blue whales sampled in three locations: 
# 1) Monterey Bay, CA, USA ('CA'),
# 2) Corcovado Gulf, Chile ('Chile'),  
# 3) South Taranaki Bight, New Zealand ('NZ').

# Data from: 
# Barlow et al. (2024). Shaped by their environment: variation in 
# blue whale morphology across three productive coastal ecosystems. 
# Integrative Organismal Biology, obad039, https://doi.org/10.1093/iob/obad039


# Using multiple imputations to account for uncertainty in the ANOVA
# BAI ~ Population


# Load packages
library(tidyverse)
library(car)
library(multcomp)
library(coda)
library(ggpubr)
library(mvtnorm)

#
# import and set up data
#


bw_data <- read.csv("Multiple_Imputations/data/Barlow_etal_2024_Blue_whale_body_area_index.csv") %>% 
  rename(Population = Location)
nrow(bw_data)  # n = 53

View(bw_data)

#
# Plot body condition estimates with uncertainty for each individual for each group.
#

## Boxplots
bw_bai_boxplots <- bw_data %>% ggplot() + theme_bw() +
  geom_boxplot(aes(x = Population, y = BAI.mean), alpha = 0.6) +
  geom_pointrange(aes(x = Population, y = BAI.mean, ymin = BAI.lower, ymax = BAI.upper), position = "jitter", color = "skyblue2") +
  ylab("BAI") + xlab("")
bw_bai_boxplots
ggsave(file.path("Multiple_Imputations/figures/BW_comparison_BAI.png"), width = 5, height = 5)

# But how do we incorporate the uncertainty associated with each individual in an ANOVA to compare
# body condition across populations, BAI ~ Location?

# We will use Multiple Imputations to incorporate uncertainty associated with BAI into an ANOVA, BAI ~ Location
# Uncertainty:
# 1. We don't know the true BAI, but have an estimate of BAI with uncertainty
# 2. Uncertainty in ANOVA incorporate multivariate normal distribution acounts 

#
# MI  model
#

tick = proc.time()[3] # to track how long it takes to run based on n replicates
mi_anova = replicate(n = 1e4, expr = { # set n replicates
  
  # save data frame
  dat = bw_data
  
  # For each individual, draw a random sample from a normal distribution that is described by
  # that individual's mean and sd of their posterior distribution
  dat$BAI = rnorm(n = nrow(dat), mean = dat$BAI.mean, sd = sqrt(dat$BAI.var))
  
  # fit model
  fit.mod <- lm(BAI ~ Population-1, data = dat) 
  
 
  # sample coefficients and fitted values from a multivariate normal distribution to quantify their uncertainty 
  matrix_coef <- mvtnorm::rmvnorm(n = 1, mean = coef(fit.mod), sigma = vcov(fit.mod))

  # Set up ANOVA to analyze results from fit model and set as Type II, due to unbalanced samples sizes across locations. 
  bai_anova <- Anova(fit.mod, type = 2)
  
  # select outputs to save
  c(# Summary outputs from ANOVA:
    SS_population = bai_anova$`Sum Sq`[1], # sum of squares population
    SS_residuals = bai_anova$`Sum Sq`[2],  # sum of squares residuals
    DF_between = bai_anova$DF[1],          # degrees of freedom between
    DF_within = bai_anova$DF[2],           # degrees of freedom between
    F_value = bai_anova$`F value`[1],      # f-value
    p_value = bai_anova$`Pr(>F)`[1],       # p-value
    
    # Coefficients
    CA_coef = matrix_coef[1], 
    Chile_coef = matrix_coef[2], 
    NZ_coef = matrix_coef[3], 
    
    # Multiple comparisons 
    "CA - Chile" = matrix_coef[1] - matrix_coef[2], 
    "NZ - CA" = matrix_coef[3] - matrix_coef[1], 
    "NZ - Chile" = matrix_coef[3] - matrix_coef[2])
})
tock = proc.time()[3]
tock - tick     

#
# Save outputs as a tibble
#
mi_outputs <- tibble(label = rownames(mi_anova), 
                       means = round(rowMeans(mi_anova), 4),
                           HPD_lwr = round(HPDinterval(mcmc(t(mi_anova)))[,1],4),
                           HPD_upr = round(HPDinterval(mcmc(t(mi_anova)))[,2],4))

View(mi_outputs)

# Can save outputs
write.csv(mi_outputs, file.path("Multiple_Imputations/outputs/BW_comparison_outputs.csv"))

#
# plot results
#
# an overlap with 0 indicates no statistical difference
MI_bw_bai_output <- mi_outputs %>%
  filter(label == "CA - Chile" | label == "NZ - CA" | label == "NZ - Chile") %>%
  ggplot(aes(x = label, y = means, ymin = HPD_lwr, ymax = HPD_upr)) + geom_pointrange() + 
  geom_abline(intercept = 0, slope = 0, lty = 2) + xlab("") + ylab("Difference in BAI") + 
  theme_bw()
MI_bw_bai_output
ggsave(file.path("Multiple_Imputations/figures/BW_mult_comparison_BAI.png"), width = 5, height = 5)


# Plot Boxplots with MC results
ggarrange(bw_bai_boxplots, MI_bw_bai_output)
ggsave(file.path("Multiple_Imputations/figures/Blue_whale_results.png"), width = 10, height = 5)



# {!} Challenge {!}
#
# Set up an ANOVA comparing TL across these blue whale populations. 
