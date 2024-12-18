# Allometry_example_Bayesian_vs_multiple_imputations
# KC Bierlich & Josh Hewitt
# Dec. 18, 2024

#
# Multiple Imputations vs. full Bayesian approach
#

# Here we will compare the relationship between rostrum-blowhole (RB) and total length (TL) of Antarctic minke whales using a linear model.  
# 
# RB ~ TL      
# 
# We will compare the results using two methods for including measurement uncertainty into the analysis:  
# 1) full Bayesian model  (single-stage )
# 2) multiple imputations (two-stage)  
# 
# Will we use the RB and TL measurement outputs with uncertainty from Bierlich et al., 2021 MEPS, 
# which uses a Bayesian model to estimate photogrammetric measurements and associated uncertainty.
# We will then calculate the linear relationship of RB ~ TL with associated uncertainty 
# using multiple imputations and compare these results with the full Bayesian linear model results 
# from  Bierlich et al., 2021. 
# 
# 
# Bierlich KC, Schick RS, Hewitt J, Dale J, Goldbogen JA, Friedlaender AS, 
# Johnston DW (2021) Bayesian approach for predicting photogrammetric uncertainty 
# in morphometric measurements derived from drones. Mar Ecol Prog Ser 673:193-210. 
# https://doi.org/10.3354/meps13814


# load packages
library(tidyverse)
library(coda)

# import data
meps <- read.csv("data/Bierlich_etal_2021_AMW_outputs_Scenario2_Model2.csv")


#
# Multiple Imputations
#
# (run this chunk, lines 42-89)
tick = proc.time()[3] # to track how long it takes to run based on n replicates
# save data frame
dat = meps

multi_imp2 = replicate(n = 1e4, expr = { # set n replicates
  
  
  
  # For each individual, draw a random sample from a normal distribution that is described by
  # that individual's mean and sd of their posterior distribution
  dat$RB = rnorm(n = nrow(dat), mean = dat$RB_mean, sd = dat$RB_sd)
  dat$TL = rnorm(n = nrow(dat), mean = dat$TL_mean, sd = dat$TL_sd)
  
  # fit model
  fit.mod <- lm(RB ~ TL, data = dat) 
  
  # sample coefficients and fitted values from a multivariate normal distribution to quantify their uncertainty 
  matrix_coef <- mvtnorm::rmvnorm(n = 1, mean = coef(fit.mod), sigma = vcov(fit.mod))
  
  
  # save fitted samples to estimate regression
  fitted_lm_samples <- predict(
    object = fit.mod, 
    se.fit = TRUE, 
    # Build a data frame 
    newdata = expand.grid(TL = seq(from = min(dat$TL_mean), to = max(dat$TL_mean, length.out = 5)))
  )
  
  # sample fitted values from a normal distribution
  fitted_lm_samples <- rnorm(
    n = length(fitted_lm_samples$fit), 
    mean = fitted_lm_samples$fit,
    sd = fitted_lm_samples$se.fit)
  
  
  # select outputs to save
  c(# Coefficients
    matrix_coef_int = matrix_coef[1], 
    matrix_coef_tl = matrix_coef[2], 
    
    # means
    meanTL = mean(dat$TL),
    meanRB = mean(dat$RB), 
    fitted_lm_samples
  )
})
tock = proc.time()[3]
tock - tick


#
# Save coefficients
#
mi_outputs2 <- tibble(label = rownames(multi_imp2), 
                      means = round(rowMeans(multi_imp2), 4),
                      HPD_lwr = round(HPDinterval(mcmc(t(multi_imp2)))[,1],4),
                      HPD_upr = round(HPDinterval(mcmc(t(multi_imp2)))[,2],4))
mi_outputs2[1:4,]


## save the expected vales of linear model
expected2 <- data.frame(exp_TL = seq(from = min(dat$TL_mean), to = max(dat$TL_mean, length.out = 5)), 
                        exp_RB = as.list(mi_outputs2[5:9,2]),
                        as.list(mi_outputs2[5:9,3]),
                        as.list(mi_outputs2[5:9,4]))


#
# AMW Bayesian outputs
#

bayes_lm_sum <- read.csv("data/Bierlich_etal_2021_Bayesian_lm_sum.csv")

# filter out betas (for ease of plotting below)
bayes_lm_sum_expTL <- bayes_lm_sum %>% filter(label != "beta1" & label != "beta2")

#
# Plot results
#
meps %>%
  ggplot() + theme_bw() + 
  geom_pointrange(aes(x = TL_mean, y = RB_mean, ymin = RB_lower, ymax = RB_upper)) + 
  geom_point(aes(x = TL_mean, y = RB_mean)) +
  geom_errorbarh(aes(xmin = TL_lower, xmax = TL_upper, y = RB_mean)) + 
  xlab("predicted TL (m)") + ylab("predicted RB (m)") + 
  #xlim(4, 11) + ylim(0.25, 2) +
  geom_line(data = expected2, aes(x = exp_TL, y = means, color = "blue")) +
  geom_ribbon(data = expected2, aes(x = exp_TL, ymin = HPD_lwr, ymax = HPD_upr), fill = "blue", alpha = 0.2) +
  
  #geom_abline(aes(slope = 0.175, intercept = -0.223, color = "red")) +
  geom_line(data = bayes_lm_sum_expTL, 
            aes(x = expectedTL, y = mean, color = "red")) +
  geom_ribbon(data = bayes_lm_sum_expTL, 
              aes(x = expectedTL, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
  
  scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
  scale_colour_manual(name = 'model', 
                      values =c('blue'='blue','red'='red'), 
                      labels = c('multiple imputations','full Bayesian'))
# to save
#ggsave("figures/Allometry_example_Bayesian_vs_MI.png", height = 5, width = 7)


# {!} Both models yield similar results.
