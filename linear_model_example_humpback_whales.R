# Multiple Imputations Linear Model
#
# KC Bierlich 
# Center of Drone Excellence (CODEX), Marine Mammal Institute (MMI), Oregon State University
# 241025

# Incorporating uncertainty into a linear model


# Example: Exploring how body condition, here measured as body area index (BAI), 
# increases over the foraging season for humpback whales along the Western
# Antarctic Peninsula across different reproductive classes.

# We will use a linear model and incorporate uncertainty associated BAI by using the 
# mean and sd of the posterior distribution for BAI for each individual.

# BAI ~ Reproductive Class * Seasonal Duration


# Data from:
# Bierlich et al. (2022). Seasonal gain in body condition of foraging humpback 
# whales along the Western Antarctic Peninsula. Frontiers in Marine Science. 9, 1–16. 
# doi:10.3389/fmars.2022.1036860 https://www.frontiersin.org/articles/10.3389/fmars.2022.1036860/full  


# load packages
library(tidyverse)
library(coda)
library(mvtnorm)


# Import data
humpback_data <- read.csv("Multiple_Imputations/data/Bierlich_etal_2022_WAP_humpback_whale_body_condition_data.csv") # n = 203


#
# MI model
#

tick = proc.time()[3] # to track how long it takes to run based on n replicates
mi_bai = replicate(n = 1e4, expr = {
  
  dat = humpback_data
  
  # for each individual, draw a random sample from a normal distribution described by
  # that individual's mean and sd of the posterior distribution
  dat$BAI = rnorm(n = nrow(dat), mean = dat$BAI.mean, sd = sqrt(dat$BAI.var))
  
  # fit linear model
  fit.bai <- lm(BAI ~ class*seasonal_duration, data = dat)
  
  # sample coefficients and fitted values from a multivariate normal distribution to quantify their uncertainty 
  matrix_bai_coefs <- mvtnorm::rmvnorm(n = 1, mean = coef(fit.bai), sigma = vcov(fit.bai))
  
  # save fitted samples to estimate regression
  fitted_bai_samples <- predict(
    object = fit.bai, 
    se.fit = TRUE, 
    # Build a data frame with all of the combinations of predictors you want, here seasonal duration (~every 2 weeks) and class.
    newdata = expand.grid(seasonal_duration = seq(from = 20, to = 225, length.out = 10), # this is about every two weeks.
                          class = c('Calf', 'Juvenile Female', 'Juvenile Unknown', 'Lactating Female', 'Mature Unknown', 'Non-pregnant Female'))
    )
  
  # sample fitted values from a normal distribution
  fitted_bai_samples <- rnorm(
    n = length(fitted_bai_samples$fit), 
    mean = fitted_bai_samples$fit,
    sd = fitted_bai_samples$se.fit)
  
  # select outputs to save
  c(matrix_bai_coefs, fitted_bai_samples)
})

tock = proc.time()[3]
tock - tick 


#
# build dataframe
#

length(rowMeans(mi_bai)) # 72
12+(6*10) # 12 coefs + (6 classes, 10 sampling points)  = 72

mi_bai_output_df <- data.frame(label = NA,
                               means = rowMeans(mi_bai), 
                               HPD_lwr = HPDinterval(mcmc(t(mi_bai)))[,1],
                               HPD_upr = HPDinterval(mcmc(t(mi_bai)))[,2])

View(mi_bai_output_df)

# assign labels to linear model coefficient outputs
fit.lm.bai <- lm(BAI.mean ~ class*seasonal_duration, data = humpback_data)
mi_bai_output_df$label[1:12] <- names(fit.lm.bai$coefficients)


# make the linear model coefficients one dataframe and the fitted values another one
mi_bai_coefs <- mi_bai_output_df %>% filter(!is.na(label))   # coefs
mi_bai_fitted_samples <- mi_bai_output_df %>% filter(is.na(label)) # fitted values

## build up the fitted samples dataframe 
mi_bai_fitted_samples_df <- data.frame(seasonal_duration = rep(seq(from = 20, to = 225, length.out = 10), 6),
                                       class = rep(c('Calf', 'Juvenile Female', 'Juvenile Unknown',
                                                     'Lactating Female', 'Mature Unknown', 'Non-pregnant Female'), each = 10),
                                       mi_bai_fitted_samples)

#
# Plot predicted values
#


## filter out predicted BAI that extends beyond sampling days for each class
# here are the min and max days in season for each class
humpback_data %>% group_by(class) %>% 
  summarise(min = min(seasonal_duration), max = max(seasonal_duration))

mi_fit_filt <- mi_bai_fitted_samples_df %>% 
  filter(!(class == "Calf" & (seasonal_duration < 40 | seasonal_duration > 157)) &
           !(class == "Juvenile Female" & (seasonal_duration > 138)) &
           !(class == "Juvenile Unknown" & (seasonal_duration > 157)) &
           !(class == "Lactating Female" & (seasonal_duration < 40 | seasonal_duration > 160))
         ) 

  
# Plot results
          
ggplot() + 
  geom_pointrange(data = humpback_data,
                  aes(x = seasonal_duration, y = BAI.mean, ymin = BAI.lower, ymax = BAI.upper), 
                  color = "skyblue4", alpha = 0.5) + 
  geom_ribbon(data = mi_fit_filt,
              aes(x = seasonal_duration, ymin = HPD_lwr, ymax = HPD_upr), 
              fill = "skyblue4", alpha = .15) + 
  geom_line(data = mi_fit_filt,
            aes(x = seasonal_duration, y = means), color = "skyblue4", size = 1) +
  facet_wrap(~class) +
  scale_x_continuous(limits = c(20, 230)) +
  xlab("Days since Nov. 1") +
  ylab("Predicted BAI") + 
  theme_bw()
ggsave(file.path("Multiple_Imputations/figures/Humpback_BAI_results.png"), width = 8, height = 7)



