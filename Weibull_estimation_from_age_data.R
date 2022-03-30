### Information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022 (not provided on github)

# The purpose of this code is to get the parameters for the Weibull prior to use in the age estimation for the other elephant population

### Set up ####
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(survival)
library(brms)
library(tidybayes)
library(ggthemes)
library(dplyr)

### Import nodes data ####
setwd("~/Downloads/ElephantAges")
all_nodes <- readxl::read_excel('data_raw/Raw_ATE_AllElephants_Lee220118.xlsx')
boxplot(as.numeric(all_nodes$Age_death) ~ all_nodes$Death_accuracy) # adds NA values for anything classed as "living"
str(all_nodes)

### calculate oldest age each observed ####
# create dataframe that only considers node ages
ages <- all_nodes[,c(1:3,15,16,18:20)] %>% 
  janitor::clean_names()
ages$id_no <- ifelse(ages$sex == 'unknown', paste('U', ages$id, sep = ''),
                     ifelse(ages$sex == 'Female', paste('F', ages$id, sep = ''),
                            ifelse(ages$sex == 'Male', paste('M', ages$id, sep = ''),
                                   'check')))

# generate numeric variable for age in 2021 of living elephants or age at death
ages$age_death_num <- as.numeric(ifelse(ages$age_death == 'living', 2021 - ages$birth_year, ages$age_death))
ages$age_death_round <- round(ages$age_death_num, 0)
table(ages$age_death_round, ages$sex)
hist(ages$age_death_round, breaks = 75)                               # overall age distribution
hist(ages$age_death_round[which(ages$sex == 'Male')], breaks = 70)    # male age distribution
hist(ages$age_death_round[which(ages$sex == 'Female')], breaks = 75)  # female age distribution
ages <- ages[,c(9,1,3,8,10)]

# create new 
ages$censor <- ifelse(ages$age_death == 'living', 'TRUE', 'FALSE')    # censor = TRUE when elephant is still alive

# clean up
age_cens <- ages[!is.na(ages$age_death_num),c('age_death_num','censor','sex')] # single data frame for age at death/2021
rm(ages, all_nodes)

### estimate Weibull distribution ####
colnames(age_cens)[1] <- 'age'
age_cens$age_non0 <- age_cens$age+0.01 # wouldn't run when I allowed age = 0, so increased all values by 0.01 to make it run

age_cens$censor <- as.integer(!as.logical(age_cens$censor)) # Format needed for brms censoring

# Examine default brms priors (male elephants)
bfit_m_default <- age_cens %>% filter(sex == 'Male') %>%
  brm(age_non0 | cens(censor) ~ 1,
      data = ., family = "weibull", 
      chains = 4, cores = 4)
prior_summary(bfit_m_default)
# Defaults
#intercept - student_t(3, 2.3, 2.5)
#shape - gamma(0.01, 0.01)

# Prior predictive check
num_samples <- 500

# Default priors
prior_intercept <- rstudent_t(num_samples, 3, 2.3, 2.5)
prior_shape <- rgamma(num_samples, 0.01, 0.01)

# Priors used
prior_intercept <- rstudent_t(num_samples, 3, 5, 5)
prior_shape <- rexp(num_samples, 0.5)
prior_scale <- exp(prior_intercept) / (gamma(1 + 1 / prior_shape))

hist(prior_shape)
hist(prior_scale[prior_scale<100])

probs <- dweibull(1:100, shape = prior_shape[1], scale = prior_scale[1])
probs <- probs / sum(probs)
plot(probs, type="l", lwd=0.5)
for( i in 2:num_samples) {
  probs <- dweibull(1:100, shape = prior_shape[i], scale = prior_scale[i])
  probs <- probs / sum(probs)
  lines(probs, lwd=0.5)
}


# Run model for male elephants with slightly better priors
bfit_m <- age_cens %>% filter(sex == 'Male') %>%
  brm(age_non0 | cens(censor) ~ 1,
      prior = c(
        prior(student_t(3, 5, 5), class = Intercept),
        prior(exponential(0.5), class = shape)
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)

plot(bfit_m)

# extract scale and shape from posterior draws
# to convert to scale we need to both undo the link function by taking the exponent
# and then refer to the brms documentation to understand how the mean relates to the scale
bfit_m_draws <- as_draws_df(bfit_m) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_m_draws$shape)
est_scale <- mean(bfit_m_draws$scale)

# Let's do a posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs) 

# Lets simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale))

# Shape = 1.29, scale = 29.5 (so we will use shape = 1.3 and scale = 30)
