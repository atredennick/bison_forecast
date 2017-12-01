################################################################################
##  bison_forecast.R: R script to fit a population growth model for YNP Bison,
##  forecast 7 new years.
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: December 1, 2017
################################################################################

##  Clear everything...
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location
root <- "~/Repos/bison_forecast/"
setwd(paste0(root,"code/"))



####
####  LOAD LIBRARIES ----
####
library(tidyverse) # Data science functions
library(dplyr)     # Data wrangling
library(rjags)     # Fitting Bayesian models with JAGS
library(coda)      # MCMC summaries
# library(devtools) # For installing packages from GitHub
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR) # MCMC manipulation (by M. Dietze)



####
####  LOAD DATA ----------------------------------------------------------------
####
snow_ynp  <- read.csv("../data/west_yellowstone_snotel_summary.csv", row.names = 1) 
bison_raw <- read.csv("../data/YNP_bison_population_size.csv")
bison_dat <- bison_raw %>% 
  dplyr::select(-source) %>%     # drop the source column
  mutate(set = ifelse(year < 2011, "training", "validation")) %>% # make new column for data splits
  left_join(snow_ynp, by="year") # merge in SNOTEL data



####
####  JAGS State-Space Model ---------------------------------------------------
####
r_mu_prior <- log(1.11) # lambda = 1.11 in Hobbs et al. 2015
r_sd_prior <- sd(log(rnorm(100000,1.11,0.024))) # sd_lambda = 0.024 in Hobbs et al. 2015

my_model <- "  
  model{

    #### Variance Priors
    sigma_proc ~ dgamma(0.01,0.01)
    tau_proc <- 1/sigma_proc^2
    
    #### Fixed Effects Priors
    r  ~ dnorm(0.1, 1/0.02^2) # intrinsic growth rate, informed prior
    b  ~ dnorm(0,0.0001)      # strength of density dependence (r/K)
    b1 ~ dnorm(0,0.0001)      # effect of snow
    
    #### Initial Conditions
    z[1] ~ dnorm(Nobs[1], tau_obs[1]) # varies around observed abundance at t = 1
    
    #### Process Model
    for(t in 2:npreds){
      mu[t] <- max( 1, log( z[t-1]*exp(r + b*z[t-1] + b1*x[t]) ) )
      z[t] ~ dlnorm(mu[t], tau_proc)
    }
    
    #### Data Model
    for(j in 2:n){
      Nobs[j] ~ dnorm(z[j], tau_obs[j])
    }
  
  }"



####
####  Fit Bison Forecasting Model ----------------------------------------------
####
##  For years without observation error, set to max observed standard deviation
##  TODO: Impute in the model?
na_sds                       <- which(is.na(bison_dat$count.sd)==T)
bison_dat[na_sds,"count.sd"] <- max(bison_dat$count.sd, na.rm=T)

##  Split into training and validation sets
training_dat   <- filter(bison_dat, set == "training")
validation_dat <- filter(bison_dat, set == "validation")

##  Set up SWE knowns (2011-2017), relative to scaling of observations
swe_mean     <- mean(training_dat$mean_snow_water_equiv_mm)
swe_sd       <- sd(training_dat$mean_snow_water_equiv_mm)
forecast_swe <- snow_ynp %>%
  filter(year %in% validation_dat$year) %>%
  pull(mean_snow_water_equiv_mm)
scl_fut_swe  <- (forecast_swe - swe_mean) / swe_sd

##  Set initial values for unkown parameters
inits <- list(
  list(sigma_proc = 0.01,
       r = 0.05,
       b = -0.001,
       b1 = -0.5),
  list(sigma_proc = 0.3,
       r = 0.4,
       b = -0.1,
       b1 = -0.01),
  list(sigma_proc = 0.1,
       r = 0.7,
       b = -0.00001,
       b1 = -0.2)
)



####
####  FIT AND FORECAST WITH KNOWN SWE
####
##  Prepare data list
mydat <- list(Nobs    = training_dat$count.mean, # mean counts
              n       = nrow(training_dat), # number of observations
              tau_obs = 1/training_dat$count.sd^2, # transform s.d. to precision
              x       = c(as.numeric(scale(training_dat$mean_snow_water_equiv_mm)),scl_fut_swe), # snow depth, plus forecast years
              npreds  = nrow(training_dat)+nrow(validation_dat)) # number of total predictions (obs + forecast)

##  Random variables to collect
out_variables <- c("r", "b", "b1", "sigma_proc", "z")

##  Send to JAGS
mc3     <- jags.model(file = textConnection(my_model), 
                      data = mydat, 
                      n.chains = length(inits), 
                      n.adapt = 50000, 
                      inits = inits) 
           update(mc3, n.iter = 100000) 
mc3.out <- coda.samples(model=mc3, 
                        variable.names=out_variables, 
                        n.iter=100000) 

##  Split MCMC output for file constraints
saveRDS(mc3.out[[1]],"../results/swe_est_posteriors_chain1.RDS")
saveRDS(mc3.out[[2]],"../results/swe_est_posteriors_chain2.RDS")
saveRDS(mc3.out[[3]],"../results/swe_est_posteriors_chain3.RDS")



####
####  FIT AND FORECAST ASSUMING AVG SWE
####
scl_fut_swe[] <- 0 # average is 0 by definition
##  Prepare data list
mydat <- list(Nobs    = training_dat$count.mean, # mean counts
              n       = nrow(training_dat), # number of observations
              tau_obs = 1/training_dat$count.sd^2, # transform s.d. to precision
              x       = c(as.numeric(scale(training_dat$mean_snow_water_equiv_mm)),scl_fut_swe), # snow depth, plus forecast years
              npreds  = nrow(training_dat)+nrow(validation_dat)) # number of total predictions (obs + forecast)

##  Random variables to collect
out_variables <- c("r", "b", "b1", "sigma_proc", "z")

##  Send to JAGS
mc3     <- jags.model(file=textConnection(my_model), data=mydat, n.chains=3, n.adapt = 50000) 
update(mc3, n.iter = 100000) 
mc3.out <- coda.samples(model=mc3, 
                        variable.names=out_variables, 
                        n.iter=100000) 

##  Split MCMC output for file constraints
saveRDS(mc3.out[[1]],"../results/swe_avg_posteriors_chain1.RDS")
saveRDS(mc3.out[[2]],"../results/swe_avg_posteriors_chain2.RDS")
saveRDS(mc3.out[[3]],"../results/swe_avg_posteriors_chain3.RDS")

