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
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse) # Data science functions
library(dplyr)     # Data wrangling
library(rjags)     # Fitting Bayesian models with JAGS
library(coda)      # MCMC summaries
library(ggmcmc)    # MCMC-to-dataframe functions
# library(devtools) # For installing packages from GitHub
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR) # MCMC manipulation (by M. Dietze; modified by A. Tredennick)



####
####  LOAD DATA ----------------------------------------------------------------
####
snow_ynp  <- read.csv("../data/west_yellowstone_snotel_summary.csv", row.names = 1)
bison_raw <- read.csv("../data/YNP_bison_population_size.csv")
bison_dat <- bison_raw %>% 
  dplyr::select(-ends_with("source")) %>%     # drop the source column
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
    tau_proc   <- 1/sigma_proc^2
    eta        ~ dunif(0, 50)
    
    #### Fixed Effects Priors
    r  ~ dnorm(0.1, 1/0.02^2) # intrinsic growth rate, informed prior
    b  ~ dnorm(0,0.0001)      # strength of density dependence
    b1 ~ dnorm(0,0.0001)      # effect of snow
    
    #### Initial Conditions
    z[1]    ~ dnorm(Nobs[1], tau_obs[1]) # varies around observed abundance at t = 1
    zlog[1] <- log(z[1]) # set first zlog
    
    #### Process Model
    for(t in 2:npreds){
      # Gompertz growth, on log scale
      mu[t]   <- zlog[t-1] + r + b*zlog[t-1] + b1*x[t]
      zlog[t] ~ dnorm(mu[t], tau_proc)
      z[t]    <- exp(zlog[t]) # back transform to arithmetic scale
    }
    
    #### Data Model
    for(j in 2:n){
      p[j]     <- eta/(eta + z[j]) # calculate NB centrality parameter
      Nobs[j]  ~ dnegbin(p[j], eta) # NB likelihood
    }
    
    ####  Derived Quantities for Model Evaluation
    for(i in 1:n){
      # For autocorrelation test
      epsilon.obs[i] <- Nobs[i] - z[i]

      # Simulate new data
      p2[i]        <- eta/(eta + z[i])
      Nnew [i]     ~ dnegbin(p2[i], eta)
      sqerr[i]     <- (Nobs[i] - z[i])^2
      sqerr_new[i] <- (Nnew[i] - z[i])^2
    }
    fit     <- sum(sqerr[])
    fit.new <- sum(sqerr_new[])
    pvalue  <- step(fit.new-fit)

  }"



####
####  Fit Bison Forecasting Model ----------------------------------------------
####
##  For years without observation error, set to max observed standard deviation
##  TODO: Impute sigma_obs in the model?
na_sds                       <- which(is.na(bison_dat$count.sd)==T)
bison_dat[na_sds,"count.sd"] <- max(bison_dat$count.sd, na.rm=T)

##  Split into training and validation sets
training_dat   <- filter(bison_dat, set == "training")
validation_dat <- filter(bison_dat, set == "validation")

##  Set up SWE knowns (2011-2017), relative to scaling of observations
swe_mean     <- mean(training_dat$accum_snow_water_equiv_mm)
swe_sd       <- sd(training_dat$accum_snow_water_equiv_mm)
forecast_swe <- snow_ynp %>%
  filter(year %in% validation_dat$year) %>%
  pull(accum_snow_water_equiv_mm)
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
mydat <- list(Nobs    = round(training_dat$count.mean), # mean counts
              n       = nrow(training_dat), # number of observations
              tau_obs = 1/training_dat$count.sd^2, # transform s.d. to precision
              x       = c(as.numeric(scale(training_dat$accum_snow_water_equiv_mm)),scl_fut_swe), # snow depth, plus forecast years
              npreds  = nrow(training_dat)+nrow(validation_dat)) # number of total predictions (obs + forecast)

##  Random variables to collect
out_variables <- c("r", "b", "b1", "eta", "sigma_proc", "z", "pvalue", "fit", "fit.new")

##  Send to JAGS
mc3     <- jags.model(file = textConnection(my_model), 
                      data = mydat, 
                      n.chains = length(inits), 
                      n.adapt = 5000, 
                      inits = inits) 
           update(mc3, n.iter = 10000) 
mc3.out <- coda.samples(model=mc3, 
                        variable.names=out_variables, 
                        n.iter=10000) 

# ggs(mc3.out) %>%
#   filter(Parameter %in% c("r")) %>%
#   ggplot(aes(x=Iteration,y = value, color = as.factor(Chain)))+
#   geom_line()
# 
# ggs(mc3.out) %>%
#   filter(Parameter %in% c("fit", "fit.new")) %>%
#   spread(Parameter, value) %>%
#   ggplot(aes(x=fit, y=fit.new))+
#    geom_point()+
#    geom_abline(aes(intercept=0, slope=1), color="red")
# 
ggs(mc3.out) %>%
  filter(Parameter %in% c("r","b","b1")) %>%
  ggplot(aes(x=value))+
  geom_histogram()+
  facet_wrap(~Parameter, scales = "free")
# 
# outstats <- summary(mc3.out)$stat
# outquant <- summary(mc3.out)$quantile
# outstats[which(rownames(outstats)=="pvalue"),"Mean"]

##  Split MCMC output for file constraints
saveRDS(mc3.out[[1]],"../results/swe_est_posteriors_chain1.RDS")
saveRDS(mc3.out[[2]],"../results/swe_est_posteriors_chain2.RDS")
saveRDS(mc3.out[[3]],"../results/swe_est_posteriors_chain3.RDS")



####
####  FIT AND FORECAST ASSUMING AVG SWE
####
scl_fut_swe[] <- 0 # average is 0 by definition
##  Prepare data list
mydat <- list(Nobs    = round(training_dat$count.mean), # mean counts
              n       = nrow(training_dat), # number of observations
              tau_obs = 1/training_dat$count.sd^2, # transform s.d. to precision
              x       = c(as.numeric(scale(training_dat$accum_snow_water_equiv_mm)),scl_fut_swe), # snow depth, plus forecast years
              npreds  = nrow(training_dat)+nrow(validation_dat)) # number of total predictions (obs + forecast)

##  Random variables to collect
out_variables <- c("r", "b", "b1", "eta", "sigma_proc", "z", "pvalue", "fit", "fit.new")

##  Send to JAGS
mc3     <- jags.model(file = textConnection(my_model), 
                      data = mydat, 
                      n.chains = 3, 
                      n.adapt = 5000)
update(mc3, n.iter = 10000)
mc3.out <- coda.samples(model = mc3,
                        variable.names = out_variables,
                        n.iter = 10000)

##  Split MCMC output for file constraints
saveRDS(mc3.out[[1]],"../results/swe_avg_posteriors_chain1.RDS")
saveRDS(mc3.out[[2]],"../results/swe_avg_posteriors_chain2.RDS")
saveRDS(mc3.out[[3]],"../results/swe_avg_posteriors_chain3.RDS")



