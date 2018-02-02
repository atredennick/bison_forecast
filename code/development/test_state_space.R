################################################################################
##  test_state_space.R: R script to simulate population dynamics from a 
##  Gompertz population model and then use a Bayesian state-space model to
##  recover model parameters and make forecasts. The main purpose is to test
##  the effect of forecast length on parameters estimates, if any.
################################################################################

rm(list = ls(all.names = TRUE))



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse) # Data science functions
library(dplyr)     # Data wrangling
library(ggthemes)  # Pleasing themese for ggplot2
library(rjags)     # Fitting Bayesian models with JAGS
library(coda)      # MCMC summaries
library(ggmcmc)    # MCMC-to-dataframe functions



####
####  GOMPERTZ POPULATION MODEL FUNCTION ---------------------------------------
####
sim_gompertz <- function(sim_time, r, b0, z0, noise = 0.2){
  z <- numeric(length(sim_time))
  z[1] <- z0
  for(t in 2:sim_time){
    z[t] <- z[t-1] + r + b0*z[t-1] + rnorm(1,0,noise)
  }
  return(rpois(n = sim_time, lambda = exp(z)))
}



####
####  SIMULATE TIME SERIES -----------------------------------------------------
####
sim_time <- 40
r        <- 0.1
b0       <- -0.01
z0       <- log(500)

z_ts <- sim_gompertz(sim_time, r, b0, z0)



####
####  DEFINE JAGS STATE-SPACE MODEL --------------------------------------------
####
my_model <- "  
  model{

    #### Variance Priors
    sigma_proc ~ dunif(0,5)
    tau_proc   <- 1/sigma_proc^2
    eta        ~ dunif(0, 50)
    error_sd ~ dunif(0,10)
    error_prec <- 1/error_sd^2
    for(i in 1:npreds){
      error[i] ~ dnorm(0, error_prec)
    }
    
    #### Fixed Effects Priors
    r  ~ dnorm(0.1, 1/0.02^2)  # intrinsic growth rate, informed prior
    b0  ~ dnorm(0,1/2^2)I(-2,2) # strength of density dependence, bounded
    
    #### Initial Conditions
    z[1]    ~ dnorm(Nobs[1], 100) # varies around observed abundance at t = 1
    zlog[1] <- log(z[1]) # set first zlog
    
    #### Process Model
    for(t in 2:npreds){
      # Gompertz growth, on log scale
      mu[t]   <- zlog[t-1] + r + b0*zlog[t-1] + error[t]
      zlog[t] ~ dnorm(mu[t], tau_proc)
      z[t]    <- exp(zlog[t]) # back transform to arithmetic scale
    }
    
    #### Data Model
    for(j in 2:n){
      Nobs[j] ~ dpois(z[j]) # Poisson likelihood
    }

}"


####
####  FIT BAYESIAN STATE-SPACE MODEL -------------------------------------------
####
##  Set initial values for unkown parameters
inits <- list(
  list(sigma_proc = 0.01,
       r = 0.05,
       b0 = -0.001),
  list(sigma_proc = 0.3,
       r = 0.4,
       b0 = -0.1),
  list(sigma_proc = 0.1,
       r = 0.7,
       b0 = -0.00001)
)

##  Set up data list for JAGS
# nforecasts <- 30
out_df <- {}
for(nforecasts in c(0, 10, 20, 30, 40)){
  mydat <- list(Nobs    = z_ts, # mean counts
                n       = length(z_ts), # number of observations
                npreds  = length(z_ts) + nforecasts) # number of total predictions (obs + forecast)
  
  ##  Random variables to collect
  out_variables <- c("r", "b0", "z")
  
  ##  Send to JAGS
  mc3     <- jags.model(file = textConnection(my_model), 
                        data = mydat, 
                        n.chains = length(inits), 
                        n.adapt = 5000, 
                        inits = inits) 
  update(mc3, n.iter = 5000) 
  mc3.out <- coda.samples(model=mc3, 
                          variable.names=out_variables, 
                          n.iter=5000) 
  
  quants <- summary(mc3.out)$quantiles
  beta_quants <- quants[grep("b0", row.names(quants)), c(1,3,5)]
  
  beta_df <- data.frame(forecasts = nforecasts,
                        low_est = beta_quants[1],
                        med_est = beta_quants[2],
                        hi_est = beta_quants[3])
  
  out_df <- rbind(out_df, beta_df)
}



####
####  PLOT PARAMETER ESTIMATES -------------------------------------------------
####
ggplot(out_df, aes(x = forecasts, y = med_est))+
  geom_point()+
  geom_errorbar(aes(ymin = low_est, ymax = hi_est), width = 0.5)+
  geom_hline(aes(yintercept = b0), col = "red")+
  ylab("Density-dependence estimate")+
  xlab("Number of forecasts made")


####
####  PLOT THE FIT AND FORECAST ------------------------------------------------
####
# par(mfrow = c(1,3))
# zquants <- quants[grep("z", row.names(quants)),]
# plot(zquants[,3], type = "l")
# lines(zquants[,1], lty = 2)
# lines(zquants[,5], lty = 2)
# points(z_ts, pch=19)
# 
# rparam <- ggs(mc3.out) %>% filter(Parameter %in% c("r"))
# bparam <- ggs(mc3.out) %>% filter(Parameter %in% c("b0"))
# hist(rparam$value)
# hist(bparam$value)
