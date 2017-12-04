################################################################################
##  partition_forecast.R: R script to plot observations and posterior
##  predictions from the bison population model and then partition the forecast
##  variance following Dietze 2017, Ecological Applications
##  http://onlinelibrary.wiley.com/doi/10.1002/eap.1589/full
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: December, 4 2017
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
library(ggthemes)  # Pleasing themes for ggplot2
library(cowplot)   # Combining ggplots
library(rjags)     # Fitting Bayesian models with JAGS
library(coda)      # MCMC summaries
# library(devtools) # For installing packages from GitHub
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR) # MCMC manipulation (by M. Dietze)
source("./utilities/plotting_theme.R")



####
####  LOAD DATA AND MCMCs ------------------------------------------------------
####
##  Data
snow_ynp  <- read.csv("../data/west_yellowstone_snotel_summary.csv", row.names = 1) 
bison_raw <- read.csv("../data/YNP_bison_population_size.csv")
bison_dat <- bison_raw %>% 
  dplyr::select(-source) %>%     # drop the source column
  mutate(set = ifelse(year < 2011, "training", "validation")) %>% # make new column for data splits
  left_join(snow_ynp, by="year") # merge in SNOTEL data

##  MCMC
swe_est_posts <- as.mcmc.list(list(readRDS("../results/swe_est_posteriors_chain1.RDS"),
                                   readRDS("../results/swe_est_posteriors_chain2.RDS"),
                                   readRDS("../results/swe_est_posteriors_chain3.RDS")))
swe_avg_posts <- as.mcmc.list(list(readRDS("../results/swe_avg_posteriors_chain1.RDS"),
                                   readRDS("../results/swe_avg_posteriors_chain2.RDS"),
                                   readRDS("../results/swe_avg_posteriors_chain3.RDS")))



####
####  REFORMAT KNOWN SWE MCMC RESULTS ------------------------------------------
####
## Split output
out          <- list(params=NULL, predict=NULL)
mfit         <- as.matrix(swe_est_posts,chains=TRUE)
pred.cols    <- union(grep("z[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

## Collate predictions
predictions        <- rbind(fitted_model$predict[[1]],
                            fitted_model$predict[[2]],
                            fitted_model$predict[[3]])
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})
prediction_df      <- data.frame(year = bison_dat$year,
                                 set = bison_dat$set,
                                 observation = bison_dat$count.mean,
                                 upper_observation = bison_dat$count.mean+bison_dat$count.sd,
                                 lower_observation = bison_dat$count.mean-bison_dat$count.sd,
                                 median_prediction = median_predictions,
                                 upper_prediction = upper_predictions,
                                 lower_prediction = lower_predictions)



####
####  PLOT DATA AND POSTERIOR PREDICTIONS --------------------------------------
####
pred_color <- "#CF4C26"
obs_color  <- "#278DAF"
calibration_plot <- ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction),
              fill=pred_color, 
              color=NA, 
              alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_errorbar(aes(ymin=lower_observation, ymax=upper_observation), 
                width=0.5, 
                color=obs_color, 
                size=0.2)+
  geom_point(aes(y=observation), color=obs_color, size=0.5)+
  geom_vline(aes(xintercept=2010), linetype=2,color="grey55")+
  ylab("Number of bison")+
  xlab("Year")+
  my_theme



####
####  PARTITION FORECAST UNCERTAINTY -------------------------------------------
####
##  Function for the ecological process (Gompertz population growth)
iterate_process <- function(Nnow, xnow, r, b, b1, sd_proc) { 
  mu <- log(Nnow) + r + b*log(Nnow) + b1*xnow # determinstic process; log scale
  zlog <- rnorm(length(mu), mu, sd_proc) # stochastic process; log scale
  N <- exp(zlog) # back transform to arithmetic scale
}

##  Set up SWE knowns (2011-2017), relative to scaling of observations
training_dat <- filter(bison_dat, set == "training")
validation_dat <- filter(bison_dat, set == "validation")
swe_mean     <- mean(training_dat$accum_snow_water_equiv_mm)
swe_sd       <- sd(training_dat$accum_snow_water_equiv_mm)
forecast_swe <- snow_ynp %>%
  filter(year %in% validation_dat$year) %>%
  pull(accum_snow_water_equiv_mm)
scl_fut_swe  <- (forecast_swe - swe_mean) / swe_sd


##  Initial condition uncertainty: make forecasts from all MCMC iterations of
##  the final year, but use mean parameter values and no process error.
forecast_steps <- nrow(validation_dat)
num_iters      <- 1000
x              <- scl_fut_swe
z              <- sample(predictions[,nrow(training_dat)], num_iters, replace = TRUE)
param_summary  <- summary(fitted_model$params)$quantile
b              <- param_summary[1,3]
b1             <- param_summary[2,3]
r              <- param_summary[7,3]
sd_proc        <- param_summary[8,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[t], r = r, b = b, b1 = b1, sd_proc = 0)
  forecasts[,t] <- z
}
varI <- apply(forecasts,2,var)


##  Initial conditions and parameter uncertainty
x              <- scl_fut_swe
z              <- sample(predictions[,nrow(training_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
r              <- params[sample_params,"r"]
sd_proc        <- param_summary[3,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[t], r = r, b = b, b1 = b1, sd_proc = 0)
  forecasts[,t] <- z
}
varIP <- apply(forecasts,2,var)


##  Initial conditions, parameter, and driver uncertainty
xfun           <- function(x){rnorm(num_iters,x,sigma_x)}
sigma_x        <- 0.1
x              <- lapply(scl_fut_swe, FUN = xfun)
z              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
r              <- params[sample_params,"r"]
sd_proc        <- params[sample_params,"sigma_proc"]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[[t]], r = r, b = b, b1 = b1, sd_proc = 0)
  forecasts[,t] <- z
}
varIPD <- apply(forecasts,2,var)


##  Initial conditions, parameter, driver, and process uncertainty
sigma_x        <- c(0.1,0.2,0.4,0.8,1.6,2.4,4.8)
x              <- matrix(data = NA, ncol = length(scl_fut_swe), nrow = num_iters)
for(i in 1:ncol(x)){
  x[,i] <- rnorm(num_iters, scl_fut_swe[i], sigma_x[i])
}
z              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
r              <- params[sample_params,"r"]
sd_proc        <- params[sample_params,"sigma_proc"]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[[t]], r = r, b = b, b1 = b1, sd_proc = sd_proc)
  forecasts[,t] <- z
}
varIPDE <- apply(forecasts,2,var)


V.pred.sim     <- rbind(varIPDE,varIPD,varIP,varI)
V.pred.sim.rel <- apply(V.pred.sim,2,function(x) {x/max(x)})



####
####  PLOT THE FORECASTING UNCERTAINTY PARTITION -------------------------------
####
var_rel_preds <- as.data.frame(t(V.pred.sim.rel*100))
var_rel_preds$x <- 1:nrow(var_rel_preds)
my_cols <- c("#0A4D5B", "#139AB8", "#39B181","grey")
variance_plot <- ggplot(data=var_rel_preds, aes(x=x))+
  geom_ribbon(aes(ymin=0, ymax=varI), fill=my_cols[1])+
  geom_ribbon(aes(ymin=varI, ymax=varIP), fill=my_cols[2])+
  geom_ribbon(aes(ymin=varIP, ymax=varIPD), fill=my_cols[3])+
  geom_ribbon(aes(ymin=varIPD, ymax=varIPDE), fill=my_cols[4])+
  ylab("Percent of uncertainty")+
  xlab("Forecast steps")+
  scale_x_continuous(breaks=seq(1,forecast_steps,by=1), 
                     labels=paste(seq(1,forecast_steps,by=1), "yrs"))+
  scale_y_continuous(labels=paste0(seq(0,100,25),"%"))+
  my_theme



####
####  COMBINE PLOTS AND SAVE ---------------------------------------------------
####
suppressWarnings( # ignore warning about missing values, we know they are missing
  plot_grid(calibration_plot, variance_plot, nrow = 2, labels = "AUTO")
)
ggsave(filename = "../figures/bison_combined.png", 
       width = 4, 
       height = 6, 
       units = "in", 
       dpi =120)

