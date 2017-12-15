################################################################################
##  bison_neural_net.R: R script to fit a neural network model to the bison
##  data and make forecasts. Just for fun, for now.
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: December 15, 2017
################################################################################

##  Clear everything
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
library(forecast)  # Has neural network functions, and other forecasting fxns
source("./utilities/plotting_theme.R") # Source my plotting theme



####
####  LOAD DATA ----------------------------------------------------------------
####
snow_ynp  <- read.csv("../data/west_yellowstone_snotel_summary.csv", row.names = 1) 
bison_raw <- read.csv("../data/YNP_bison_population_size.csv")
bison_dat <- bison_raw %>% 
  dplyr::select(-count.source, -removal.source) %>%     # drop the source column
  mutate(set = ifelse(year < 2011, "training", "validation")) %>% # make new column for data splits
  left_join(snow_ynp, by="year") # merge in SNOTEL data



####
####  FIT AND FORECAST ---------------------------------------------------------
####
bison_training <- filter(bison_dat, set == "training")
y <- pull(bison_training, count.mean)
y_ts <- ts(data = y, 
           start = min(bison_training$year), 
           end = max(bison_training$year), 
           frequency = 1)

X <- cbind(pull(bison_training, accum_snow_water_equiv_mm),
           pull(bison_training, wint.removal))
X[is.na(X[,2]==TRUE),2] <- 0

##  Fit the model to training data
fit <- nnetar(y_ts, xreg = X)

##  Forecast with known SWE
bison_testing <- filter(bison_dat, set == "validation")
Xcast <- cbind(pull(bison_testing, accum_snow_water_equiv_mm),
               pull(bison_testing, wint.removal))
fcast <- forecast(fit, xreg = Xcast, PI = TRUE)
plot(fcast)
lines(bison_testing$year, bison_testing$count.mean, col = "red")
