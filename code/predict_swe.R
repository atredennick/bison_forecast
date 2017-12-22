################################################################################
##  predict_swe.R: R script to predict SWE at West Yellowstone SNOTEL using the
##  algorithm from Tercek et al. 2016, PLoS One.
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##                Mike Tercek
##  Date created: December 1, 2017
################################################################################

##  Clear everything...
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location
root      <- "~/Repos/bison_forecast/"
cmip_dir  <- paste0(root, "data/tercek_data/bias_corrected_daily_t_and_p")
calib_dir <- paste0(root, "data/tercek_data/")
setwd(paste0(root,"code/"))



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)



####
####  DEFINE FUNCTIONS FOR PREDICTING SWE --------------------------------------
####
sim_snow_for_one_water_year <- function(tavg, precip, melt_factor, melt_thresh_temperature){
  todays_snow <- 0 # start on October 1
  snow_vector <- numeric(length(tavg))
  for(i in 1:length(tavg)){
    if(is.na(tavg[i]) == FALSE){
      todays_snow <- melt_one_day(todays_snow,tavg[i],melt_factor, melt_thresh_temperature)
      todays_snow <- accum_one_day(todays_snow,tavg[i],precip[i],melt_factor, melt_thresh_temperature)
      snow_list[i] <- todays_snow 
    }
  }
  return(snow_vector)
}

melt_one_day <- function(start_swe, tavg, melt_factor, melt_thresh_temperature){
  if(tavg <= melt_thresh_temperature){
    return(start_swe)
  }else{
    swe_delta <-  (tavg - melt_thresh_temperature) * melt_factor
    end_swe <- max(0, start_swe - swe_delta)
    return(end_swe)
  }
}

accum_one_day <- function(start_swe, tavg, precip, precip_fraction, melt_thresh_temperature){
  if(tavg <= melt_thresh_temperature){
    rain <- 0.0
  }else if(tavg > (melt_thresh_temperature + 6)){
    rain <- 1.0
  }else{
    rain <- precip_fraction * (tavg - melt_thresh_temperature)
  }
  
  snow_increment <- max(0, (1.0 - rain) * precip)
  end_swe <- start_swe + snow_increment
  return(end_swe)
}



####
####  READ IN CMIP5 WEATHER PROJECTIONS ----------------------------------------
####
weather_files <- as.data.frame(list.files(cmip_dir)) %>%
  dplyr::rename(fname = `list.files(cmip_dir)`) %>%
  dplyr::mutate(fnamesep = fname) %>%
  separate(fnamesep, into = c("source", "variable", "rcp", "model", "biassep"), sep = "_") %>%
  separate(biassep, into = c("bias", "filetype"), sep = "[.]")



