################################################################################
##  format_weather_projections: R script to collate, reformat, and summarize
##  CMIP5 weather projections by model, run, and scenario.
##
##  ----------------------------------------------------------------------------
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: January 2, 2018
################################################################################

##  Clear the workspace
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # only for RStudio



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(dplyr)
library(stringr)



####
####  READ IN DATA AND ADD ROW COLUMN INFORMATION ------------------------------
####
col_names <- c("year",
               "month",
               "day",
               as.character(as.data.frame(read.delim("../data/CMIP_YNP/bcca5/COLS_pr.txt"))[,1])
               )
raw_precip_projs <- read_csv("../data/CMIP_YNP/bcca5/pr.csv", col_names = FALSE)

