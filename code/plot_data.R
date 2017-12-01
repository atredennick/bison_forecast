################################################################################
##  plot_data.R: R script to plot the time series of YNP bison counts and
##  the West Yellowstone SNOTEL soil water equivalent (annual mean).
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: December 1, 2017
################################################################################

##  Clear everything
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location...only for RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



####
####  LOAD LIBRARIES ----
####
library(tidyverse) # Data science functions
library(dplyr)     # Data wrangling
library(ggthemes)  # Pleasing themes for ggplot2
library(cowplot)   # Combining ggplots
source("./utilities/plotting_theme.R") # Source my plotting theme



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
####  PLOT BISON AND SNOW DATA -------------------------------------------------
####
plot_data <- bison_dat %>%
  dplyr::select(year, set, count.mean, count.sd, mean_snow_water_equiv_mm) %>%
  dplyr::rename(avg_swe = mean_snow_water_equiv_mm)

docolor  <- "#278DAF"
altcolor <- "#CF4C26"
bison_plot <- ggplot(plot_data, aes(x = year, y = count.mean, color = set))+
  geom_line(alpha = 0.6)+
  geom_point(size=1.5)+
  geom_errorbar(aes(ymin = count.mean-count.sd, ymax = count.mean+count.sd), width=0.5, size=0.5)+
  scale_color_manual(values = c(docolor, altcolor))+
  scale_y_continuous(breaks = seq(0,6000,1000))+
  ylab("Number of bison")+
  xlab("Year")+
  my_theme+
  guides(color = FALSE)

snow_plot <- ggplot(plot_data, aes(x = year, y = avg_swe, color = set))+
  geom_line(alpha = 0.6)+
  geom_point(size=1.5)+
  scale_color_manual(values = c(docolor, altcolor))+
  ylab("Mean SWE (mm)")+
  xlab("Year")+
  my_theme+
  guides(color = FALSE)

the_plots <- list(bison_plot, snow_plot)
suppressWarnings( # Ignore warning about 6 NA rows for errorbars where sd not reported
  plot_grid(plotlist = the_plots, labels = "AUTO", ncol = length(the_plots))
)
ggsave(filename = "../figures/bison_data_plots.png", height = 3, width = 8, units = "in", dpi = 120)



####
####  OLD
####
# bison_growth_data <- bison_dat %>%
#   dplyr::select(year, count.mean) %>%
#   mutate(id = 1) %>% # constant id to work with ave()
#   mutate(growth_rate = ave(count.mean, id, FUN=function(x) c(0, diff(log(x)))))

# bison_growth <- ggplot(bison_growth_data, aes(x = year, y = growth_rate))+
#   geom_line(color = docolor, alpha = 0.6)+
#   geom_point(size=1.5, color = docolor)+
#   ylab("Population growth rate (r)")+
#   xlab("Year")+
#   my_theme
