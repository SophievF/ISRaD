# Explore 14C profiles in ISRaD #
# Relationship between 14C and SOC with Alox/Feox #
# Sophie von Fromm #
# 15/08/2022 #

## Depth corrected values ##

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(mpspline2)

# Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-08-12"))

# Filter data for mspline function
lyr_mpspline <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
  #only use studies that have Alox/Feox
  drop_na(lyr_al_ox|lyr_fe_ox)