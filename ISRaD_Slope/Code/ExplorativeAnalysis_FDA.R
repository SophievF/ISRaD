# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 22/06/2022 #

## Functional data analysis (FDA) ##
#https://www.psych.mcgill.ca/misc/fda/ex-goods-a1.html

library(tidyverse)
library(fda)

#Load filtered lyr data
lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-06-21"))

##Depth~14C; SOC~14C; built fda object with both options; compare to mspline