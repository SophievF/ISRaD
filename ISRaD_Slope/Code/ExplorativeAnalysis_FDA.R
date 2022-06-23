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

## Functional PCA ##
#https://cran.r-project.org/web/packages/fdapace/index.html

library(fdapace)

lyr_list <- lyr_data %>%
  #remove studies that have multiple depth layers for now
  filter(entry_name != "Drake_2019" ,
         entry_name != "Richer_1999",
         entry_name != "Giardina_2014") %>% 
  #overlapping depth layers; not enough depth layers
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
  dplyr::select(id, depth, lyr_14c) %>% 
  group_by(id) %>% 
  summarise(across(, list))

FPCAsparse <- FPCA(lyr_list$lyr_14c, lyr_list$depth, list(plot = TRUE))

par(mfrow=c(1,2))
CreatePathPlot(FPCAsparse, subset = c(3,5,135), main = 'K = 11', pch = 4); grid()
CreatePathPlot(FPCAsparse, subset = c(3,5,135), K = 3, main = 'K = 3', pch = 4) ; grid()

par(mfrow=c(1,1))
CreateOutliersPlot(FPCAsparse, optns = list(K = 3, variant = 'KDE'))

CreateFuncBoxPlot(FPCAsparse, xlab = 'CORG', ylab = 'Delta14C', 
                  optns = list(K =3, variant='pointwise'))



