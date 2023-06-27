# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 10/02/2023 #

## Code and functions for modelling (SoilR) provided by Carlos Sierra

devtools::install_github('MPIBGC-TEE/SoilR-exp/pkg')

# library(tidyverse)
# library(ggpubr)
# library(RColorBrewer)
library(SoilR)

### Load filtered and splined lyr data
# lyr_msp_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-02-06.csv") %>% 
#   dplyr::select(entry_name, pro_name, pro_long, pro_lat, lyr_name, id,
#                 lyr_14c_msp, CORG_msp, pro_AI, pro_MAT_mod, 
#                 pro_GPP_Fluxcom_2001_2012_gC_m2d1, lyr_clay_mod,
#                 ClimateZoneAnd, MineralType, UD, LD, lyr_obs_date_y)
# 
# skimr::skim(lyr_msp_data)
# 
# lyr_msp_data %>% 
#   distinct(id, .keep_all = TRUE) %>% 
#   count(ClimateZoneAnd)
# 
# lyr_msp_data %>% 
#   distinct(id, .keep_all = TRUE) %>% 
#   count(MineralType)

### Set-up a depth-resolved model with SoilR (Carlos) ###
source("./Code/SoilR/VTLM.R")
source("./Code/SoilR/collect_by_layer.R")
source("./Code/SoilR/collect_by_pool.R")

## One-pool-model

# Step 1. Prepare what you need for a SoilR pool model
years <- seq(-48050,2019, by = 0.5)
LitterInput <- 700 
k1 <- 1/1000
initialC <- 0 # You can put any number here. They will be ignored because steady-state will be assumed
initialF14 <- 0 # Same here. For the depth simulation the age of the steady-state carbon will be used to set initial F14C
AtmF14 <- bind.C14curves(prebomb=IntCal13,postbomb=Hua2021$NHZone2,time.scale="AD")[,1:2] # Atmospheric radiocarbon curve for Northern Hemisphere zone 2 in Delta14C

# Step 2. Create SoilR pool model
OnePoolModel <- OnepModel14(t = years, 
                            k = k1, 
                            C0 = initialC,
                            F0_Delta14C = initialF14, 
                            In = LitterInput,
                            inputFc = AtmF14)

# Step 3. Create Vertical Transport model
depthLayers <- c(1:10)
rootInputs <- rep(10, times = length(depthLayers))
downwardTransferRate <- 0.005
upwardTransferRate <- 0

OnepDepth <- VTLM(Model = OnePoolModel, 
                  lyrs = depthLayers, 
                  latIn = rootInputs, 
                  d = downwardTransferRate, 
                  u = upwardTransferRate)

# Step 4. Solve the vertical transport model and aggregate by pool and layer
Cdt <- getC(OnepDepth)
C14dt <- getF14(OnepDepth)

# not needed for one pool model
# L14t <- collect_by_layer(Cdt, C14dt, nlayer = 10, npool = 1)
# P14t <- collect_by_pool(Cdt, C14dt, nlayer=10, npool = 1)

# Step 5. Plot the results
plot(AtmF14[5135:5694,], type="l", lty = 1)
matlines(OnepDepth@times, C14dt, lty = 1, col = 2:11) # Find a better color palette

matplot(OnepDepth@times, Cdt, type = "l", lty = 1, col = 2:11)

plot(tail(Cdt, 1), -depthLayers, xlab = "C stock", ylab = "Depth")
plot(tail(C14dt , 1), -depthLayers, xlab = "Delta14C", ylab = "Depth")

plot(tail(Cdt, 1), tail(C14dt, 1), xlab = "C stock", ylab = "Delta14C")


