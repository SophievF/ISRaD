### Download and access SOilGrids data
## https://rdrr.io/cran/geodata/man/soil_grids.html
## Sophie von Fromm
## 22-06-2023

library(geodata)
library(raster)
library(tidyverse)
library(tmap)
library(sf)

# https://geodata.ucdavis.edu/geodata/soil/soilgrids/
soil_world(var = "sand", depth = 5, stat = "mean", path = "./Data/SoilGrids/")

## Load raster files: sand and silt
sand_5_dir <- "./Data/SoilGrids/soil_world/sand_0-5cm_mean_30s.tif" 
sand_5 <- raster::raster(sand_5_dir)

silt_5_dir <- "./Data/SoilGrids/soil_world/silt_0-5cm_mean_30s.tif" 
silt_5 <- raster::raster(silt_5_dir)

# calculate clay content: clay = 100 - sand - silt
clay_5 <- 100 - sand_5 - silt_5
plot(clay_5)
writeRaster(clay_5, "./Data/SoilGrids/soil_world/clay_0-5cm_mean_30s.tif",
            overwrite = TRUE)

### Write David Rossiter about tif files ###

## Calculate SOC content / m² for each climate zone; try with 0-5 cm first
# unit: 1 dg/kg = 0.1 g/kg; 1g/kg = 0.1 %; 1dg/kg = 0.01 %
soc_5_dir <- "./Data/SoilGrids/soil_world/soc_0-5cm_mean_30s.tif" 
soc_5 <- raster::raster(soc_5_dir)  

plot(soc_5)

## Load climate zones
KG_p_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/Beck_KG_V1/Beck_KG_V1_present_0p0083.tif"
ClimateZone <- raster::raster(KG_p_dir)

plot(ClimateZone)

#Match extent
bb <- c(-180, 180, -60, 90)
ClimateZone <- setExtent(ClimateZone, bb, keepres = TRUE)
ClimateZone[ClimateZone == 0] <- NA

soc_climate <- raster::stack(soc_5, ClimateZone)

stackSave(soc_climate, "./Data/SoilGrids/soc_climate.tif")

## select random datapoints based on world map 
data("World")
world_map <- World %>%
  st_transform(4326) %>% 
  dplyr::filter(continent != "Seven seas (open ocean)",
                continent != "Antarctica")

world_map[c(2,4,5)]

set.seed(42)
sf_use_s2(FALSE)

rdm_p <- st_sample(world_map[2,4,5], size = 10)

plot(rdm_p)

# soc_climate_df <- sampleRandom(soc_climate, size=100, xy=TRUE, na.rm = FALSE) %>% 
#   as.data.frame() %>% 
#   rename(ClimateZone = Beck_KG_V1_present_0p0083,
#          Latitude = x,
#          Longitude = y)
# 
# soc_climate_df


# soc_climate_mean <- calc(soc_climate, fun = mean)


head(as.data.frame(soc_climate))


# data("World")
# world_map <- World %>%
#   st_transform(4326) %>%
#   dplyr::filter(continent != "Seven seas (open ocean)",
#                 continent != "Antarctica")
# 
# #Crop climate zones
# ClimateZone_crop <- raster::crop(ClimateZone, world_map)
# ClimateZone_global <- raster::mask(ClimateZone_crop, world_map)
# 
# soc_5_crop <- raster::crop(soc_5, world_map)
# soc_5_global <- raster::mask(soc_5_crop, world_map)
# 
# soc_climate <- raster::stack(soc_5_global, ClimateZone_global)
# 
# stackSave(soc_climate, "./Data/SoilGrids/soc_climate.tif")

