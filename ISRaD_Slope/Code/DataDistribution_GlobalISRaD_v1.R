# Data distribution: Global vs ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 08/03/2023 #

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(raster)
library(sf)
library(tmap)
library(terra)

# Load filtered database data
lyr_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-03-09.csv") %>% 
  #remove layers that have CORG > 20
  filter(CORG_msp <= 20) %>% 
  mutate(MineralGroupsNew = case_when(
    MineralType == "low-activity clay" ~ "low-activity clay",
    MineralType == "amorphous" ~ "amorphous",
    pro_usda_soil_order == "Mollisols" ~ "high-activity clay",
    pro_usda_soil_order == "Spodosols" ~ "high-activity clay",
    pro_usda_soil_order == "Vertisols" ~ "high-activity clay",
    TRUE ~ "primary mineral"
  ))

# Boundaries World
data("World")
world_map <- World %>%
  st_transform(4326) %>% 
  dplyr::filter(continent != "Seven seas (open ocean)",
                continent != "Antarctica")

plot(world_map$geometry)

### Data distribution
## Load and crop data

# MAT
MAT_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_01.tif"
MAT <- raster(MAT_dir)
MAT_crop <- raster::crop(MAT, world_map)
MAT_global <- raster::mask(MAT_crop, world_map)

# MAP
MAP_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_12.tif"
MAP <- raster(MAP_dir)

# PET
PET_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/global-et0_annual.tif/et0_yr/et0_yr.tif"
PET <- raster(PET_dir)

## Climate zone
KG_p_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/Beck_KG_V1/Beck_KG_V1_present_0p0083.tif"
KG_p_raster <- raster::raster(KG_p_dir)

# calculate area
x <- terra::aggregate(rast(KG_p_dir), 100)
#not sure where the / 10000 comes from; https://gis.stackexchange.com/questions/433375/calculate-area-for-raster-in-r
a <- terra::cellSize(x, unit = "km") / 10000
b <- terra::resample(a, rast(KG_p_dir))

minmax(a)

z <- terra::zonal(b, rast(KG_p_dir), sum, na.rm = TRUE)

z %>%
  tibble() %>%
  rename(KG_p_present = Beck_KG_V1_present_0p0083) %>%
  filter(KG_p_present != 0) %>%
  mutate(area_km2 = area/1000000) %>%
  mutate(rel_dis = area_km2/sum(area_km2)*100)

climate_legend <- read.csv("D:/Sophie/PhD/AfSIS_GlobalData/Beck_KG_V1/KG_present_legend.csv")

ClimateZone_area <- z %>%
  tibble() %>%
  rename(pro_KG_present = Beck_KG_V1_present_0p0083) %>%
  filter(pro_KG_present != 0) %>%
  left_join(climate_legend) %>%
  mutate(ClimateZone = case_when(
    grepl("A", pro_KG_present_short) ~ "tropical",
    grepl("B", pro_KG_present_short) ~ "arid",
    grepl("C", pro_KG_present_short) ~ "warm temperate",
    grepl("D", pro_KG_present_short) ~ "cold temperate",
    grepl("E", pro_KG_present_short) ~ "tundra/polar"
  )) %>%
  dplyr::select(ClimateZone, area) %>%
  group_by(ClimateZone) %>%
  summarise(area_km2_sum = sum(area)) %>% 
  mutate(area_prop = area_km2_sum / sum(area_km2_sum))

x_grp <- x
x_grp_crop <- raster::crop(x_grp, world_map)
x_grp_global <- raster::mask(x_grp_crop, world_map)

x_grp_global[x_grp_global == 0] <- NA
#tropical
x_grp_global[x_grp_global <= 3] <- 1
#arid
x_grp_global[x_grp_global > 3 & x_grp_global <= 7] <- 2
#warm temperate
x_grp_global[x_grp_global > 7 & x_grp_global <= 16] <- 3
#cold temperate
x_grp_global[x_grp_global > 16 & x_grp_global <= 28] <- 4
#polar/tundra 
x_grp_global[x_grp_global > 28 & x_grp_global <= 30] <- 5

tm_shape(x_grp_global) +
  tm_raster(title = "Climate groups", style = "cat", 
            palette = "Paired", labels = c("tropical", "arid", "warm temperate",
                                             "cold temperate", "polar/tundra")) +
  tm_layout()

# # SOC (SoilGrids) v2
# # unit: 1g/kg = 0.1 %
# soc_5_dir <- "./Data/SoilGrids/soil_world/soc_0-5cm_mean_30s.tif" 
# soc_5_raster <- raster::raster(soc_5_dir)  
# 
# soc_15_dir <- "./Data/SoilGrids/soil_world/soc_5-15cm_mean_30s.tif" 
# soc_15_raster <- raster::raster(soc_15_dir) 
# 
# soc_30_dir <- "./Data/SoilGrids/soil_world/soc_15-30cm_mean_30s.tif" 
# soc_30_raster <- raster::raster(soc_30_dir) 
# 
# soc_60_dir <- "./Data/SoilGrids/soil_world/soc_30-60cm_mean_30s.tif" 
# soc_60_raster <- raster::raster(soc_60_dir) 

# SoilGrids v1
# SOC stocks (0-1m); unit: t/ha
socs_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/SoilGrids_v1/OCSTHA_M_100cm_250m_ll.tif"
socs_raster <- raster::raster(socs_dir)

## Randomly selected 10,00 sampling points globally
set.seed(42)
sf_use_s2(FALSE)

rdm_p <- sampleRandom(MAT_global, size = 10000, sp = TRUE)

plot(rdm_p)

##Extract data from global data points
# MAP
MAP_rdm <- raster::extract(MAP, rdm_p, df = TRUE) %>% 
  rename(MAP = wc2.0_bio_30s_12)

# PET
PET_rdm <- raster::extract(PET, rdm_p, df = TRUE) %>% 
  rename(PET = OBJECTID)

# Climate zones
CZ_rdm <- raster::extract(KG_p_raster, rdm_p, df = TRUE) %>% 
  rename(pro_KG_present = Beck_KG_V1_present_0p0083)

#Load climate legend
climate_legend <- read.csv("D:/Sophie/PhD/AfSIS_GlobalData/Beck_KG_V1/KG_present_legend.csv")

CZ_rdm_clean <- CZ_rdm %>% 
  left_join(climate_legend) %>% 
  mutate(ClimateZone = case_when(
    grepl("A", pro_KG_present_short) ~ "tropical",
    grepl("B", pro_KG_present_short) ~ "arid",
    grepl("C", pro_KG_present_short) ~ "warm temperate",
    grepl("D", pro_KG_present_short) ~ "cold temperate",
    grepl("E", pro_KG_present_short) ~ "tundra/polar"
  )) %>% 
  dplyr::select(ClimateZone)

# ## SoilGrids data v2
# # unit:  1g/kg = 0.1 %
# # SOC: 0-5 cm
# soc_5_rdm <- raster::extract(soc_5_raster, rdm_p, df = TRUE) 
# 
# # SOC: 5-15 cm
# soc_15_rdm <- raster::extract(soc_15_raster, rdm_p, df = TRUE) 
# 
# # SOC: 15-30 cm
# soc_30_rdm <- raster::extract(soc_30_raster, rdm_p, df = TRUE) 
# 
# # SOC: 30-60 cm
# soc_60_rdm <- raster::extract(soc_60_raster, rdm_p, df = TRUE)

## SoilGrids data v1
# SOC stocks (0-1m); unit: t/ha
socs_rdm <- raster::extract(socs_raster, rdm_p, df = TRUE)

##Merge data
rdm_data <- as.data.frame(rdm_p) %>% 
  rename(MAT = wc2.0_bio_30s_01,
         Longitude = x,
         Latitude = y) %>% 
  cbind(MAP_rdm[-1], PET_rdm[-1], CZ_rdm_clean, socs_rdm[-1]) %>% 
  mutate(dataset = "global") %>% 
  tibble() %>% 
  rename(SOC_t_ha_0_1m = OCSTHA_M_100cm_250m_ll)

head(rdm_data)

write.csv(rdm_data, file = paste0("./Data/Global_rdm_data_distribution_", 
                                     Sys.Date(), ".csv"),
          row.names = FALSE)

rdm_data <- read_csv("./Data/Global_rdm_data_distribution_2023-06-27.csv")

lyr_mod <- lyr_data %>% 
  dplyr::select(ClimateZone, pro_MAT_mod, pro_MAP_mod, pro_PET_mm_yr_mod,
                pro_long, pro_lat) %>% 
  rename(MAT = pro_MAT_mod,
         MAP = pro_MAP_mod,
         PET = pro_PET_mm_yr_mod,
         Longitude = pro_long,
         Latitude = pro_lat) %>% 
  mutate(dataset = "ISRaD")

climate_all <- rdm_data %>% 
  dplyr::select(c(1:6,8)) %>% 
  rbind(lyr_mod) %>% 
  rename("MAT [°C]" = MAT,
         "MAP [mm]" = MAP,
         "PET [mm]" = PET) %>% 
  tibble() %>% 
  pivot_longer(c("MAT [°C]", "MAP [mm]", "PET [mm]"), 
               names_to = "names", values_to = "values")

## Plot data
climate_all$ClimateZone <- factor(climate_all$ClimateZone,
                                  levels = c("tundra/polar", "cold temperate",
                                             "warm temperate", "tropical", "arid"))

climate_all %>% 
  ggplot(aes(fill = dataset, x = ClimateZone, y = values)) +
  geom_boxplot() +
  # geom_violin() +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.position = "top",
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~names, scales = "free", nrow = 1) +
  scale_x_discrete(labels = c("tundra/\npolar", "cold\ntemperate",
                              "warm\ntemperate", "tropical", "arid"))
ggsave(file = paste0("./Figure/ISRaD_msp_Global_climate_dis_box_", Sys.Date(),
                     ".jpeg"), width = 10, height = 6)

## Calculate mean for SOC for each climate zone
rdm_data %>% 
  group_by(ClimateZone) %>% 
  skimr::skim_without_charts(SOC_t_ha_0_1m)

rdm_data$ClimateZone <- factor(rdm_data$ClimateZone,
                               levels = c("tundra/polar", "cold temperate",
                                          "warm temperate", "tropical", "arid"))
# rdm_data %>% 
#   left_join(ClimateZone_area) %>% 
#   pivot_longer(cols = c(soc_0.5cm, soc_5.15cm, soc_15.30cm, soc_30.60cm),
#                names_to = "names", values_to = "values") %>% 
#   mutate(ClimateZone = factor(ClimateZone, 
#                               levels = c("tundra/polar", "cold temperate",
#                                          "warm temperate", "tropical", "arid"))) %>% 
#   mutate(names = factor(names, 
#                         levels = c("soc_0.5cm", "soc_5.15cm",
#                                    "soc_15.30cm", "soc_30.60cm"))) %>% 
#   ggplot(aes(fill = names, x = ClimateZone, y = values*area_prop)) +
#   geom_boxplot(notch = TRUE) +
#   theme_bw(base_size = 16) +
#   theme(axis.text = element_text(color = "black"),
#         panel.grid = element_blank(),
#         axis.title.x = element_blank(),
#         legend.position = "top",
#         strip.background = element_rect(fill = "white"),
#         axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_y_continuous("Soil organic carbon content [g/kg] / area [km²]", expand = c(0,0),
#                      limits = c(0,100)) +
#   facet_wrap(~names, scales = "free", nrow = 1) +
#   scale_x_discrete(labels = c("tundra/\npolar", "cold\ntemperate",
#                               "warm\ntemperate", "tropical", "arid"))
# ggsave(file = paste0("./Figure/SoilGrids_SOC_area_ClimateZone_box_", Sys.Date(),
#                      ".jpeg"), width = 14, height = 6)

rdm_data %>%
  left_join(ClimateZone_area) %>% 
  group_by(ClimateZone) %>% 
  summarise(mean_SOC_t_ha = mean(SOC_t_ha_0_1m, na.rm = TRUE),
            area_km2_sum = mean(area_km2_sum),
            area_prop = mean(area_prop)) %>% 
  mutate(amount_C_Pg = (mean_SOC_t_ha*100) * area_km2_sum*10^-9,
         prop_c = mean_SOC_t_ha * area_prop)
  

rdm_data %>%
  left_join(ClimateZone_area) %>%
  mutate(ClimateZone = factor(ClimateZone,
                              levels = c("tundra/polar", "cold temperate",
                                         "warm temperate", "tropical", "arid"))) %>%
  ggplot(aes(fill = ClimateZone, y = SOC_t_ha_0_1m, x = ClimateZone)) +
  geom_boxplot(notch = TRUE) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  scale_y_continuous("SOC stocks [t/ha]", expand = c(0,0), limits = c(0,3200),
                     breaks = seq(0,3000,500)) 
ggsave(file = paste0("./Figure/SoilGrids1_SOCs_ClimateZone_box_", Sys.Date(),
                     ".jpeg"), width = 8, height = 8)

