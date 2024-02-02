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

## Randomly selected 10,00 sampling points globally
set.seed(42)
sf_use_s2(FALSE)

rdm_p <- sampleRandom(MAT_global, size = 10000, sp = TRUE)

plot(rdm_p)

##Extract data from global data points
# MAP
MAP_rdm <- raster::extract(MAP, rdm_p, df = TRUE) %>% 
  dplyr::rename(MAP = wc2.0_bio_30s_12)

# PET
PET_rdm <- raster::extract(PET, rdm_p, df = TRUE) %>% 
  dplyr::rename(PET = OBJECTID)

# Climate zones
CZ_rdm <- raster::extract(KG_p_raster, rdm_p, df = TRUE) %>% 
  dplyr::rename(pro_KG_present = Beck_KG_V1_present_0p0083)

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

##Merge data
rdm_data <- as.data.frame(rdm_p) %>% 
  dplyr::rename(MAT = wc2.0_bio_30s_01,
                Longitude = x,
                Latitude = y) %>% 
  cbind(MAP_rdm[-1], PET_rdm[-1], CZ_rdm_clean) %>% 
  mutate(dataset = "global") %>% 
  tibble() 

head(rdm_data)

lyr_mod <- lyr_data %>% 
  dplyr::select(ClimateZone, pro_MAT_mod, pro_MAP_mod, pro_PET_mm_yr_mod,
                pro_long, pro_lat) %>% 
  dplyr::rename(MAT = pro_MAT_mod,
                MAP = pro_MAP_mod,
                PET = pro_PET_mm_yr_mod,
                Longitude = pro_long,
                Latitude = pro_lat) %>% 
  mutate(dataset = "ISRaD")

climate_all <- rdm_data %>% 
  rbind(lyr_mod) %>% 
  dplyr::rename("MAT [°C]" = MAT,
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
ggsave(file = paste0("./Figure/Figure_A1_", Sys.Date(),
                     ".jpeg"), width = 10, height = 6)


