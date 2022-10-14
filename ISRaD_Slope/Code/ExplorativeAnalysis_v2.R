# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 01/06/2022 #

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-10-05"))

lyr_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

lyr_data <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(ClimateZone = case_when(
    entry_name == "Gentsch_2018" ~ "polar",
    pro_usda_soil_order == "Gelisols" ~ "polar",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold",
    str_detect(pro_KG_present_long, "Polar") ~ "polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  mutate(ClimateZoneAnd = case_when(
    entry_name == "Gentsch_2018" ~ "polar",
    pro_usda_soil_order == "Gelisols" ~ "polar",
    pro_usda_soil_order == "Andisols" ~ "andisols",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold",
    str_detect(pro_KG_present_long, "Polar") ~ "polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  #remove for now: need to fix depth
  filter(entry_name != "Fernandez_1993a")

lyr_data %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

names(lyr_data)

## Check climate
lyr_data %>% 
  filter(pro_country == "Russia") %>% 
  count(entry_name, ClimateZone, pro_usda_soil_order)

lyr_data %>% 
  filter(pro_usda_soil_order == "Gelisols") %>% 
  count(entry_name, ClimateZone, pro_country, pro_usda_soil_order)

## Mapping sampling locations ##

world <- map_data("world") %>% 
  filter(region != "Antarctica")

library(raster)
climate_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/Beck_KG_V1/Beck_KG_V1_present_0p083.tif"
climate_raster <- raster::raster(climate_dir)

plot(climate_raster)

#reclassify climate raster
recal <- c(0,3,1, 3,7,2, 7,16,3, 16,28,4, 28,30,5)
recal_mat <- matrix(recal, ncol = 3, byrow = TRUE)

climate_grp <- reclassify(climate_raster, recal_mat)
plot(climate_grp)

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
climate_pts <- rasterToPoints(climate_grp, spatial = TRUE)
# Then to a 'conventional' dataframe
climate_df  <- data.frame(climate_pts) %>%
  drop_na() %>% 
  rename(ClimateZone = Beck_KG_V1_present_0p083) %>% 
  filter(ClimateZone != 0) %>% 
  mutate(ClimateZone = case_when(
    ClimateZone == 1 ~ "tropical",
    ClimateZone == 2 ~ "arid",
    ClimateZone == 3 ~ "temperate",
    ClimateZone == 4 ~ "cold",
    ClimateZone == 5 ~ "polar"
  ))

rm(climate_pts)

summary(climate_df)

climate_df$ClimateZone <- factor(climate_df$ClimateZone,
                                 levels = c("polar", "cold", "arid",
                                            "temperate", "tropical"))

climate_color <- c("#e66101", "#fdb863", "#f6e8c3", "#b2abd2", "#5e3c99")

ggplot() +
  geom_raster(data = climate_df, 
              aes(x = x, y = y, fill = ClimateZone)) +
  geom_point(data = lyr_data,
             aes(x = pro_long, y = pro_lat),
             fill = "#252525", size = 1.5, shape = 21, color = "white") +
  theme_bw(base_size = 12) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = c(0.11,0.29),
        legend.title = element_text(size = 10)) +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), expand = c(0,0),
                     breaks = c(-100,0,100), limits = c(-170,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), expand = c(0,0),
                     breaks = c(-50,0,50), limits = c(-55,80)) +
  scale_fill_manual("Climate zones", values = climate_color)

ggsave(file = paste0("./Figure/ISRaD_14C_Climate_map_", Sys.Date(),
                     ".jpeg"), width = 6, height = 3.5)
  
ggplot() +  
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "#d9d9d9")  +
  geom_point(data = lyr_data, 
             aes(x = pro_long, y = pro_lat),
             fill = "#116656", size = 2, shape = 21, color = "white") +
  theme_bw(base_size = 16) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")) +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), 
                     breaks = c(-100,0,100), limits = c(-160,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), 
                     breaks = c(-50,0,50), limits = c(-55,80))
ggsave(file = paste0("./Figure/ISRaD_14C_map_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

lyr_data %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

### Apply mspline function

## mspline 14C
lyr_data_mpspline_14c <- lyr_data %>% 
  filter(entry_name != "Fernandez_1993a") %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 0.5)

## mspline CORG
lyr_data_mpspline_c <- lyr_data %>% 
  filter(entry_name != "Fernandez_1993a") %>%  
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 60, lam = 0.5)

## 14C and SOC
mspline_14c_c <- lyr_data_mpspline_14c$est_1cm %>% 
  rename(lyr_14c_msp = SPLINED_VALUE) %>% 
  full_join(lyr_data_mpspline_c$est_1cm %>% 
              rename(CORG_msp = SPLINED_VALUE)) %>% 
  filter(LD < 101) %>% 
  tibble()

mspline_14c_c_all <- mspline_14c_c %>%
  dplyr::left_join(lyr_data %>% 
                     filter(entry_name != "Fernandez_1993a") %>%  
                     distinct(id, .keep_all = TRUE), 
                   by = "id") %>% 
  group_by(id) %>% 
  arrange(UD) %>% 
  ungroup() 

## Data exploration ##

# Data distribution
lyr_data %>% 
  group_by(id) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  ggplot(aes(x = pro_usda_soil_order, fill = site_name)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) +
  facet_wrap(~ClimateZone) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  scale_x_discrete("") +
  scale_y_continuous("Number of profiles", expand = c(0,0), limits = c(0,65),
                     breaks = seq(0,65,20))
ggsave(file = paste0("./Figure/ISRaD_climate_soiltype_dis_usda_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

lyr_data %>% 
  group_by(ClimateZone) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

## Raw data points: SOC vs 14C
lyr_data %>% 
  filter(depth < 101) %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.005,60)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,25)) +
  coord_cartesian(ylim = c(-1000,400)) 

ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

lyr_data %>% 
  filter(depth < 101) %>% 
  filter(ClimateZone == "polar") %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", size = 1),
        strip.background = element_rect(fill = "#e66101"),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  scale_x_continuous("", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,25)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_polar_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)

lyr_data %>% 
  filter(depth < 101) %>% 
  filter(ClimateZone == "cold") %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", size = 1),
        strip.background = element_rect(fill = "#fdb863"),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,25)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_cold_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)

lyr_data %>% 
  filter(depth < 101) %>% 
  filter(ClimateZone == "arid") %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", size = 1),
        strip.background = element_rect(fill = "#f6e8c3"),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,25)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_arid_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)

lyr_data %>% 
  filter(depth < 101) %>% 
  filter(ClimateZone == "temperate") %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", size = 1),
        strip.background = element_rect(fill = "#b2abd2"),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,25), breaks = c(0,10,20)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_temp_", Sys.Date(),
                     ".jpeg"), width = 5, height = 3.9)

lyr_data %>% 
  filter(depth < 101) %>% 
  filter(ClimateZone == "tropical") %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", size = 1),
        strip.background = element_rect(fill = "#5e3c99"),
        strip.text = element_text(face = "bold", color = "white"),
        legend.position = "none") +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,25)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_trop_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)

## Raw data + averaged data: SOC and climate
climate_all <- mspline_14c_c_all %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

ggplot() + 
  geom_hex(data = lyr_data %>% 
             filter(depth < 101), 
           aes(x = CORG, y = lyr_14c), alpha = 0.7) +
  geom_path(data = climate_all %>% 
              filter(n > 4 & n_rel > 33),
            aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(data = climate_all %>% 
                  filter(n > 4 & n_rel > 33),
                aes(ymin = lci_14c, ymax = uci_14c, x = median_c, color = ClimateZoneAnd), 
                alpha = 0.3) +
  geom_errorbarh(data = climate_all %>% 
                  filter(n > 4 & n_rel > 33),
                 aes(xmin = lci_c, xmax = uci_c, y = median_14c, color = ClimateZoneAnd), 
                 alpha = 0.3) +
  facet_wrap(~ClimateZoneAnd) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.005,60)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1) +
  coord_cartesian(ylim = c(-1000,400)) 


### Add HLZ data
library(sf)
library(tmap)
HLZ_directory <- "D:/Sophie/PhD/AfSIS_GlobalData/HLZ/holdrid/holdrid.shp"

HLZ <- sf::st_read(HLZ_directory)
st_crs(HLZ) <- 4326

lyr_data_sf <- sf::st_as_sf(lyr_data, 
                            coords = c("pro_long", "pro_lat"), 
                            crs = 4326)

tmap_mode("view")
tm_shape(HLZ, projection = 4326) +
  tm_polygons(col = "DESC") +
  tm_shape(lyr_data_sf) +
  tm_dots(popup.vars = "id")

lyr_data_HLZ_sf <- sf::st_join(lyr_data_sf, HLZ, 
                               left = TRUE)

lyr_data_HLZ_NA <- lyr_data_HLZ_sf %>% 
  tibble() %>% 
  dplyr::select(-c(geometry, AREA:FREQUENCY, SYMBOL)) %>% 
  rename(HLZ_Zone = DESC)

names(lyr_data_HLZ_NA)

lyr_data_HLZ_NA %>% 
  filter(is.na(HLZ_Zone)) %>% 
  count(id) %>% view()

# lyr_data_HLZ_NA %>% 
#   filter(entry_name == "Kramer_2012") %>% 
#   count(id, HLZ_Zone)
# 
# lyr_data_HLZ_NA %>% 
#   filter(entry_name == "Torn_1997") %>% 
#   filter(is.na(HLZ_Zone)) %>% 
#   count(id)
# 
# tmap_mode("view")
# tm_shape(HLZ, projection = 4326) +
#   tm_polygons(col = "DESC") +
#   tm_shape(lyr_data_sf[(lyr_data_sf$entry_name == "Crow_2015"|
#                           lyr_data_sf$entry_name == "Cusack_2012"|
#                           lyr_data_sf$entry_name == "Kramer_2012"|
#                           lyr_data_sf$entry_name == "Torn_1997"),]) +
#   tm_dots(popup.vars = "id")

#Manually assign missing HLZ
lyr_data_HLZ <- lyr_data_HLZ_NA %>% 
  mutate(HLZ_Zone = ifelse(entry_name == "Basile_Doelsch_2005", 
                            "Subtropical dry forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Heckman_2018_CA_Mollisol_SCT2", 
                           "Warm temperate dry forest", HLZ_Zone)) %>%
  mutate(HLZ_Zone = ifelse(id == "Heckman_2018_MI_Spodosol_MiB1", 
                           "Cool temperate moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Lassey_1996_Westland_Hokotika", 
                           "Cool temperate wet forest", HLZ_Zone)) %>% 
  #Could also be "Cool temperate wet forest"
  mutate(HLZ_Zone = ifelse(id == "Lassey_1996_Hawera_Whareroa road", 
                           "Warm temperate moist forest", HLZ_Zone)) %>%
  #Could also be "Cool temperate steppe"
  mutate(HLZ_Zone = ifelse(grepl("Lawrence_2021_Santa Cruz_SC", id),
                           "Warm temperate dry forest", HLZ_Zone)) %>%
  mutate(HLZ_Zone = ifelse(grepl("McFarlane_2013_MI-Coarse UMBS", id),
                           "Cool temperate moist forest", HLZ_Zone)) %>%
  mutate(HLZ_Zone = ifelse(grepl("Sanaiotti_2002_Alter do Chão_Alter do Chão", id),
                           "Subtropical moist forest", HLZ_Zone)) %>%
  #Could also be "Subtropical moist forest" or "Subtropical dry forest"
  mutate(HLZ_Zone = ifelse(entry_name == "Schwartz_1992", 
                           "Tropical dry forest", HLZ_Zone)) %>%
  #Most of Hawaii has missing data
  mutate(HLZ_Zone = ifelse(grepl("Crow_2015_AND_EUC_Composite", id),
                           "Warm temperate wet forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(grepl("Crow_2015_MOL_AG_Composite", id), 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(grepl("Crow_2015_OX_AG_Composite", id), 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(entry_name == "Cusack_2012", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(grepl("Grant_2022_Kohala", id), 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Kramer_2012_Pololu_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_1997_Amalu-precipitation_Amalu-precipitation_profile_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_1997_Kohala (150ky)_Kohala (150ky)_profile_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_1997_Kokee (4.1my)_Kokee (4.1my)_profile_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_1997_Kolekole (1.4my)_Kolekole (1.4my)_profile_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_2005_Kohala_Kohala_150", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_2005_Kokee, Kauai_Kokee_Kauai_4100", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_2005_Kolekole, Molokai_Kolekole_Molokai_1400", 
                           "Subtropical moist forest", HLZ_Zone)) 

lyr_data_HLZ %>% 
  filter(is.na(HLZ_Zone)) %>% 
  count(id) %>% view()

lyr_data_HLZ %>% 
  count(HLZ_Zone)

lyr_data_HLZ %>% 
  filter(depth <= 200) %>%
  ggplot(aes(x = depth, y = lyr_14c, fill = pro_BIO12_mmyr_WC2.1)) + 
  geom_point(shape = 21, size = 4, alpha = 0.7) +
  facet_wrap(~HLZ_Zone) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Delta14C") +
  scale_fill_viridis_c("MAP [mm]", trans = "log10", direction = -1) +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2))

lyr_data_HLZ %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) +
  geom_point(aes(color = pro_BIO12_mmyr_WC2.1), size = 3) +
  # geom_line(aes(group = id, color = pro_BIO12_mmyr_WC2.1), orientation = "y", 
  #           alpha = 0.7) +
  facet_wrap(~HLZ_Zone) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  scale_x_continuous("SOC [%]", trans = "log10") +
  # geom_smooth(orientation = "y", method = "gam") +
  scale_color_viridis_c("MAP [mm]", trans = "log10", direction = -1) +
  guides(color = guide_colorbar(barheight = 10, frame.colour = "black", 
                                ticks.linewidth = 2))

