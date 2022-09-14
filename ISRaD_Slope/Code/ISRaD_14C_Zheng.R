# Explore 14C profiles in ISRaD #
# Relationship between 14C and SOC #
# Sophie von Fromm #
# 15/08/2022 #

# Extract data from Zheng and compare with ISRaD profiles #
library(tidyverse)
library(ncdf4)
library(raster)

# Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-09-13"))

# Filter data for mspline function
lyr_mpspline <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  ))

### Global 14C data (Zheng et al. 2020)
# Data source: https://zenodo.org/record/3823612#.X6lGvVNKjUI

lyr_sf <- lyr_mpspline %>% 
  dplyr::select(id, entry_name, site_name, pro_long, pro_lat) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  sf::st_as_sf(coords = c("pro_long", "pro_lat"), crs = 4326)

Global14C_raster <- raster::brick("D:/Sophie/PhD/AfSIS_GlobalData/ZhengGlobal14C/global_delta_14C.nc")
plot(Global14C_raster$X0)

# Select random points globally
rdm_p <- sampleRandom(Global14C_raster, size = 1000, sp = TRUE)

plot(rdm_p)

Global14C <- raster::extract(Global14C_raster, rdm_p, df = TRUE)

rdm_14c <- cbind(coordinates(rdm_p), Global14C) %>% 
  pivot_longer(cols = c(X0:X99), names_to = "depth_X", values_to = "lyr_14c") %>%
  separate(depth_X, into = c("X", "depth"),
           sep = "(?<=[X])(?=[0-9])") %>% 
  mutate(depth = as.numeric(depth)) %>% 
  dplyr::select(-X) %>% 
  rename(Longitude = x, Latitude = y)

summary(rdm_14c)

world <- map_data("world") %>% 
  filter(region != "Antarctica")

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgrey")  +
  geom_point(data = rdm_14c, 
             aes(x = Longitude, y = Latitude),
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

rdm_14c %>% 
  ggplot(aes(x = lyr_14c, y = depth, group = factor(ID), 
             color = factor(ID))) +
  geom_line(orientation = "y", size = 1) +
  theme_bw(base_size = 16) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous(position = "top") +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "magma")

# Extract other global products
# Load present climate zones

KG_p_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/Beck_KG_V1/Beck_KG_V1_present_0p0083.tif"
KG_p_raster <- raster::raster(KG_p_dir)
pro_KG_present <- raster::extract(KG_p_raster, cbind(rdm_14c$Longitude,
                                                     rdm_14c$Latitude))
KG_p_legend <- read_csv("D:/Sophie/PhD/AfSIS_GlobalData/Beck_KG_V1/KG_present_legend.csv")

KG_p_long <- cbind(KG_p, KG_p_legend)

#Merge both data sets
rdm_14c_KG <- cbind(rdm_14c, pro_KG_present) %>% 
  tibble() %>% 
  left_join(KG_p_legend, by = "pro_KG_present") %>% 
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  ))

summary(rdm_14c_KG)

global_rdm_14c_climate <- rdm_14c_KG %>% 
  drop_na(lyr_14c) %>%
  group_by(ClimateZone, depth) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2]) %>% 
  ungroup()

global_rdm_14c_climate %>%  
  dplyr::select(-c(ID, lyr_14c)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(depth) %>% 
  ggplot(aes(x = median_14c, y = depth, color = ClimateZone)) +
  geom_path() +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClimateZone),
              alpha = 0.3) +
  theme_classic(base_size = 16) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), position = "top",
                     expand = c(0,0), limits = c(-750,250)) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        legend.position = c(0.2,0.8))
ggsave(file = paste0("./Figure/Global_rdm_14C_depth_climate_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

## Extract data for each profile in ISRaD
Global14C <- raster::extract(Global14C_raster, lyr_sf,
                             #add buffer ~ 10000 m and return median value
                             #one grid ~ 55x55km; 55+1/2diagonal; 55 + 55*sqrt(2)/2
                             buffer = 100000, fun = median)

summary(Global14C)

# Combine both datasets
lyr_14c_global <- cbind(lyr_mpspline %>% 
                          distinct(id, .keep_all = TRUE), Global14C) %>% 
  tibble() %>% 
  dplyr::select(id, entry_name, site_name, pro_name, pro_long, pro_lat, c(X0:X99),
                pro_usda_soil_order, ClimateZone) %>% 
  pivot_longer(cols = c(X0:X99), names_to = "depth_X", values_to = "lyr_14c") %>%
  separate(depth_X, into = c("X", "depth"),
           sep = "(?<=[X])(?=[0-9])") %>% 
  mutate(depth = as.numeric(depth)) %>% 
  dplyr::select(-X)

summary(lyr_14c_global)

lyr_14c_global %>% 
  filter(is.na(lyr_14c)) %>% 
  count(id) #4 profiles with NA

lyr_14c_global %>% 
  drop_na(lyr_14c) %>% 
  ggplot(aes(x = lyr_14c, y = depth, group = id, color = id)) +
  geom_line(orientation = "y", size = 1) +
  theme_bw(base_size = 16) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous(position = "top") +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "magma")

global_14c_climate <- lyr_14c_global %>% 
  drop_na(lyr_14c) %>%
  group_by(ClimateZone, depth) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup()

global_14c_climate %>%  
  dplyr::select(-c(id, lyr_14c)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(depth) %>% 
  ggplot(aes(x = median_14c, y = depth, color = ClimateZone)) +
  geom_path() +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClimateZone),
              alpha = 0.3) +
  theme_classic(base_size = 16) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), position = "top",
                     expand = c(0,0), limits = c(-750,250)) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        legend.position = c(0.2,0.8))
ggsave(file = paste0("./Figure/Global_14C_depth_climate_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

ggplot() +
  geom_path(data = global_14c_climate %>%  
              dplyr::select(-c(id, lyr_14c)) %>% 
              distinct(median_14c, .keep_all = TRUE) %>%
              arrange(depth), 
            aes(x = median_14c, y = depth), color = "blue") +
  geom_ribbon(data = global_14c_climate %>%  
                dplyr::select(-c(id, lyr_14c)) %>% 
                distinct(median_14c, .keep_all = TRUE) %>%
                arrange(depth),
              aes(xmin = lci_14c, xmax = uci_14c, x = median_14c, y = depth),
              alpha = 0.3,  fill = "blue") +
  #Data comes from ExplorativeAnalysis_Depth.R
  geom_path(data = climate_14c %>% 
              arrange(UD) %>% 
              filter(n > 4 & n_rel > 60),
            aes(x = median_pseudo, y = UD), color = "red") +
  geom_ribbon(data = climate_14c %>% 
                arrange(UD) %>% 
                filter(n > 4 & n_rel > 60),
              aes(xmin = lci_median, xmax = uci_median, x = median_pseudo,
                  y = UD), alpha = 0.3, fill = "red") +
  theme_classic(base_size = 16) +
  facet_wrap(~ClimateZone) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), position = "top",
                     expand = c(0,0), limits = c(-750,250)) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2))
ggsave(file = paste0("./Figure/Global_ISRaD_14C_depth_climate_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

ggplot() +
  geom_path(data = global_rdm_14c_climate %>%  
              dplyr::select(-c(ID, lyr_14c)) %>% 
              distinct(median_14c, .keep_all = TRUE) %>%
              arrange(depth), 
            aes(x = median_14c, y = depth), color = "blue") +
  geom_ribbon(data = global_rdm_14c_climate %>%  
                dplyr::select(-c(ID, lyr_14c)) %>% 
                distinct(median_14c, .keep_all = TRUE) %>%
                arrange(depth),
              aes(xmin = lci_14c, xmax = uci_14c, x = median_14c, y = depth),
              alpha = 0.3,  fill = "blue") +
  #Data comes from ExplorativeAnalysis_Depth.R
  geom_path(data = climate_14c %>% 
              arrange(UD) %>% 
              filter(n > 4 & n_rel > 60),
            aes(x = median_pseudo, y = UD), color = "red") +
  geom_ribbon(data = climate_14c %>% 
                arrange(UD) %>% 
                filter(n > 4 & n_rel > 60),
              aes(xmin = lci_median, xmax = uci_median, x = median_pseudo,
                  y = UD), alpha = 0.3, fill = "red") +
  theme_classic(base_size = 16) +
  facet_wrap(~ClimateZone) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), position = "top",
                     expand = c(0,0), limits = c(-750,250)) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2))
ggsave(file = paste0("./Figure/Global_ISRaD_rdm_14C_depth_climate_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)
