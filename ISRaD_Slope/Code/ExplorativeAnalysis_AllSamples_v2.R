# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC: all profiles #
# Sophie von Fromm #
# 15/09/2022 #

library(tidyverse)
library(ggpubr)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-10-17"))

lyr_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

lyr_data <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  group_by(id) %>%
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
  )) 

lyr_data %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id),
            n_samples = n())

lyr_data$ClimateZone <- factor(lyr_data$ClimateZone,
                               levels = c("polar", "cold", "arid",
                                          "temperate", "tropical"))

lyr_data$ClimateZoneAnd <- factor(lyr_data$ClimateZoneAnd,
                                  levels = c("polar", "cold", "arid",
                                             "temperate", "tropical",
                                             "andisols"))

## Mapping sampling locations ##

world <- map_data("world") %>% 
  filter(region != "Antarctica")

library(raster)
climate_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/Beck_KG_V1/Beck_KG_V1_present_0p083.tif"
climate_raster <- raster::raster(climate_dir)

#reclassify climate raster
recal <- c(0,3,1, 3,7,2, 7,16,3, 16,28,4, 28,30,5)
recal_mat <- matrix(recal, ncol = 3, byrow = TRUE)

climate_grp <- reclassify(climate_raster, recal_mat)

## convert to a df for plotting in two steps,
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

ggsave(file = paste0("./Figure/ISRaD_14C_Climate_map_all_", Sys.Date(),
                     ".jpeg"), width = 6, height = 3.5)

### Density distribution
p1 <- lyr_data %>% 
  filter(depth < 101) %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(trans = "log10", direction = -1, limits = c(1,340)) +
  coord_cartesian(ylim = c(-1000,400))

p2 <- lyr_data %>% 
  filter(depth < 101) %>% 
  ggplot(aes(x = depth, y = lyr_14c)) + 
  geom_hex(binwidth = c(3,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(-5,110)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,400)) +
  scale_fill_viridis_c(direction = -1, trans = "log10", limits = c(1,340)) +
  coord_cartesian(xlim = c(0,101))

ggarrange(p1, p2, common.legend = TRUE)

ggsave(file = paste0("./Figure/ISRaD_14C_SOC_depth_hex_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

## Raw data points: SOC vs 14C
lyr_data %>%
  filter(depth < 101) %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  geom_point(data = lyr_data %>% 
               filter(depth < 101) %>% 
               group_by(ClimateZoneAnd) %>% 
               summarise(m_c = mean(CORG),
                         m_c14 = mean(lyr_14c)),
             aes(x = m_c, y = m_c14), color = "red", size = 3) +
  facet_wrap(~ClimateZoneAnd) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1) 
  coord_cartesian(ylim = c(-1000,400))

lyr_data %>%
  filter(depth < 101) %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1) +
  # geom_smooth(orientation = "y", method = "loess", color = "blue") +
  # geom_smooth(method = "loess", color = "red") +
  coord_cartesian(ylim = c(-1000,400)) 
  
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_all_", Sys.Date(),
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
  scale_fill_viridis_c(direction = -1, limits = c(0,35), breaks = c(0,10,20,30)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_polar_all_", Sys.Date(),
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
  scale_fill_viridis_c(direction = -1, limits = c(0,35), breaks = c(0,10,20,30)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_cold_all_", Sys.Date(),
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
  scale_fill_viridis_c(direction = -1, limits = c(0,35), breaks = c(0,10,20,30)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_arid_all_", Sys.Date(),
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
  scale_fill_viridis_c(direction = -1, limits = c(0,35), breaks = c(0,10,20,30)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_temp_all_", Sys.Date(),
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
  scale_fill_viridis_c(direction = -1, limits = c(0,35), breaks = c(0,10,20,30)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_trop_all_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)


