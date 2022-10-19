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
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-10-17"))

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

lyr_data$ClimateZone <- factor(lyr_data$ClimateZone,
                               levels = c("polar", "cold", "arid",
                                          "temperate", "tropical"))

lyr_data$ClimateZoneAnd <- factor(lyr_data$ClimateZoneAnd,
                                  levels = c("polar", "cold", "arid",
                                             "temperate", "tropical",
                                             "andisols"))

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
  scale_y_continuous("Number of profiles", expand = c(0,0), limits = c(0,60),
                     breaks = seq(0,60,20))
ggsave(file = paste0("./Figure/ISRaD_climate_soiltype_dis_usda_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

lyr_data %>% 
  group_by(ClimateZoneAnd) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

## Raw data points (from mpsline): SOC vs 14C
mspline_14c_c_all %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.005,60)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,260)) +
  coord_cartesian(ylim = c(-1000,400)) 

ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

mspline_14c_c_all %>% 
  filter(ClimateZone == "polar") %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) + 
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
  scale_fill_viridis_c(direction = -1,limits = c(0,260)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_polar_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)

mspline_14c_c_all %>% 
  filter(ClimateZone == "cold") %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) + 
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
  scale_fill_viridis_c(direction = -1, limits = c(0,260)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_cold_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)

mspline_14c_c_all %>%
  filter(ClimateZone == "arid") %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) + 
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
  scale_fill_viridis_c(direction = -1, limits = c(0,260)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_arid_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)

mspline_14c_c_all %>%
  filter(ClimateZone == "temperate") %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) + 
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
  scale_fill_viridis_c(direction = -1, limits = c(0,260)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_temp_", Sys.Date(),
                     ".jpeg"), width = 5, height = 3.9)

mspline_14c_c_all %>% 
  filter(ClimateZone == "tropical") %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) + 
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
  scale_fill_viridis_c(direction = -1, limits = c(0,260)) +
  coord_cartesian(ylim = c(-1000,400)) 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_trop_", Sys.Date(),
                     ".jpeg"), width = 4, height = 3.9)

## Raw data + averaged data: SOC and climate
climate_all_and <- mspline_14c_c_all %>%
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

climate_all <- mspline_14c_c_all %>%
  group_by(ClimateZone, UD) %>% 
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

mspline_14c_c_all %>% 
  ggplot() + 
  geom_hex(aes(x = CORG_msp, y = lyr_14c_msp), alpha = 0.7) +
  geom_path(data = climate_all %>%
              filter(n > 4 & n_rel > 33),
            aes(x = median_c, y = median_14c, color = ClimateZone), size = 2) +
  geom_errorbar(data = climate_all %>%
                  filter(n > 4 & n_rel > 33),
                aes(ymin = lci_14c, ymax = uci_14c, x = median_c, color = ClimateZone),
                alpha = 0.3) +
  geom_errorbarh(data = climate_all %>%
                  filter(n > 4 & n_rel > 33),
                 aes(xmin = lci_c, xmax = uci_c, y = median_14c, color = ClimateZone),
                 alpha = 0.3) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8,0.2), 
        legend.box = "horizontal") +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1) +
  coord_cartesian(ylim = c(-1000,400)) 

ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_fit", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

climate_all_and$ClimateZoneAnd <- factor(climate_all_and$ClimateZoneAnd,
                                  levels = c("andisols", "polar", "cold", 
                                             "temperate", "arid",
                                             "tropical"))

climate_all_and %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.15,0.25)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Climate zones")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_avg_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

climate_soil_all <- mspline_14c_c_all %>%
  group_by(ClimateZone, pro_usda_soil_order, UD) %>% 
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

climate_soil_all %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZone), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZone), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZone), alpha = 0.3) +
  facet_wrap(~pro_usda_soil_order) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Climate zones")

## Arid soils
#Krull_2005: "medium heavy clay" ~ 55% clay content
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$entry_name == "Krull_2005"), 55)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$entry_name == "Krull_2005"), 55)

#Harden_1987_PM13&14_PM14-2: based on close by profile: ~ 40% clay content
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Harden_1987_PM13&14_PM14-2"), 40)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Harden_1987_PM13&14_PM14-2"), 40)

#Harden_1987_PM22_PM22: "silty loam": ~ 13% clay content
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Harden_1987_PM22_PM22"), 13)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Harden_1987_PM22_PM22"), 13)

#Khomo_2017_MG-550-C_MG-550-C1: based on other profiles/depth layers: ~ 15%
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Khomo_2017_MG-550-C_MG-550-C1"), 15)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Khomo_2017_MG-550-C_MG-550-C1"), 15)

#Krull_2003_CHI:-26.72,150.6_CHI:-26.72,150.6_346: based on reported values in publication: ~ 60%
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Krull_2003_CHI:-26.72,150.6_CHI:-26.72,150.6_346"), 65)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Krull_2003_CHI:-26.72,150.6_CHI:-26.72,150.6_346"), 65)

#Krull_2003_PG:-27.45,150.52_PG:-27.45,150.52_347: based on reported values in publication: ~ 60%
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Krull_2003_PG:-27.45,150.52_PG:-27.45,150.52_347"), 60)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Krull_2003_PG:-27.45,150.52_PG:-27.45,150.52_347"), 60)

#Leavitt_2007_MTS_MTS: "loam": ~ 17%
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Leavitt_2007_MTS_MTS"), 17)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Leavitt_2007_MTS_MTS"), 17)

#Leavitt_2007_BLS_BLS: "clay loam": ~ 35%
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Leavitt_2007_BLS_BLS"), 35)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Leavitt_2007_BLS_BLS"), 35)

#Leavitt_2007_COS_COS: "silty loam": ~ 13%
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Leavitt_2007_COS_COS"), 13)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Leavitt_2007_COS_COS"), 13)

#Leavitt_2007_DHS_DHS: "fine sandy loam": ~ 10%
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "Leavitt_2007_DHS_DHS"), 10)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "Leavitt_2007_DHS_DHS"), 10)

#McClaran_2000_Empire-Cienega:31.7528,-110.6278_Empire-Cienega:31.7528,-110.6278_556: "sandy-loam to loam": ~ 13%
mspline_14c_c_all$lyr_clay_tot_psa <- replace(mspline_14c_c_all$lyr_clay_tot_psa,
                                              which(mspline_14c_c_all$id == "McClaran_2000_Empire-Cienega:31.7528,-110.6278_Empire-Cienega:31.7528,-110.6278_556"), 13)
climate_soil_all$lyr_clay_tot_psa <- replace(climate_soil_all$lyr_clay_tot_psa,
                                             which(climate_soil_all$id == "McClaran_2000_Empire-Cienega:31.7528,-110.6278_Empire-Cienega:31.7528,-110.6278_556"), 13)


mspline_14c_c_all %>%
  filter(ClimateZone == "arid") %>% 
  dplyr::select(entry_name, id, lyr_clay_tot_psa) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  filter(is.na(lyr_clay_tot_psa))

mspline_14c_c_all %>% 
  filter(ClimateZone == "arid") %>% 
  mutate(clay_group = case_when(
    lyr_clay_tot_psa < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
  )) %>% 
  group_by(clay_group, UD) %>% 
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
  mutate(n_rel = n * 100 / max(n)) %>% 
  ggplot() +
  geom_path(aes(x = median_c, y = median_14c, color = clay_group), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c,
                    color = clay_group), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c,
                     color = clay_group), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.08,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Clay content")

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_arid_clay_avg_", Sys.Date(),
                     ".jpeg"), width = 9, height = 5)

mspline_14c_c_all %>% 
  filter(ClimateZone == "arid") %>% 
  mutate(clay_group = case_when(
    lyr_clay_tot_psa < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
  )) %>%
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, group = id, color = clay_group)) +
  geom_path(size = 2) + 
  facet_wrap(~pro_KG_present_long) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  scale_color_discrete("Clay content")

mspline_14c_c_all %>% 
  filter(ClimateZone == "arid") %>% 
  mutate(clay_group = case_when(
    lyr_clay_tot_psa < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
  )) %>% 
  group_by(clay_group) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(pro_name),
            n_profiles = n_distinct(id))

mspline_14c_c_all %>% 
  filter(ClimateZone == "arid") %>% 
  mutate(clay_group = case_when(
    lyr_clay_tot_psa < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
  )) %>%
  group_by(clay_group) %>%
  dplyr::select(lyr_clay_tot_psa) %>% 
  skimr::skim()

# Polar climates
mspline_14c_c_all %>% 
  filter(ClimateZone == "polar") %>% 
  group_by(pro_KG_present_long, UD) %>% 
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
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() +
  geom_path(aes(x = median_c, y = median_14c, color = pro_KG_present_long), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c,
                    color = pro_KG_present_long), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c,
                     color = pro_KG_present_long), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.08,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Climate")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_polar_avg_climate_", Sys.Date(),
                     ".jpeg"), width = 9, height = 5)


climate_all %>%
  filter(ClimateZone == "polar") %>%
  # filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZone), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZone), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZone), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Climate zones")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_polar_avg_", Sys.Date(),
                     ".jpeg"), width = 9, height = 5)

mspline_14c_c_all %>% 
  filter(ClimateZoneAnd == "polar") %>% 
  group_by(pro_KG_present_long) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(pro_name),
            n_profiles = n_distinct(id))

# Tropical climates
climate_all_and %>%
  filter(ClimateZoneAnd == "tropical") %>%
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZone), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZone), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZone), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Climate zones")

mspline_14c_c_all %>% 
  filter(ClimateZoneAnd == "tropical") %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, group = id, color = pro_clay_100_cm_SG)) +
  geom_path(size = 2) +
  # facet_wrap(~pro_wrb_soil_order) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  scale_color_viridis_c()

climate_soil_all %>% 
  filter(ClimateZone == "tropical") %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = pro_usda_soil_order), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = pro_usda_soil_order), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = pro_usda_soil_order), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8,0.2)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Soil type")

mspline_14c_c_all %>% 
  filter(ClimateZoneAnd == "tropical") %>%
  group_by(pro_KG_present_long, UD) %>% 
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
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>%
  ggplot() +
  geom_path(aes(x = median_c, y = median_14c, 
                color = factor(pro_KG_present_long, 
                               levels = c("Tropical, rainforest",
                                          "Tropical, monsoon",
                                          "Tropical, savannah"))), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c,
                    color = pro_KG_present_long), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c,
                     color = pro_KG_present_long), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8,0.2)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Climate zone")

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_tropical_climate_avg_", Sys.Date(),
                     ".jpeg"), width = 10, height = 5)

mspline_14c_c_all %>% 
  filter(ClimateZoneAnd == "tropical") %>%
  mutate(clay_sg = cut_interval(pro_clay_0_cm_SG, 2)) %>% 
  mutate(clay_group = case_when(
    # lyr_clay_tot_psa < 36 ~ "< 36 %",
    # lyr_clay_tot_psa > 36 ~ "> 36 %",
    clay_sg == "[16,35]" ~ "< 36 %",
    clay_sg == "(35,54]" ~ "> 36 %"
  )) %>% 
  group_by(clay_group, UD) %>% 
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
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>%
  ggplot() +
  geom_path(aes(x = median_c, y = median_14c, color = clay_group), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c,
                    color = clay_group), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c,
                     color = clay_group), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Clay content")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_tropical_clay_avg_", Sys.Date(),
                     ".jpeg"), width = 9, height = 5)

mspline_14c_c_all %>% 
  filter(ClimateZoneAnd == "tropical") %>% 
  mutate(clay_group = cut_interval(pro_clay_0_cm_SG, 2)) %>% 
  group_by(clay_group) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(pro_name),
            n_profiles = n_distinct(id))
  
  
# Temperate climates
climate_all_and %>%
  filter(ClimateZoneAnd == "temperate") %>%
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZone), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZone), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZone), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Climate zones")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_temp_usda_avg_", Sys.Date(),
                    ".jpeg"), width = 9, height = 5)

climate_soil_all %>%
  filter(ClimateZone == "temperate") %>%
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = pro_usda_soil_order), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = pro_usda_soil_order), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = pro_usda_soil_order), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,60)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Soil type")

mspline_14c_c_all %>% 
  filter(ClimateZoneAnd == "temperate") %>%
  mutate(wrb_group = case_when(
    pro_wrb_soil_order == "Acrisols" ~ "a",
    pro_wrb_soil_order == "Ferralsols" ~ "a",
    pro_wrb_soil_order == "Nitisols" ~ "a",
    pro_wrb_soil_order == "Alisols" ~ "a",
    pro_wrb_soil_order == "Vertisols" ~ "b",
    pro_wrb_soil_order == "Andosols" ~ "b",
    TRUE ~ "c"
  )) %>% 
  group_by(pro_wrb_soil_order, UD) %>% 
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
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>%
  ggplot() +
  geom_path(aes(x = median_c, y = median_14c, color = pro_wrb_soil_order), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c,
                    color = pro_wrb_soil_order), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c,
                     color = pro_wrb_soil_order), alpha = 0.3) +
  facet_wrap(~wrb_group) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.9, 0.25)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,30)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,200)) +
  scale_color_discrete("Soil type")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_temp_wrb_avg_", Sys.Date(),
                    ".jpeg"), width = 10, height = 5)

mspline_14c_c_all %>% 
  filter(ClimateZoneAnd == "temperate") %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, group = id, color = pro_BIO12_mmyr_WC2.1)) +
  geom_path(size = 2) +
  # facet_wrap(~pro_usda_soil_order) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  scale_color_viridis_c(direction = -1)
