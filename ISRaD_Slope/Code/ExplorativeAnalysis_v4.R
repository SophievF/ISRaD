# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 09/01/2023 #

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2023-02-08"))

# lyr_all %>% 
#   summarise(n_studies = n_distinct(entry_name),
#             n_sites = n_distinct(site_name),
#             n_profiles = n_distinct(id))
# 
# lyr_data <- lyr_all %>% 
#   # filter(lyr_obs_date_y > 1959) %>% 
#   group_by(id) %>%
#   #Filter for studies that have more than 2 depth layers
#   filter(n() > 2) %>%
#   arrange(depth, .by_group = TRUE) %>% 
#   ungroup() %>% 
#   mutate(ClimateZone = case_when(
#     entry_name == "Gentsch_2018" ~ "tundra/polar",
#     pro_usda_soil_order == "Gelisols" ~ "tundra/polar",
#     pro_usda_soil_order == "Aridisols" ~ "arid",
#     str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
#     str_detect(pro_KG_present_long, "Temperate") ~ "warm temperate",
#     str_detect(pro_KG_present_long, "Cold") ~ "cold temperate",
#     str_detect(pro_KG_present_long, "Polar") ~ "tundra/polar",
#     str_detect(pro_KG_present_long, "Arid") ~ "arid",
#   )) %>% 
#   mutate(ClimateZoneAnd = case_when(
#     entry_name == "Gentsch_2018" ~ "tundra/polar",
#     pro_usda_soil_order == "Gelisols" ~ "tundra/polar",
#     pro_usda_soil_order == "Andisols" ~ "volcanic soils",
#     pro_usda_soil_order == "Aridisols" ~ "arid",
#     str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
#     str_detect(pro_KG_present_long, "Temperate") ~ "warm temperate",
#     str_detect(pro_KG_present_long, "Cold") ~ "cold temperate",
#     str_detect(pro_KG_present_long, "Polar") ~ "tundra/polar",
#     str_detect(pro_KG_present_long, "Arid") ~ "arid",
#   )) %>% 
#   mutate(MineralType = case_when(
#     pro_usda_soil_order == "Andisols" ~ "amorphous",
#     pro_usda_soil_order == "Oxisols" ~ "low-activity clay",
#     pro_usda_soil_order == "Ultisols" ~ "low-activity clay",
#     TRUE ~ "high-activity clay"
#   )) %>% 
#   #remover overlapping depth layer in Fernandez_1993
#   filter(lyr_name != "B_61.1-122.2") %>% 
#   mutate(pro_MAT_mod = case_when(
#     is.na(pro_MAT) ~ pro_BIO1_C_WC2.1,
#     TRUE ~ pro_MAT
#   )) %>% 
#   mutate(pro_MAP_mod = case_when(
#     is.na(pro_MAP) ~ pro_BIO12_mmyr_WC2.1,
#     TRUE ~ pro_MAP
#   ))

lyr_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-03-09.csv")

lyr_data %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

names(lyr_data)

lyr_data$ClimateZone <- factor(lyr_data$ClimateZone,
                               levels = c("tundra/polar", "cold temperate", "arid",
                                          "warm temperate", "tropical"))

lyr_data$ClimateZoneAnd <- factor(lyr_data$ClimateZoneAnd,
                                  levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                             "warm temperate", "arid","tropical"))

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

# plot(climate_raster)

#reclassify climate raster
recal <- c(0,3,1, 3,7,2, 7,16,3, 16,28,4, 28,30,5)
recal_mat <- matrix(recal, ncol = 3, byrow = TRUE)

climate_grp <- reclassify(climate_raster, recal_mat)
# plot(climate_grp)

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
    ClimateZone == 3 ~ "warm temperate",
    ClimateZone == 4 ~ "cold temperate",
    ClimateZone == 5 ~ "tundra/polar"
  ))

rm(climate_pts)
summary(climate_df)

climate_df$ClimateZone <- factor(climate_df$ClimateZone,
                                 levels = c("tundra/polar", "cold temperate", "arid",
                                            "warm temperate", "tropical"))

climate_color <- c("#e66101", "#fdb863", "#f6e8c3", "#b2abd2", "#5e3c99")
# color_palette(palette = "Purples")

ggplot() +
  geom_raster(data = climate_df, 
              aes(x = x, y = y, fill = ClimateZone)) +
  geom_point(data = lyr_data,
             aes(x = pro_long, y = pro_lat),
             fill = "white", size = 1.5, shape = 21, color = "black") +
  theme_bw(base_size = 16) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = c(0.11,0.2)) +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), expand = c(0,0),
                     breaks = c(-100,0,100), limits = c(-170,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), expand = c(0,0),
                     breaks = c(-50,0,50), limits = c(-55,80)) +
  # scale_fill_brewer("Climate zones", palette = "Purples") +
  scale_fill_manual("Climate zones",
                     values = c("#8c96c6", "#bcbddc", "#9e9ac8", "#756bb1", "#54278f"))
  
ggsave(file = paste0("./Figure/ISRaD_14C_Climate_map_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)
  
ggplot() +  
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "#d9d9d9", fill = "#d9d9d9") +
  geom_point(data = lyr_data, 
             aes(x = pro_long, y = pro_lat),
             fill = "purple", size = 2, shape = 21, color = "white") +
  theme_bw(base_size = 16) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")) +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), 
                     breaks = c(-100,0,100), limits = c(-170,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), 
                     breaks = c(-50,0,50), limits = c(-55,80))
ggsave(file = paste0("./Figure/ISRaD_14C_map_", Sys.Date(),
                     ".jpeg"), width = 8, height = 4)

lyr_data %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

# Data example
p1 <- ggplot() +
  geom_line(data = lyr_data %>% 
              filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
            aes(y = UD, x = lyr_14c_msp, group = id, color = id), orientation = "y",
            size = 1.5) +
  geom_pointrange(data = lyr_all %>% 
                    filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
                  aes(x = lyr_14c, y = depth, ymin = lyr_top, ymax = lyr_bot,
                      color = id), size = 1) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.29,0.79),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)) +
  scale_y_reverse("Depth [cm]", expand = c(0,0)) +
  scale_x_continuous(expression(paste(Delta^14,"C [‰]")), 
                     expand = c(0,0), limits = c(-1000,200),
                     position = "top") +
  scale_color_discrete(paste0("Example profile:"), 
                       label = c("P1", "P2")) +
  coord_cartesian(ylim = c(100,0))

p2 <- ggplot() +
  geom_line(data = lyr_data %>% 
              filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
            aes(y = UD, x = CORG_msp, group = id, color = id), orientation = "y",
            size = 1.5) +
  geom_pointrange(data = lyr_all %>% 
                    filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
                  aes(x = CORG, y = depth, ymin = lyr_top, ymax = lyr_bot,
                      color = id), size = 1) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_y_reverse("", expand = c(0,0)) +
  scale_x_continuous("SOC [wt-%]; log-scaled", expand = c(0,0), position = "top", 
                     trans = "log10") +
  scale_color_discrete(paste0("Example profile:"), 
                       label = c("P1", "P2")) +
  coord_cartesian(ylim = c(100,0))

ggarrange(p1, p2)
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_depth_example_", Sys.Date(),
                     ".jpeg"), width = 10, height = 5)

lyr_all %>% 
  filter(lyr_bot > 100) %>% 
  summarise(n_id = n_distinct(id),
            n_site = n_distinct(site_name))



### Data distribution
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
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, trans = "log10", limits = c(1,340)) +
  coord_cartesian(xlim = c(0,100), ylim = c(-1000,400))

ggarrange(p2, p1, common.legend = TRUE)

lyr_data %>% 
  filter(depth < 101) %>% 
  ggplot(aes(x = lyr_13c, y = lyr_14c)) + 
  geom_hex(binwidth = c(1,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous(expression(paste(delta^13, "C [‰]"))) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(trans = "log10", direction = -1, limits = c(1,340)) +
  coord_cartesian(ylim = c(-1000,400))

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

lyr_data %>% 
  group_by(id) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  ggplot(aes(x = MineralType)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) +
  facet_wrap(~ClimateZone) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  scale_x_discrete("") +
  scale_y_continuous("Number of profiles", expand = c(0,0), limits = c(0,170),
                     breaks = seq(0,155,50))

lyr_data %>% 
  group_by(id) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  ggplot(aes(x = pro_land_cover)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) +
  facet_wrap(~ClimateZone) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  scale_x_discrete("") +
  scale_y_continuous("Number of profiles", expand = c(0,0), limits = c(0,170),
                     breaks = seq(0,155,50))

##GAP-FILL LAND COVER DATA

lyr_data %>% 
  ggplot(aes(x = depth)) +
  stat_ecdf(color = "red") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Cummulative density [%]", expand = c(0,0)) +
  scale_x_continuous("Mid-sampling depth [cm]", expand = c(0,0),
                     limits = c(0,600)) +
  geom_vline(xintercept = 100, linetype = "dashed") +
  geom_hline(yintercept = 0.835, linetype = "dashed")

lyr_data %>% 
  ggplot(aes(x = lyr_top)) +
  stat_ecdf(color = "red") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Cummulative density [%]", expand = c(0,0)) +
  scale_x_continuous("Upper-sampling depth [cm]", expand = c(0,0),
                     limits = c(0,600)) +
  geom_vline(xintercept = 100, linetype = "dashed") +
  geom_hline(yintercept = 0.85, linetype = "dashed")

lyr_data %>% 
  ggplot(aes(x = lyr_bot)) +
  stat_ecdf(color = "red") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Cummulative density [%]", expand = c(0,0)) +
  scale_x_continuous("Lower-sampling depth [cm]", expand = c(0,0),
                     limits = c(0,600)) +
  geom_vline(xintercept = 100, linetype = "dashed") +
  geom_hline(yintercept = 0.785, linetype = "dashed")


### Apply mspline function
summary(lyr_data$depth)
summary(lyr_data$lyr_14c)
summary(lyr_data$CORG)

## mspline 14C
lyr_data_mpspline_14c <- lyr_data %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, vhigh = 350, lam = 0.5)

## mspline CORG
lyr_data_mpspline_c <- lyr_data %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.007, vhigh = 60, lam = 0.5)

## 14C and SOC
mspline_14c_c <- lyr_data_mpspline_14c$est_1cm %>% 
  rename(lyr_14c_msp = SPLINED_VALUE) %>% 
  full_join(lyr_data_mpspline_c$est_1cm %>% 
              rename(CORG_msp = SPLINED_VALUE)) %>% 
  filter(UD < 101) %>% 
  tibble()

mspline_14c_c_all <- mspline_14c_c %>%
  dplyr::left_join(lyr_data %>% 
                     distinct(id, .keep_all = TRUE), 
                   by = "id") %>% 
  group_by(id) %>% 
  arrange(UD) %>% 
  ungroup() 

summary(mspline_14c_c$lyr_14c_msp)
summary(mspline_14c_c$CORG_msp)
summary(mspline_14c_c$UD)
summary(mspline_14c_c$LD)

mspline_14c_c_all %>% 
  group_by(ClimateZoneAnd) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

mspline_14c_c_all %>% 
  group_by(MineralType) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

### Data exploration ###
mspline_14c_c_all %>% 
  group_by(UD, ClimateZoneAnd) %>% 
  mutate(median_14c = median(lyr_14c_msp),
         median_CORG = median(CORG_msp),
         n = n()) %>% 
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(n_rel > 33) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_CORG, median_14c, n_rel) %>% 
  arrange(desc(median_14c)) %>% 
  mutate(diff_c = rev(cumsum(rev(median_CORG*100/sum(median_CORG))))) %>% 
  mutate(cum_c = cumsum(median_CORG*100/sum(median_CORG))) %>%
  mutate_at("cum_c", ~replace(., which.min(.), 0)) %>% 
  ggplot(aes(x = median_14c, y = cum_c, color = ClimateZoneAnd)) +
  geom_path(size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125),
                     breaks = c(-1000,-750,-500,-250,0,125)) +
  scale_y_continuous("Cummulative SOC content [%]", expand = c(0,0),
                     trans = "reverse", limits = c(100,0)) +
  scale_color_discrete("Climate zone")
ggsave(file = paste0("./Figure/ISRaD_14C_cumC_climate_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)


mspline_14c_c_all %>% 
  group_by(UD, MineralType) %>% 
  mutate(median_14c = median(lyr_14c_msp),
         median_CORG = median(CORG_msp),
         n = n()) %>% 
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>% 
  filter(n_rel > 33) %>% 
  dplyr::select(MineralType, UD, median_CORG, median_14c, n_rel) %>% 
  arrange(desc(median_14c)) %>% 
  mutate(diff_c = rev(cumsum(rev(median_CORG*100/sum(median_CORG))))) %>% 
  mutate(cum_c = cumsum(median_CORG*100/sum(median_CORG))) %>%
  mutate_at("cum_c", ~replace(., which.min(.), 0)) %>%
  ggplot(aes(x = median_14c, y = cum_c, color = MineralType)) +
  geom_path(size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), 
                     breaks = c(-1000,-750,-500,-250,0,125)) +
  scale_y_continuous("Cummulative SOC content [%]", expand = c(0,0),
                     trans = "reverse", limits = c(100,0)) +
  scale_color_discrete("Mineral type")
ggsave(file = paste0("./Figure/ISRaD_14C_cumC_mineral_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

mspline_14c_c_all %>% 
  pivot_longer(cols = c(lyr_14c_msp, lyr_14c), names_to = "data_14c",
               values_to = "values_14c") %>% 
  ggplot() +
  stat_ecdf(aes(x = values_14c, color = data_14c), pad = FALSE) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Cummulative probaility [%]", expand = c(0,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     trans = "reverse")

mspline_14c_c_all %>% 
  pivot_longer(cols = c(CORG_msp, CORG), names_to = "data_CORG",
               values_to = "values_CORG") %>% 
  ggplot() +
  stat_ecdf(aes(x = values_CORG, color = data_CORG), pad = FALSE) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Cummulative propability [%]", expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", expand = c(0,0))

mspline_14c_c_all %>% 
  ggplot(aes(y = UD, x = lyr_14c_msp, group = id, color = ClimateZoneAnd)) + 
  geom_path() +
  facet_wrap(~ClimateZoneAnd) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(),
        legend.position = "none",
        panel.spacing = unit(1, "cm")) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,300)) 
# ggsave(file = paste0("./Figure/ISRaD_14C_depth_climate_", Sys.Date(),
#                      ".jpeg"), width = 12, height = 6)

mspline_14c_c_all %>% 
  ggplot(aes(y = UD, x = CORG_msp, group = id, color = ClimateZoneAnd)) + 
  geom_path() +
  facet_wrap(~ClimateZoneAnd) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(),
        legend.position = "none",
        panel.spacing = unit(1, "cm")) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous("SOC [wt-%]", expand = c(0,0), trans = "log10")


mspline_14c_c_all %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) + 
  # geom_path(aes(group = id)) +
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZoneAnd) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.005,60)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,260)) +
  coord_cartesian(ylim = c(-1000,400)) 
# ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_", Sys.Date(),
#                      ".jpeg"), width = 12, height = 6)


## Splined and averaged data
climate_all_and <- lyr_data %>%
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


climate_all_and$ClimateZoneAnd <- factor(climate_all_and$ClimateZoneAnd,
                                  levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                             "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_and %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)

climate_all_and %>% 
  filter(n > 4 & n_rel > 33) %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::summarise(UD_min = min(UD),
                   UD_max = max(UD),
                   SOC_max = max(median_c),
                   SOC_min = min(median_c),
                   C14_max = max(median_14c),
                   C14_min = min(median_14c)) %>% 
  mutate(SOC_diff = SOC_max - SOC_min,
         C14_diff = C14_max - C14_min)


c1 <- climate_all_and %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_", Sys.Date(),
                     ".jpeg"), width = 8, height = 5)

c1_1 <- climate_all_and %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = UD, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = UD, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.15),
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Depth [cm]", limits = c(0,100), expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 

ggarrange(c1, c1_1, common.legend = TRUE)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_depth_Climate_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

mineral_type <- lyr_data %>%
  group_by(MineralType, UD) %>% 
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

depth_sum_mineral <- mineral_type %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(MineralType, UD, median_c, median_14c) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)

mineral_type %>% 
  filter(n > 4 & n_rel > 33) %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::summarise(UD_min = min(UD),
                   UD_max = max(UD),
                   SOC_max = max(median_c),
                   SOC_min = min(median_c),
                   C14_max = max(median_14c),
                   C14_min = min(median_14c)) %>% 
  mutate(SOC_diff = SOC_max - SOC_min,
         C14_diff = C14_max - C14_min)


m1 <- mineral_type %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = MineralType), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralType), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralType), alpha = 0.3) +
  geom_point(data = depth_sum_mineral, aes(x = median_c, y = median_14c), shape = 21,
             size = 2, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        legend.position = c(0.3,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,35), expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 

plot(m1)

# ggarrange(c1,c2, widths = c(1.2,1))
# ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_Climate_ClayType_avg_", Sys.Date(),
#                      ".jpeg"), width = 20, height = 11.5)

m1_1 <- mineral_type %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = UD, y = median_14c, color = MineralType), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = UD, 
                     color = MineralType), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Depth [cm]", limits = c(0,100), expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 

ggarrange(m1, m1_1, common.legend = TRUE)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_depth_mineral_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

## Plot with same number of profiles in each group
lyr_data %>% 
  group_by(ClimateZoneAnd) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

# set.seed(42)
# set.seed(123)
# set.seed(456)
set.seed(789)
rdm_id <- mspline_14c_c_all %>% 
  distinct(id, .keep_all = TRUE) %>% 
  group_by(ClimateZoneAnd) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  sample_n(28) %>% 
  ungroup()

climate_rdm <- mspline_14c_c_all %>% 
  filter(id %in% rdm_id$id) %>% 
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

climate_rdm$ClimateZoneAnd <- factor(climate_rdm$ClimateZoneAnd,
                                     levels = c("tundra/polar", "cold temperate", 
                                                "warm temperate", "arid", "tropical"))

depth_sum <- climate_rdm %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)


c1_rdm <- climate_rdm %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 

# plot(c1_rdm)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_Climate_rdm_4_", Sys.Date(),
                     ".jpeg"), width = 8, height = 5)
