# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 01/06/2022 #

# library(ISRaD)
library(tidyverse)
library(ggpubr)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-09-13"))

lyr_all %>% 
  count(entry_name)

lyr_data <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  # filter(lyr_bot <= 100) %>% 
  # filter(lyr_top <= 100) %>%
  group_by(id) %>%
  # filter(min(lyr_top) == 0) %>% 
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
  )) %>% 
  ungroup()

lyr_data %>% 
  count(entry_name) %>% view()

names(lyr_data)

## Mapping sampling locations ##
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)

world <- map_data("world") %>% 
  filter(region != "Antarctica")

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgrey")  +
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
ggsave(file = paste0("./Figure/ISRaD_14C_depth_map_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

lyr_data %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

lyr_data %>% 
  group_by(ClimateZone, pro_usda_soil_order) %>% 
  summarise(n_profiles = n_distinct(id))

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
ggsave(file = paste0("./Figure/ISRaD_climate_soiltype_dis_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)


lyr_data %>% 
  group_by(id) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  ggplot(aes(x = pro_usda_soil_order, fill = site_name)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete("") +
  scale_y_continuous("Number of profiles", expand = c(0,0), limits = c(0,80))

lyr_data %>% 
  group_by(id) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  ggplot(aes(x = pro_wrb_soil_order, fill = site_name)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete("") +
  scale_y_continuous("Number of profiles", expand = c(0,0), limits = c(0,100))

#lyr_14c
plotly::ggplotly(
  lyr_data %>% 
    ggplot(aes(x = depth, y = lyr_14c, group = entry_name)) + 
    geom_point(aes(), size = 3, shape = 21) +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_x_continuous("Depth [cm]", expand = c(0.01,0.01)) +
    scale_y_continuous() 
)

# Density distribution
p0 <- lyr_data %>% 
  ggplot(aes(x = depth, y = lyr_14c)) + 
  geom_hex(color = NA, binwidth = c(10,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous(limits = c(-1005,305)) +
  scale_fill_viridis_c(trans = "log10", limits = c(1,350))

p1 <- lyr_data %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(color = NA, binwidth = c(0.1,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(limits = c(-1005,305)) +
  scale_fill_viridis_c(trans = "log10", limits = c(1,350))

p2 <- lyr_data %>% 
  ggplot(aes(x = CORG, y = lyr_dd14c)) + 
  geom_hex(color = NA, binwidth = c(0.1,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(limits = c(-1505,305)) +
  scale_fill_viridis_c(trans = "log10")

ggarrange(p0, p1, common.legend = TRUE)
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_depth_hex_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

# Colored by sampling year
p1 <- lyr_data %>% 
  mutate(sampl_yr = cut(lyr_obs_date_y,
                        breaks = c(1899,1960,1984,1995,1999,2009,2012,2018))) %>% 
  ggplot(aes(x = depth, y = lyr_14c,
             fill = sampl_yr)) + 
  geom_point(aes(group = entry_name),
             size = 5, alpha = 0.8, shape = 21) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,205)) +
  scale_y_continuous(limits = c(-1005,305)) +
  scale_fill_viridis_d()
ggsave(file = paste0("./Figure/ISRaD_14C_depth_yr_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

p2 <- lyr_data %>% 
  mutate(sampl_yr = cut(lyr_obs_date_y,
                        breaks = c(1899,1960,1984,1995,1999,2009,2012,2022))) %>% 
  ggplot(aes(y = depth, x = lyr_dd14c,
             fill = sampl_yr)) + 
  geom_point(aes(group = entry_name),
             size = 5, alpha = 0.8, shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_reverse("Depth [cm]") +
  scale_x_continuous(limits = c(-1505,305)) +
  scale_fill_viridis_d()

ggarrange(p1, p2, common.legend = TRUE)

lyr_data %>% 
  ggplot(aes(x = depth, y = CORG, z = lyr_14c)) + 
  stat_summary_hex(color = NA, binwidth = c(5,0.1),
                   fun = ~median(.x)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("SOC [wt-%]", trans = "log10") +
  scale_fill_viridis_c("Delta14C", limits = c(-1005,255),
                       option = "A") +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2))

lyr_data %>% 
  ggplot(aes(x = depth, y = lyr_14c)) +
  geom_line(aes(group = id), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(log(x)),
              fill = "lightblue") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,205)) +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250)) 
ggsave(file = paste0("./Figure/ISRaD_14C_depth_profile_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

lyr_data %>%  
  ggplot(aes(y = lyr_14c, x = CORG)) +
  geom_path(aes(group = id), alpha = 0.5) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250))
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_profile_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

lyr_data %>% 
  ggplot(aes(y = lyr_14c, x = CORG, color = pro_BIO12_mmyr_WC2.1)) +
  geom_point(size = 5) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_color_viridis_c("MAP [mm]", trans = "log10", direction = -1)

lyr_data %>% 
  ggplot(aes(y = lyr_14c, x = CORG, color = pro_BIO1_C_WC2.1)) +
  geom_point(size = 5) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_color_viridis_c("MAT [C]")

# Manually assign climate zone for Czimczik_Unpublished
lyr_data %>% 
  filter(is.na(pro_KG_present_long)) %>% 
  count(entry_name)

lyr_data_KG <- lyr_data %>% 
  mutate(pro_KG_present_reclas = case_when(
    is.na(pro_KG_present_long) ~ "Polar, tundra",
    TRUE ~ pro_KG_present_long
  ))

lyr_data_KG %>% 
  filter(is.na(pro_KG_present_reclas)) %>% 
  count(entry_name)

lyr_data_KG %>% 
  filter(depth <= 200) %>% 
  ggplot(aes(x = depth, y = lyr_14c, fill = pro_BIO12_mmyr_WC2.1)) + 
  geom_point(shape = 21, size = 4, alpha = 0.7) +
  facet_wrap(~pro_KG_present_reclas) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Delat14C") +
  scale_fill_viridis_c("MAP [mm]", trans = "log10", direction = -1) +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2))

lyr_data_KG %>% 
  ggplot(aes(x = CORG, y = lyr_14c, fill = pro_BIO12_mmyr_WC2.1)) + 
  geom_point(shape = 21, size = 4, alpha = 0.7) +
  facet_wrap(~pro_KG_present_reclas) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  scale_x_continuous("SOC", trans = "log10") +
  scale_y_continuous("Delat14C") +
  scale_fill_viridis_c("MAP [mm]", trans = "log10", direction = -1) +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2))

lyr_data %>% 
  filter(is.na(pro_usda_soil_order)) %>% 
  count(entry_name)

lyr_data %>% 
  drop_na(pro_usda_soil_order) %>% 
  filter(depth <= 200) %>%
  group_by(id) %>%
  # # Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  ungroup() %>%
  filter(pro_usda_soil_order != "Aridisols",
         pro_usda_soil_order != "Histosols") %>%
  #reclassify soil type Schuur_2001: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Schuur_2001" & pro_usda_soil_order == "Inceptisols",
  #                                      "Andisols")) %>% 
  # #reclassify soil type Guillet_1988: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Guillet_1988",
  #                                      "Andisols")) %>% 
  # #reclassify soil type Torn_1997: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Torn_1997",
  #                                      "Andisols")) %>% 
  # #reclassify soil type Kramer_2012: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Kramer_2012",
  #                                      "Andisols")) %>%
  ggplot(aes(x = depth, y = lyr_14c, fill = pro_BIO12_mmyr_WC2.1)) + 
  geom_point(shape = 21, size = 4, alpha = 0.7) +
  facet_wrap(~pro_usda_soil_order) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Delat14C") +
  scale_fill_viridis_c("MAP [mm]", trans = "log10", direction = -1) +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2)) +
  geom_smooth()
ggsave(file = paste0("./Figure/ISRaD_14C_depth_soiltype_MAP_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

plotly::ggplotly(lyr_data %>% 
  drop_na(pro_usda_soil_order) %>%
  # Filter for studies that have more than 2 depth layers
  filter(pro_usda_soil_order != "Aridisols",
         pro_usda_soil_order != "Histosols") %>%
  #reclassify soil type Schuur_2001: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Schuur_2001" & pro_usda_soil_order == "Inceptisols",
  #                                      "Andisols")) %>% 
  # #reclassify soil type Guillet_1988: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Guillet_1988",
  #                                      "Andisols")) %>% 
  # #reclassify soil type Torn_1997: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Torn_1997",
  #                                      "Andisols")) %>% 
  # #reclassify soil type Kramer_2012: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Kramer_2012",
  #                                      "Andisols")) %>% 
  ggplot(aes(x = CORG, y = lyr_14c, group = entry_name)) +
  geom_point(aes(color = pro_BIO12_mmyr_WC2.1), size = 3, alpha = 0.7) +
  facet_wrap(~pro_usda_soil_order) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  scale_x_continuous("SOC [%]", trans = "log10") +
  scale_color_viridis_c("MAP [mm]", trans = "log10", direction = -1) +
  guides(color = guide_colorbar(barheight = 10, frame.colour = "black", 
                                ticks.linewidth = 2))
)
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_soiltype_MAP_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

lyr_data %>% 
  drop_na(pro_usda_soil_order) %>%
  # group_by(id) %>%
  # # Filter for studies that have more than 2 depth layers
  # filter(n() > 2) %>%
  # ungroup() %>% 
  filter(pro_usda_soil_order != "Aridisols",
         pro_usda_soil_order != "Histosols") %>%
  #reclassify soil type Schuur_2001: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Schuur_2001" & pro_usda_soil_order == "Inceptisols",
                                       "Andisols")) %>%
  #reclassify soil type Guillet_1988: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Guillet_1988",
                                       "Andisols")) %>%
  #reclassify soil type Torn_1997: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Torn_1997",
                                       "Andisols")) %>%
  ggplot(aes(x = CORG, y = lyr_14c)) +
  # geom_point(aes(color = pro_BIO12_mmyr_WC2.1), size = 3) +
  geom_line(aes(group = id, color = pro_BIO12_mmyr_WC2.1), orientation = "x",
            size = 1) +
  facet_wrap(~pro_usda_soil_order) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  scale_x_continuous("SOC [%]", trans = "log10") +
  # geom_smooth(orientation = "y", method = "gam") +
  scale_color_viridis_c("MAP [mm]", trans = "log10", direction = -1) +
  guides(color = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2))
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_soiltype_MAP_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

### Add WRB data
library(raster)
WRB_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/soilOrder_250_x_x_m__WRB.tif"
WRB_raster <- raster::raster(WRB_dir)
WRB_soil_type <- raster::extract(WRB_raster, cbind(lyr_data$pro_long,
                                                   lyr_data$pro_lat))


summary(WRB_soil_type)

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
  mutate(HLZ_Zone = ifelse(id == "Crow_2015_AND_EUC_Composite", 
                           "Warm temperate wet forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Crow_2015_MOL_AG_Composite", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Crow_2015_OX_AG_Composite", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(entry_name == "Cusack_2012", 
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
                           "Subtropical moist forest", HLZ_Zone)) 

lyr_data_HLZ %>% 
  filter(is.na(HLZ_Zone)) %>% 
  count(id)

lyr_data_HLZ %>% 
  count(HLZ_Zone)

lyr_data_HLZ %>% 
  filter(depth <= 200) %>%
  group_by(id) %>%
  # # Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  ungroup() %>%
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

