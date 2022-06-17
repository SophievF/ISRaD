# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 01/06/2022 #

library(ISRaD)
library(tidyverse)
library(ggpubr)

#Load filtered lyr data
lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-06-16"))

lyr_data %>% 
  count(entry_name)

lyr_data %>% 
  count(id)

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
             color = "#4D36C6", shape = 1, size = 3) +
  theme_bw(base_size = 14) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")) +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), 
                     breaks = c(-100,0,100), limits = c(-160,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), 
                     breaks = c(-50,0,50), limits = c(-55,80))
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_map_", Sys.Date(),
                     ".jpeg"), width = 10, height = 5)

# Data exploration ##

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

#lyr_dd14c
plotly::ggplotly(
  lyr_data %>% 
    ggplot(aes(x = depth, y = lyr_dd14c, group = entry_name)) + 
    geom_point(size = 3, shape = 21) +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_x_continuous("Depth [cm]", expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(-1505,305)) 
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
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_depth_hex_", format(Sys.time(), "%Y%m%d"),
                     ".jpeg"), width = 11, height = 6)

# Colored by sampling year
p1 <- lyr_data %>% 
  mutate(sampl_yr = cut(lyr_obs_date_y,
                        breaks = c(1899,1960,1984,1995,1999,2009,2012,2018))) %>% 
  ggplot(aes(y = depth, x = lyr_14c,
             fill = sampl_yr)) + 
  geom_point(aes(group = entry_name),
             size = 5, alpha = 0.8, shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_reverse("Depth [cm]") +
  scale_x_continuous(limits = c(-1505,305)) +
  scale_fill_viridis_d()

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
  drop_na(pro_usda_soil_order) %>% 
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  ungroup() %>% 
  ggplot(aes(x = depth, y = lyr_14c)) +
  geom_line(aes(group = id), alpha = 0.5) +
  # geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(log(x)),
              fill = "lightblue") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,205)) +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250)) 

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
  filter(depth <= 100) %>%
  group_by(id) %>%
  # # Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  ungroup() %>%
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
  #reclassify soil type Kramer_2012: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Kramer_2012",
                                       "Andisols")) %>%
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
  #reclassify soil type Kramer_2012: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Kramer_2012",
                                       "Andisols")) %>% 
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
  filter(CORG > 0) %>% 
  drop_na(pro_usda_soil_order) %>%
  group_by(id) %>%
  # Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  ungroup() %>% 
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
  #reclassify soil type Kramer_2012: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Kramer_2012",
                                       "Andisols")) %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) +
  geom_point(aes(color = pro_BIO12_mmyr_WC2.1), size = 1, alpha = 0.7) +
  geom_line(aes(group = id, color = pro_BIO12_mmyr_WC2.1), orientation = "y", 
            alpha = 0.7) +
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


