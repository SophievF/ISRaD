# Data analysis: Profile data ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 08/03/2023 #

library(tidyverse)
library(ggpubr)
library(scales)
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(iml)
library(mlr3spatial)
library(mlr3spatiotempcv) 
library(RColorBrewer)

# Load filtered and splined data
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

# Check data
lyr_data %>% 
  dplyr::summarise(n_studies = n_distinct(entry_name),
                   n_sites = n_distinct(site_name),
                   n_profiles = n_distinct(id))

lyr_data %>% 
  dplyr::select(lyr_14c_msp, CORG_msp, UD, pro_MAT_mod, pro_AI, lyr_clay_mod) %>% 
  cor()

skimr::skim_without_charts(lyr_data$lyr_obs_date_y)

lyr_data %>% 
  ggplot(aes(x = lyr_obs_date_y)) +
  stat_ecdf() +
  theme_bw()

# Convert characters into factors
lyr_data$ClimateZone <- factor(lyr_data$ClimateZone,
                               levels = c("tundra/polar", "cold temperate",
                                          "warm temperate", "tropical", "arid"))

lyr_data$ClimateZoneAnd <- factor(lyr_data$ClimateZoneAnd,
                                  levels = c("volcanic soils", "tundra/polar", 
                                             "cold temperate", "warm temperate", 
                                             "tropical", "arid"))

lyr_data$MineralGroupsNew <- factor(lyr_data$MineralGroupsNew,
                                    levels = c("amorphous", 
                                               "primary mineral",
                                               "high-activity clay",
                                               "low-activity clay"))

lyr_data %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::summarise(n_studies = n_distinct(entry_name),
                   n_sites = n_distinct(site_name),
                   n_profiles = n_distinct(id),
                   n_countries = n_distinct(pro_country))

lyr_data %>% 
  group_by(MineralGroupsNew) %>% 
  dplyr::summarise(n_studies = n_distinct(entry_name),
                   n_sites = n_distinct(site_name),
                   n_profiles = n_distinct(id),
                   n_countries = n_distinct(pro_country))

lyr_data %>% 
  group_by(ClimateZoneAnd) %>% 
  skimr::skim_without_charts(pro_MAT_mod, pro_AI, lyr_clay_mod)

lyr_data %>% 
  group_by(MineralGroupsNew) %>% 
  skimr::skim_without_charts(pro_MAT_mod, pro_AI, lyr_clay_mod)

## Define color code 
# cold, warm, tropical, arid
color_climate_wo_polar <- c("#bdd7e7", "#7fbf7b", "#1b7837", "#dfc27d")
# polar, cold, warm, tropical, arid
color_climate <- c("#3182bd", "#bdd7e7", "#7fbf7b", "#1b7837", "#dfc27d")
# amorphous, polar, cold, warm, tropical, arid
color_climate_w_amorph <- c("#c083be", "#3182bd", "#bdd7e7", "#7fbf7b", 
                            "#1b7837", "#dfc27d")
# amorphous, primary mineral, high-activity clay, low-activity clay
color_mineral <- c("#c083be", "#bf812d", "#bdbdbd", "#fc9272")

## Mapping sampling locations ##
world <- map_data("world") %>% 
  filter(region != "Antarctica")

set.seed(123)
map_loc <- ggplot() +  
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "#d9d9d9", fill = "#d9d9d9") +
  geom_jitter(data = lyr_data %>% 
                distinct(id, .keep_all = TRUE), 
              aes(x = pro_long, y = pro_lat, fill = ClimateZone),
              size = 2, shape = 21, color = "white",
              width = 1, height = 1) +
  theme_bw(base_size = 14) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.margin = margin(c(t = 0, r = 0, l = 0, b = 0))) +
  coord_sf() +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), 
                     breaks = c(-100,0,100), limits = c(-170,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), 
                     breaks = c(-50,0,50), limits = c(-55,80)) +
  scale_fill_manual("Grouping", values = color_climate)

### Data distribution
amorph <- "#c994c7"
prim <- c("#f6e8c3", "#dfc27d", "#bf812d", "#8c510a", "#6c3e06")
hac <- brewer.pal(3, "Greys")
lac <- c("#fee0d2", "#fc9272")

MineralGroup_color <- c(amorph, prim, hac, lac)

bar_dis <- lyr_data %>% 
  group_by(id) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  mutate(pro_usda_soil_order = factor(pro_usda_soil_order,
                                      levels = c("Andisols", "Alfisols", 
                                                 "Aridisols",  "Entisols", 
                                                 "Gelisols", "Inceptisols",
                                                 "Mollisols", "Spodosols", 
                                                 "Vertisols", "Oxisols", 
                                                 "Ultisols", ""))) %>% 
  ggplot(aes(x = MineralGroupsNew, fill = pro_usda_soil_order)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) +
  facet_wrap(~ClimateZone, nrow = 3) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 5, l = 5, r = 5, b = 5),
        legend.position = c(0.75,0.05)) +
  scale_x_discrete("", labels = c("amorphous", "primary\nmineral",
                                  "high-activity\nclay", "low-activity\nclay")) +
  scale_y_continuous("Number of profiles", expand = c(0,0), limits = c(0,100),
                     breaks = seq(0,100,20)) +
  scale_fill_manual("Soil types", values = MineralGroup_color) +
  guides(fill = guide_legend(ncol = 2))
# ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_Mineral_bar_", Sys.Date(),
#                      ".jpeg"), width = 11, height = 7)

g_bar_dis <- ggplot_gtable(ggplot_build(bar_dis))
strip_top <- which(grepl('strip-t', g_bar_dis$layout$name))
strip_top <- strip_top[-2]
fills <- c("#dfc27d", "#7fbf7b", "#1b7837", "#3182bd", "#bdd7e7")
k <- 1
for(i in strip_top){
  j <- which(grepl('rect', g_bar_dis$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_bar_dis$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- 1 + k
}

# Figure 1 - manually modified for publication
ggarrange(map_loc, g_bar_dis, widths = c(1.5,1))
ggsave(file = paste0("./Figure/Figure_1_", 
                     Sys.Date(),".jpeg"), width = 15, height = 9)

## Boxplots
lyr_data %>% 
  group_by(ClimateZoneAnd) %>% 
  skimr::skim_without_charts(pro_MAT_mod, pro_AI, lyr_clay_mod)

lyr_data %>% 
  group_by(MineralGroupsNew) %>% 
  skimr::skim_without_charts(pro_MAT_mod, pro_AI, lyr_clay_mod)

dis_clim <- lyr_data %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(ClimateZoneAnd, pro_MAT_mod, pro_MAP_mod, pro_PET_mm_yr_mod,
                lyr_clay_mod) %>% 
  rename("Clay content [%]" = lyr_clay_mod,
         "PET [mm]" = pro_PET_mm_yr_mod,
         "MAP [mm]" = pro_MAP_mod,
         "MAT [°C]" = pro_MAT_mod) %>% 
  pivot_longer(!ClimateZoneAnd, names_to = "names", values_to = "values") %>% 
  # factor(x = character(names), 
  #        levels = c("MAT [°C]", "MAP [mm]", "PET [mm]", "Clay content [%]")) %>% 
  group_by(ClimateZoneAnd) %>% 
  ggplot(aes(x = ClimateZoneAnd, y = values, fill = ClimateZoneAnd)) +
  geom_boxplot() +
  # geom_violin() +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        strip.background = element_rect(fill = "white")) +
  facet_wrap(~names, scales = "free", nrow = 1) +
  scale_fill_manual("Climate groups", values = color_climate)

dis_min <- lyr_data %>% 
  dplyr::select(MineralGroupsNew, pro_MAT_mod, pro_MAP_mod, pro_PET_mm_yr_mod,
                lyr_clay_mod) %>% 
  rename("Clay content [%]" = lyr_clay_mod,
         "PET [mm]" = pro_PET_mm_yr_mod,
         "MAP [mm]" = pro_MAP_mod,
         "MAT [°C]" = pro_MAT_mod) %>% 
  pivot_longer(!MineralGroupsNew, names_to = "names", values_to = "values") %>% 
  group_by(MineralGroupsNew) %>% 
  ggplot(aes(x = MineralGroupsNew, y = values, fill = MineralGroupsNew)) +
  geom_boxplot() +
  # geom_violin() +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        strip.background = element_rect(fill = "white")) +
  facet_wrap(~names, scales = "free", nrow = 1) +
  scale_fill_manual("Mineral groups", values = color_mineral)

# Figure 2
ggarrange(dis_clim, dis_min, nrow = 2, labels = c("a)", "b)"), vjust = c(5,5))
ggsave(file = paste0("./Figure/Figure_2_", 
                     Sys.Date(),".jpeg"), width = 9, height = 6.5, dpi = 500)

