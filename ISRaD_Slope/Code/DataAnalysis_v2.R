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

# Load filtered database data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2023-02-08"))

lyr_all %>% 
  skimr::skim_without_charts(lyr_14c, CORG)

# Load filtered and splined data
lyr_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-03-09.csv") %>% 
  #remove layers that have CORG > 20
  filter(CORG_msp <= 20)

# Check data
lyr_data %>% 
  dplyr::summarise(n_studies = n_distinct(entry_name),
                   n_sites = n_distinct(site_name),
                   n_profiles = n_distinct(id))

lyr_data %>% 
  dplyr::select(lyr_14c_msp, CORG_msp, UD, pro_MAT_mod, pro_AI, lyr_clay_mod) %>% 
  cor()

# Convert characters into factors
lyr_data$ClimateZone <- factor(lyr_data$ClimateZone,
                               levels = c("tundra/polar", "cold temperate",
                                          "warm temperate", "tropical", "arid"))

lyr_data$ClimateZoneAnd <- factor(lyr_data$ClimateZoneAnd,
                                  levels = c("volcanic soils", "tundra/polar", 
                                             "cold temperate", "warm temperate", 
                                             "tropical", "arid"))

lyr_data %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::summarise(n_studies = n_distinct(entry_name),
                   n_sites = n_distinct(site_name),
                   n_profiles = n_distinct(id),
                   n_countries = n_distinct(pro_country))

lyr_data %>% 
  group_by(MineralType) %>% 
  dplyr::summarise(n_studies = n_distinct(entry_name),
                   n_sites = n_distinct(site_name),
                   n_profiles = n_distinct(id),
                   n_countries = n_distinct(pro_country))

## Define color code 
# cold, warm, tropical, arid
color_climate_wo_polar <- c("#762a83", "#7fbf7b", "#1b7837", "#dfc27d")
# polar, cold, warm, tropical, arid
color_climate <- c("#af8dc3", "#762a83", "#7fbf7b", "#1b7837", "#dfc27d")
# amorphous, polar, cold, warm, tropical, arid
color_climate_w_amorph <- c("#225ea8", "#af8dc3", "#762a83", "#7fbf7b", "#1b7837", 
                            "#dfc27d")
# amorphous, high-activity clay, low-activity clays
color_mineral <- c("#225ea8", "#41b6c4", "#a1dab4")

## Mapping sampling locations ##
world <- map_data("world") %>% 
  filter(region != "Antarctica")

set.seed(123)
ggplot() +  
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "#d9d9d9", fill = "#d9d9d9") +
  geom_jitter(data = lyr_data %>% 
                distinct(id, .keep_all = TRUE), 
              aes(x = pro_long, y = pro_lat, fill = ClimateZoneAnd),
              size = 2, shape = 21, color = "white",
              width = 1, height = 1) +
  theme_bw(base_size = 14) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = "top") +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), 
                     breaks = c(-100,0,100), limits = c(-170,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), 
                     breaks = c(-50,0,50), limits = c(-55,80)) +
  scale_fill_manual("Grouping", values = color_climate_w_amorph)

ggsave(file = paste0("./Figure/ISRaD_14C_map_", Sys.Date(),
                     ".jpeg"), width = 8, height = 5)

## Spline profile example ##
# Data example
p1 <- ggplot() +
  geom_line(data = lyr_data %>% 
              filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
            aes(y = UD, x = lyr_14c_msp, group = id, color = id), orientation = "y",
            linewidth = 1.5) +
  geom_pointrange(data = lyr_all %>% 
                    filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
                  aes(x = lyr_14c, y = depth, ymin = lyr_top, ymax = lyr_bot,
                      color = id), size = 1) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank()) +
  scale_y_reverse("Depth [cm]", expand = c(0,0)) +
  scale_x_continuous(expression(paste(Delta^14,"C [‰]")), 
                     expand = c(0,0), limits = c(-1000,200),
                     position = "top") +
  scale_color_discrete(paste0("Example profile:"), 
                       label = c("P1", "P2")) +
  coord_cartesian(ylim = c(100,0))

# p2 <- ggplot() +
#   geom_line(data = lyr_data %>% 
#               filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
#             aes(y = UD, x = CORG_msp, group = id, color = id), orientation = "y",
#             linewidth = 1.5) +
#   geom_pointrange(data = lyr_all %>% 
#                     filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
#                   aes(x = CORG, y = depth, ymin = lyr_top, ymax = lyr_bot,
#                       color = id), size = 1) +
#   theme_bw(base_size = 16) +
#   theme(axis.text = element_text(color = "black"),
#         panel.grid.minor = element_blank()) +
#   scale_y_reverse("", expand = c(0,0)) +
#   scale_x_continuous("SOC [wt-%]; log-scaled", expand = c(0,0), position = "top",
#                      trans = "log10") +
#   scale_color_discrete(paste0("Example profile:"), 
#                        label = c("P1", "P2")) +
#   annotation_logticks(sides = "t", 
#                       short = unit(1.5,"mm"),
#                       mid = unit(3,"mm"),
#                       long = unit(4,"mm"))

p2 <- ggplot() +
  geom_line(data = lyr_data %>% 
              filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
            aes(y = UD, x = CORG_msp, group = id, color = id), orientation = "y",
            linewidth = 1.5) +
  geom_pointrange(data = lyr_all %>% 
                    filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
                  aes(x = CORG, y = depth, ymin = lyr_top, ymax = lyr_bot,
                      color = id), size = 1) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank()) +
  scale_y_reverse("", expand = c(0,0)) +
  scale_x_continuous("SOC [wt-%]; log-scaled", expand = c(0,0), position = "top",
                     breaks = c(0.1,1,5,10),
                     labels = c(0.1,1,5,10)) +
  coord_trans(x = "log10", ylim = c(100,0), xlim = c(0.05,10)) +
  scale_color_discrete(paste0("Example profile:"), 
                       label = c("P1", "P2")) +
  annotation_logticks(sides = "t", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm"))

ggarrange(p1, p2, common.legend = TRUE)
# ggsave(file = paste0("./Figure/ISRaD_14C_SOC_depth_example_", Sys.Date(),
#                      ".jpeg"), width = 10, height = 5)

### Plot profile data ###
## Climate data
climate_14c_c <- lyr_data %>%
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(median_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
                lci_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
                uci_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
                median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
                lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
                uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
                n = n(),
                n_site = n_distinct(site_name)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  dplyr::mutate(n_rel = n * 100 / max(n))

climate_14c_c %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::summarise(median_14c_max = max(median_14c),
                   median_14c_min = min(median_14c),
                   median_c_max = max(median_c),
                   median_c_min = min(median_c)) %>% 
  rowwise() %>% 
  mutate(diff_14c = abs(median_14c_max - median_14c_min),
         diff_c = abs(median_c_max - median_c_min),
         rel_change_c = 100-(median_c_min * 100/median_c_max))

depth_sum_c <- climate_14c_c %>% 
  dplyr::filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
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
                  UD == 100)

depth_sum_c %>% 
  group_by(ClimateZoneAnd) %>% 
  filter(UD == max(UD))

c1_14c_c <- climate_14c_c %>% 
  dplyr::filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), linewidth = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  geom_point(data = depth_sum_c, aes(x = median_c, y = median_14c, 
                                     shape = as.factor(UD)), 
             size = 2, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 25, b = 10)) +
  scale_x_continuous("Soil organic carbon [wt-%]; log-scaled", expand = c(0,0), 
                     breaks = c(0.1,0.5,1,5),
                     labels = c(0.1,0.5,1,5)) +
  coord_trans(x = "log10", xlim = c(0.1,6)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_manual("Climate grouping", values = color_climate) +
  guides(shape = "none")

c1_14c <- climate_14c_c %>% 
  dplyr::filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(y = UD, x = median_14c, color = ClimateZoneAnd), linewidth = 2) +
  geom_errorbarh(aes(xmin = lci_14c, xmax = uci_14c, y = UD, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.8),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 15, b = 5, t = 10, l = 5),
        axis.title.y = element_blank()) +
  scale_y_reverse(expand = c(0,0), limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_manual("Climate grouping", values = color_climate)

c1_c <- climate_14c_c %>% 
  dplyr::filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(y = UD, x = median_c, color = ClimateZoneAnd), linewidth = 2) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = UD, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8,0.8),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 15, b = 10, t = 5, l = 5),
        axis.title.y = element_blank()) +
  scale_y_reverse(expand = c(0,0), limits = c(100,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", expand = c(0,0),
                     limits = c(0,6)) +
  scale_color_manual("Climate grouping", values = color_climate)

ggarrange(annotate_figure(ggarrange(c1_14c, c1_c, nrow = 2, legend = "none",
                                    labels = c("a)", "b)"), vjust = c(2.2,1.8),
                                    hjust = c(-3,-2.8)),
                          left = text_grob("Depth [cm]", rot = 90, size = 16)), 
          c1_14c_c, common.legend = TRUE, widths = c(0.8,1), labels = c("", "c)"),
          vjust = c(1.5,2.2), hjust = c(-0.5,-5.8))
ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14C_climateAnd_depth_", Sys.Date(),
                     ".jpeg"), width = 10, height = 6)

## Mineral data
mineral_14c_c <- lyr_data %>%
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

mineral_14c_c %>% 
  group_by(MineralType) %>% 
  dplyr::summarise(median_14c_max = max(median_14c),
                   median_14c_min = min(median_14c),
                   median_c_max = max(median_c),
                   median_c_min = min(median_c)) %>% 
  rowwise() %>% 
  mutate(diff_14c = median_14c_max - median_14c_min,
         diff_c = median_c_max - median_c_min)

depth_sum_m <- mineral_14c_c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
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
                  UD == 100)

depth_sum_m %>% 
  group_by(MineralType) %>% 
  filter(UD == max(UD))


m1_14c_c <- mineral_14c_c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = MineralType), linewidth = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralType), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralType), alpha = 0.3) +
  geom_point(data = depth_sum_m, aes(x = median_c, y = median_14c,
                                     shape = as.factor(UD)), 
             size = 2, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 25, b = 10)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_x_continuous("Soil organic carbon [wt-%]; log-scaled", expand = c(0,0), 
                     breaks = c(0.1,0.5,1,5,10),
                     labels = c(0.1,0.5,1,5,10)) +
  coord_trans(x = "log10", xlim = c(0.1,16)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_manual("Mineral grouping", values = color_mineral) +
  guides(shape = "none")

m1_14c <- mineral_14c_c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(y = UD, x = median_14c, color = MineralType), linewidth = 2) +
  geom_errorbarh(aes(xmin = lci_14c, xmax = uci_14c, y = UD, 
                     color = MineralType), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.8),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 15, b = 5, t = 10, l = 5),
        axis.title.y = element_blank()) +
  scale_y_reverse(expand = c(0,0), limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_manual("Mineral grouping", values = color_mineral) 

m1_c <- mineral_14c_c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(y = UD, x = median_c, color = MineralType), linewidth = 2) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = UD, 
                     color = MineralType), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8,0.5),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 15, b = 10, t = 5, l = 5),
        axis.title.y = element_blank()) +
  scale_y_reverse(expand = c(0,0), limits = c(100,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", expand = c(0,0),
                     limits = c(0,16)) +
  scale_color_manual("Mineral grouping", values = color_mineral) 

ggarrange(annotate_figure(ggarrange(m1_14c, m1_c, nrow = 2, legend = "none",
                                    labels = c("a)", "b)"), vjust = c(2.2,1.8),
                                    hjust = c(-3,-2.8)),
                          left = text_grob("Depth [cm]", rot = 90, size = 16)), 
          m1_14c_c, common.legend = TRUE, widths = c(0.8,1), labels = c("", "c)"),
          vjust = c(1.5,2.2), hjust = c(-0.5,-5.8))
ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14C_mineral_depth_", Sys.Date(),
                     ".jpeg"), width = 10, height = 6)

### Random Forest Analysis ###
set.seed(42)
lyr_data_sf <- sf::st_as_sf(lyr_data, coords = c("pro_long", "pro_lat"), crs = 4326)

# Filter data (only keep variables of interest)
rf_data_14c <- lyr_data_sf %>% 
  arrange(id) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

rf_data_c <- lyr_data_sf %>% 
  arrange(id) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

## Set-up random forest
task_rf_14c <- as_task_regr_st(x = rf_data_14c %>%
                              dplyr::select(-entry_name),
                            target = "lyr_14c_msp")

# task_rf_14c <- as_task_regr(x = rf_data_14c %>%
#                               dplyr::select(-entry_name),
#                             target = "lyr_14c_msp")

lrn_rf_14c <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000)

# Add id as group for CV (same id kept together)
task_rf_14c$set_col_roles("id", roles = "group")
print(task_rf_14c)

task_rf_c <- as_task_regr_st(x = rf_data_c %>% 
                            dplyr::select(-entry_name), 
                          target = "CORG_msp")

lrn_rf_c <- lrn("regr.ranger", importance = "permutation",
                num.trees = 1000)

# Add id as group for CV (same id kept together)
task_rf_c$set_col_roles("id", roles = "group")
print(task_rf_c)

# cross-validation
resampling_14c <- rsmp("cv", folds = 10)
resampling_14c$instantiate(task_rf_14c)

# autoplot(resampling_14c, task_rf_14c, fold_id = c(1,2,5,10)) +
#   theme_bw(base_size = 14) +
#   theme(axis.text = element_blank())
# 
# ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14C_RF_cv_", Sys.Date(),
#                      ".jpeg"), width = 12, height = 8)

resampling_c <- rsmp("cv", folds = 10)
resampling_c$instantiate(task_rf_c)

# autoplot(resampling_c, task_rf_c, fold_id = c(1,2,5,10))

## Train model & check performance
rr_14c <- mlr3::resample(task = task_rf_14c, learner = lrn_rf_14c, 
                         resampling = resampling_14c, store_models = TRUE)

rr_14c$aggregate(measures = msrs(c("regr.rmse", "regr.mse", "regr.rsq")))
rr_score_14c_ind <- rr_14c$score(measures = msrs(c("regr.rmse", "regr.mse", 
                                                   "regr.rsq")))

rr_14c_measure <- data.frame(regr_rmse = rr_score_14c_ind$regr.rmse,
                             regr_mse = rr_score_14c_ind$regr.mse,
                             regr_rsq = rr_score_14c_ind$regr.rsq)

write_csv(rr_14c_measure, file = paste0("./Data/ISRaD_RF_wo_andipol_14c_measures_", 
                                        Sys.Date(), ".csv"))

rr_14c_measure <- read_csv("./Data/ISRaD_RF_wo_andipol_14c_measures_2023-03-20.csv")


rm(rr_score_14c_ind)

rr_pred_14c <- rr_14c$prediction(predict_sets = "test")
rr_14c_prediction <- data.frame(truth = rr_pred_14c$truth,
                                response = rr_pred_14c$response)

write_csv(rr_14c_prediction, file = paste0("./Data/ISRaD_RF_wo_andipol_14c_prediction_",
                                           Sys.Date(), ".csv"))

rr_14c_prediction <- read_csv("./Data/ISRaD_RF_wo_andipol_14c_prediction_2023-03-20.csv")


rm(rr_pred_14c)

vi_14c_all <- lapply(rr_14c$learners, function(x) x$model$variable.importance) 

vi_14c_sum <- vi_14c_all %>%   
  plyr::ldply() %>% 
  mutate(row_sum = rowSums(pick(where(is.numeric)))) %>% 
  mutate(UD = UD/row_sum*100,
         lyr_clay_mod = lyr_clay_mod/row_sum*100,
         pro_AI = pro_AI/row_sum*100,
         pro_MAT_mod = pro_MAT_mod/row_sum*100) %>% 
  dplyr::select(-row_sum) %>% 
  dplyr::summarise(across(everything(), list(median = median, mad = mad), 
                          .names = "{.col}.{.fn}")) %>% 
  pivot_longer(everything(), names_to = "names", values_to = "values") %>% 
  separate_wider_delim(names, ".", names = c("predictor", "measure")) %>% 
  pivot_wider(names_from = measure, values_from = values)

write_csv(vi_14c_sum, file = paste0("./Data/ISRaD_RF_wo_andipol_14c_VaribleImportance_", 
                                  Sys.Date(), ".csv"))

vi_14c_sum <- read_csv("./Data/ISRaD_RF_wo_andipol_14c_VaribleImportance_2023-03-20.csv")

rm(vi_14c_all)
rm(rr_14c)

p_rf_14c <- rr_14c_prediction %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous(expression(paste("Observed ", Delta^14,"C [‰]")), 
                     limits = c(-1000,500), breaks = seq(-1000,500,250)) +
  scale_x_continuous(expression(paste("Predicted ", Delta^14,"C [‰]")), 
                     limits = c(-1000,300), breaks = seq(-1000,250,250)) +
  geom_text(aes(x = -720, y = 430), 
            label = paste("R² = ", round(median(rr_14c_measure$regr_rsq), 2), 
                          " ± ", round(mad(rr_14c_measure$regr_rsq), 2))) +
  geom_text(aes(x = -720, y = 370), 
            label = paste("MSE = ", round(median(rr_14c_measure$regr_mse), 0), 
                          " ± ", round(mad(rr_14c_measure$regr_mse), 0))) +
  geom_text(aes(x = -720, y = 310), 
            label = paste("RMSE = ", round(median(rr_14c_measure$regr_rmse), 0), 
                          " ± ", round(mad(rr_14c_measure$regr_rmse), 0), "‰"))

rr_c <- mlr3::resample(task = task_rf_c, learner = lrn_rf_c, 
                       resampling = resampling_c, store_models = TRUE)

rr_c$aggregate(measures = msrs(c("regr.rmse", "regr.mse", "regr.rsq")))

rr_score_c_ind <- rr_c$score(measures = msrs(c("regr.rmse", "regr.mse", "regr.rsq")))
rr_c_measure <- data.frame(regr_rmse = rr_score_c_ind$regr.rmse,
                           regr_mse = rr_score_c_ind$regr.mse,
                           regr_rsq = rr_score_c_ind$regr.rsq)

write_csv(rr_c_measure, file = paste0("./Data/ISRaD_RF_wo_andipol_c_measures_", 
                                      Sys.Date(), ".csv"))

rr_c_measure <- read_csv("./Data/ISRaD_RF_wo_andipol_c_measures_2023-03-20.csv")

rm(rr_score_c_ind)

rr_pred_c <- rr_c$prediction(predict_sets = "test")
rr_c_prediction <- data.frame(truth = rr_pred_c$truth,
                              response = rr_pred_c$response)

write_csv(rr_c_prediction, file = paste0("./Data/ISRaD_RF_wo_andipol_c_prediction_", 
                                         Sys.Date(), ".csv"))

rr_c_prediction <- read_csv("./Data/ISRaD_RF_wo_andipol_c_prediction_2023-03-20.csv")

rm(rr_pred_c)

vi_c_all <- lapply(rr_c$learners, function(x) x$model$variable.importance)

vi_c_sum <- vi_c_all %>%   
  plyr::ldply() %>% 
  mutate(row_sum = rowSums(pick(where(is.numeric)))) %>% 
  mutate(UD = UD/row_sum*100,
         lyr_clay_mod = lyr_clay_mod/row_sum*100,
         pro_AI = pro_AI/row_sum*100,
         pro_MAT_mod = pro_MAT_mod/row_sum*100) %>% 
  dplyr::select(-row_sum) %>% 
  dplyr::summarise(across(everything(), list(median = median, mad = mad), 
                          .names = "{.col}.{.fn}")) %>% 
  pivot_longer(everything(), names_to = "names", values_to = "values") %>% 
  separate_wider_delim(names, ".", names = c("predictor", "measure")) %>% 
  pivot_wider(names_from = measure, values_from = values)

write_csv(vi_c_sum, file = paste0("./Data/ISRaD_RF_wo_andipol_c_VaribleImportance_", 
                                         Sys.Date(), ".csv"))

vi_c_sum <- read_csv("./Data/ISRaD_RF_wo_andipol_c_VaribleImportance_2023-03-20.csv")

rm(vi_c_all)
rm(rr_c)

p_rf_c <- rr_c_prediction %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed SOC [wt-%]", limits = c(0,25)) +
  scale_x_continuous("Predicted SOC [wt-%]", limits = c(0,20)) +
  geom_text(aes(x = 7, y = 24), 
            label = paste("R² = ", round(median(rr_c_measure$regr_rsq), 2), 
                          " ± ", round(mad(rr_c_measure$regr_rsq), 2))) +
  geom_text(aes(x = 7, y = 23), 
            label = paste("MSE = ", round(median(rr_c_measure$regr_mse), 2), 
                          " ± ", round(mad(rr_c_measure$regr_mse), 2))) +
  geom_text(aes(x = 7, y = 22), 
            label = paste("RMSE = ", round(median(rr_c_measure$regr_rmse), 2), 
                          " ± ", round(mad(rr_c_measure$regr_rmse), 2), "wt-%"))
  
ggarrange(p_rf_14c, p_rf_c)
ggsave(file = paste0("./Figure/ISRaD_msp_RF_Perform_14C_SOC_", Sys.Date(),
                     ".jpeg"), width = 10, height = 5)

#Variable importance
vi_14c <- vi_14c_sum %>% 
  ggplot() +
  geom_bar(aes(x = reorder(predictor, -median), y = median, fill = predictor),
           stat = "identity") +
  geom_errorbar(aes(ymin = median-mad, ymax = median+mad, 
                    x = reorder(predictor, -median)), width = 0.2) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_y_continuous("Relative importance [%]", expand = c(0,0),
                     limits = c(0,45)) +
  scale_x_discrete("", labels = c("Depth", "MAT", "PET/MAP", "Clay content")) +
  scale_fill_manual(values = c("#a6611a", "#80cdc1", "#018571", "#dfc27d"))
  
vi_c <- vi_c_sum  %>% 
  ggplot() +
  geom_bar(aes(x = reorder(predictor, -median), y = median, fill = predictor),
           stat = "identity") +
  geom_errorbar(aes(ymin = median-mad, ymax = median+mad, 
                    x = reorder(predictor, -median)), width = 0.2) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(0,45)) +
  scale_x_discrete("", labels = c("PET/MAP", "MAT", "Depth", "Clay content")) +
  scale_fill_manual(values = c("#a6611a", "#80cdc1", "#018571", "#dfc27d"))
  
ggarrange(vi_14c, vi_c, labels = c("a) Radiocarbon model",  "b) SOC model"), 
          vjust = 2.5, hjust = c(-0.35,-0.6))
ggsave(file = paste0("./Figure/ISRaD_msp_vi_14c_SOC_", Sys.Date(),
                     ".jpeg"), width = 10, height = 5)

#Partial Dependence plots
#not working with spatial data
rf_data_14c_pdp <- lyr_data %>% 
  arrange(id) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

task_rf_14c_pdp <- as_task_regr(x = rf_data_14c_pdp %>%
                                  dplyr::select(-entry_name, -id),
                                target = "lyr_14c_msp")

lrn_rf_14c_pdp <- lrn("regr.ranger", importance = "permutation",
                      num.trees = 1000) 

lrn_rf_14c_pdp$train(task_rf_14c_pdp)

model_14c <- Predictor$new(lrn_rf_14c_pdp, data = rf_data_14c_pdp %>%  
                             dplyr::select(-id, -entry_name))

effect_14c_ice_MAT <- FeatureEffects$new(model_14c, method = "ice",
                                         features = "pro_MAT_mod")

ice_mat <- lyr_data %>% 
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>%
  rownames_to_column(var = ".id") %>%
  mutate(.id = as.integer(.id)) %>%
  full_join(effect_14c_ice_MAT$effects$pro_MAT_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_MAT_mod) %>%
  dplyr::mutate(pred_mean_14c = median(.value)) %>%
  ungroup() %>%
  ggplot() +
  geom_path(aes(x = pro_MAT_mod, y = pred_mean_14c, group = .id, color = ClimateZoneAnd),
            linewidth = 1) +
  geom_rug(data = lyr_data %>% 
             filter(ClimateZoneAnd != "volcanic soils") %>%
             filter(ClimateZoneAnd != "tundra/polar") %>%
             distinct(pro_MAT_mod, .keep_all = TRUE), 
           aes(x = pro_MAT_mod, color = ClimateZoneAnd), sides = "b") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("MAT [°C]") +
  scale_y_continuous(limits = c(-520,100), expand = c(0,0)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

# AI_bin <- lyr_data %>%
#   arrange(id) %>%
#   dplyr::filter(ClimateZoneAnd != "volcanic soils") %>%
#   dplyr::filter(ClimateZoneAnd != "tundra/polar") %>%
#   distinct(id, .keep_all = TRUE) %>%
#   dplyr::mutate(bin_AI = cut_number(pro_AI, 20)) %>%
#   group_by(bin_AI) %>%
#   dplyr::mutate(bin_mean_AI = mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "",
#                                                             as.character(bin_AI)), ","))))) %>%
#   ungroup() %>%
#   distinct(bin_mean_AI, .keep_all = TRUE) %>%
#   arrange(bin_mean_AI) %>%
#   dplyr::select(bin_AI, bin_mean_AI)


# effect_14c_ice_AI <- FeatureEffect$new(model_14c, method = "ice",
#                                         feature = "pro_AI", 
#                                         grid.points = as.double(AI_bin$bin_mean_AI))

effect_14c_ice_AI <- FeatureEffects$new(model_14c, method = "ice",
                                        features = "pro_AI")

ice_AI <- lyr_data %>%
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>%
  rownames_to_column(var = ".id") %>%
  mutate(.id = as.integer(.id)) %>%
  full_join(effect_14c_ice_AI$effects$pro_AI$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_AI) %>%
  dplyr::mutate(pred_mean_14c = median(.value)) %>%
  ungroup() %>%
  ggplot() +
  geom_path(aes(x = pro_AI, y = pred_mean_14c, group = .id, color = ClimateZoneAnd),
            linewidth = 1) +
  geom_rug(data = lyr_data %>% 
             filter(ClimateZoneAnd != "volcanic soils") %>%
             filter(ClimateZoneAnd != "tundra/polar") %>%
             distinct(pro_AI, .keep_all = TRUE), 
           aes(x = pro_AI, color = ClimateZoneAnd), sides = "b") +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("PET/MAP") +
  scale_y_continuous(limits = c(-520,100), expand = c(0,0)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_wo_AndiPolar_ice_AI_5GridsMan_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

#UD
effect_14c_ice_UD <- FeatureEffects$new(model_14c, method = "ice", 
                                        features = "UD")

ice_UD <- lyr_data %>% 
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_14c_ice_UD$effects$UD$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(pred_mean_14c = median(.value)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_path(aes(x = UD, y = pred_mean_14c, group = .id, color = ClimateZoneAnd),
            linewidth = 1) +
  geom_rug(data = lyr_data %>% 
             filter(ClimateZoneAnd != "volcanic soils") %>%
             filter(ClimateZoneAnd != "tundra/polar") %>%
             distinct(UD, .keep_all = TRUE), 
           aes(x = UD, color = ClimateZoneAnd), sides = "b") +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous(limits = c(-520,100), expand = c(0,0)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

#clay
effect_14c_ice_clay <- FeatureEffects$new(model_14c, method = "ice", 
                                          features = "lyr_clay_mod")

ice_clay <- lyr_data %>% 
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_14c_ice_clay$effects$lyr_clay_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::mutate(pred_mean_14c = median(.value)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_path(aes(x = lyr_clay_mod, y = pred_mean_14c, group = .id, color = ClimateZoneAnd),
            linewidth = 1) +
  geom_rug(data = lyr_data %>% 
             filter(ClimateZoneAnd != "volcanic soils") %>%
             filter(ClimateZoneAnd != "tundra/polar") %>%
             distinct(lyr_clay_mod, .keep_all = TRUE), 
           aes(x = lyr_clay_mod, color = ClimateZoneAnd), sides = "b") +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("Clay content [%]") +
  scale_y_continuous(limits = c(-520,100), expand = c(0,0)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

# ggarrange(ice_UD, ice_clay, ice_mat, ice_AI, common.legend = TRUE)

annotate_figure(ggarrange(ice_UD, ice_clay, ice_mat, ice_AI, common.legend = TRUE),
                left = text_grob(expression(paste("Median predicted ", Delta^14, "C [‰]")), 
                                 rot = 90, size = 14))

ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_wo_AndiPolar_ice_climate_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

#SOC
rf_data_c_pdp <- lyr_data %>% 
  arrange(id) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

task_rf_c_pdp <- as_task_regr(x = rf_data_c_pdp %>%
                                dplyr::select(-entry_name, -id),
                              target = "CORG_msp")

lrn_rf_c_pdp <- lrn("regr.ranger", importance = "permutation",
                      num.trees = 1000) 

lrn_rf_c_pdp$train(task_rf_c_pdp)

model_c <- Predictor$new(lrn_rf_c_pdp, data = rf_data_c_pdp %>%  
                             dplyr::select(-id, -entry_name))

effect_c_ice_MAT <- FeatureEffects$new(model_c, method = "ice",
                                       features = "pro_MAT_mod")

ice_mat_c <- lyr_data %>%
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>%
  rownames_to_column(var = ".id") %>%
  mutate(.id = as.integer(.id)) %>%
  full_join(effect_c_ice_MAT$effects$pro_MAT_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_MAT_mod) %>%
  dplyr::mutate(pred_mean_c = median(.value)) %>%
  ungroup() %>%
  ggplot() +
  geom_path(aes(x = pro_MAT_mod, y = pred_mean_c, group = .id, color = ClimateZoneAnd),
            linewidth = 1) +
  geom_rug(data = lyr_data %>% 
             filter(ClimateZoneAnd != "volcanic soils") %>%
             filter(ClimateZoneAnd != "tundra/polar") %>%
             distinct(pro_MAT_mod, .keep_all = TRUE), 
           aes(x = pro_MAT_mod, color = ClimateZoneAnd), sides = "b") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("MAT [°C]") +
  scale_y_continuous(limits = c(-0.1,9), expand = c(0,0)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

effect_c_ice_AI <- FeatureEffects$new(model_c, method = "ice",
                                      features = "pro_AI")

ice_AI_c <- lyr_data %>%
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>%
  rownames_to_column(var = ".id") %>%
  mutate(.id = as.integer(.id)) %>%
  full_join(effect_c_ice_AI$effects$pro_AI$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_AI) %>%
  dplyr::mutate(pred_mean_c = median(.value)) %>%
  ungroup() %>%
  ggplot() +
  geom_path(aes(x = pro_AI, y = pred_mean_c, group = .id, color = ClimateZoneAnd),
            linewidth = 1) +
  geom_rug(data = lyr_data %>% 
             filter(ClimateZoneAnd != "volcanic soils") %>%
             filter(ClimateZoneAnd != "tundra/polar") %>%
             distinct(pro_AI, .keep_all = TRUE), 
           aes(x = pro_AI, color = ClimateZoneAnd), sides = "b") +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("PET/MAP") +
  scale_y_continuous(limits = c(-0.1,9), expand = c(0,0)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

#UD
effect_c_ice_UD <- FeatureEffects$new(model_c, method = "ice", 
                                      features = "UD")

ice_UD_c <- lyr_data %>% 
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_c_ice_UD$effects$UD$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(pred_mean_c = median(.value)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_path(aes(x = UD, y = pred_mean_c, group = .id, color = ClimateZoneAnd),
            linewidth = 1) +
  geom_rug(data = lyr_data %>% 
             filter(ClimateZoneAnd != "volcanic soils") %>%
             filter(ClimateZoneAnd != "tundra/polar") %>%
             distinct(UD, .keep_all = TRUE), 
           aes(x = UD, color = ClimateZoneAnd), sides = "b") +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous(limits = c(-0.1,9), expand = c(0,0)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

#clay
effect_c_ice_clay <- FeatureEffects$new(model_c, method = "ice", 
                                        features = "lyr_clay_mod")

ice_clay_c <- lyr_data %>% 
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_c_ice_clay$effects$lyr_clay_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::mutate(pred_mean_c = median(.value)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_path(aes(x = lyr_clay_mod, y = pred_mean_c, group = .id, color = ClimateZoneAnd),
            linewidth = 1) +
  geom_rug(data = lyr_data %>% 
             filter(ClimateZoneAnd != "volcanic soils") %>%
             filter(ClimateZoneAnd != "tundra/polar") %>%
             distinct(lyr_clay_mod, .keep_all = TRUE), 
           aes(x = lyr_clay_mod, color = ClimateZoneAnd), sides = "b") +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("Clay content [%]") +
  scale_y_continuous(limits = c(-0.1,9), expand = c(0,0)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

annotate_figure(ggarrange(ice_UD_c, ice_clay_c, ice_mat_c, ice_AI_c, 
                          common.legend = TRUE),
                left = text_grob("Median predicted SOC [wt-%]", rot = 90, 
                                 size = 14))

ggsave(file = paste0("./Figure/ISRaD_msp_predSOC_AndiPolar_ice_climate_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

#Merge both PDP's
mat_ice <- lyr_data %>% 
  arrange(id) %>%
  dplyr::filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  dplyr::mutate(.id = as.integer(.id)) %>% 
  left_join(effect_14c_ice_MAT$effects$pro_MAT_mod$results %>% 
              dplyr::rename(pred_14c = .value), multiple = "all") %>% 
  left_join(effect_c_ice_MAT$effects$pro_MAT_mod$results %>% 
              dplyr::rename(pred_c = .value), multiple = "all") %>% 
  dplyr::group_by(ClimateZoneAnd, pro_MAT_mod) %>% 
  dplyr::mutate(pred_median_c = median(pred_c),
                pred_median_14c = median(pred_14c)) %>% 
  ungroup() %>% 
  arrange(pro_MAT_mod) %>% 
  dplyr::select(ClimateZoneAnd, pro_MAT_mod, pred_median_c, pred_median_14c) %>% 
  distinct(.keep_all = TRUE)
  
write_csv(mat_ice, file = paste0("./Data/ISRaD_RF_wo_andipol_ice_mat_14c_c_", 
                                 Sys.Date(), ".csv"))

mat_ice <- read_csv("./Data/ISRaD_RF_wo_andipol_ice_mat_14c_c_2023-03-15.csv")

ai_ice <- lyr_data %>% 
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(effect_14c_ice_AI$effects$pro_AI$results %>% 
              dplyr::rename(pred_14c = .value), multiple = "all") %>% 
  left_join(effect_c_ice_AI$effects$pro_AI$results %>% 
              dplyr::rename(pred_c = .value), multiple = "all") %>% 
  dplyr::group_by(ClimateZoneAnd, pro_AI) %>% 
  dplyr::mutate(pred_median_c = median(pred_c),
                pred_median_14c = median(pred_14c)) %>% 
  ungroup() %>% 
  arrange(pro_AI) %>% 
  dplyr::select(ClimateZoneAnd, pro_AI, pred_median_c, pred_median_14c) %>% 
  distinct(.keep_all = TRUE)

write_csv(ai_ice, file = paste0("./Data/ISRaD_RF_wo_andipol_ice_ai_14c_c_", 
                                Sys.Date(), ".csv"))

ai_ice <- read_csv("./Data/ISRaD_RF_wo_andipol_ice_ai_14c_c_2023-03-15.csv")

ud_ice <- lyr_data %>% 
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(effect_14c_ice_UD$effects$UD$results %>% 
              dplyr::rename(pred_14c = .value), multiple = "all") %>% 
  left_join(effect_c_ice_UD$effects$UD$results %>% 
              dplyr::rename(pred_c = .value), multiple = "all") %>% 
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(pred_median_c = median(pred_c),
                pred_median_14c = median(pred_14c)) %>% 
  ungroup() %>% 
  arrange(UD) %>% 
  dplyr::select(ClimateZoneAnd, UD, pred_median_c, pred_median_14c) %>% 
  distinct(.keep_all = TRUE)

write_csv(ud_ice, file = paste0("./Data/ISRaD_RF_wo_andipol_ice_ud_14c_c_", 
                                Sys.Date(), ".csv"))

ud_ice <- read_csv("./Data/ISRaD_RF_wo_andipol_ice_ud_14c_c_2023-03-15.csv")

clay_ice <- lyr_data %>% 
  arrange(id) %>%
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(effect_14c_ice_clay$effects$lyr_clay_mod$results %>% 
              dplyr::rename(pred_14c = .value), multiple = "all") %>% 
  left_join(effect_c_ice_clay$effects$lyr_clay_mod$results %>% 
              dplyr::rename(pred_c = .value), multiple = "all") %>% 
  dplyr::group_by(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::mutate(pred_median_c = median(pred_c),
                pred_median_14c = median(pred_14c)) %>% 
  ungroup() %>% 
  arrange(lyr_clay_mod) %>% 
  dplyr::select(ClimateZoneAnd, lyr_clay_mod, pred_median_c, pred_median_14c) %>% 
  distinct(.keep_all = TRUE)

write_csv(clay_ice, file = paste0("./Data/ISRaD_RF_wo_andipol_ice_clay_14c_c_", 
                                  Sys.Date(), ".csv"))

clay_ice <- read_csv("./Data/ISRaD_RF_wo_andipol_ice_clay_14c_c_2023-03-15.csv")

##PDP plots
mat_ice_pdp <- mat_ice %>% 
  rename(pred_value = pro_MAT_mod) %>% 
  mutate(predictor = "MAT [°C]")
ai_ice_pdp <- ai_ice %>% 
  rename(pred_value = pro_AI) %>% 
  mutate(predictor = "PET/MAP")
ud_ice_pdp <- ud_ice %>% 
  rename(pred_value = UD) %>% 
  mutate(predictor = "Depth [cm]")

pred_pdp <- rbind(mat_ice_pdp, ai_ice_pdp, ud_ice_pdp) %>% 
  tibble() %>% 
  mutate(ClimateZoneAnd = factor(ClimateZoneAnd, levels = c("cold temperate",
                                                            "warm temperate",
                                                            "tropical", "arid")))

lyr_data_rug <- lyr_data %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  group_by(ClimateZoneAnd) %>% 
  dplyr::distinct(ClimateZoneAnd, pro_MAT_mod, pro_AI) %>% 
  rename('PET/MAP' = pro_AI,
         'MAT [°C]' = pro_MAT_mod) %>% 
  pivot_longer(!ClimateZoneAnd, values_to = "pred_value", names_to = "predictor")

# Group by ClimateZoneAnd and then Summarise 95% of data
quantiles_climate <- lyr_data %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::distinct(ClimateZoneAnd, pro_MAT_mod, pro_AI) %>% 
  dplyr::summarise(P5_MAT = quantile(pro_MAT_mod, 0.025),
                   P95_MAT = quantile(pro_MAT_mod, 0.975),
                   P5_PET_MAP = quantile(pro_AI, 0.025),
                   P95_PET_MAP = quantile(pro_AI, 0.975))

quantiles_depth <- lyr_data %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::summarise(P5_depth = quantile(UD, 0.025),
                   P95_depth = quantile(UD, 0.975))

quantiles <- quantiles_climate %>% 
  left_join(quantiles_depth)

#Split the data based on the predictor so that we can create new vectors
#that can be interpolated onto
pdp_split <- split(pred_pdp, f = pred_pdp$predictor)

#create the higher resolution vectors to interpolate onto

#one for MAT
pred_value_high_res_mat <- seq(min(pdp_split$`MAT [°C]`$pred_value),
                               max(pdp_split$`MAT [°C]`$pred_value),
                               length.out = 200)

#one for PET_MAP
pred_value_high_res_pet <- seq(min(pdp_split$`PET/MAP`$pred_value),
                               max(pdp_split$`PET/MAP`$pred_value),
                               length.out = 200) #the higher this is the great the resolution

pred_value_high_res_depth <- seq(min(pdp_split$`Depth [cm]`$pred_value),
                                 max(pdp_split$`Depth [cm]`$pred_value),
                                 length.out = 200)

#Split the data again by ClimateZoneAnd
pdp_split_mat <- split(pdp_split$`MAT [°C]`, f = pdp_split$`MAT [°C]`$ClimateZoneAnd)
pdp_split_pet <- split(pdp_split$`PET/MAP`, f = pdp_split$`PET/MAP`$ClimateZoneAnd)
pdp_split_depth <- split(pdp_split$`Depth [cm]`, f = pdp_split$`Depth [cm]`$ClimateZoneAnd)

#Create a function that will compute a natural spline of the data based on the
#higher resolution vectors
spline_data_14c <- function(data, new_x) {
  
  splined_data <- approx(x = data$pred_value,
                         y = data$pred_median_14c,
                         method = "linear",
                         xout = new_x) |> data.frame()
  
  names(splined_data) <- c("pred_value", "pred_median_14c")
  
  splined_data$ClimateZoneAnd <- unique(data$ClimateZoneAnd)
  
  splined_data$predictor <- unique(data$predictor)
  
  return(splined_data)
  
}

spline_data_c <- function(data, new_x) {
  
  splined_data <- approx(x = data$pred_value,
                         y = data$pred_median_c,
                         method = "linear",
                         xout = new_x) |> data.frame()
  
  names(splined_data) <- c("pred_value", "pred_median_c")
  
  splined_data$ClimateZoneAnd <- unique(data$ClimateZoneAnd)
  
  splined_data$predictor <- unique(data$predictor)
  
  return(splined_data)
  
}

#Apply the spline function to the MAT data
pdp_split_mat_splined_14c <- lapply(pdp_split_mat, spline_data_14c,
                                    new_x = pred_value_high_res_mat)

pdp_split_mat_splined_c <- lapply(pdp_split_mat, spline_data_c,
                                  new_x = pred_value_high_res_mat)

#Apply the spline function to the PET_MAP data
pdp_split_pet_splined_14c <- lapply(pdp_split_pet, spline_data_14c,
                                    new_x = pred_value_high_res_pet)

pdp_split_pet_splined_c <- lapply(pdp_split_pet, spline_data_c,
                                  new_x = pred_value_high_res_pet)

pdp_split_depth_splined_14c <- lapply(pdp_split_depth, spline_data_14c,
                                      new_x = pred_value_high_res_depth)

pdp_split_depth_splined_c <- lapply(pdp_split_depth, spline_data_c,
                                    new_x = pred_value_high_res_depth)

#Bind them all together
pdp_splined_14c <- do.call(rbind, c(pdp_split_mat_splined_14c,
                                    pdp_split_pet_splined_14c,
                                    pdp_split_depth_splined_14c))

pdp_splined_c <- do.call(rbind, c(pdp_split_mat_splined_c,
                                  pdp_split_pet_splined_c,
                                  pdp_split_depth_splined_c))

pdp_splined <- pdp_splined_14c %>% 
  left_join(pdp_splined_c) %>% 
  tibble()

#Split them again by ClimateZoneAnd so that we can run the outlier loop
pdp_splined_split <- split(pdp_splined, f = pdp_splined$ClimateZoneAnd)

#Just a quick check that I'm dealing with data in the right order because
#I'm about to right a clunky for loop that uses data in the pdp list and in the
#quantiles data table.
identical(as.character(quantiles$ClimateZoneAnd), names(pdp_splined_split))

#Run the loop that assigns outliers
for (i in 1:length(pdp_splined_split )) {
  
  pdp_splined_split[[i]] <- pdp_splined_split[[i]] %>% 
    mutate(outlier = case_when(
    predictor == "MAT [°C]" & 
      (pred_value < quantiles$P5_MAT[i] | pred_value > quantiles$P95_MAT[i]) ~ TRUE,
    predictor == "PET/MAP" & 
      (pred_value < quantiles$P5_PET_MAP[i] | pred_value > quantiles$P95_PET_MAP[i]) ~ TRUE,
    predictor == "Depth [cm]" & 
      (pred_value < quantiles$P5_depth[i] | pred_value > quantiles$P95_depth[i]) ~ TRUE,
    .default = FALSE
  ))
  
}

pdp_splined <- do.call(rbind, pdp_splined_split) %>% 
  tibble()

pred_pdp_14c <- pdp_splined %>% 
  ggplot(aes(x = pred_value, color = ClimateZoneAnd)) +
  geom_line(aes(y = pred_median_14c), linewidth = 1.5, alpha = 0.3,
            linetype = "dashed") +
  geom_line(data = pdp_splined[-which(pdp_splined$outlier == TRUE), ], 
            aes(y = pred_median_14c),
            linewidth = 1.5) +
  geom_rug(data = lyr_data_rug, sides = "b") +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        plot.margin = margin(r = 5, l = 5, t = 5, b = 5),
        axis.title.x = element_blank()) +
  facet_wrap(~predictor, scales = "free_x") +
  scale_x_continuous("") +
  scale_y_continuous(expression(paste("Predicted  ", Delta^14, "C [‰]")),
                     expand = c(0,0), limits = c(-500,100)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

pred_pdp_c <- pdp_splined %>% 
  ggplot(aes(x = pred_value, color = ClimateZoneAnd)) +
  geom_line(aes(y = pred_median_c), linewidth = 1.5, alpha = 0.3, 
            linetype = "dashed") +
  geom_line(data = pdp_splined[-which(pdp_splined$outlier == TRUE), ], 
            aes(y = pred_median_c),
            linewidth = 1.5) +
  geom_rug(data = lyr_data_rug, sides = "b") +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(r = 5, l = 20, t = 5, b = 10),
        axis.title.x = element_blank()) +
  facet_wrap(~predictor, scales = "free_x") +
  scale_x_continuous("") +
  scale_y_continuous("Predicted SOC [wt-%]",
                     expand = c(0,0), limits = c(0.2,10),
                     breaks = c(0.5,1,5,10), labels = c(0.5,1,5,10)) +
  coord_trans(y = "log10") +
  annotation_logticks(sides = "l", scaled = FALSE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

ggarrange(pred_pdp_14c, pred_pdp_c, common.legend = TRUE, nrow = 2)
ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14c_AndiPolar_ice_climate_", 
                     Sys.Date(), ".jpeg"), width = 10, height = 7)


### Point PDP plots (95th data range)
## Clay content
# select 95% of the data for each climate group
clay_p_sum <- lyr_data %>% 
  dplyr::filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::filter(ClimateZoneAnd != "tundra/polar") %>%
  group_by(ClimateZoneAnd) %>% 
  dplyr::select(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::summarise(P5 = quantile(lyr_clay_mod, 0.025),
                   P95 = quantile(lyr_clay_mod, 0.975))

clay_p <- clay_ice %>%
  left_join(clay_p_sum) %>% 
  dplyr::mutate(
    lyr_clay_mod_p = case_when(
      lyr_clay_mod < P5 ~ NA,
      lyr_clay_mod > P95 ~ NA,
      TRUE ~ lyr_clay_mod),
    pred_median_14c_p = case_when(
      lyr_clay_mod < P5 ~ NA,
      lyr_clay_mod > P95 ~ NA,
      TRUE ~ pred_median_14c),
    pred_median_c_p = case_when(
      lyr_clay_mod < P5 ~ NA,
      lyr_clay_mod > P95 ~ NA,
      TRUE ~ pred_median_c))

clay_p$ClimateZoneAnd <- factor(clay_p$ClimateZoneAnd,
                                levels = c("cold temperate", "warm temperate",
                                           "tropical", "arid"))

clay_14c_p <- clay_p %>%
  ggplot(aes(color = ClimateZoneAnd)) +
  geom_path(aes(x = lyr_clay_mod, y = pred_median_14c), 
            linewidth = 1, alpha = 0.3, linetype = "dashed") +
  geom_path(aes(x = lyr_clay_mod_p, y = pred_median_14c_p), 
            linewidth = 1) +
  geom_point(aes(x = lyr_clay_mod_p, y = pred_median_14c_p, size = lyr_clay_mod_p)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        plot.margin = margin(r = 5, l = 15, t = 10)) +
  scale_y_continuous(expression(paste("Median predicted ", Delta^14,"C [‰]")), 
                     limits = c(-300,0), expand = c(0,0)) +
  scale_x_continuous("", expand = c(0.01,0.01)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar) +
  guides(color = guide_legend(override.aes = list(size = 2)))

clay_c_p <- clay_p %>% 
  ggplot(aes(color = ClimateZoneAnd)) +
  geom_path(aes(x = lyr_clay_mod, y = pred_median_c), 
            linewidth = 1, alpha = 0.3, linetype = "dashed") +
  geom_path(aes(x = lyr_clay_mod_p, y = pred_median_c_p), 
            linewidth = 1) +
  geom_point(aes(x = lyr_clay_mod_p, y = pred_median_c_p, size = lyr_clay_mod_p)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        plot.margin = margin(r = 5, l = 25, b = 10)) +
  scale_y_continuous("Median predicted SOC [wt-%]", 
                     limits = c(0,2), expand = c(0,0)) +
  scale_x_continuous("Clay content [%]", expand = c(0.01,0.01)) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar) +
  guides(color = guide_legend(override.aes = list(size = 2)))

clay_p %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::slice(which.min(lyr_clay_mod_p),
               which.max(lyr_clay_mod_p))

clay_14c_c_p <- clay_p %>% 
  ggplot() +
  geom_point(aes(y = pred_median_14c_p, x = pred_median_c_p, 
                 size = lyr_clay_mod_p, color = ClimateZoneAnd)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(r = 15, t = 10, b = 10, l = 5)) +
  scale_y_continuous(expression(paste("Median predicted ", Delta^14,"C [‰]")), 
                     limits = c(-300,0), expand = c(0,0)) +
  scale_x_continuous("Median predicted SOC [wt-%]", expand = c(0,0), 
                     limits = c(0.3,2), breaks = c(0.5,1,2), labels = c(0.5,1,2)) +
  coord_trans(x = "log10") +
  annotation_logticks(sides = "b", scaled = FALSE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  geom_segment(x = 0.35, y = -95, xend = 0.67, yend = -195,
               color = color_climate_wo_polar[4],
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(x = 0.73, y = -75, xend = 0.78, yend = -120,
               color = color_climate_wo_polar[3],
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(x = 0.86, y = -120, xend = 1.3, yend = -200,
               color = color_climate_wo_polar[2],
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(x = 1.1, y = -128, xend = 1.5, yend = -235,
               color = color_climate_wo_polar[1],
               arrow = arrow(length = unit(0.5, "cm"))) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar) + 
  scale_size_continuous("Clay content [%]") +
  guides(color = guide_legend(override.aes = list(size = 2)))
  
ggarrange(ggarrange(clay_14c_p, clay_c_p, nrow = 2, legend = "none"),
          clay_14c_c_p, common.legend = TRUE, widths = c(0.8,1))
ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14c_AndiPolar_ice_climate_clay_p_", 
                     Sys.Date(), ".jpeg"), width = 11, height = 6)

## AI
# select 95% of the data for each climate group
ai_p_sum <- lyr_data %>% 
  distinct(id, .keep_all = TRUE) %>% 
  dplyr::filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::filter(ClimateZoneAnd != "tundra/polar") %>%
  group_by(ClimateZoneAnd) %>% 
  dplyr::select(ClimateZoneAnd, pro_AI) %>% 
  dplyr::summarise(P5 = quantile(pro_AI, 0.05),
                   P95 = quantile(pro_AI, 0.95))

ai_p <- ai_ice %>%
  left_join(ai_p_sum) %>% 
  dplyr::mutate(
    pro_AI_p = case_when(
      pro_AI < P5 ~ NA,
      pro_AI > P95 ~ NA,
      TRUE ~ pro_AI),
    pred_median_14c_p = case_when(
      pro_AI < P5 ~ NA,
      pro_AI > P95 ~ NA,
      TRUE ~ pred_median_14c),
    pred_median_c_p = case_when(
      pro_AI < P5 ~ NA,
      pro_AI > P95 ~ NA,
      TRUE ~ pred_median_c))

ai_p %>% 
  dplyr::filter(ClimateZoneAnd == "tropical") %>% 
  approxfun(x = pro_AI, y = pred_median_14c)

ai_14c_p <- ai_p %>%
  ggplot(aes(color = ClimateZoneAnd, x = pro_AI, y = pred_median_14c)) +
  geom_path(linewidth = 1, alpha = 0.3, linetype = "dashed") +
  geom_path(aes(x = pro_AI_p, y = pred_median_14c_p), 
            linewidth = 1) +
  geom_point(aes(x = pro_AI_p, y = pred_median_14c_p, size = pro_AI_p)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        plot.margin = margin(r = 5, l = 15, t = 10)) +
  scale_y_continuous(expression(paste("Median predicted ", Delta^14,"C [‰]")), 
                     limits = c(-500,0), expand = c(0,0)) +
  scale_x_continuous("") +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

ai_c_p <- ai_p %>% 
  ggplot(aes(color = ClimateZoneAnd)) +
  geom_path(aes(x = pro_AI, y = pred_median_c), 
            linewidth = 1, alpha = 0.3, linetype = "dashed") +
  geom_path(aes(x = pro_AI_p, y = pred_median_c_p), 
            linewidth = 1) +
  geom_point(aes(x = pro_AI_p, y = pred_median_c_p, size = pro_AI_p)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        plot.margin = margin(r = 5, l = 25, b = 10)) +
  scale_y_continuous("Median predicted SOC [wt-%]", 
                     limits = c(0,9), expand = c(0,0)) +
  scale_x_continuous("PET/MAP") +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

ai_14c_c_p <- ai_p %>% 
  ggplot() +
  geom_point(aes(y = pred_median_14c_p, x = pred_median_c_p, 
                 size = pro_AI_p, color = ClimateZoneAnd)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(r = 15, t = 10, b = 10, l = 5)) +
  scale_y_continuous(expression(paste("Median predicted ", Delta^14,"C [‰]")), 
                     limits = c(-500,0), expand = c(0,0)) +
  scale_x_continuous("Median predicted SOC [wt-%]", expand = c(0,0), 
                     limits = c(0.3,9), breaks = c(0.5,1,5), labels = c(0.5,1,5)) +
  coord_trans(x = "log10") +
  annotation_logticks(sides = "b", scaled = FALSE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar) + 
  scale_size_continuous("PET/MAP") 

ggarrange(ggarrange(ai_14c_p, ai_c_p, nrow = 2, legend = "none"),
          ai_14c_c_p, common.legend = TRUE, widths = c(0.8,1))
ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14c_AndiPolar_ice_climate_ai_p_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

## MAT
# select 95% of the data for each climate group
mat_p_sum <- lyr_data %>% 
  distinct(id, .keep_all = TRUE) %>% 
  dplyr::filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::filter(ClimateZoneAnd != "tundra/polar") %>%
  group_by(ClimateZoneAnd) %>% 
  dplyr::select(ClimateZoneAnd, pro_MAT_mod) %>% 
  dplyr::summarise(P5 = quantile(pro_MAT_mod, 0.025),
                   P95 = quantile(pro_MAT_mod, 0.975))

mat_p <- mat_ice %>%
  left_join(mat_p_sum) %>% 
  dplyr::mutate(
    pro_MAT_mod_p = case_when(
      pro_MAT_mod < P5 ~ NA,
      pro_MAT_mod > P95 ~ NA,
      TRUE ~ pro_MAT_mod),
    pred_median_14c_p = case_when(
      pro_MAT_mod < P5 ~ NA,
      pro_MAT_mod > P95 ~ NA,
      TRUE ~ pred_median_14c),
    pred_median_c_p = case_when(
      pro_MAT_mod < P5 ~ NA,
      pro_MAT_mod > P95 ~ NA,
      TRUE ~ pred_median_c))

mat_14c_p <- mat_p %>%
  ggplot(aes(color = ClimateZoneAnd, x = pro_MAT_mod, y = pred_median_14c)) +
  geom_path(linewidth = 1, alpha = 0.3, linetype = "dashed") +
  geom_path(aes(x = pro_MAT_mod_p, y = pred_median_14c_p), 
            linewidth = 1) +
  geom_point(aes(x = pro_MAT_mod_p, y = pred_median_14c_p, size = pro_MAT_mod_p)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        plot.margin = margin(r = 5, l = 15, t = 10)) +
  scale_y_continuous(expression(paste("Median predicted ", Delta^14,"C [‰]")), 
                     limits = c(-500,0), expand = c(0,0)) +
  scale_x_continuous("") +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

mat_c_p <- mat_p %>% 
  ggplot(aes(color = ClimateZoneAnd)) +
  geom_path(aes(x = pro_MAT_mod, y = pred_median_c), 
            linewidth = 1, alpha = 0.3, linetype = "dashed") +
  geom_path(aes(x = pro_MAT_mod_p, y = pred_median_c_p), 
            linewidth = 1) +
  geom_point(aes(x = pro_MAT_mod_p, y = pred_median_c_p, size = pro_MAT_mod_p)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        plot.margin = margin(r = 5, l = 25, b = 10)) +
  scale_y_continuous("Median predicted SOC [wt-%]", 
                     limits = c(0,9), expand = c(0,0)) +
  scale_x_continuous("MAT [°C]") +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar)

mat_14c_c_p <- mat_p %>% 
  ggplot() +
  geom_point(aes(y = pred_median_14c_p, x = pred_median_c_p, 
                 size = pro_MAT_mod_p, color = ClimateZoneAnd)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(r = 15, t = 10, b = 10, l = 5)) +
  scale_y_continuous(expression(paste("Median predicted ", Delta^14,"C [‰]")), 
                     limits = c(-300,0), expand = c(0,0)) +
  scale_x_continuous("Median predicted SOC [wt-%]", expand = c(0,0), 
                     limits = c(0.3,2), breaks = c(0.5,1,2), labels = c(0.5,1,2)) +
  coord_trans(x = "log10") +
  annotation_logticks(sides = "b", scaled = FALSE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar) + 
  scale_size_continuous("MAT [°C]") 

ggarrange(ggarrange(mat_14c_p, mat_c_p, nrow = 2, legend = "none"),
          mat_14c_c_p, common.legend = TRUE, widths = c(0.8,1))
ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14c_AndiPolar_ice_climate_mat_p_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)
