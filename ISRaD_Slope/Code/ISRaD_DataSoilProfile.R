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
        legend.position = c(0.23,0.15),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 5, l = 5),
        axis.title.x = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]; log-scaled", expand = c(0,0), 
                     breaks = c(0.1,0.5,1,5,10),
                     labels = c(0.1,0.5,1,5,10)) +
  coord_trans(x = "log10", xlim = c(0.1,16)) +
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
        plot.margin = margin(r = 5, b = 5, t = 5, l = 10)) +
  scale_y_reverse("Depth [cm]",expand = c(0,0), limits = c(100,0)) +
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
        plot.margin = margin(r = 10, b = 5, t = 5, l = 10),
        axis.title.y = element_blank()) +
  scale_y_reverse("Depth [cm]", expand = c(0,0), limits = c(100,0)) +
  scale_x_continuous("SOC [wt-%]", expand = c(0,0), 
                     breaks = c(0.1,0.5,1,5,10),
                     labels = c(0.1,0.5,1,5,10)) +
  coord_trans(x = "log10", xlim = c(0.1,16)) + 
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_color_manual("Climate grouping", values = color_climate)

ggarrange(annotate_figure(ggarrange(c1_14c, c1_c, nrow = 2, legend = "none",
                                    labels = c("a)", "b)"), vjust = c(2.2,1.8),
                                    hjust = c(-3,-2.8)),
                          left = text_grob("Depth [cm]", rot = 90, size = 16)),
          c1_14c_c, common.legend = TRUE, widths = c(0.8,1), labels = c("", "c)"),
          vjust = c(1.5,2.2), hjust = c(-0.5,-5.8))
# ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14C_climateAnd_depth_", Sys.Date(),
#                      ".jpeg"), width = 10, height = 6)

## Mineral data
mineral_14c_c <- lyr_data %>%
  group_by(MineralGroupsNew, UD) %>% 
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
  group_by(MineralGroupsNew) %>% 
  dplyr::summarise(median_14c_max = max(median_14c),
                   median_14c_min = min(median_14c),
                   median_c_max = max(median_c),
                   median_c_min = min(median_c)) %>% 
  rowwise() %>% 
  mutate(diff_14c = median_14c_max - median_14c_min,
         diff_c = median_c_max - median_c_min)

depth_sum_m <- mineral_14c_c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  dplyr::select(MineralGroupsNew, UD, median_c, median_14c) %>% 
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
  group_by(MineralGroupsNew) %>% 
  filter(UD == max(UD))


m1_14c_c <- mineral_14c_c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = MineralGroupsNew), linewidth = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralGroupsNew), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralGroupsNew), alpha = 0.3) +
  geom_point(data = depth_sum_m, aes(x = median_c, y = median_14c,
                                     shape = as.factor(UD)), 
             size = 2, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.3,0.14),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 5, r = 5),
        axis.title.x = element_blank()) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_x_continuous("Soil organic carbon [wt-%]; log-scaled", expand = c(0,0), 
                     breaks = c(0.1,0.5,1,5,10),
                     labels = c(0.1,0.5,1,5,10)) +
  coord_trans(x = "log10", xlim = c(0.1,16)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_manual("Dominating mineral type", values = color_mineral) +
  guides(shape = "none")

m1_14c <- mineral_14c_c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(y = UD, x = median_14c, color = MineralGroupsNew), linewidth = 2) +
  geom_errorbarh(aes(xmin = lci_14c, xmax = uci_14c, y = UD, 
                     color = MineralGroupsNew), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.8),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 5, b = 5, t = 5, l = 10)) +
  scale_y_reverse("Depth [cm]", expand = c(0,0), limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_manual("Dominating mineral type", values = color_mineral) 

m1_c <- mineral_14c_c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(y = UD, x = median_c, color = MineralGroupsNew), linewidth = 2) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = UD, 
                     color = MineralGroupsNew), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8,0.5),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 10, b = 5, t = 5, l = 10),
        axis.title.y = element_blank()) +
  scale_y_reverse(expand = c(0,0), limits = c(100,0)) +
  scale_x_continuous("SOC [wt-%]", expand = c(0,0), 
                     breaks = c(0.1,0.5,1,5,10),
                     labels = c(0.1,0.5,1,5,10)) +
  coord_trans(x = "log10", xlim = c(0.1,16)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_color_manual("Dominating mineral type", values = color_mineral) 

## Figure 3
clim_min <- annotate_figure(ggarrange(c1_14c_c, m1_14c_c, labels = c("a)", "b)")),
                            bottom = text_grob("Soil organic carbon [wt-%]; log-scaled",
                                               size = 16))
ggsave(file = paste0("./Figure/Figure_3_a_b_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

climate_plot <- ggarrange(c1_14c, c1_c, ncol = 2, legend = "none",
                          labels = c("c)", "d)"))

mineral_plot <- ggarrange(m1_14c, m1_c, ncol = 2, legend = "none",
                          labels = c("e)", "f)"))

ggarrange(climate_plot, mineral_plot, nrow = 2, legend = "none")
ggsave(file = paste0("./Figure/Figure_3_c_f_", Sys.Date(),
                     ".jpeg"), width = 6, height = 6)

