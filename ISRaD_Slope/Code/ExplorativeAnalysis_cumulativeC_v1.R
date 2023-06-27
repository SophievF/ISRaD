# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 09/01/2023 #

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mpspline2)

### Cumulative analysis ###

# Load filtered database data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2023-02-08"))

lyr_all %>% 
  skimr::skim_without_charts(lyr_14c, CORG)

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
color_climate_wo_polar <- c("#9ecae1", "#7fbf7b", "#1b7837", "#dfc27d")
# polar, cold, warm, tropical, arid
color_climate <- c("#3182bd", "#9ecae1", "#7fbf7b", "#1b7837", "#dfc27d")
# amorphous, polar, cold, warm, tropical, arid
color_climate_w_amorph <- c("#c083be", "#3182bd", "#9ecae1", "#7fbf7b", "#1b7837", 
                            "#dfc27d")
# amorphous, primary mineral, high-activity clay, low-activity clay
color_mineral <- c("#c083be", "#bf812d", "#7fcdbb", "#1b7837")


## Cumulative C - all profiles
cum_all <- lyr_data %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  group_by(UD) %>% 
  summarise(median_c = median(CORG_msp),
            median_14c = median(lyr_14c_msp)) %>% 
  mutate(rel_c = median_c*100/sum(median_c)) %>% 
  mutate(cum_c = cumsum(median_c*100/sum(median_c))) %>% 
  mutate(prc_rank = percent_rank(desc(median_c))*100) %>% 
  mutate(cum_c = case_when(
    UD == 1 ~ 0,
    TRUE ~ cum_c
  ))

lyr_data %>% 
  group_by(id) %>%
  mutate(cum_c_each = cumsum(CORG_msp*100/sum(CORG_msp))) %>% 
  ungroup() %>% 
  group_by(UD) %>% 
  summarise(cum_c = median(cum_c_each),
            median_14c = median(lyr_14c_msp)) %>% 
  mutate(cum_c = case_when(
    UD == 1 ~ 0,
    TRUE ~ cum_c
  )) 

depth_all_cum <- cum_all %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 89|
                  UD == 99)

cum_all %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_14c, y = cum_c)) +
  geom_path(linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = depth_all_cum,
             aes(y = cum_c, x = median_14c), size = 3,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"), 
        plot.margin = margin(t = 10, b = 5, l = 5, r = 15)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125),
                     breaks = c(-1000,-750,-500,-250,0,125)) +
  scale_y_continuous("Cumulative SOC content [%]", expand = c(0,0),
                     trans = "reverse", limits = c(100,0))
ggsave(file = paste0("./Figure/ISRaD_msp_cum_C_all_", Sys.Date(),
                     ".jpeg"), width = 6, height = 4)

## Cumulative C - by climate
climate_median <- lyr_data %>%
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

climate_cum <- climate_median %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c, n_rel, n) %>% 
  arrange(UD) %>% 
  #calculate based on sum for each climate zone  
  filter(n > 4 & n_rel > 33) %>%
  mutate(cum_c = cumsum(median_c*100/sum(median_c))) %>% 
  mutate(prc_rank = percent_rank(desc(median_c))*100) %>% 
  mutate(cum_c = case_when(
    UD == 1 ~ 0,
    TRUE ~ cum_c
  ))

depth_climate_cum <- climate_cum %>% 
  filter(n > 4 & n_rel > 33) %>% 
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

climate_cum %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_14c, y = cum_c, color = ClimateZoneAnd)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_path(linewidth = 1.5) +
  geom_point(data = depth_climate_cum,
             aes(y = cum_c, x = median_14c), size = 3,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"), 
        plot.margin = margin(t = 10, b = 5, l = 5, r = 15),
        legend.position = c(0.15,0.8)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125),
                     breaks = c(-1000,-750,-500,-250,0,125)) +
  scale_y_continuous("Cumulative SOC content [%]", expand = c(0,0),
                     trans = "reverse", limits = c(100,0)) +
  scale_color_manual("Grouping", values = color_climate_w_amorph)
ggsave(file = paste0("./Figure/ISRaD_msp_cum_C_climate_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)

climate_cum %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ungroup() %>% 
  summarise(sum_c = sum(median_c)) 

climate_cum %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(sum_c = sum(median_c)) %>% 
  mutate(rel_c = sum_c*100/sum(sum_c))

climate_cum_weighted <- climate_cum %>% 
  filter(n > 4 & n_rel > 33) %>% 
  mutate(cum_c_weighted = cumsum(median_c*100/1497)) %>% 
  mutate_at("cum_c_weighted", ~replace(., which(UD == 1), 0))

depth_climate_cum_weighted <- climate_cum_weighted %>% 
  filter(n > 4 & n_rel > 33) %>% 
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


climate_cum_weighted %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_14c, y = cum_c_weighted, color = ClimateZoneAnd)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_path(linewidth = 1.5) +
  geom_point(data = depth_climate_cum_weighted, 
             aes(y = cum_c_weighted, x = median_14c), size = 3,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"), 
        plot.margin = margin(t = 10, b = 5, l = 5, r = 15),
        legend.position = c(0.8,0.2)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125),
                     breaks = c(-1000,-750,-500,-250,0,125)) +
  scale_y_continuous("Cumulative SOC content [%]", expand = c(0,0),
                     trans = "reverse", limits = c(60,0)) +
  scale_color_manual("Grouping", values = color_climate_w_amorph)
ggsave(file = paste0("./Figure/ISRaD_msp_cum_C_climate_weighted_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)  




  
climate_cum %>% 
  # filter(ClimateZoneAnd != "volcanic soils",
  #        ClimateZoneAnd != "tundra/polar") %>% 
  ungroup() %>% 
  summarise(sum_c = sum(median_c)) 

climate_cum %>% 
  filter(ClimateZoneAnd != "volcanic soils",
         ClimateZoneAnd != "tundra/polar") %>% 
  summarise(sum_c = sum(median_c)) %>% 
  mutate(rel_c = sum_c*100/sum(sum_c))

climate_cum %>% 
  # filter(ClimateZoneAnd != "volcanic soils",
  #        ClimateZoneAnd != "tundra/polar") %>% 
  mutate(cum_c_weighted = cumsum(median_c*100/1497)) %>% 
  mutate_at("cum_c_weighted", ~replace(., which(UD == 1), 0)) %>% view()
arrange(UD) %>% 
  ggplot(aes(x = median_14c, y = cum_c_weighted, color = ClimateZoneAnd)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_path(linewidth = 1.5) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"), 
        plot.margin = margin(t = 10, b = 5, l = 5, r = 15),
        legend.position = c(0.15,0.8)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-500,125),
                     breaks = c(-500,-250,0,125)) +
  scale_y_continuous("Cumulative SOC content [%]", expand = c(0,0),
                     trans = "reverse", limits = c(40,0)) +
  scale_color_manual("Grouping", values = color_climate_w_amorph)


