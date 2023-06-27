# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 01/06/2022 #

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mpspline2)

source("./Code/function_mpspline_mod.R")

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-10-21"))

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
    entry_name == "Gentsch_2018" ~ "tundra/polar",
    pro_usda_soil_order == "Gelisols" ~ "tundra/polar",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "warm temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold temperate",
    str_detect(pro_KG_present_long, "Polar") ~ "tundra/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  mutate(ClimateZoneAnd = case_when(
    entry_name == "Gentsch_2018" ~ "tundra/polar",
    pro_usda_soil_order == "Gelisols" ~ "tundra/polar",
    pro_usda_soil_order == "Andisols" ~ "volcanic soils",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "warm temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold temperate",
    str_detect(pro_KG_present_long, "Polar") ~ "tundra/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  #remove for now: need to fix depth
  filter(entry_name != "Fernandez_1993a")

lyr_data %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>%
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

### Apply mspline function

## mspline 14C
lyr_data_mpspline_14c <- lyr_data %>% 
  # filter(pro_name == "HAK1468_1.3") %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 0.5,
                   d = c(0, 5, 15, 30, 60, 100, 200))

## mspline CORG
lyr_data_mpspline_c <- lyr_data %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 60, lam = 0.5,
                d = c(0, 5, 15, 30, 60, 100, 200))

## 14C and SOC
mspline_14c_c <- lyr_data_mpspline_14c$est_1cm %>% 
  rename(lyr_14c_msp = SPLINED_VALUE) %>% 
  full_join(lyr_data_mpspline_c$est_1cm %>% 
              rename(CORG_msp = SPLINED_VALUE)) %>% 
  # filter(LD < 102) %>% 
  tibble()

mspline_14c_c_all <- mspline_14c_c %>%
  dplyr::left_join(lyr_data %>% 
                     distinct(id, .keep_all = TRUE), 
                   by = "id") %>% 
  group_by(id) %>% 
  arrange(UD) %>% 
  ungroup() 

### Plotting

## 14C and depth

# spline data
p1 <- mspline_14c_c_all %>% 
  # filter(LD < 101) %>% 
  group_by(UD) %>% 
  summarize(median = median(lyr_14c_msp, na.rm = TRUE),
            mad = mad(lyr_14c_msp, na.rm = TRUE),
            n = n()) %>% 
  ggplot(aes(x = median, y = UD)) +
  geom_ribbon(aes(xmin = median - mad, xmax = median + mad), fill = "darkgrey") +
  geom_path(color= "red", size = 2) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Delta 14C", expand = c(0,0), limits = c(-1000,200),
                     position = "top") +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0),
                     limits = c(701,0))

# raw data
p2 <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>%
  dplyr::select(id, depth, lyr_14c) %>% 
  mutate(depth_bin = cut(depth, breaks = seq(0,110,10))) %>% 
  extract(depth_bin, c("top", "bot"), "(-*[0-9]+\\.*[0-9]*),(-*[0-9]+\\.*[0-9]*)", 
          convert = TRUE) %>% 
  mutate(depth = (bot-top)/2+top) %>%
  group_by(top) %>% 
  summarize(median = median(lyr_14c),
            mad = mad(lyr_14c)) %>%
  drop_na() %>% 
  ggplot(aes(x = median, y = top)) +
  geom_ribbon(aes(xmin = median - mad, xmax = median + mad), fill = "darkgrey") +
  geom_path(color= "red", size = 2) +
  geom_point() +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Delta 14C", expand = c(0,0), limits = c(-600,200),
                     position = "top") +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0))

ggarrange(p1, p2)
ggsave(file = paste0("./Figure/ISRaD_14C_depth_mspline_raw_avg_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)


## 14C and SOC

# spline data
p3 <- mspline_14c_c_all %>% 
  group_by(UD) %>% 
  summarize(median_14c = median(lyr_14c_msp),
            mad_14c = mad(lyr_14c_msp),
            median_CORG = median(CORG_msp),
            mad_CORG = mad(CORG_msp),
            n = n(),
            n_site = n_distinct(site_name),
            lower_CORG = median_CORG - mad_CORG,
            upper_CORG = median_CORG + mad_CORG) %>% 
  mutate(lower_CORG = if_else(lower_CORG <= 0.01, 0.01, lower_CORG)) %>% 
  ungroup() %>% 
  mutate(n_rel = n/max(n)) %>% 
  arrange(desc(median_14c)) %>% 
  filter(n_site > 4) %>%
  filter(n_rel > 1/9) %>%
  ggplot(aes(y = median_14c, x = median_CORG)) +
  geom_errorbarh(aes(xmin = lower_CORG, xmax = upper_CORG), color = "grey") +
  geom_errorbar(aes(ymin = median_14c - mad_14c, ymax = median_14c + mad_14c), color = "grey") +
  geom_path(color= "red", size = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Delta 14C", expand = c(0,0), limits = c(-1000,200)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", expand = c(0,0), limits = c(0.01, 7))

# raw data
p4 <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>%
  mutate(depth_bin = cut_number(lyr_14c, 20)) %>% 
  extract(depth_bin, c("bot_14c", "top_14c"), "(-*[0-9]+\\.*[0-9]*),(-*[0-9]+\\.*[0-9]*)", 
          convert = TRUE) %>% 
  mutate(mean_14c = (bot_14c-top_14c)/2+top_14c) %>%
  group_by(mean_14c, top_14c, bot_14c) %>% 
  summarize(median_CORG = median(CORG),
            mad_CORG = mad(CORG), 
            n = n(),
            n_site = n_distinct(site_name),
            lower_CORG = median_CORG - mad_CORG,
            upper_CORG = median_CORG + mad_CORG) %>%
  mutate(lower_CORG = if_else(lower_CORG <= 0.01, 0.01, lower_CORG)) %>% 
  ungroup() %>% 
  mutate(n_rel = n/max(n)) %>% 
  arrange(desc(mean_14c)) %>% 
  filter(n_site > 4) %>% 
  ggplot(aes(x = median_CORG, y = mean_14c)) +
  geom_errorbarh(aes(xmin = lower_CORG, xmax = upper_CORG), color = "grey",
                 height = 0) +
  geom_errorbar(aes(ymin = bot_14c, ymax = top_14c), color = "grey",
                width = 0) +
  geom_path(color = "red", size = 2) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Delta 14C", expand = c(0,0), limits = c(-1000,305)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", expand = c(0,0), limits = c(0.01, 7))

ggarrange(p3, p4)
ggsave(file = paste0("./Figure/ISRaD_14C_C_mspline_raw_avg_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>%
  mutate(depth_bin = cut_interval(lyr_14c, 10)) %>% 
  filter(depth_bin == "[-981,-852]") %>% count(entry_name, pro_country)
