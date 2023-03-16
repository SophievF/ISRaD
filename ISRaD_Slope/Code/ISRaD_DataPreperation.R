# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 25/01/2023 #

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2023-02-08"))

lyr_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

lyr_data <- lyr_all %>% 
  # filter(CORG <= 20) %>% 
  # filter(lyr_obs_date_y > 1959) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  dplyr::arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(ClimateZone = case_when(
    entry_name == "Gentsch_2018" ~ "tundra/polar",
    pro_usda_soil_order == "Gelisols" ~ "tundra/polar",
    pro_usda_soil_order == "Aridisols" ~ "arid",
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
    pro_usda_soil_order == "Aridisols" ~ "arid",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "warm temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold temperate",
    str_detect(pro_KG_present_long, "Polar") ~ "tundra/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  mutate(MineralType = case_when(
    pro_usda_soil_order == "Andisols" ~ "amorphous",
    pro_usda_soil_order == "Oxisols" ~ "low-activity clay",
    pro_usda_soil_order == "Ultisols" ~ "low-activity clay",
    TRUE ~ "high-activity clay"
  )) %>% 
  mutate(pro_AI = pro_PET_mm_yr_mod/pro_MAP_mod) %>% 
  #remover overlapping depth layer in Fernandez_1993
  filter(lyr_name != "B_61.1-122.2") %>% 
  #remove Scharpenseel_1973a_Romania_Romania:44.67,22.33: only starts at 90 cm
  filter(id != "Scharpenseel_1973a_Romania_Romania:44.67,22.33") %>% 
  #remove Scharpensel_1973a site Riohlon St: buried soil (really old: -600 at 20 cm)
  filter(site_name != "Riohlon St")

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

## mspline 14C
lyr_data_mpspline_14c <- lyr_data %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, vhigh = 350, lam = 0.5)

## mspline CORG
lyr_data_mpspline_c <- lyr_data %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.005, vhigh = 60, lam = 0.5)

## mspline clay content
lyr_data_mpspline_clay <- lyr_data %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_clay_tot_psa) %>% 
  mpspline_tidy(vlow = 0, vhigh = 100, lam = 0.5)

## 14C and SOC and clay
mspline_14c_c_clay <- lyr_data_mpspline_14c$est_1cm %>% 
  dplyr::rename(lyr_14c_msp = SPLINED_VALUE) %>% 
  full_join(lyr_data_mpspline_c$est_1cm %>% 
              dplyr::rename(CORG_msp = SPLINED_VALUE)) %>% 
  left_join(lyr_data_mpspline_clay$est_1cm %>%
              dplyr::rename(clay_msp = SPLINED_VALUE)) %>%
  filter(UD < 101) %>% 
  tibble() %>% 
  #some profiles only have a few NA's - replace with value below/above
  group_by(id) %>%
  tidyr::fill(clay_msp, .direction = "downup") %>%
  ungroup()

summary(mspline_14c_c_clay)

mspline_14c_c_clay$clay_msp <- as.numeric(mspline_14c_c_clay$clay_msp)
lyr_data$pro_clay_10_cm_SG <- as.numeric(lyr_data$pro_clay_10_cm_SG)
lyr_data$pro_clay_30_cm_SG <- as.numeric(lyr_data$pro_clay_30_cm_SG)
lyr_data$pro_clay_60_cm_SG <- as.numeric(lyr_data$pro_clay_60_cm_SG)
lyr_data$pro_clay_100_cm_SG <- as.numeric(lyr_data$pro_clay_100_cm_SG)

mspline_14c_c_all <- mspline_14c_c_clay %>%
  dplyr::left_join(lyr_data %>% 
                     distinct(id, .keep_all = TRUE), 
                   by = "id") %>% 
  group_by(id) %>% 
  dplyr::arrange(UD) %>% 
  ungroup() %>% 
  #gap-fill missing clay data
  mutate(lyr_clay_mod = case_when(
    is.na(clay_msp) & 
      LD <= 10 ~ pro_clay_10_cm_SG,
    is.na(clay_msp) &
      LD > 10 & LD <= 30 ~ pro_clay_30_cm_SG,
    is.na(clay_msp) &
      LD > 30 & LD <= 60 ~ pro_clay_60_cm_SG,
    is.na(clay_msp) &
      LD > 60 & LD <= 101 ~ pro_clay_100_cm_SG,
    TRUE ~ clay_msp
  ))

summary(mspline_14c_c_all$lyr_clay_mod)
summary(mspline_14c_c_all$clay_msp)
summary(mspline_14c_c_all$pro_clay_10_cm_SG)

write_csv(mspline_14c_c_all, file = paste0("./Data/ISRaD_flat_splined_filled_", 
                                           Sys.Date(), ".csv"))
