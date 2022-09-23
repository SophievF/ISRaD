# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 23/09/2022 #

library(tidyverse)
library(ggpubr)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-09-22"))

lyr_all %>% 
  count(entry_name)

names(lyr_all)

lyr_mpspline <- lyr_all %>% 
  #reclassify soil type Schuur_2001: all Andisols
  # mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
  #                                      entry_name == "Schuur_2001" & 
  #                                        pro_usda_soil_order == "Inceptisols",
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
  filter(lyr_obs_date_y > 1959) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(ClimateZone = case_when(
    pro_usda_soil_order == "Andisols" ~ "andisols",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold",
    str_detect(pro_KG_present_long, "Polar") ~ "polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  #remove for now: need to fix depth
  filter(entry_name != "Fernandez_1993a")

lyr_mpspline$ClimateZone <- factor(lyr_mpspline$ClimateZone,
                                   levels = c("andisols", "polar", "cold",
                                              "temperate", "arid", "tropical"))

summary(lyr_mpspline$CORG)
summary(lyr_mpspline$lyr_14c)
summary(lyr_mpspline$ClimateZone)

lyr_mpspline %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

### Apply mspline function

## mspline 14C
lyr_data_mpspline_14c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 0.5)

## mspline CORG
lyr_data_mpspline_c <- lyr_mpspline %>% 
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
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id, .keep_all = TRUE), 
                   by = "id") %>% 
  group_by(id) %>% 
  arrange(UD) %>% 
  ungroup()

mspline_14c_c_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

## All splined profiles
plotly::ggplotly(
mspline_14c_c_all %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) +
  geom_path(aes(group = id, color = ClimateZone)) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  scale_x_continuous("Soil organic carbon [%]", trans = "log10")
)

## Climate Zones and SOC
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

p_climate_all <- climate_all %>%
  filter(n > 4 & n_rel > 60)  %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.3) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.3),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")

p_climate_all +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path()

mspline_14c_c_all %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

## Look into mechanisms: filter study with clay content, alox, feox, ...
data_alox <- lyr_mpspline %>% 
  drop_na(lyr_al_ox)

data_alox %>% 
  group_by(ClimateZone) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

data_alox %>% 
  distinct(entry_name) 

data_alox %>% 
  count(entry_name, ClimateZone) 

data_alox %>% 
  group_by(id) %>% 
  arrange(lyr_top) %>% 
  ungroup() %>% 
  ggplot(aes(x = CORG, y = lyr_14c, fill = lyr_al_ox)) +
  geom_path(aes(group = id)) +
  geom_point(size = 3, shape = 21) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_fill_viridis_c(trans = "log10")
