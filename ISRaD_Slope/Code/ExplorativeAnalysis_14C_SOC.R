# Explore 14C profiles in ISRaD #
# Relationship between 14C and SOC #
# Sophie von Fromm #
# 15/08/2022 #

## Depth corrected values ##

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(mpspline2)

# Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-09-13"))

# Filter data for mspline function
lyr_mpspline <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  group_by(id) %>%
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
  ))

### Apply mspline function

## mspline 14C
lyr_data_mpspline_14c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 1)

## mspline CORG
lyr_data_mpspline_c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 60, lam = 1)

## 14C and SOC
mspline_14c_c <- lyr_data_mpspline_14c$est_1cm %>% 
  rename(lyr_14c = SPLINED_VALUE) %>% 
  full_join(lyr_data_mpspline_c$est_1cm %>% 
              rename(CORG = SPLINED_VALUE)) %>% 
  filter(LD < 101) %>% 
  tibble()

### Controls on relationship between SOC and 14C

## Climate Zones
climate_soil <- mspline_14c_c %>%
  tibble() %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, ClimateZone), 
                   by = "id") %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup()

# climate_soil %>% 
#   dplyr::select(entry_name, site_name, ClimateZone, median_14c, n, n_site, UD) %>% 
#   group_by(ClimateZone, UD) %>% 
#   mutate(n_rel = n/max(n)*100)

p_climate <- climate_soil %>%
  filter(n_site > 4) %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.3) +
  theme_classic(base_size = 16) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        legend.position = c(0.1,0.8))

p_climate +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path()
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_climate_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

# 14C and depth by climate zone
climate_soil %>%
  filter(n_site > 4) %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_14c, y = UD, color = ClimateZone)) +
  geom_path() +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClimateZone),
              alpha = 0.3) +
  geom_path(aes(x = n), linetype = "dashed") +
  theme_classic(base_size = 16) +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), position = "top",
                     expand = c(0,0), limits = c(-550,250)) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        legend.position = c(0.2,0.8))
ggsave(file = paste0("./Figure/ISRaD_14C_mspline_depth_climate_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

## USDA soil type
usda_soil <- mspline_14c_c %>%
  tibble() %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order), 
                   by = "id") %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup()  

p_usda <- usda_soil %>%
  filter(n_site > 4) %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c), color = "#fee0d2") +
  facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank())

p_usda +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c), color = "#fee0d2") +
  geom_path()
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_sum_soilt_1m_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

usda_soil %>% 
  filter(n_site > 4) %>%
  group_by(pro_usda_soil_order) %>% 
  count(UD) %>% 
  dplyr::filter(UD == 1| UD == 60| UD == 90 | UD == 97| UD == 99) %>% 
  view()

## Climate and clay type
cc_soil <- mspline_14c_c %>%
  tibble() %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order,
                                   ClimateZone), 
                   by = "id") %>% 
  mutate(clay_type = case_when(
    pro_usda_soil_order == "Oxisols" ~ "low-activity clays",
    pro_usda_soil_order == "Ultisols" ~ "low-activity clays",
    TRUE ~ "high-activity clays"
  )) %>% 
  group_by(ClimateZone, clay_type, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup()  

p_cc <- cc_soil %>%
  filter(n_site > 4) %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_type)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c), alpha = 0.3) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank())

p_cc +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c), alpha = 0.3) +
  geom_path()
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_clay_climate_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)


## WRB soil type
wrb_soil <- mspline_14c_c %>%
  tibble() %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_wrb_soil_order), 
                   by = "id") %>% 
  group_by(pro_wrb_soil_order, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name))  %>% 
  ungroup()

p_wrb <- wrb_soil %>% 
  # filter(n_site > 4) %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c), color = "#fee0d2") +
  facet_wrap(~pro_wrb_soil_order) +
  theme_classic(base_size = 16) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank())

p_wrb +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c), color = "#fee0d2") +
  geom_path()
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_wrb_all_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

wrb_soil %>% 
  filter(n_site > 4) %>% 
  group_by(pro_wrb_soil_order, UD) %>% 
  summarise(n = n_distinct(id)) %>% 
  filter(UD == 1| UD == 80| UD == 97) %>% view()

