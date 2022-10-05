# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 23/09/2022 #

library(tidyverse)
library(ggpubr)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-10-05"))

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
    entry_name == "Gentsch_2018" ~ "polar",
    pro_usda_soil_order == "Gelisols" ~ "polar",
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
  # group_by(ClimateZone) %>% 
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
  group_by(ClimateZone) %>%
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

## All splined profiles
fig_int <- plotly::ggplotly(
mspline_14c_c_all %>% 
  group_by(ClimateZone) %>% 
  arrange(UD, .by_group = TRUE) %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  # filter(n > 4 & n_rel > 60) %>% 
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

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_raw_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

htmlwidgets::saveWidget(
  widget = fig_int, #the plotly object
  file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_raw_", Sys.Date(),
                ".html"), #the path & file name
  selfcontained = TRUE #creates a single html file
)

max_depth <- mspline_14c_c_all %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  summarise(max_UD = max(UD))

fig_int <- plotly::ggplotly(
mspline_14c_c_all %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  # filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD, x = lyr_14c_msp)) +
  geom_path(aes(group = id, color = ClimateZone)) +
  geom_hline(data = max_depth, aes(yintercept = max_UD), linetype = "dashed") +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none") +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), position = "top") +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0))
)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_climate_raw_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

htmlwidgets::saveWidget(
  widget = fig_int, #the plotly object
  file = paste0("./Figure/ISRaD_msp_14C_depth_climate_raw_", Sys.Date(),
                ".html"), #the path & file name
  selfcontained = TRUE #creates a single html file
)

## Climate Zones
mx_sum <- mspline_14c_c_all %>% 
  group_by(ClimateZone, UD) %>%
  mutate(n = n()) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>%
  filter(n > 4 & n_rel > 60) %>%
  summarise(n_max = max(UD))

mspline_14c_c_all %>% 
  ggplot(aes(y = UD, x = lyr_14c_msp, color = ClimateZone)) +
  geom_path(aes(group = id)) +
  geom_hline(data = mx_sum, aes(yintercept = n_max), linetype = "dashed") +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none") +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     position = "top") +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0),
                     limits = c(100,0))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_climate_raw_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

mspline_14c_c_all %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, color = ClimateZone)) +
  geom_path(aes(group = id)) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [%]", trans = "log10")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_raw_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

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

cd1 <- climate_all %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_14c, color = ClimateZone), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClimateZone, 
                  color = ClimateZone), alpha = 0.5) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,110),
                     position = "top", breaks = seq(-1000,100,200)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))

cd2 <- climate_all %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = n, color = ClimateZone), orientation = "y") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = "none") +
  scale_x_continuous("Number of profiles", expand = c(0,0), limits = c(0,210),
                     position = "top", breaks = seq(0,200,50)) +
  scale_y_reverse("", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))

ggarrange(cd1, cd2, widths = c(0.7,0.3))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_climate_avg_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

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
        legend.position = c(0.1,0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")

p_climate_all +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path(size = 1)

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_avg_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

mspline_14c_c_all %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

### Soil type
fig_int <- plotly::ggplotly(
  mspline_14c_c_all %>% 
    group_by(pro_usda_soil_order) %>% 
    arrange(UD, .by_group = TRUE) %>% 
    group_by(pro_usda_soil_order, UD) %>% 
    mutate(n = n()) %>%
    ungroup(UD) %>% 
    mutate(n_rel = n * 100 / max(n)) %>% 
    # filter(n > 4 & n_rel > 60) %>% 
    ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) +
    geom_path(aes(group = id, color = pro_usda_soil_order)) +
    facet_wrap(~pro_usda_soil_order) +
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

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_soiltype_raw_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

htmlwidgets::saveWidget(
  widget = fig_int, #the plotly object
  file = paste0("./Figure/ISRaD_msp_14C_SOC_soiltype_raw_", Sys.Date(),
                ".html"), #the path & file name
  selfcontained = TRUE #creates a single html file
)

max_depth <- mspline_14c_c_all %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  summarise(max_UD = max(UD))

fig_int <- plotly::ggplotly(
  mspline_14c_c_all %>% 
    group_by(pro_usda_soil_order, UD) %>% 
    mutate(n = n()) %>%
    ungroup(UD) %>% 
    mutate(n_rel = n * 100 / max(n)) %>% 
    # filter(n > 4 & n_rel > 60) %>% 
    ggplot(aes(y = UD, x = lyr_14c_msp)) +
    geom_path(aes(group = id, color = pro_usda_soil_order)) +
    geom_hline(data = max_depth, aes(yintercept = max_UD), linetype = "dashed") +
    facet_wrap(~pro_usda_soil_order) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                          size = 0.3),
          panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                          size = 0.2),
          panel.spacing.x = unit(2, "lines"),
          legend.position = "none") +
    scale_x_continuous(expression(paste(Delta^14, "C [‰]")), position = "top") +
    scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                       limits = c(100,0))
)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_soiltype_raw_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

htmlwidgets::saveWidget(
  widget = fig_int, #the plotly object
  file = paste0("./Figure/ISRaD_msp_14C_depth_climate_raw_", Sys.Date(),
                ".html"), #the path & file name
  selfcontained = TRUE #creates a single html file
)

## Soil types and SOC (USDA)
soil_all <- mspline_14c_c_all %>%
  group_by(pro_usda_soil_order, UD) %>% 
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


cd1 <- soil_all %>% 
  filter(pro_usda_soil_order == "Vertisols"| 
           n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_14c, color = pro_usda_soil_order), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = pro_usda_soil_order, 
                  color = pro_usda_soil_order), alpha = 0.5) +
  # facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.7),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,110),
                     position = "top", breaks = seq(-1000,100,200)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))

cd2 <- soil_all %>% 
  filter(pro_usda_soil_order == "Vertisols"| 
           n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = n, color = pro_usda_soil_order), orientation = "y") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = "none") +
  scale_x_continuous("Number of profiles", expand = c(0,0), limits = c(0,110),
                     position = "top", breaks = seq(0,100,50)) +
  scale_y_reverse("", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))

ggarrange(cd1, cd2, widths = c(0.7,0.3))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_soiltype_avg_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

p_soil_all <- soil_all %>%
  filter(pro_usda_soil_order == "Vertisols"| 
           n > 4 & n_rel > 60)  %>% 
  ggplot(aes(x = median_c, y = median_14c, color = pro_usda_soil_order)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = pro_usda_soil_order),
                alpha = 0.3) +
  facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = "none",
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")

p_soil_all +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = pro_usda_soil_order),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_soiltype_avg_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

soil_all %>%
  filter(pro_usda_soil_order == "Vertisols"| 
           n > 4 & n_rel > 60)  %>% 
  filter(pro_usda_soil_order == "Andisols"|
           pro_usda_soil_order == "Gelisols"|
           pro_usda_soil_order == "Vertisols") %>% 
  ggplot(aes(x = median_c, y = median_14c, color = pro_usda_soil_order)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = pro_usda_soil_order),
                alpha = 0.4) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = pro_usda_soil_order),
                 alpha = 0.4) +
  geom_path() +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_soiltype_AndGelVer_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

soil_all %>%
  filter(n > 4 & n_rel > 60)  %>% 
  filter(pro_usda_soil_order != "Andisols",
         pro_usda_soil_order != "Gelisols",
         pro_usda_soil_order != "Vertisols",
         pro_usda_soil_order != "Entisols",
         pro_usda_soil_order != "Spodosols") %>% 
  ggplot(aes(x = median_c, y = median_14c, color = pro_usda_soil_order)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = pro_usda_soil_order),
                alpha = 0.4) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = pro_usda_soil_order),
                 alpha = 0.4) +
  geom_path(size = 1) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_soiltype_AlfIncepMolOxiUlti_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

mspline_14c_c_all %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

## Soil type and climate
soil_climate <- mspline_14c_c_all %>%
  group_by(ClimateZone, pro_usda_soil_order, UD) %>% 
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


soil_climate %>%
  filter(n > 4 & n_rel > 60)  %>% 
  filter(ClimateZone == "cold"|
           ClimateZone == "temperate") %>% 
  filter(pro_usda_soil_order == "Entisols") %>%
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.4) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path() +
  facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_Ent_TempCold_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

mspline_14c_c_all %>%
  filter(ClimateZone == "cold"|
           ClimateZone == "temperate") %>% 
  filter(pro_usda_soil_order == "Entisols") %>%
  group_by(ClimateZone, pro_country) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

mspline_14c_c_all %>%
  filter(ClimateZone == "tropical"|
           ClimateZone == "temperate") %>% 
  filter(pro_usda_soil_order == "Oxisols"|
           pro_usda_soil_order == "Ultisols"|
           pro_usda_soil_order == "Inceptisols") %>%
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
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.4) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path() +
  # facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_TropTemp_IncepOxUlt_avg_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

soil_climate %>%
  filter(n > 4 & n_rel > 60)  %>% 
  filter(ClimateZone == "tropical"|
           ClimateZone == "temperate") %>% 
  filter(pro_usda_soil_order == "Oxisols"|
           pro_usda_soil_order == "Ultisols"|
           pro_usda_soil_order == "Inceptisols") %>%
  ggplot(aes(x = UD, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.4) +
  # geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
  #                alpha = 0.4) +
  geom_path() +
  facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_TropTemp_IncepOxUlt_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)



mspline_14c_c_all %>%
  filter(ClimateZone == "tropical"|
           ClimateZone == "temperate") %>% 
  filter(pro_usda_soil_order == "Oxisols"|
           pro_usda_soil_order == "Ultisols"|
           pro_usda_soil_order == "Inceptisols") %>%
  group_by(ClimateZone, pro_usda_soil_order) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_countries = n_distinct(pro_country),
            n_profiles = n_distinct(id))

## Soil types and SOC (WRB)
soil_all_wrb <- mspline_14c_c_all %>%
  group_by(pro_wrb_soil_order, UD) %>% 
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

cd1 <- soil_all_wrb %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_14c, color = pro_wrb_soil_order), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = pro_usda_soil_order, 
                  color = pro_usda_soil_order), alpha = 0.5) +
  # facet_wrap(~pro_wrb_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.7),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,110),
                     position = "top", breaks = seq(-1000,100,200)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))

cd2 <- soil_all_wrb %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = n, color = pro_wrb_soil_order), orientation = "y") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = "none") +
  scale_x_continuous("Number of profiles", expand = c(0,0), limits = c(0,110),
                     position = "top", breaks = seq(0,100,50)) +
  scale_y_reverse("", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))

ggarrange(cd1, cd2, widths = c(0.7,0.3))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_soiltype_wrb_avg_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

p_soil_all_wrb <- soil_all_wrb %>%
  filter(n > 4 & n_rel > 60)  %>% 
  ggplot(aes(x = median_c, y = median_14c, color = pro_wrb_soil_order)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = pro_wrb_soil_order),
                alpha = 0.3) +
  facet_wrap(~pro_wrb_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = "none",
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10")

p_soil_all_wrb +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = pro_wrb_soil_order),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_soiltype_wrb_avg_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

soil_all_wrb %>%
  filter(n > 4 & n_rel > 60)  %>% 
  filter(pro_wrb_soil_order == "Andosols"|
           pro_wrb_soil_order == "Cryosols"|
           pro_wrb_soil_order == "Vertisols") %>% 
  ggplot(aes(x = median_c, y = median_14c, color = pro_wrb_soil_order)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = pro_wrb_soil_order),
                alpha = 0.3) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = pro_wrb_soil_order),
                 alpha = 0.4) +
  geom_path()


mspline_14c_c_all %>% 
  group_by(pro_wrb_soil_order, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id)) 

### Holdridge Life Zones
## Add HLZ data
library(sf)
library(tmap)
HLZ_directory <- "D:/Sophie/PhD/AfSIS_GlobalData/HLZ/holdrid/holdrid.shp"

HLZ <- sf::st_read(HLZ_directory)
st_crs(HLZ) <- 4326

mspline_all_sf <- sf::st_as_sf(mspline_14c_c_all, 
                               coords = c("pro_long", "pro_lat"), 
                               crs = 4326)

tmap_mode("view")
tm_shape(HLZ, projection = 4326) +
  tm_polygons(col = "DESC") +
  tm_shape(lyr_data_sf) +
  tm_dots(popup.vars = "id")

lyr_data_HLZ_sf <- sf::st_join(mspline_all_sf, HLZ, 
                               left = TRUE)

lyr_data_HLZ_NA <- lyr_data_HLZ_sf %>% 
  tibble() %>% 
  dplyr::select(-c(geometry, AREA:FREQUENCY, SYMBOL)) %>% 
  rename(HLZ_Zone = DESC)

names(lyr_data_HLZ_NA)

lyr_data_HLZ_NA %>% 
  filter(is.na(HLZ_Zone)) %>% 
  count(id) %>% view()

#Manually assign missing HLZ
mspline_HLZ <- lyr_data_HLZ_NA %>% 
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
  mutate(HLZ_Zone = ifelse(entry_name == "Quero_2022",
                           "Warm temperate dry forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(grepl("Sanaiotti_2002_Alter do Chão_Alter do Chão", id),
                           "Subtropical moist forest", HLZ_Zone)) %>%
  #Could also be "Subtropical moist forest" or "Subtropical dry forest"
  mutate(HLZ_Zone = ifelse(entry_name == "Schwartz_1992", 
                           "Tropical dry forest", HLZ_Zone)) %>%
  #Most of Hawaii has missing data
  mutate(HLZ_Zone = ifelse(grepl("Grant_2022_Kohala", id), 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_1997_Amalu-precipitation_Amalu-precipitation_profile_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_1997_Kohala (150ky)_Kohala (150ky)_profile_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_1997_Kokee (4.1my)_Kokee (4.1my)_profile_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_1997_Kolekole (1.4my)_Kolekole (1.4my)_profile_1", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_2005_Kohala_Kohala_150", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_2005_Kokee, Kauai_Kokee_Kauai_4100", 
                           "Subtropical moist forest", HLZ_Zone)) %>% 
  mutate(HLZ_Zone = ifelse(id == "Torn_2005_Kolekole, Molokai_Kolekole_Molokai_1400", 
                           "Subtropical moist forest", HLZ_Zone)) 

mspline_HLZ %>% 
  filter(is.na(HLZ_Zone)) %>% 
  count(id)

mspline_HLZ %>% 
  group_by(HLZ_Zone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  # filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD, x = lyr_14c_msp)) +
  geom_path(aes(group = id, color = HLZ_Zone)) +
  facet_wrap(~HLZ_Zone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none") +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), position = "top") +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0), 
                     limits = c(100,0))

mspline_HLZ %>% 
  group_by(HLZ_Zone) %>% 
  arrange(UD, .by_group = TRUE) %>% 
  group_by(HLZ_Zone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) +
  geom_path(aes(group = id, color = entry_name)) +
  facet_wrap(~HLZ_Zone) +
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



### Look into mechanisms: filter study with clay content, alox, feox, ...
## Linear mixed model


## Clay content
data_clay <- lyr_mpspline %>% 
  drop_na(lyr_clay_tot_psa)

data_clay %>% 
  group_by(ClimateZone) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

data_clay %>%  
  group_by(id) %>% 
  arrange(lyr_top) %>% 
  ungroup() %>% 
  ggplot(aes(x = CORG, y = lyr_14c, fill = lyr_clay_tot_psa)) +
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
  scale_fill_viridis_c()

## Alox
data_alox <- lyr_mpspline %>% 
  drop_na(lyr_al_ox)

data_alox %>% 
  group_by(ClimateZone) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

data_alox %>% 
  distinct(entry_name) %>% 
  view()

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
