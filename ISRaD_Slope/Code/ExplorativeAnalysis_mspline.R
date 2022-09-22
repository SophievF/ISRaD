# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 15/09/2022 #

library(tidyverse)
library(ggpubr)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-09-22"))

lyr_all %>% 
  count(entry_name)

names(lyr_all)

lyr_mpspline <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  # filter(lyr_bot <= 100) %>% 
  # filter(lyr_top <= 100) %>%
  group_by(id) %>%
  # filter(min(lyr_top) == 0) %>% 
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
  )) %>% 
  #remove for now: need to fix depth
  filter(entry_name != "Fernandez_1993a")

summary(lyr_mpspline$CORG)
summary(lyr_mpspline$lyr_14c)

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
                   by = "id") 

mspline_14c_c_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

mspline_14c_c_all %>% 
  ggplot(aes(y = lyr_14c_msp, x = UD)) +
  geom_line(aes(group = id), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous(expression(paste("fitted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,350)) +
  scale_x_continuous("Depth [cm]", limits = c(0,105),
                     expand = c(0,0)) +
  geom_smooth(method = "gam", formula = y ~ s(log(x)))

mspline_14c_c_all %>% 
  ggplot(aes(y = CORG_msp, x = UD)) +
  geom_line(aes(group = id), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("fitted soil organic carbon [wt-%]", expand = c(0,0), limits = c(0,55)) +
  scale_x_continuous("Depth [cm]", limits = c(0,105), expand = c(0,0)) +
  geom_smooth(method = "gam", formula = y ~ s(log(x)))

## "Raw" profiles
mspline_14c_c_all %>% 
  filter(entry_name == "Lassey_1996") %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) +
  geom_path(aes(group = id)) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [%]", trans = "log10")


mspline_14c_c_all %>% 
  mutate(ClimateZone = case_when(
    pro_usda_soil_order == "Andisols" ~ "andisols",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) +
  geom_path(aes(group = id)) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0)) +
  scale_x_continuous("Soil organic carbon [%]", trans = "log10")


## Check Andisols
mspline_14c_c_all %>% 
  filter(entry_name == "Giardina_2014") %>% 
  ggplot() +
  geom_line(aes(y = UD, x = lyr_14c_msp, group = id, color = plot_name), orientation = "y") +
  facet_wrap(~plot_name) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,300),
                     position = "top", breaks = seq(-1000,200,200)) 

mspline_14c_c_all %>% 
  filter(entry_name == "Grant_2022") %>% 
  ggplot() +
  geom_line(aes(y = UD, x = lyr_14c_msp, group = id, color = pro_land_cover), 
            orientation = "y") +
  facet_wrap(~pro_land_cover) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,300),
                     position = "top", breaks = seq(-1000,200,200)) 

mspline_14c_c_all %>% 
  filter(entry_name == "Giardina_2014"|
           entry_name == "Grant_2022") %>% 
  group_by(entry_name, UD) %>% 
  mutate(median_pseudo = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_median = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2]) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_pseudo, color = entry_name), orientation = "y") +
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median, fill = entry_name, 
                 color = entry_name), alpha = 0.5) +
  facet_wrap(~entry_name) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25)) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,300),
                     position = "top", breaks = seq(-1000,200,200)) 
  

## Climate and depth
climate_14c_depth <-  mspline_14c_c_all %>% 
  dplyr::filter(LD != 101) %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(median_pseudo = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_median = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_pseudo, .keep_all = TRUE) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  ungroup()

cd1 <- climate_14c_depth %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_pseudo, color = ClimateZone), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median, fill = ClimateZone, 
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
                     limits = c(-450,110),
                     position = "top", breaks = seq(-400,100,100)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25)) +
  scale_fill_viridis_d(direction = -1) +
  scale_color_viridis_d(direction = -1)

cd2 <- climate_14c_depth %>% 
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
                  breaks = seq(0,100,25)) +
  scale_color_viridis_d(direction = -1) 

ggarrange(cd1, cd2, widths = c(0.7,0.3))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_climate_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

# Climate and depth w/o andisols
climate_14c_depth_and <- mspline_14c_c_all %>% 
  mutate(ClimateZone = case_when(
    pro_usda_soil_order == "Andisols" ~ "andisols",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  dplyr::filter(LD != 101) %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(median_pseudo = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_median = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_pseudo, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  ungroup()

climate_14c_depth_and %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_pseudo, color = ClimateZone), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median, fill = ClimateZone, 
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
                     limits = c(-650,110),
                     position = "top", breaks = seq(-600,100,100)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25)) +
  scale_fill_manual(values = c("#fb6a4a", "#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF")) +
  scale_color_manual(values = c("#fb6a4a", "#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_climate_and_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

# Average andisols
# Calculate mean for each study: Giardina_2014 & Grant_2022
andi <- mspline_14c_c_all %>% 
  filter(entry_name == "Giardina_2014"|
           entry_name == "Grant_2022") %>% 
  group_by(entry_name, pro_land_cover, UD) %>% 
  mutate(lyr_14c_msp = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         CORG_msp = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate) %>% 
  distinct(lyr_14c_msp, .keep_all = TRUE) %>% 
  ungroup()

climate_14c_depth_avg <- mspline_14c_c_all %>% 
  filter(entry_name != "Giardina_2014",
           entry_name != "Grant_2022") %>% 
  add_row(andi) %>% 
  dplyr::filter(LD != 101) %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(median_pseudo = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_median = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_pseudo, .keep_all = TRUE) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  ungroup()

cd1 <- climate_14c_depth_avg %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_pseudo, color = ClimateZone), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median, fill = ClimateZone, 
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
                     limits = c(-450,110),
                     position = "top", breaks = seq(-400,100,100)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25)) +
  scale_fill_viridis_d(direction = -1) +
  scale_color_viridis_d(direction = -1)

cd2 <- climate_14c_depth_avg %>% 
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
                  breaks = seq(0,100,25)) +
  scale_color_viridis_d(direction = -1) 

ggarrange(cd1, cd2, widths = c(0.7,0.3))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_climate_avg_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

## Climate Zones and SOC
# All profiles
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
  dplyr::select(-c(lyr_14c_msp, CORG_msp)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.3) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-450,110), breaks = seq(-400,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_fill_viridis_d(direction = -1) +
  scale_color_viridis_d(direction = -1)

p_climate_all +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_all_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

# Averaged for Giardina_2014 and Grant_2022
climate_avg <- mspline_14c_c_all %>%
  filter(entry_name != "Giardina_2014",
         entry_name != "Grant_2022") %>% 
  add_row(andi) %>% 
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

p_climate <- climate_avg %>%
  filter(n > 4 & n_rel > 60)  %>% 
  dplyr::select(-c(lyr_14c_msp, CORG_msp)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.3) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-450,110), breaks = seq(-400,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_fill_viridis_d(direction = -1) +
  scale_color_viridis_d(direction = -1)

p_climate +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_avg_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

# w/o Andisols
climate_and <- mspline_14c_c_all %>%
  mutate(ClimateZone = case_when(
    pro_usda_soil_order == "Andisols" ~ "andisols",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
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

p_climate_and <- climate_and %>%
  filter(n > 4 & n_rel > 60)  %>% 
  dplyr::select(-c(lyr_14c_msp, CORG_msp)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.3) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-650,110), breaks = seq(-600,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_fill_manual(values = c("#fb6a4a", "#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF")) +
  scale_color_manual(values = c("#fb6a4a", "#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"))

p_climate_and +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_and_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

### Clay types

## Averaged andisols
clay_avg <- mspline_14c_c_all %>%
  filter(entry_name != "Giardina_2014",
         entry_name != "Grant_2022") %>%
  add_row(andi) %>%
  mutate(ClayType = case_when(
    pro_usda_soil_order == "Oxisols" ~ "low-activity",
    pro_usda_soil_order == "Ultisols" ~ "low-activity",
    TRUE ~ "high-activity"
  )) %>% 
  group_by(ClayType, UD) %>% 
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

## with depth
clay_avg %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_14c, color = ClayType), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClayType, 
                  color = ClayType), alpha = 0.5) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-450,110),
                     position = "top", breaks = seq(-400,100,100)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))   

p_clay <- clay_avg %>% 
  filter(n > 4 & n_rel > 60)  %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClayType)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClayType),
                alpha = 0.3) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-450,110), breaks = seq(-400,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") 

p_clay +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClayType),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_clay_avg_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

## w/o andisols
clay_and <- mspline_14c_c_all %>%
  # filter(entry_name != "Giardina_2014",
  #        entry_name != "Grant_2022") %>%
  # add_row(andi) %>%
  mutate(ClayType = case_when(
    pro_usda_soil_order == "Andisols" ~ "amorphous",
    pro_usda_soil_order == "Oxisols" ~ "low-activity",
    pro_usda_soil_order == "Ultisols" ~ "low-activity",
    TRUE ~ "high-activity"
  )) %>% 
  group_by(ClayType, UD) %>% 
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

## with depth
clay_and %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_14c, color = ClayType), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClayType, 
                  color = ClayType), alpha = 0.5) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-650,110),
                     position = "top", breaks = seq(-600,100,100)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))   

p_clay_and <- clay_and %>% 
  filter(n > 4 & n_rel > 60)  %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClayType)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClayType),
                alpha = 0.3) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-650,110), breaks = seq(-600,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") 

p_clay_and +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClayType),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_clay_and_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

### Climate zones and clay types

## Averaged andisols
climate_clay_avg <- mspline_14c_c_all %>%
  filter(entry_name != "Giardina_2014",
         entry_name != "Grant_2022") %>%
  add_row(andi) %>%
  mutate(ClayType = case_when(
    pro_usda_soil_order == "Oxisols" ~ "low-activity",
    pro_usda_soil_order == "Ultisols" ~ "low-activity",
    TRUE ~ "high-activity"
  )) %>% 
  group_by(ClimateZone, UD, ClayType) %>% 
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

## with depth
climate_clay_avg %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_14c, color = ClayType), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClayType, 
                  color = ClayType), alpha = 0.5) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-450,110),
                     position = "top", breaks = seq(-400,100,100)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))   

ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_climate_clay_avg_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

climate_clay_avg %>% 
  filter(ClimateZone == "temperate" & ClayType == "high-activity"|
         ClimateZone == "cold/polar" & ClayType == "high-activity"|
         ClimateZone == "temperate" & ClayType == "low-activity"|
         ClimateZone == "tropical" & ClayType == "low-activity") %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_14c, color = ClimateZone), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClimateZone, 
                  color = ClimateZone), alpha = 0.5) +
  facet_wrap(~ClayType, nrow = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-450,110),
                     position = "top", breaks = seq(-400,100,100)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25)) +
  scale_fill_manual(values = c("#35B779FF", "#31688EFF", "#440154FF")) +
  scale_color_manual(values = c("#35B779FF", "#31688EFF", "#440154FF"))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_depth_climate_clay_avg_red_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

p_climate_clay <- climate_clay_avg %>% 
  filter(n > 4 & n_rel > 60)  %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClayType)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClayType),
                alpha = 0.3) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-650,110), breaks = seq(-600,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") 

p_climate_clay +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClayType),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_clay_avg_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

climate_clay_avg %>% 
  filter(ClimateZone == "temperate" & ClayType == "high-activity"|
           ClimateZone == "cold/polar" & ClayType == "high-activity"|
           ClimateZone == "temperate" & ClayType == "low-activity"|
           ClimateZone == "tropical" & ClayType == "low-activity") %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path() +
  facet_wrap(~ClayType, nrow = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-450,110), breaks = seq(-400,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_fill_manual(values = c("#35B779FF", "#31688EFF", "#440154FF")) +
  scale_color_manual(values = c("#35B779FF", "#31688EFF", "#440154FF"))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_clay_avg_red", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

# w/o andisols

climate_clay_and <- mspline_14c_c_all %>%
  mutate(ClayType = case_when(
    pro_usda_soil_order == "Andisols" ~ "amorphous",
    pro_usda_soil_order == "Oxisols" ~ "low-activity",
    pro_usda_soil_order == "Ultisols" ~ "low-activity",
    TRUE ~ "high-activity"
  )) %>% 
  group_by(ClayType, ClimateZone, UD) %>% 
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

## with depth
climate_clay_and %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_14c, color = ClayType), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_14c, xmax = uci_14c, fill = ClayType, 
                  color = ClayType), alpha = 0.5) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-650,110),
                     position = "top", breaks = seq(-600,100,100)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25))   

p_climate_clay_and <- climate_clay_and %>% 
  filter(n > 4 & n_rel > 60)  %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClayType)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClayType),
                alpha = 0.3) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-650,110), breaks = seq(-600,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") 

p_climate_clay_and +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClayType),
                 alpha = 0.4) +
  geom_path()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_clay_and_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

climate_clay_and %>% 
  filter(ClimateZone == "temperate" & ClayType == "high-activity"|
           ClimateZone == "cold/polar" & ClayType == "high-activity"|
           ClimateZone == "temperate" & ClayType == "low-activity"|
           ClimateZone == "tropical" & ClayType == "low-activity") %>% 
  filter(n > 4 & n_rel > 60) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, color = ClimateZone),
                alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, color = ClimateZone),
                 alpha = 0.4) +
  geom_path() +
  facet_wrap(~ClayType, nrow = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-450,110), breaks = seq(-400,100,100)) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_fill_manual(values = c("#35B779FF", "#31688EFF", "#440154FF")) +
  scale_color_manual(values = c("#35B779FF", "#31688EFF", "#440154FF"))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_clay_and_red", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

