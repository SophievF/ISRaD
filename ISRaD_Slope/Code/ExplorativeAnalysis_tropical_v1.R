# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 17/01/2023 #

#### TROPICAL SOILS ONLY ####

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mpspline2)

#Load filtered lyr data
mspline_14c_c_all <- read_csv("./Data/ISRaD_flat_splined_filled_2023-02-06.csv")

trop_msp <- mspline_14c_c_all %>% 
  filter(ClimateZone == "tropical") %>% 
  filter(pro_usda_soil_order != "Andisols")

trop_msp %>% 
  arrange(UD) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, group = id)) +
  geom_path() +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400))

trop_sum <- trop_msp %>% 
  group_by(UD) %>% 
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

trop_sum %>% 
  filter(n > 4 & n_rel > 33) %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c)) +
  geom_path(size = 1) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c), 
                color = "red", alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c), 
                 color = "red", alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10), 
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_trop_only_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

## DATA exploration ##
trop_msp %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id),
            n_countries = n_distinct(pro_country))

trop_msp %>% 
  distinct(id, .keep_all = TRUE) %>% 
  count(pro_country)

trop_msp %>% 
  group_by(pro_usda_soil_order) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

trop_msp %>% 
  group_by(pro_land_cover) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

trop_msp %>% 
  group_by(pro_KG_present_long) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

#MAP, MAT, clay content
trop_msp %>% 
  distinct(id, .keep_all = TRUE) %>%
  skimr::skim(pro_MAT_mod, pro_MAP_mod, lyr_clay_mod)

trop_msp %>% 
  filter(lyr_clay_mod == 0) %>% 
  count(entry_name, pro_name)

### FURTHER data exploration ###
trop_mineral <- trop_msp %>% 
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
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = MineralType)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralType), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralType), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.16),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10), 
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Mineral type")
plot(trop_mineral)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_trop_mineral_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

trop_msp %>% 
  group_by(MineralType, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(max = max(LD))

trop_climate <- trop_msp %>% 
  group_by(pro_KG_present_long, UD) %>% 
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
  filter(n > 4 & n_rel > 33) %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = pro_KG_present_long)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = pro_KG_present_long), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = pro_KG_present_long), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.7, 0.16),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10), 
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Climate")
plot(trop_climate)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_trop_climate_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

trop_clay %>% 
  group_by(pro_KG_present_long, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(max = max(LD))

trop_lc <- trop_msp %>% 
  #number of profiles likely to low ~ 6
  filter(pro_land_cover != "shrubland") %>% 
  group_by(pro_land_cover, UD) %>% 
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
  filter(n > 4 & n_rel > 33) %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = pro_land_cover)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = pro_land_cover), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = pro_land_cover), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.16),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10), 
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Land cover")
plot(trop_lc)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_trop_lc_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

trop_msp %>% 
  group_by(pro_land_cover, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(max = max(LD))

trop_clay_grp <- trop_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 40 ~ "< 40 %",
    TRUE ~ "> 40 %",
  )) %>% 
  group_by(clay_group, UD) %>% 
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
  filter(n > 4 & n_rel > 33) %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.16),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10), 
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Clay content")
plot(trop_clay_grp)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_trop_clay_40_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

trop_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 40 ~ "< 40 %",
    TRUE ~ "> 40 %",
  )) %>% 
  group_by(clay_group, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  summarise(n_profile = n_distinct(id)) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(max = max(LD)) %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod, lyr_clay_mod)





