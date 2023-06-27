# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 21/02/2023 #

#### COLD TEMPERATE SOILS ONLY ####

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mpspline2)

lyr_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-02-08.csv")

ctemp_msp <- lyr_data %>% 
  filter(ClimateZone == "cold temperate") %>% 
  filter(pro_usda_soil_order != "Andisols")

ctemp_msp %>% 
  arrange(UD) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, group = id)) +
  geom_path() +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.005,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400))

ctemp_sum <- ctemp_msp %>% 
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

ctemp_sum %>% 
  filter(n > 4 & n_rel > 33) %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c)) +
  geom_path() +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c), 
                color = "red", alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c), 
                 color = "red", alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.1,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400))

## DATA exploration ##
ctemp_msp %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id),
            n_countries = n_distinct(pro_country))

ctemp_msp %>% 
  distinct(id, .keep_all = TRUE) %>% 
  count(pro_country)

ctemp_msp %>% 
  group_by(pro_usda_soil_order) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

ctemp_msp %>% 
  group_by(pro_land_cover) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

ctemp_msp %>% 
  group_by(pro_KG_present_long) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

#MAP, MAT, clay content
ctemp_msp %>% 
  distinct(id, .keep_all = TRUE) %>%
  skimr::skim(pro_MAT_mod, pro_MAP_mod, pro_AI, lyr_clay_mod,
              pro_GPP_Fluxcom_2001_2012_gC_m2d1)

ctemp_soil <- ctemp_msp %>% 
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
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = pro_usda_soil_order)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = pro_usda_soil_order), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = pro_usda_soil_order), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.4, 0.01),
        legend.justification = c("right", "bottom"),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Soil type")
plot(ctemp_soil)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_ctemp_soil_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

ctemp_mineral <- ctemp_msp %>% 
  group_by(MineralType_new, UD) %>% 
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
  ggplot(aes(x = median_c, y = median_14c, color = MineralType_new)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralType_new), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralType_new), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.5, 0.01),
        legend.justification = c("right", "bottom"),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Soil grouping")
plot(ctemp_mineral)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_ctemp_mineral_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

ctemp_climate <- ctemp_msp %>% 
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
        legend.position = c(0.95, 0.01),
        legend.justification = c("right", "bottom"),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Climate zones")
plot(ctemp_climate)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_ctempclimate_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

ctemp_msp %>% 
  group_by(pro_KG_present_long) %>% 
  # summarise(n_profile = n_distinct(id)) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  # summarise(max = max(LD)) %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod, pro_AI, lyr_clay_mod)

ctemp_climate1 <- ctemp_msp %>% 
  mutate(ClimateGroup = case_when(
    grepl("no dry", pro_KG_present_long) ~ "no dry season",
    TRUE ~ "dry season"
  )) %>% 
  group_by(ClimateGroup, UD) %>% 
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
  ggplot(aes(x = median_c, y = median_14c, color = ClimateGroup)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateGroup), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateGroup), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.95, 0.01),
        legend.justification = c("right", "bottom"),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Seasonality")
plot(wtemp_climate1)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_wtempclimate_seas_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

ctemp_clay_grp <- ctemp_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 30 ~ "< 30 %",
    TRUE ~ "> 30 %",
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
        legend.position = c(0.27, 0.25),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,30), 
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Clay content")
plot(ctemp_clay_grp)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_ctemp_clay_30_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

ctemp_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 30 ~ "< 30 %",
    TRUE ~ "> 30 %",
  )) %>% 
  group_by(clay_group, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  # summarise(n_profile = n_distinct(id)) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  # summarise(max = max(LD)) %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod)

ctemp_lc <- ctemp_msp %>% 
  filter(!is.na(pro_land_cover)) %>% 
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
        legend.position = c(0.27, 0.25),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,30), 
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Land cover")
plot(wtemp_lc)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_wtemp_lc_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

wtemp_msp %>% 
  group_by(pro_land_cover, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  # summarise(n_profile = n_distinct(id)) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(max = max(LD)) %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod)
