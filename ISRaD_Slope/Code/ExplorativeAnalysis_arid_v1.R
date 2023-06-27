# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 16/01/2023 #

#### ARID SOILS ONLY ####

library(tidyverse)
library(ggpubr)
library(RColorBrewer)

mspline_14c_c_all <- read_csv("./Data/ISRaD_flat_splined_filled_2023-02-06.csv")

arid_msp <- mspline_14c_c_all %>% 
  filter(ClimateZone == "arid") %>% 
  filter(pro_usda_soil_order != "Andisols")

arid_msp %>% 
  arrange(UD) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, group = id)) +
  geom_path() +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.008,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400))
  
arid_sum <- arid_msp %>% 
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

arid_sum %>% 
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
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_arid_only_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)


## DATA exploration ##
arid_msp %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id),
            n_countries = n_distinct(pro_country))

arid_msp %>% 
  distinct(id, .keep_all = TRUE) %>% 
  count(pro_country)

arid_msp %>% 
  group_by(pro_usda_soil_order) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

arid_msp %>% 
  group_by(pro_land_cover) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

# cultivated based on google earth
arid_msp %>% 
  group_by(pro_land_cover) %>% 
  filter(is.na(pro_land_cover)) %>% 
  count(entry_name, site_name, id)

arid_clt <- arid_msp %>% 
  mutate(pro_land_cover_mod = case_when(
    is.na(pro_land_cover) ~ "cultivated",
    TRUE ~ pro_land_cover
  ))

arid_clt %>% 
  group_by(pro_land_cover_mod) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

arid_msp %>% 
  group_by(pro_KG_present_long) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

#MAP, MAT, clay content
arid_msp %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod, lyr_clay_mod, pro_PET_mm_yr_mod, 
              pro_GPP_Fluxcom_2001_2012_gC_m2d1, pro_AI)

arid_msp %>% 
  distinct(id, .keep_all = TRUE) %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod, lyr_clay_mod, pro_PET_mm_yr_mod, 
              pro_GPP_Fluxcom_2001_2012_gC_m2d1, pro_AI)

### FURTHER data exploration ###
arid_soil <- arid_msp %>% 
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
        legend.position = c(0.8, 0.16),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10),
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Soil type")
plot(arid_soil)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_arid_soil_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

arid_msp %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(max = max(LD))

arid_lc <- arid_clt %>% 
  group_by(pro_land_cover_mod, UD) %>% 
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
  ggplot(aes(x = median_c, y = median_14c, color = pro_land_cover_mod)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = pro_land_cover_mod), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = pro_land_cover_mod), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.7, 0.15),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10),
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Land cover")
plot(arid_lc)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_arid_lc_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

arid_clt %>% 
  group_by(pro_land_cover_mod, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(max = max(LD))

arid_climate <- arid_msp %>% 
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
plot(arid_climate)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_arid_climate_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

arid_msp %>% 
  group_by(pro_KG_present_long, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  summarise(max = max(LD))

arid_climate_mod <- arid_msp %>% 
  mutate(climate_mod = case_when(
    grepl("Temperate", pro_KG_present_long) ~ "cold",
    grepl("hot", pro_KG_present_long) ~ "hot",
    grepl("cold", pro_KG_present_long) ~ "cold",
  )) %>% 
  group_by(climate_mod, UD) %>% 
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
  ggplot(aes(x = median_c, y = median_14c, color = climate_mod)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = climate_mod), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = climate_mod), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10), 
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_manual("Climate", values = c("#00BFC4", "#F8766D"))
plot(arid_climate_mod)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_arid_climate_mod_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

arid_msp %>% 
  mutate(climate_mod = case_when(
    grepl("Temperate", pro_KG_present_long) ~ "cold",
    grepl("hot", pro_KG_present_long) ~ "hot",
    grepl("cold", pro_KG_present_long) ~ "cold",
  )) %>% 
  group_by(climate_mod, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  # summarise(max = max(LD)) %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod, pro_AI, lyr_clay_mod)

arid_climate_mod_2 <- arid_msp %>% 
  mutate(climate_mod = case_when(
    grepl("desert", pro_KG_present_long) ~ "desert",
    grepl("steppe", pro_KG_present_long) ~ "steppe",
    TRUE ~ "steppe"
  )) %>% 
  group_by(climate_mod, UD) %>% 
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
  ggplot(aes(x = median_c, y = median_14c, color = climate_mod)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = climate_mod), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = climate_mod), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10),
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Climate")
plot(arid_climate_mod_2)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_arid_climate_mod_2_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

arid_msp %>% 
  mutate(climate_mod = case_when(
    grepl("desert", pro_KG_present_long) ~ "desert",
    grepl("steppe", pro_KG_present_long) ~ "steppe",
    TRUE ~ "steppe"
  )) %>% 
  group_by(climate_mod, UD) %>%
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  # summarise(max = max(LD)) %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod, pro_AI, lyr_clay_mod)

arid_clay_grp <- arid_msp %>% 
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
        legend.position = c(0.8, 0.15),
        legend.background = element_blank()) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.05,10),
                     expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_discrete("Clay content")
plot(arid_clay_grp)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_arid_clay_30_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

arid_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 5 ~ "< 5 %",
    TRUE ~ "> 5 %",
  )) %>% 
  group_by(clay_group) %>% 
  # summarise(n_profile = n_distinct(id)) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 33) %>% 
  # summarise(max = max(LD)) %>% 
  skimr::skim(pro_MAT_mod, pro_MAP_mod, pro_AI, lyr_clay_mod)


