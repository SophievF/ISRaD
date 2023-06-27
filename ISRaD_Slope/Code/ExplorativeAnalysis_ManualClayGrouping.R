# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 16/01/2023 #

library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# Load filtered database data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2023-02-08"))

lyr_all %>% 
  skimr::skim_without_charts(lyr_14c, CORG)

# Load filtered and splined data
lyr_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-03-09.csv") %>% 
  #remove layers that have CORG > 20
  filter(CORG_msp <= 20)

## Define color code 
#cold, warm, tropical, arid
color_climate_wo_polar <- c("#762a83", "#7fbf7b", "#1b7837", "#dfc27d")
# color_climate_wo_polar <- c("#a6dba0", "#c2a5cf", "#7b3294", "#dfc27d")
#polar, cold, warm, tropical, arid
color_climate <- c("#af8dc3", "#762a83", "#7fbf7b", "#1b7837", "#dfc27d")
# color_climate <- c("#008837", "#a6dba0", "#c2a5cf", "#7b3294", "#dfc27d")
#amorphous, high-activity clay, low-activity clays
color_mineral <- c("#225ea8", "#41b6c4", "#a1dab4")

#### ARID SOILS ONLY ####

arid_msp <- lyr_data %>% 
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
  geom_path(linewidth = 1) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c), 
                    color = "red", alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c), 
                     color = "red", alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  coord_trans(x = "log10", xlim = c(0.05,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous("SOC [wt-%]", limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
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
  skimr::skim_without_charts(pro_MAT_mod, pro_MAP_mod, lyr_clay_mod, pro_PET_mm_yr_mod, 
                             pro_AI)

arid_msp %>% 
  distinct(id, .keep_all = TRUE) %>% 
  skimr::skim_without_charts(pro_MAT_mod, pro_MAP_mod, lyr_clay_mod, pro_PET_mm_yr_mod, 
              pro_AI)

### Manual grouping: clay content ###
#10%, 20%, 30%
arid_clay_bins <- arid_msp %>% 
  mutate(clay_group_10 = case_when(
    lyr_clay_mod < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
  )) %>% 
  mutate(clay_group_20 = case_when(
    lyr_clay_mod < 20 ~ "< 20 %",
    TRUE ~ "> 20 %",
  )) %>% 
  mutate(clay_group_30 = case_when(
    lyr_clay_mod < 30 ~ "< 30 %",
    TRUE ~ "> 30 %",
  )) %>% 
  dplyr::select(id, site_name, UD, clay_group_10, clay_group_20, clay_group_30, 
                lyr_14c_msp, CORG_msp)

arid_clay_grp_10 <- arid_clay_bins %>% 
  group_by(clay_group_10, UD) %>% 
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
  arrange(UD) 

pc10 <- arid_clay_grp_10 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_10)) +
  geom_path(linewidth = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_10), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_10), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

arid_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
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

arid_clay_grp_20 <- arid_clay_bins %>% 
  group_by(clay_group_20, UD) %>% 
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
  arrange(UD) 

pc20 <- arid_clay_grp_20 %>%  
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_20)) +
  geom_path(linewidth = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_20), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_20), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous("", expand = c(0,0), limits = c(-1000,125), 
                     breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

arid_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 20 ~ "< 20 %",
    TRUE ~ "> 20 %",
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

arid_clay_grp_30 <- arid_clay_bins %>% 
  group_by(clay_group_30, UD) %>% 
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
  arrange(UD) 

pc30 <- arid_clay_grp_30 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_30)) +
  geom_path(linewidth = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_30), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_30), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

arid_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 30 ~ "< 30 %",
    TRUE ~ "> 30 %",
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

annotate_figure(ggpubr::ggarrange(pc10, pc20, pc30, nrow = 1),
                bottom = text_grob("Soil organic carbon [wt-%]; log-scaled", 
                                   size = 16))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_arid_clay_bin_", Sys.Date(),
                     ".jpeg"), width = 12, height = 5)

arid_res <- arid_clay_grp_20 %>% 
  dplyr::select(UD, clay_group_20, median_14c, median_c) %>% 
  ungroup() %>% 
  pivot_wider(names_from = clay_group_20, values_from = c(median_14c, median_c)) %>% 
  dplyr::mutate(res_14c = `median_14c_> 20 %` - `median_14c_< 20 %`,
                res_c = `median_c_> 20 %` - `median_c_< 20 %`) %>% 
  mutate(Group = "arid") %>% 
  dplyr::select(UD, res_14c, res_c, Group)

#### WARM TEMPERATE SOILS ONLY ####
wtemp_msp <- lyr_data %>% 
  filter(ClimateZone == "warm temperate") %>% 
  filter(pro_usda_soil_order != "Andisols")

wtemp_msp %>% 
  arrange(UD) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, group = id)) +
  geom_path() +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", limits = c(0.005,57)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400))

wtemp_sum <- wtemp_msp %>% 
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

wtemp_sum %>% 
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
  coord_trans(x = "log10", xlim = c(0.05,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous("SOC [wt-%]", limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400))

## DATA exploration ##
wtemp_msp %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id),
            n_countries = n_distinct(pro_country))

wtemp_msp %>% 
  distinct(id, .keep_all = TRUE) %>% 
  count(pro_country)

wtemp_msp %>% 
  group_by(pro_usda_soil_order) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

wtemp_msp %>% 
  group_by(pro_land_cover) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

wtemp_msp %>% 
  group_by(pro_KG_present_long) %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

#MAP, MAT, clay content
wtemp_msp %>% 
  distinct(id, .keep_all = TRUE) %>%
  skimr::skim_without_charts(pro_MAT_mod, pro_MAP_mod, pro_AI, lyr_clay_mod)

### Manual grouping: clay content ###
#10%, 20%, 30%
wtemp_clay_bins <- wtemp_msp %>% 
  mutate(clay_group_10 = case_when(
    lyr_clay_mod < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
  )) %>% 
  mutate(clay_group_20 = case_when(
    lyr_clay_mod < 20 ~ "< 20 %",
    TRUE ~ "> 20 %",
  )) %>% 
  mutate(clay_group_30 = case_when(
    lyr_clay_mod < 30 ~ "< 30 %",
    TRUE ~ "> 30 %",
  )) %>% 
  dplyr::select(id, site_name, UD, clay_group_10, clay_group_20, clay_group_30, 
                lyr_14c_msp, CORG_msp)

wtemp_clay_grp_10 <- wtemp_clay_bins %>% 
  group_by(clay_group_10, UD) %>% 
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
  arrange(UD) 

pc10_wtemp <- wtemp_clay_grp_10 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_10)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_10), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_10), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

wtemp_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
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

wtemp_clay_grp_20 <- wtemp_clay_bins %>% 
  group_by(clay_group_20, UD) %>% 
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
  arrange(UD) 

pc20_wtemp <- wtemp_clay_grp_20 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_20)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_20), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_20), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

wtemp_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 20 ~ "< 20 %",
    TRUE ~ "> 20 %",
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

wtemp_clay_grp_30 <- wtemp_clay_bins %>% 
  group_by(clay_group_30, UD) %>% 
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
  arrange(UD) 

pc30_wtemp <- wtemp_clay_grp_30 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_30)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_30), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_30), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

wtemp_msp %>% 
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

annotate_figure(ggpubr::ggarrange(pc10_wtemp, pc20_wtemp, pc30_wtemp, nrow = 1),
                bottom = text_grob("Soil organic carbon [wt-%]; log-scaled", 
                                   size = 16))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_wtemp_clay_bin_", Sys.Date(),
                     ".jpeg"), width = 12, height = 5)

wtemp_res <- wtemp_clay_grp_20 %>% 
  dplyr::select(UD, clay_group_20, median_14c, median_c) %>% 
  ungroup() %>% 
  pivot_wider(names_from = clay_group_20, values_from = c(median_14c, median_c)) %>% 
  dplyr::mutate(res_14c = `median_14c_> 20 %` - `median_14c_< 20 %`,
                res_c = `median_c_> 20 %` - `median_c_< 20 %`) %>% 
  mutate(Group = "warm temperate") %>% 
  dplyr::select(UD, res_14c, res_c, Group)

#### COLD TEMPERATE SOILS ONLY ####
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
  coord_trans(x = "log10", xlim = c(0.05,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous("SOC [wt-%]", limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
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
  skimr::skim_without_charts(pro_MAT_mod, pro_MAP_mod, pro_AI, lyr_clay_mod)

### Manual grouping: clay content ###
#10%, 20%, 30%
ctemp_clay_bins <- ctemp_msp %>% 
  mutate(clay_group_10 = case_when(
    lyr_clay_mod < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
  )) %>% 
  mutate(clay_group_20 = case_when(
    lyr_clay_mod < 20 ~ "< 20 %",
    TRUE ~ "> 20 %",
  )) %>% 
  mutate(clay_group_30 = case_when(
    lyr_clay_mod < 30 ~ "< 30 %",
    TRUE ~ "> 30 %",
  )) %>% 
  dplyr::select(id, site_name, UD, clay_group_10, clay_group_20, clay_group_30, 
                lyr_14c_msp, CORG_msp)

ctemp_clay_grp_10 <- ctemp_clay_bins %>% 
  group_by(clay_group_10, UD) %>% 
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
  arrange(UD) 

pc10_ctemp <- ctemp_clay_grp_10 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_10)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_10), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_10), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

ctemp_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 10 ~ "< 10 %",
    TRUE ~ "> 10 %",
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

ctemp_clay_grp_20 <- ctemp_clay_bins %>% 
  group_by(clay_group_20, UD) %>% 
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
  arrange(UD) 

pc20_ctemp <- ctemp_clay_grp_20 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_20)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_20), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_20), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

ctemp_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 20 ~ "< 20 %",
    TRUE ~ "> 20 %",
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

ctemp_clay_grp_30 <- ctemp_clay_bins %>% 
  group_by(clay_group_30, UD) %>% 
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
  arrange(UD) 

pc30_ctemp <- ctemp_clay_grp_30 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_30)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_30), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_30), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

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

annotate_figure(ggpubr::ggarrange(pc10_ctemp, pc20_ctemp, pc30_ctemp, nrow = 1),
                bottom = text_grob("Soil organic carbon [wt-%]; log-scaled", 
                                   size = 16))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_ctemp_clay_bin_", Sys.Date(),
                     ".jpeg"), width = 12, height = 5)

ctemp_res <- ctemp_clay_grp_20 %>% 
  dplyr::select(UD, clay_group_20, median_14c, median_c) %>% 
  ungroup() %>% 
  pivot_wider(names_from = clay_group_20, values_from = c(median_14c, median_c)) %>% 
  dplyr::mutate(res_14c = `median_14c_> 20 %` - `median_14c_< 20 %`,
                res_c = `median_c_> 20 %` - `median_c_< 20 %`) %>% 
  mutate(Group = "cold temperate") %>% 
  dplyr::select(UD, res_14c, res_c, Group)

#### TROPICAL SOILS ONLY ####
trop_msp <- lyr_data %>% 
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
  coord_trans(x = "log10", xlim = c(0.05,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous("SOC [wt-%]", limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400))

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
  skimr::skim_without_charts(pro_MAT_mod, pro_MAP_mod, lyr_clay_mod)

### Manual grouping: clay content ###
#only 30 & 40% possible
trop_clay_bins <- trop_msp %>% 
  mutate(clay_group_30 = case_when(
    lyr_clay_mod < 30 ~ "< 30 %",
    TRUE ~ "> 30 %",
  )) %>% 
  mutate(clay_group_40 = case_when(
    lyr_clay_mod < 40 ~ "< 40 %",
    TRUE ~ "> 40 %",
  )) %>% 
  dplyr::select(id, site_name, UD, clay_group_30, clay_group_40, 
                lyr_14c_msp, CORG_msp)

trop_clay_grp_30 <- trop_clay_bins %>% 
  group_by(clay_group_30, UD) %>% 
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
  arrange(UD) 

pc30_trop <- trop_clay_grp_30 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_30)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_30), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_30), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

trop_msp %>% 
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

trop_clay_grp_40 <- trop_clay_bins %>% 
  group_by(clay_group_40, UD) %>% 
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
  arrange(UD) 

pc40_trop <- trop_clay_grp_40 %>% 
  ggplot(aes(x = median_c, y = median_14c, color = clay_group_40)) +
  geom_path(size = 1.5) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = clay_group_40), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = clay_group_40), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank()) +
  coord_trans(x = "log10", xlim = c(0.08,10)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) + 
  scale_x_continuous(limits = c(0.05,10), 
                     expand = c(0,0), breaks = c(0.1, 0.5,1,5,10),
                     labels = c(0.1, 0.5,1,5,10)) +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Clay content")

trop_msp %>% 
  mutate(clay_group = case_when(
    lyr_clay_mod < 40 ~ "< 40 %",
    TRUE ~ "> 40 %",
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

annotate_figure(ggpubr::ggarrange(pc30_trop, pc40_trop,nrow = 1),
                bottom = text_grob("Soil organic carbon [wt-%]; log-scaled", 
                                   size = 16))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_trop_clay_bin_", Sys.Date(),
                     ".jpeg"), width = 12, height = 5)

trop_res <- trop_clay_grp_30 %>% 
  dplyr::select(UD, clay_group_30, median_14c, median_c) %>% 
  ungroup() %>% 
  pivot_wider(names_from = clay_group_30, values_from = c(median_14c, median_c)) %>% 
  dplyr::mutate(res_14c = `median_14c_> 30 %` - `median_14c_< 30 %`,
                res_c = `median_c_> 30 %` - `median_c_< 30 %`) %>% 
  mutate(Group = "tropical") %>% 
  dplyr::select(UD, res_14c, res_c, Group)

climate_res <- rbind(arid_res, wtemp_res, ctemp_res, trop_res) %>% 
  tibble() 

climate_res_depth <- climate_res %>% 
  mutate(layer = case_when(
    UD == 1 ~ 0,
    UD == 10 ~ 10,
    UD == 20 ~ 20,
    UD == 30 ~ 30,
    UD == 40 ~ 40,
    UD == 50 ~ 50,
    UD == 60 ~ 60,
    UD == 70 ~ 70,
    UD == 80 ~ 80,
    UD == 90 ~ 90,
    UD == 100 ~ 100,
    TRUE ~ NA
  )) %>% 
  drop_na()

climate_res %>% 
  mutate(Group = factor(Group, levels = c("cold temperate", "warm temperate",
                                         "tropical", "arid"))) %>% 
  ggplot(aes(x = res_c, y = res_14c)) +
  geom_path(aes(color = Group), linewidth = 1.5) +
  geom_point(data = climate_res_depth, aes(fill = as.factor(layer)), size = 2,
             shape = 21, color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual("Climate grouping", values = color_climate_wo_polar) +
  scale_x_continuous("Residuals SOC: high - low clay", expand = c(0,0),
                     limits = c(-2,2)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-300,100)) +
  scale_fill_discrete("Depth layers")
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_clay_bin_res_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)
  
annotate_figure(ggarrange(pc30_ctemp, pc30_wtemp, pc30_trop, pc30),
                bottom = text_grob("Soil organic carbon [wt-%]; log-scaled", 
                                   size = 16))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_all_clay_bin_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)
