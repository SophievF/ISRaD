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
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-08-12"))

# Filter data for mspline function
lyr_mpspline <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup()

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
  group_by(pro_wrb_soil_order, UD) %>% 
  summarise(n = n_distinct(id)) %>% 
  filter(UD == 1| UD == 97) %>% view()

