# Data analysis: Profile data ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 08/03/2023 #

## Derivatives of carbon with depth (based on Carlos' idea)

library(tidyverse)

# Load splined and filtered data
lyr_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-02-08.csv")


lyr_data_deri <- lyr_data %>% 
  dplyr::select(id, UD, LD, CORG_msp, lyr_14c_msp, ClimateZoneAnd, MineralType,
                site_name) %>% 
  group_by(id) %>% 
  mutate(CORG_d1 = lead(CORG_msp, 1) - CORG_msp,
         CORG_d2 = lead(CORG_d1, 1) - CORG_d1) %>% 
  ungroup()

# Example
lyr_data_deri %>% 
  filter(id == "Baisden_2007_China hat_China hat_122") %>% 
  pivot_longer(cols = c(CORG_msp, CORG_d1, CORG_d2), 
               names_to = "derivative", values_to = "CORG") %>% 
  arrange(UD) %>% 
  ggplot(aes(x = CORG, y = UD, color = derivative)) +
  geom_path() +
  scale_y_reverse("Depth", limits = c(100,0), expand = c(0,0)) +
  theme_bw()

lyr_data_deri %>% 
  filter(id == "Baisden_2007_China hat_China hat_122") %>% 
  pivot_longer(cols = c(CORG_msp, CORG_d1, CORG_d2), 
               names_to = "derivative", values_to = "CORG") %>% 
  arrange(UD) %>% 
  ggplot(aes(x = CORG, y = lyr_14c_msp, color = derivative)) +
  geom_path() +
  theme_bw()

## Climate grouping
# by depth
climate_depth <- lyr_data_deri %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(Derivative_1 = median(CORG_d1, na.rm = TRUE),
         Derivative_2 = median(CORG_d2, na.rm = TRUE),
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_depth$ClimateZoneAnd <- factor(climate_depth$ClimateZoneAnd,
                                       levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                  "warm temperate", "arid", "tropical"))

c_depth_1 <- climate_depth %>% 
  drop_na() %>% 
  filter(n > 4 & n_rel > 33 & n_site > 3) %>% 
  ggplot() + 
  geom_path(aes(x = Derivative_1, y = UD, color = ClimateZoneAnd), linewidth = 1.5) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        legend.background = element_blank(), 
        plot.margin = margin(l = 25, t = 25, b = 25)) +
  scale_x_continuous("Change in soil organic carbon [gC/gSoil/cm]", 
                     limits = c(-0.6, 0.05), expand = c(0,0),
                     breaks = seq(-0.6,0,0.2)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0)) +
  scale_color_discrete("Profile grouping")

c_depth_2 <- climate_depth %>% 
  drop_na() %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(x = Derivative_2, y = UD, color = ClimateZoneAnd), linewidth = 1.5) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        legend.background = element_blank(), 
        plot.margin = margin(r = 25, t = 25, b = 25)) +
  scale_x_continuous("Change in soil organic carbon [gC/gSoil/cm²]", 
                      expand = c(0,0), limits = c(-0.09, 0.03)) +
  scale_y_reverse("", limits = c(100,0), expand = c(0,0)) +
  scale_color_discrete("Profile grouping")

ggpubr::ggarrange(c_depth_1, c_depth_2, common.legend = TRUE, 
                  labels = c("a) 1st derivative", "b) 2nd derivative"), 
                  font.label = list(face = "plain"), hjust = c(-0.8, -0.5))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_depth_deriv_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)


# by 14C
climate_14c <- lyr_data_deri %>%
  drop_na() %>% 
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = median(lyr_14c_msp),
         Derivative_1 = median(CORG_d1, na.rm = TRUE),
         Derivative_2 = median(CORG_d2, na.rm = TRUE),
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_14c$ClimateZoneAnd <- factor(climate_14c$ClimateZoneAnd,
                                   levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                              "warm temperate", "arid", "tropical"))

depth_sum <- climate_14c %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 98)

c1 <- climate_14c %>% 
  drop_na() %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(x = Derivative_1, y = median_14c, color = ClimateZoneAnd), linewidth = 1.5) +
  geom_point(data = depth_sum, aes(x = Derivative_1, y = median_14c), size = 1.5,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        legend.background = element_blank(), 
        plot.margin = margin(l = 25, t = 25, b = 25)) +
  scale_x_continuous("Change in soil organic carbon [gC/gSoil/cm]", 
                     limits = c(-0.6, 0.05), expand = c(0,0),
                     breaks = seq(-0.6,0,0.2)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Profile grouping")

c2 <- climate_14c %>% 
  drop_na() %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(x = Derivative_2, y = median_14c, color = ClimateZoneAnd), linewidth = 1.5) +
  geom_point(data = depth_sum, aes(x = Derivative_2, y = median_14c), size = 1.5,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        legend.background = element_blank(), 
        plot.margin = margin(r = 25, t = 25, b = 25)) +
  scale_x_continuous("Change in soil organic carbon [gC/gSoil/cm²]", 
                     expand = c(0,0), limits = c(-0.09, 0.03)) +
  scale_y_continuous("", expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) +
  scale_color_discrete("Profile grouping")


ggpubr::ggarrange(c1, c2, common.legend = TRUE, 
                  labels = c("a) 1st derivative", "b) 2nd derivative"), 
                  font.label = list(face = "plain"), hjust = c(-1, -0.6))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_deriv_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)
