# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC: all profiles #
# Sophie von Fromm #
# 15/09/2022 #

library(tidyverse)
library(ggpubr)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-10-05"))

lyr_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

lyr_data <- lyr_all %>% 
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
  ungroup() 

lyr_data$ClimateZone <- factor(lyr_data$ClimateZone,
                               levels = c("andisols", "polar", "cold",
                                          "temperate", "arid", "tropical"))

### Density distribution
p1 <- lyr_data %>% 
  filter(depth < 101) %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(trans = "log10", direction = -1, limits = c(1,340)) +
  coord_cartesian(ylim = c(-1000,400))

p2 <- lyr_data %>% 
  filter(depth < 101) %>% 
  ggplot(aes(x = depth, y = lyr_14c)) + 
  geom_hex(binwidth = c(3,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(-5,110)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,400)) +
  scale_fill_viridis_c(direction = -1, trans = "log10", limits = c(1,340)) +
  coord_cartesian(xlim = c(0,101))

ggarrange(p1, p2, common.legend = TRUE)

ggsave(file = paste0("./Figure/ISRaD_14C_SOC_depth_hex_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

lyr_data %>% 
  filter(depth < 101) %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1) +
  coord_cartesian(ylim = c(-1000,400)) 

ggsave(file = paste0("./Figure/ISRaD_14C_SOC_climate_hex_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

lyr_data %>% 
  filter(depth < 101) %>% 
  ggplot(aes(x = depth, y = lyr_14c)) + 
  geom_hex(binwidth = c(2,50)) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(-5,510)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1300,400)) +
  scale_fill_viridis_c(direction = -1, trans = "log10", limits = c(1,340)) +
  coord_cartesian(xlim = c(0,100), ylim = c(-1000,400)) 

ggsave(file = paste0("./Figure/ISRaD_14C_depth_climate_hex_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

lyr_data %>% 
  filter(depth < 101) %>% 
  filter(pro_usda_soil_order != "Aridisols",
         pro_usda_soil_order != "Histosols") %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(binwidth = c(0.1,50)) +
  facet_wrap(~pro_usda_soil_order) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_fill_viridis_c(direction = -1) +
  coord_cartesian(ylim = c(-1000,400)) 

ggsave(file = paste0("./Figure/ISRaD_14C_SOC_soiltype_hex_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

lyr_data %>%
  filter(depth < 101) %>% 
  filter(pro_usda_soil_order != "Aridisols",
         pro_usda_soil_order != "Histosols") %>% 
  ggplot(aes(x = depth, y = lyr_14c)) + 
  geom_hex(binwidth = c(2,50)) +
  facet_wrap(~pro_usda_soil_order) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(-5,510)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1300,400)) +
  scale_fill_viridis_c(direction = -1, trans = "log10", limits = c(1,340)) +
  coord_cartesian(xlim = c(0,100), ylim = c(-1000,400)) 

ggsave(file = paste0("./Figure/ISRaD_14C_depth_soiltype_hex_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

lyr_data %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  drop_na(lyr_al_ox) %>% 
  ggplot(aes(x = CORG, y = lyr_14c, color = lyr_al_ox)) + 
  geom_point(size = 3) +
  facet_wrap(~ClimateZone) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1005,400)) +
  scale_color_viridis_c(trans = "log", option = "D") 

## Data distribution with raw data
library(nlme)

lyr_data_101 <- lyr_data %>% 
  filter(lyr_bot <= 101)

lm_all <- lm(lyr_14c ~ log10(CORG)*ClimateZone, data = lyr_data)
plot(lm_all)
summary(lm_all)

lyr_data$predicted <- predict(lm_all)


lyr_data %>% 
  filter((lyr_bot-lyr_top)/2+lyr_top <= 101) %>% 
  arrange(id, lyr_top) %>% 
  ggplot(aes(x = CORG, y = lyr_14c, fill = ClimateZone)) +
  geom_point(size = 3, shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  facet_wrap(~ClimateZone)

lyr_data %>% 
  filter((lyr_bot-lyr_top)/2+lyr_top <= 101) %>% 
  arrange(id, lyr_top) %>% 
  ggplot(aes(x = (lyr_bot-lyr_top)/2+lyr_top, y = lyr_14c, fill = ClimateZone)) +
  geom_point(size = 3, shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  facet_wrap(~ClimateZone) 



lyr_data %>% 
  filter(ClimateZone == "temperate"|
           ClimateZone == "cold") %>% 
  ggplot(aes(x = CORG, y = lyr_14c, fill = ClimateZone)) +
  geom_point(size = 3, shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  facet_wrap(~pro_usda_soil_order) 
  
