# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 25/01/2023 #

#### Individual Sites ####

library(tidyverse)
library(ggpubr)
library(RColorBrewer)

### Load filtered and splined lyr data
mspline_14c_c_all <- read_csv("./Data/ISRaD_flat_splined_filled_2023-02-08.csv") %>% 
  arrange(UD, id)

head(mspline_14c_c_all)
names(mspline_14c_c_all)

mspline_14c_c_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

# Load original data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2023-02-08"))

# Santa Cruz Chronosequence
mspline_14c_c_all %>% 
  filter(grepl("cruz|Cruz", site_name)) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, color = pro_name)) +
  geom_path(linewidth = 1.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [%]; log-scaled", trans = "log10", limits = c(0.1,10)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250))
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_SantaCruz_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

# Mattole Chrsonosequence
mspline_14c_c_all %>% 
  arrange(UD) %>% 
  filter(grepl("Mattole", site_name)) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, color = pro_name)) +
  geom_path(linewidth = 1.5) +
  theme_bw(base_size = 14) +
  facet_wrap(~entry_name) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [%]; log-scaled", trans = "log10", limits = c(0.1,10)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_Mattole_", Sys.Date(),
                     ".jpeg"), width = 9.5, height = 5.5)

mspline_14c_c_all %>% 
  arrange(UD) %>% 
  filter(grepl("Mattole", site_name)) %>% 
  group_by(pro_name) %>% 
  summarise(min = which.min(UD),
            max = which.max(UD))

mspline_14c_c_all %>% 
  arrange(UD) %>% 
  filter(grepl("Khomo", entry_name)) %>% 
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp, color = pro_name)) +
  geom_path(linewidth = 1.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [%]; log-scaled", trans = "log10", limits = c(0.01,10)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,150), breaks = seq(-1000,0,250))

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_Khomo_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

mspline_14c_c_all %>% 
  arrange(UD) %>% 
  filter(grepl("Khomo", entry_name)) %>% 
  group_by(pro_name) %>% 
  summarise(min = which.min(UD),
            max = which.max(UD))

lyr_all %>% 
  arrange(lyr_top) %>% 
  filter(grepl("Khomo", entry_name)) %>% 
  ggplot(aes(x = CORG, y = lyr_14c, color = pro_name)) +
  geom_path(linewidth = 1.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [%]; log-scaled", trans = "log10", limits = c(0.01,10)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,150), breaks = seq(-1000,0,250))

lyr_all %>% 
  filter(grepl("Khomo", entry_name)) %>% 
  group_by(pro_name) %>% 
  summarise(max = which.max(lyr_bot))
