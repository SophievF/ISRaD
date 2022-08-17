# Explore 14C profiles in ISRaD #
# Relationship between 14C and SOC with Alox/Feox #
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
lyr_mpspline_ox <- lyr_all %>% 
  #only use studies that have Alox/Feox
  filter(lyr_obs_date_y > 1959) %>% 
  drop_na(lyr_al_ox|lyr_fe_ox) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() 

lyr_mpspline_ox %>% 
  count(id)

lyr_mpspline_ox %>% 
  ggplot(aes(x = CORG, y = lyr_14c, color = lyr_al_ox, group = id)) + 
  geom_line(orientation = "y") +
  geom_point(size = 4) +
  facet_wrap(~entry_name) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C")), limits = c(-1000, 305)) +
  scale_color_viridis_c(expression(paste("Al"[ox]," [mg/g]")), trans = "log10") +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2))

### Apply mspline function

## mspline 14C
lyr_alox_mpspline_14c <- lyr_mpspline_ox %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 1)

## mspline CORG
lyr_alox_mpspline_c <- lyr_mpspline_ox %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 60, lam = 1)

## Alox
lyr_alox_mpspline_alox <- lyr_mpspline_ox %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_al_ox) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 60, lam = 1)

## combine all
mspline_14c_c_alox <- lyr_alox_mpspline_14c$est_1cm %>% 
  rename(lyr_14c = SPLINED_VALUE) %>% 
  full_join(lyr_alox_mpspline_c$est_1cm %>% 
              rename(CORG = SPLINED_VALUE)) %>% 
  full_join(lyr_alox_mpspline_alox$est_1cm %>% 
              rename(lyr_al_ox = SPLINED_VALUE)) %>% 
  filter(LD < 101) %>% 
  dplyr::left_join(lyr_mpspline_ox %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order,
                                   pro_KG_present_long), 
                   by = "id") %>% 
  tibble()

## Plotting
mspline_14c_c_alox %>% 
  ggplot(aes(y = UD, x = lyr_al_ox, group = id)) +
  geom_path() +
  theme_classic() +
  facet_wrap(~entry_name) +
  scale_y_continuous("Depth [cm]", trans = "reverse")

mspline_14c_c_alox %>% 
  ggplot(aes(x = CORG, y = lyr_14c, color = lyr_al_ox, group = id)) + 
  geom_path() +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill =  NA)) +
  facet_wrap(~entry_name) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C")), limits = c(-1000, 305)) +
  scale_color_viridis_c(expression(paste("Al"[ox]," [mg/g]")), trans = "log10") +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2))
