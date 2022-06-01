# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 01/06/2022 #

library(ISRaD)
library(tidyverse)
library(ggpubr)

#Load filtered lyr data
lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_", Sys.Date()))

lyr_data %>% 
  count(entry_name)

names(lyr_data)

## Mapping sampling locations ##
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)

world <- map_data("world") %>% 
  filter(region != "Antarctica")

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgrey")  +
  geom_point(data = lyr_data, 
             aes(x = pro_long, y = pro_lat),
             color = "#4D36C6", shape = 1, size = 3) +
  theme_bw(base_size = 14) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")) +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), 
                     breaks = c(-100,0,100), limits = c(-160,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), 
                     breaks = c(-50,0,50), limits = c(-55,80))
## Data exploration ##

#lyr_14c
plotly::ggplotly(
  lyr_data %>% 
    ggplot(aes(x = depth, y = lyr_14c, group = entry_name)) + 
    geom_point(size = 3, shape = 21) +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_x_continuous("Depth [cm]", expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(-1505,305)) 
)

#lyr_dd14c
plotly::ggplotly(
  lyr_data %>% 
    ggplot(aes(x = depth, y = lyr_dd14c, group = entry_name)) + 
    geom_point(size = 3, shape = 21) +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_x_continuous("Depth [cm]", expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(-1505,305)) 
)

# Density distribution
p1 <- lyr_data %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_hex(color = NA, binwidth = c(0.1,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(limits = c(-1505,305)) +
  scale_fill_viridis_c(trans = "log10")

p2 <- lyr_data %>% 
  ggplot(aes(x = CORG, y = lyr_dd14c)) + 
  geom_hex(color = NA, binwidth = c(0.1,50)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous(limits = c(-1505,305)) +
  scale_fill_viridis_c(trans = "log10")

ggarrange(p1, p2, common.legend = TRUE)

# Colored by sampling year
p1 <- lyr_data %>% 
  mutate(sampl_yr = cut(lyr_obs_date_y,
                        breaks = c(1899,1960,1984,1995,1999,2009,2012,2018))) %>% 
  ggplot(aes(y = depth, x = lyr_14c,
             fill = sampl_yr)) + 
  geom_point(aes(group = entry_name),
             size = 5, alpha = 0.8, shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_reverse("Depth [cm]", expand = c(0.01,0.01)) +
  scale_x_continuous(limits = c(-1505,305)) +
  scale_fill_viridis_d()

p2 <- lyr_data %>% 
  mutate(sampl_yr = cut(lyr_obs_date_y,
                        breaks = c(1899,1960,1984,1995,1999,2009,2012,2022))) %>% 
  ggplot(aes(y = depth, x = lyr_dd14c,
             fill = sampl_yr)) + 
  geom_point(aes(group = entry_name),
             size = 5, alpha = 0.8, shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_reverse("Depth [cm]", expand = c(0.01,0.01)) +
  scale_x_continuous(limits = c(-1505,305)) +
  scale_fill_viridis_d()

ggarrange(p1, p2, common.legend = TRUE)

lyr_data %>% 
  ggplot(aes(x = depth, y = CORG, z = lyr_14c)) + 
  stat_summary_hex(color = NA, binwidth = c(5,0.1),
                   fun = ~median(.x)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("SOC [wt-%]", trans = "log10") +
  scale_fill_viridis_c("Delta14C", limits = c(-1005,250),
                       option = "A") +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2))

