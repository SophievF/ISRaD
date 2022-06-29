# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 01/06/2022 #

## Depth corrected values ##

# library(ISRaD)
library(tidyverse)
library(ggpubr)
library(mpspline2)

#Load filtered lyr data
lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-06-29"))

lyr_data %>% 
  count(entry_name)

names(lyr_data)

lyr_mpspline <- lyr_data %>% 
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  filter(min(lyr_top) == 0) %>% 
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup()

summary(lyr_mpspline$CORG)
summary(lyr_mpspline$lyr_14c)

lyr_mpspline %>% 
  filter(CORG > 20) %>% 
  dplyr::select(entry_name, id, CORG, lyr_top, lyr_bot, lyr_name) %>% 
  view()

lyr_mpspline %>% 
  count(entry_name)

lyr_mpspline %>% 
  count(id)

## Apply mpspline function
lyr_data_mpspline_14c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 1)

lyr_data_mpspline_14c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_example_14c <- lyr_data_mpspline_14c %>% 
  map(~.x %>% 
        filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
                     id)))

ggplot() +
  geom_path(data = lyr_example_14c$est_1cm, 
            aes(x = SPLINED_VALUE, y = UD, group = id)) +
  theme_bw() +
  scale_y_reverse("Depth") +
  geom_point(data = lyr_mpspline %>% 
               filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001|Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
                            id)),
             aes(x = lyr_14c, y = depth, color = entry_name), 
             shape = 17, size = 3)

lyr_example_14c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()



lyr_data_mpspline_c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 60, lam = 1)

lyr_data_mpspline_c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_example_c <- lyr_data_mpspline_c %>% 
  map(~.x %>% 
        filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001|Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
                     id)))

ggplot() +
  geom_path(data = lyr_example_c$est_1cm, 
            aes(x = SPLINED_VALUE, y = UD, group = id)) +
  theme_bw() +
  scale_y_reverse("Depth") +
  geom_point(data = lyr_mpspline %>%
               filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001|Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
                            id)),
             aes(x = CORG, y = depth, color = entry_name), shape = 17, size = 3)

lyr_example_c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_example_14c_c <- lyr_example_14c$est_1cm %>% 
  rename(lyr_14c = SPLINED_VALUE) %>% 
  full_join(lyr_example_c$est_1cm %>% 
              rename(CORG = SPLINED_VALUE)) %>% 
  tibble()

plotly::ggplotly(
  ggplot() +
    geom_path(data = lyr_example_14c_c,
              aes(x = CORG, y = lyr_14c, color = id)) +
    theme_classic() +
    scale_x_continuous(trans = "log10") +
    geom_point(data = lyr_mpspline %>%
                 filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001|Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
                              id)),
               aes(x = CORG, y = lyr_14c, color = id), shape = 17, size = 3)
)

lyr_data_mpspline_14c$est_1cm %>% 
  filter(LD != 201) %>% 
  ggplot() +
  geom_line(aes(x = SPLINED_VALUE, y = UD, color = id),
            orientation = "y") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_reverse("Depth") 

lyr_mpspline %>% 
  ggplot(aes(x = depth, y = lyr_14c)) +
  geom_path(aes(group = id), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(log(x)),
              fill = "lightblue") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,205)) +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250))

lyr_data_mpspline_14c$est_1cm %>% 
  filter(LD != 201) %>% 
  ggplot(aes(x = UD, y = SPLINED_VALUE)) +
  geom_line(aes(group = id), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(log(x)),
              fill = "lightblue") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,205)) +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250))

lyr_data_mpspline_c$est_1cm %>% 
  filter(LD != 201) %>% 
  ggplot(aes(x = UD, y = SPLINED_VALUE)) +
  geom_path(aes(group = id), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(log(x)),
              fill = "lightblue") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,205)) 


mspline_14c_c <- lyr_data_mpspline_14c$est_1cm %>% 
  rename(lyr_14c = SPLINED_VALUE) %>% 
  full_join(lyr_data_mpspline_c$est_1cm %>% 
              rename(CORG = SPLINED_VALUE)) %>% 
  filter(LD != 201) %>% 
  tibble()

# plotly::ggplotly()

ggplot() +
  geom_path(data = mspline_14c_c,
            aes(x = CORG, y = lyr_14c, color = id)) +
  # geom_point(data = lyr_mpspline,
  #            aes(x = CORG, y = lyr_14c, color = id), shape = 17, size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(trans = "log10") 


lyr_mpspline %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) +
  geom_path(aes(color = id)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250))

mspline_14c_c %>%
  extract(col = id, into = c("entry_name"), regex = "(.*?_.*?)_",
          remove = FALSE) %>%
  ggplot() +
  geom_path(aes(x = CORG, y = lyr_14c, color = id)) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~entry_name)
  

## Cluster analysis
#https://ncss-tech.github.io/AQP/aqp/aqp-profile-dissimilarity.html

library(aqp)
library(cluster)
library(ape)
library(RColorBrewer)
library(latticeExtra)
library(plotrix)


# lyr_example_cluster <- lyr_mpspline %>% 
#   filter(grepl("Baisden_2007|Lawrence_2021", id)) %>% 
#   dplyr::select(id, lyr_top, lyr_bot, lyr_name, CORG, lyr_14c) %>% 
#   mutate(lyr_top = round(lyr_top, digits = 0),
#          lyr_bot = round(lyr_bot, digits = 0))
# 
# summary(lyr_example_cluster)

lyr_cluster <- lyr_mpspline %>% 
  dplyr::filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001|Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
                      id)) %>% 
  group_by(id) %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_name, CORG, lyr_14c) %>% 
  mutate(lyr_top = round(lyr_top, digits = 0),
         lyr_bot = round(lyr_bot, digits = 0)) %>% 
  ungroup()

summary(lyr_cluster)

mspline_cluster <- mspline_14c_c

site_cluster <- lyr_mpspline %>% 
  group_by(id) %>% 
  # filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001|Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
  #              id)) %>%
  dplyr::select(id, pro_BIO12_mmyr_WC2.1) %>%
  summarise(MAP = mean(pro_BIO12_mmyr_WC2.1)) %>% 
  ungroup()
  
# upgrade to SoilProfile Collection object

depths(lyr_cluster) <- id ~ lyr_top + lyr_bot
site(lyr_cluster) <- site_cluster

# xyplot(CORG ~ lyr_14c, groups = id, data = horizons(lyr_example_cluster),
#        auto.key = list(columns = 3, points = TRUE, lines = FALSE))

# compute between-profile dissimilarity, no depth weighting
lyr_mpspline %>% 
  filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
               id)) %>% 
  group_by(id) %>% 
  filter(lyr_bot == max(lyr_bot)) %>% 
  dplyr::select(lyr_bot) %>% 
  # view() %>% 
  summary()

d.dis <- profile_compare(lyr_cluster, 
                         vars = c("lyr_14c", "CORG"), k = 0, 
                         max_d = 100, replace_na = TRUE, add_soil_flag = TRUE)

# check total, between-profile dissimilarity, normalized to maximum
d.m <- signif(as.matrix(d.dis / max(d.dis)), 2)
print(d.m)

# group via divisive hierarchical clustering
d.diana <- diana(d.dis)

str(d.diana)
d.diana$order
d.diana$height
d.diana$dc
d.diana$diss
d.diana$order.lab
d.diana$merge

# convert classes, for better plotting
d.hclust <- as.hclust(d.diana)

d.phylo <- as.phylo(d.hclust)

plot(d.phylo, direction = "down", adj  = 0.1, srt = 0, label.offset = 0.5,
     font = 1, y.lim = c(-150, 25))

# Grouping
d.clust <- cutree(d.hclust, k = 5)

data_grp <- lyr_mpspline %>% 
  dplyr::filter(grepl("Basile_Doelsch_2005|Grant_2022|Schuur_2001|Butman_2007|Chiti_2010|Neue_1980|Castanha_2012|Koarashi_2005", 
                      id)) %>% 
  left_join(cbind.data.frame(data.frame(id = names(d.clust), group = d.clust)))

ggplot(data_grp, aes(x = CORG, y = lyr_14c, color = id)) + 
  # geom_point() + 
  geom_path() + 
  facet_wrap(~group) + 
  scale_x_continuous(trans = "log10") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

factoextra::

# return dissimilarity matrices at each depth slice
d.dis.all <- profile_compare(lyr_example_cluster, vars = c("lyr_14c", "CORG"),  
                             k = 0, max_d = 200, replace_na = TRUE, 
                             add_soil_flag = TRUE, return_depth_distances = TRUE)

# check between-profile dissimilarity, at slice 1
print(as.matrix(d.dis.all[[1]]))


