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
lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-06-24"))

lyr_data %>% 
  count(entry_name)

names(lyr_data)

lyr_mpspline <- lyr_data %>% 
  #remove studies that have multiple depth layers for now
  filter(entry_name != "Drake_2019" ,
         entry_name != "Richer_1999",
         entry_name != "Giardina_2014") %>% 
  #overlapping depth layers; not enough depth layers
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup()

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
        filter(grepl("Baisden_2007|Lawrence_2021", id)))

ggplot() +
  geom_path(data = lyr_example_14c$est_1cm, 
            aes(x = SPLINED_VALUE, y = UD, color = id)) +
  theme_bw() +
  scale_y_reverse("Depth") +
  geom_point(data = lyr_mpspline %>% 
               filter(grepl("Baisden_2007|Lawrence_2021", id)),
             aes(x = lyr_14c, y = depth, color = id), shape = 17, size = 3)

lyr_example_14c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_data_mpspline_c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 20, lam = 1)

lyr_data_mpspline_c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_example_c <- lyr_data_mpspline_c %>% 
  map(~.x %>% 
        filter(grepl("Baisden_2007|Lawrence_2021", id)))

ggplot() +
  geom_path(data = lyr_example_c$est_1cm, 
            aes(x = SPLINED_VALUE, y = UD, color = id)) +
  theme_bw() +
  scale_y_reverse("Depth") +
  geom_point(data = lyr_mpspline %>%
               filter(grepl("Baisden_2007|Lawrence_2021", id)),
             aes(x = CORG, y = depth, color = id), shape = 17, size = 3)

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
                 filter(grepl("Baisden_2007|Lawrence_2021", id)),
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

plotly::ggplotly()

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

lyr_example_cluster <- lyr_mpspline %>% 
  filter(grepl("Baisden_2007|Lawrence_2021", id)) %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_name, CORG, lyr_14c) %>% 
  mutate(lyr_top = round(lyr_top, digits = 0),
         lyr_bot = round(lyr_bot, digits = 0))

summary(lyr_example_cluster)

site_example <- lyr_mpspline %>% 
  filter(grepl("Baisden_2007|Lawrence_2021", id)) %>% 
  dplyr::select(id, pro_BIO12_mmyr_WC2.1) %>% 
  group_by(id) %>% 
  summarise(MAP = mean(pro_BIO12_mmyr_WC2.1))
  
# upgrade to SoilProfile Collection object
depths(lyr_example_cluster) <- id ~ lyr_top + lyr_bot
site(lyr_example_cluster) <- site_example

xyplot(CORG ~ lyr_14c, groups = id, data = horizons(lyr_example_cluster),
       auto.key = list(columns = 3, points = TRUE, lines = FALSE))

# compute betwee-profile dissimilarity, no depth weighting
d.dis <- profile_compare(lyr_example_cluster, 
                         vars = c("CORG", "lyr_14c"), k = 0, 
                         max_d = 150, replace_na = TRUE, add_soil_flag = TRUE)

# check total, between-profile dissimilarity, normalized to maximum
d.m <- signif(as.matrix(d.dis / max(d.dis)), 2)
print(d.m)

# group via divisive hierarchical clustering
d.diana <- diana(d.dis)

d.diana$order
d.diana$height
d.diana$dc
d.diana$diss
d.diana$order.lab
d.diana$merge

# convert classes, for better plotting
d.hclust <- as.hclust(d.diana, method = "ward.D2")

d.phylo <- as.phylo(d.hclust)

plot(d.hclust)

plot(d.phylo, direction = "down", adj  = 0.1, srt = 0, label.offset = 0.5,
     font = 1, y.lim = c(-150, 25))

# plot: 2 figures side-by-side
# layout(matrix(c(1,2), nrow = 1), widths = c(0.6, 0.4))
# par(mar = c(1,1,1,1))

# profiles
# plotSPC(lyr_example_cluster, name = "lyr_name", axis.line.offset = -1)

# annotate shallow-mod.deep break
# abline(h = 120, col = 'red', lty = 2)

# add dissimilarity matrix
# addtable2plot(0.8, 70, format(d.m, digits = 2), display.rownames = TRUE, 
#               xjust = 0, yjust = 0, cex = 0.75, title = 'Total Dissimilarity')

# plot dendrogram in next panel
# plot(d.phylo, direction = 'down', adj  = 0.5, srt = 0, 
#      label.offset = 0.5, font = 1, y.lim = c(-5, 25), cex = 0.7)


# return dissimilarity matrices at each depth slice
d.dis.all <- profile_compare(lyr_example_cluster, vars = c("lyr_14c", "CORG"),  
                             k = 0, max_d = 200, replace_na = TRUE, 
                             add_soil_flag = TRUE, return_depth_distances = TRUE)

# check between-profile dissimilarity, at slice 1
print(as.matrix(d.dis.all[[1]]))
