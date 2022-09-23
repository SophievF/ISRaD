# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 15/09/2022 #

library(tidyverse)
library(ggpubr)
library(mpspline2)

#Load filtered lyr data
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-09-22"))

lyr_all %>% 
  count(entry_name)

names(lyr_all)

lyr_mpspline <- lyr_all %>% 
  filter(lyr_obs_date_y > 1959) %>% 
  # filter(lyr_bot <= 100) %>% 
  # filter(lyr_top <= 100) %>%
  group_by(id) %>%
  # filter(min(lyr_top) == 0) %>% 
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(ClimateZone = case_when(
    pro_usda_soil_order == "Andisols" ~ "andisols",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold",
    str_detect(pro_KG_present_long, "Polar") ~ "polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  #remove for now: need to fix depth
  filter(entry_name != "Fernandez_1993a")

summary(lyr_mpspline$CORG)
summary(lyr_mpspline$lyr_14c)

lyr_mpspline %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

### Apply mspline function

## mspline 14C
lyr_data_mpspline_14c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 0.5)

## mspline CORG
lyr_data_mpspline_c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 60, lam = 0.5)

## 14C and SOC
mspline_14c_c <- lyr_data_mpspline_14c$est_1cm %>% 
  rename(lyr_14c_msp = SPLINED_VALUE) %>% 
  full_join(lyr_data_mpspline_c$est_1cm %>% 
              rename(CORG_msp = SPLINED_VALUE)) %>% 
  filter(LD < 101) %>% 
  tibble()

mspline_14c_c_all <- mspline_14c_c %>%
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id, .keep_all = TRUE), 
                   by = "id") %>% 
  group_by(id) %>% 
  arrange(UD) %>% 
  ungroup()

mspline_14c_c_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

head(mspline_14c_c_all)


mspline_14c_c_all %>%
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) 
  
mspline_14c_c_all %>% 
  dplyr::select(id, UD, lyr_14c_msp) %>% 
  pivot_wider(names_from =  id, values_from = lyr_14c_msp) %>% 
  dplyr::select(-UD) %>% 
  data.matrix() %>% 
  dist() 

### Slope analysis
library(lme4)
library(broom)

plotly::ggplotly(
  mspline_14c_c_all %>% 
    group_by(ClimateZone, UD) %>% 
    mutate(n = n()) %>%
    ungroup(UD) %>% 
    mutate(n_rel = n * 100 / max(n)) %>% 
    filter(n > 4 & n_rel > 60) %>% 
    ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) +
    geom_path(aes(group = id, color = ClimateZone)) +
    facet_wrap(~ClimateZone) +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none") +
    scale_y_continuous(expression(paste(Delta^14, "C"))) +
    scale_x_continuous("Soil organic carbon [%]", trans = "log10")
)

model_lm <- mspline_14c_c_all %>% group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>%
  group_by(id) %>% 
  do(fit = tidy(lm(lyr_14c_msp ~ log10(CORG_msp), data = .))) %>% 
  unnest(fit) %>% 
  filter(term == "log10(CORG_msp)") 

sum_data <- mspline_14c_c_all %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>%
  ungroup(UD) %>% 
  mutate(n_rel = n * 100 / max(n)) %>% 
  filter(n > 4 & n_rel > 60) %>%
  dplyr::select(id, lyr_14c_msp, CORG_msp, ClimateZone, pro_usda_soil_order) %>% 
  group_by(id, ClimateZone, pro_usda_soil_order) %>% 
  summarise(max_14c = max(lyr_14c_msp),
            min_14c = min(lyr_14c_msp),
            max_CORG = max(CORG_msp),
            min_CORG = min(CORG_msp))

model_data <- sum_data %>% 
  left_join(model_lm)

# how to interpret log10: https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
plotly::ggplotly(
  model_data %>% 
    ggplot(aes(y = estimate*log(1.01), x = (max_CORG-min_CORG), fill = ClimateZone, group = id)) +
    geom_point(size = 3, shape = 21) +
    facet_wrap(~ClimateZone) +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_y_continuous("Slope: lyr_14c ~ log10(CORG)") +
    scale_x_continuous("CORG: max-min", expand = c(0,0), limits = c(0,45))
)

model_sum <- model_data %>% 
  group_by(ClimateZone) %>% 
  summarise(median_slope = median(estimate*log(1.01)),
            mad_slope = mad(estimate*log(1.01)),
            median_CORG = median((min_CORG-max_CORG)/max_CORG*100),
            mad_CORG = mad((min_CORG-max_CORG)/max_CORG*100))

plotly::ggplotly(
  ggplot() +
    geom_point(data = model_data,
               aes(y = estimate*log(1.01), x = (min_CORG-max_CORG)/max_CORG*100,
                   fill = ClimateZone, group = id),
               size = 3, shape = 21, alpha = 0.3) +
    geom_errorbarh(data = model_sum,
                   aes(y = median_slope,
                       xmin = median_CORG-mad_CORG, xmax = median_CORG+mad_CORG)) +
    geom_errorbar(data = model_sum,
                  aes(x = median_CORG,
                      ymin = median_slope-mad_slope, ymax = median_slope+mad_slope)) +
    geom_point(data = model_sum,
               aes(x = median_CORG, y = median_slope,
                   fill = ClimateZone), size = 3, shape = 21) +
    # facet_wrap(~ClimateZone) +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_y_continuous("Slope: lyr_14c ~ log10(CORG)") +
    scale_x_continuous("rel. change in SOC")
)


### Cluster analysis
## https://ncss-tech.github.io/AQP/aqp/aqp-profile-dissimilarity.html
## https://search.r-project.org/CRAN/refmans/aqp/html/pc.html

library(aqp)
library(cluster)
library(ape)
library(RColorBrewer)
library(latticeExtra)
library(plotrix)

mspline_cluster <- mspline_14c_c_all

depths(mspline_cluster) <- id ~ UD + LD

d.1 <- profile_compare(mspline_cluster, vars = c("lyr_14c_msp", "CORG_msp"),
                       max_d = 100, k = 0, replace_na = TRUE, add_soil_flag = TRUE)
str(d.1)

labels(d.1)

d.m <- signif(as.matrix(d.1 / max(d.1)), 2)
head(d.m[1,])

d.sam <- MASS::sammon(d.1)
head(d.sam$points)

h <- diana(d.1)
p <- as.phylo(as.hclust(h))
plot(p, show.tip.label=FALSE, col=cutree(h, 5))
tiplabels(mspline_14c_c_all$pro_name, col=cutree(h, 5), bg=NA, cex=0.75)

str(d.sam)

dev.off() ; dev.new()
plot(d.sam$points, type = "p", col=cutree(h, 5))
text(d.sam$points, labels=row.names(as.data.frame(d.sam$points)), 
     cex=0.75, col=cutree(h, 5))

## Not ideal yet
# need to assign id first
groups_h <- cutree(h, k = 7)

group_id <- cbind(data.frame(groups_h), data.frame(h$order.lab)) %>% 
  rename(id = h.order.lab)

datagrp <- mspline_14c_c_all %>% 
  dplyr::select(id, UD, LD, lyr_14c_msp, CORG_msp) %>% 
  # We're joining the dataset with the group obtained from clustering
  left_join(group_id) 

ggplot(datagrp, aes(y = UD, x = lyr_14c_msp)) + 
  geom_line(aes(colour = id), orientation = "y") + 
  facet_grid(cols = vars(groups_h)) + 
  scale_y_reverse() +  
  theme(legend.position = "none")

ggplot(datagrp, aes(y = UD, x = CORG_msp)) + 
  geom_line(aes(colour = id), orientation = "y") + 
  facet_grid(cols = vars(group)) + 
  scale_y_reverse() +  
  theme(legend.position = "none")

datagrp %>% 
  arrange(UD) %>% 
  ggplot(aes(y = lyr_14c_msp, x = CORG_msp, group = id)) + 
  geom_path() + 
  facet_grid(cols = vars(group)) + 
  theme(legend.position = "none") +
  scale_x_continuous(trans = "log10")

# Based on code for gaussin process regression
hclust.method <- "ward.D2"
dist.method <- "euclidean"

clust <- d.1 %>% 
  dist(method = dist.method) %>% 
  hclust(method = hclust.method)

as.dendrogram(clust) %>% 
  dendextend::set("branches_k_color", k = 7) %>% 
  tidygraph::as_tbl_graph() %>% 
  ggraph::ggraph() + 
  ggraph::geom_edge_elbow(aes())+ 
  ggraph::geom_edge_elbow(aes(colour = col)) + 
  ggraph::geom_node_text(aes(label = label, filter = leaf), 
                         angle = 90, hjust = 1, size = 3) + 
  ggraph::theme_graph() + 
  scale_y_continuous(limits = c(-8000,NA))


NbClust::NbClust(d.1, distance = dist.method, method = hclust.method, index = "cindex" )$Best.partition
groups <- cutree(clust, 4) 

#NbClust::NbClust(dists, distance = dist.method, method = hclust.method, index = "cindex" )$Best.partition

datagrp <- mspline_14c_c_all %>% 
  dplyr::select(id, UD, LD, lyr_14c_msp, CORG_msp) %>% 
  # We're joining the dataset with the group obtained from clustering
  left_join(cbind.data.frame(data.frame(id = names(groups), group = groups))) 

ggplot(datagrp, aes(y = UD, x = lyr_14c_msp)) + 
    geom_line(aes(colour = id), orientation = "y") + 
    facet_grid(cols = vars(group)) + 
    scale_y_reverse() +  
    theme(legend.position = "none")

ggplot(datagrp, aes(y = UD, x = CORG_msp)) + 
  geom_line(aes(colour = id), orientation = "y") + 
  facet_grid(cols = vars(group)) + 
  scale_y_reverse() +  
  theme(legend.position = "none")

datagrp %>% 
  arrange(UD) %>% 
  ggplot(aes(y = lyr_14c_msp, x = CORG_msp, group = id)) + 
  geom_path() + 
  facet_grid(cols = vars(group)) + 
  theme(legend.position = "none") +
  scale_x_continuous(trans = "log10")
