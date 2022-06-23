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
lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-06-23"))

lyr_data %>% 
  count(entry_name)

names(lyr_data)

## Apply mpspline function
lyr_data_mpspline_14c <- lyr_data %>% 
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
  ungroup() %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 1)

lyr_data_mpspline_14c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_example_14c <- lyr_data_mpspline_14c %>% 
  map(~.x %>% 
        filter(grepl("Baisden_2002|Lawrence_2021", id)))

ggplot() +
  geom_path(data = lyr_example_14c$est_1cm, 
            aes(x = SPLINED_VALUE, y = UD, color = id)) +
  theme_bw() +
  scale_y_reverse("Depth") +
  geom_point(data = lyr_data %>% 
               filter(lyr_bot <= 200) %>% 
               group_by(id) %>%
               #Filter for studies that have more than 2 depth layers
               filter(n() > 2) %>%
               filter(grepl("Baisden_2002|Lawrence_2021", id)),
             aes(x = lyr_14c, y = depth, color = id), shape = 17, size = 3)

lyr_example_14c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_data_mpspline_c <- lyr_data %>%
  filter(entry_name != "Drake_2019" ,
         entry_name != "Richer_1999",
         entry_name != "Giardina_2014") %>% 
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 20, lam = 1)

lyr_data_mpspline_c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_example_c <- lyr_data_mpspline_c %>% 
  map(~.x %>% 
        filter(grepl("Baisden_2002|Lawrence_2021", id)))

ggplot() +
  geom_path(data = lyr_example_c$est_1cm, 
            aes(x = SPLINED_VALUE, y = UD, color = id)) +
  theme_bw() +
  scale_y_reverse("Depth") +
  geom_point(data = lyr_data %>% 
               filter(lyr_bot <= 200) %>% 
               group_by(id) %>%
               #Filter for studies that have more than 2 depth layers
               filter(n() > 2) %>%
               filter(grepl("Baisden_2002|Lawrence_2021", id)),
             aes(x = CORG, y = depth, color = id), shape = 17, size = 3)

lyr_example_c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_example_14c_c <- lyr_example_14c$est_1cm %>% 
  rename(lyr_14c = SPLINED_VALUE) %>% 
  full_join(lyr_example_c$est_1cm %>% 
              rename(CORG = SPLINED_VALUE)) %>% 
  tibble()

ggplot() +
  geom_path(data = lyr_example_14c_c,
            aes(x = CORG, y = lyr_14c, color = id)) +
  theme_classic() +
  scale_x_continuous() +
  geom_point(data = lyr_data %>% 
               filter(lyr_bot <= 200) %>% 
               group_by(id) %>%
               #Filter for studies that have more than 2 depth layers
               filter(n() > 2) %>%
               filter(grepl("Baisden_2002|Lawrence_2021", id)),
             aes(x = CORG, y = lyr_14c, color = id), shape = 17, size = 3)

lyr_data_mpspline_14c$est_1cm %>% 
  filter(LD != 201) %>% 
  ggplot() +
  geom_line(aes(x = SPLINED_VALUE, y = UD, color = id),
            orientation = "y") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_reverse("Depth") 

lyr_data %>% 
  filter(entry_name != "Drake_2019" ,
         entry_name != "Richer_1999",
         entry_name != "Giardina_2014") %>%
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup() %>% 
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

ggplot() +
  geom_path(data = mspline_14c_c,
            aes(x = CORG, y = lyr_14c, color = id)) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous() +
  geom_point(data = lyr_data %>% 
               filter(entry_name != "Drake_2019" ,
                      entry_name != "Richer_1999",
                      entry_name != "Giardina_2014") %>%
               filter(lyr_bot <= 200) %>% 
               group_by(id) %>%
               #Filter for studies that have more than 2 depth layers
               filter(n() > 2),
             aes(x = CORG, y = lyr_14c, color = id), shape = 17, size = 3)

mspline_14c_c %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) + 
  geom_path(aes(group = id), alpha = 0.5) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250))

lyr_data %>% 
  filter(entry_name != "Drake_2019" ,
         entry_name != "Richer_1999",
         entry_name != "Giardina_2014") %>% 
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  ungroup() %>% 
  ggplot(aes(x = CORG, y = lyr_14c)) +
  geom_path(aes(group = id), alpha = 0.5) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250))

## Cluster analysis
#https://ncss-tech.github.io/AQP/aqp/aqp-profile-dissimilarity.html

library(cluster)
library(ape)

mpspline_data <- lyr_data_mpspline$est_dcm[,c(3,4)]

d <- daisy(mpspline_data)

print(d)

d.diana <- diana(d)

plot(d.diana)

d.phylo <- as.phylo(as.hclust(d.diana))

par(mfcol=c(1,2), mar=c(4.5,4,1,1))

plot(SPLINED_VALUE ~ LD, data = mpspline_data, 
     type = "p", xlab = "Depth", ylab = "14C")
grid()
# text(mpspline_data$LD, mpspline_data$SPLINED_VALUE, row.names(mpspline_data), 
#     font = 2)

plot(d.phylo, font = 2, label.offset = 0.5, adj = 0.5, direction = "down", 
     srt = 90, y.lim = c(-1, 15))
