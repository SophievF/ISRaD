# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 01/06/2022 #

## Depth corrected values ##

library(ISRaD)
library(tidyverse)
library(ggpubr)
library(mpspline2)

#Load filtered lyr data
lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_", Sys.Date()))

lyr_data %>% 
  count(entry_name)

names(lyr_data)

## Apply mpspline function
lyr_data_mpspline <- lyr_data %>% 
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  ungroup() %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000)

head(lyr_data_mpspline$est_icm)
head(lyr_data_mpspline$est_1cm)
head(lyr_data_mpspline$est_dcm)
head(lyr_data_mpspline$tmse)

str(lyr_data_mpspline$est_1cm)
lyr_data_mpspline$est_1cm


ggplot() +
  geom_line(data = lyr_data_mpspline$est_1cm, aes(x = SPLINED_VALUE, y = UD, color = id),
            orientation = "y") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_reverse("Depth") 

lyr_data_mpspline$est_1cm %>% 
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
