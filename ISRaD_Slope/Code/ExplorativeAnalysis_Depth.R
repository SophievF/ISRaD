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
lyr_all <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-08-12"))

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
  ungroup()

summary(lyr_mpspline$CORG)
summary(lyr_mpspline$lyr_14c)

# lyr_mpspline %>% 
#   filter(CORG > 20) %>% 
#   dplyr::select(entry_name, id, CORG, lyr_top, lyr_bot, lyr_name) %>% 
#   view()

lyr_mpspline %>% 
  count(entry_name)

lyr_mpspline %>% 
  count(id)

### Apply mpspline function

## mspline 14C
lyr_data_mpspline_14c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, lyr_14c) %>% 
  mpspline_tidy(vlow = -1000, lam = 1)

lyr_data_mpspline_14c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_data_mpspline_14c$est_1cm %>% 
  filter(LD < 101) %>% 
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
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,105)) +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250))

lyr_data_mpspline_14c$est_1cm %>% 
  filter(LD < 101) %>% 
  ggplot(aes(x = UD, y = SPLINED_VALUE)) +
  geom_line(aes(group = id), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(log(x)),
              fill = "lightblue") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,105)) +
  scale_y_continuous("Delta14C", expand = c(0,0), limits = c(-1000,350),
                     breaks = seq(-1000,250,250))
ggsave(file = paste0("./Figure/ISRaD_14C_depth_mspline_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

#Example
ggplot() +
  geom_line(data = lyr_data_mpspline_14c$est_1cm %>% 
              filter(LD < 101) %>% 
              filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
            aes(y = UD, x = SPLINED_VALUE, group = id, color = id), orientation = "y",
            size = 1.5) +
  geom_pointrange(data = lyr_mpspline %>% 
                    filter(lyr_bot < 101) %>% 
                    filter(grepl("Baisden_2007_China hat|Lawrence_2021_Mattole_MT3", id)),
                  aes(x = lyr_14c, y = depth, ymin = lyr_top, ymax = lyr_bot,
                      color = id), size = 1) +
  theme_classic(base_size = 50) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_y_reverse("Depth [cm]", expand = c(0,0), limits = c(105,0)) +
  scale_x_continuous(expression(paste(Delta^14,"C [‰]")), 
                     labels = c(-1000, "", -500, "", 0, ""),
                     expand = c(0,0), limits = c(-1000,250),
                     position = "top")
ggsave(file = paste0("./Figure/ISRaD_14C_depth_example_", Sys.Date(),
                     ".jpeg"), width = 9, height = 10)



lyr_data_mpspline_14c$est_1cm %>% 
  filter(LD < 101) %>% 
  group_by(UD) %>% 
  mutate(mean = mean(SPLINED_VALUE),
         median = median(SPLINED_VALUE),
         median_pseudo = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$estimate,
         sd = sd(SPLINED_VALUE),
         mad = mad(SPLINED_VALUE),
         # lci_mean = t.test(SPLINED_VALUE, conf.level = 0.95)$conf.int[1],
         # uci_mean = t.test(SPLINED_VALUE, conf.level = 0.95)$conf.int[2],
         lci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  ungroup() %>% 
  ggplot(aes(y = UD)) +
  # geom_hline(yintercept = c(70, 100, 180), linetype = "dashed", size = 0.5,
  #            color = "darkgrey") +
  #Alternative: use LCI and UCI from wilcox test and pseudo-median
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median), 
              fill = "#fee0d2", color = "black") +
  geom_line(aes(x = median_pseudo, color = "median + 95 CI", 
                linetype = "median + 95 CI"), size = 0.7, orientation = "y") +
  # geom_line(aes(x = mean, color = "mean + sd", linetype = "mean + sd"), 
  #           size = 0.7, orientation = "y") +
  # #Alternative use LCI and UCI from t test
  # geom_ribbon(aes(xmin = mean-sd, xmax = mean+sd), alpha = 0.2,
  #             fill = "blue", color = "black") +
  geom_line(aes(x = n, color = "n", linetype = "n"),
            orientation = "y") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.9),
        legend.background = element_blank()) +
  scale_x_continuous("d14C", limits = c(-1000,505), expand = c(0,0),
                     position = "top", breaks = seq(-1000,500,250)) +
  scale_y_reverse("Depth [cm]", limits = c(105,0), expand = c(0,0),
                  breaks = seq(0,200,50)) +
  scale_color_manual(name = "", values = c("median + 95 CI" = "red",
                                           # "mean + sd" = "blue",
                                           "n" = "black")) +
  scale_linetype_manual(name = "", values = c("median + 95 CI" = "solid",
                                              # "mean + sd" = "solid",
                                              "n" = "solid")) 
ggsave(file = paste0("./Figure/ISRaD_14C_depth_mspline_sum_1_CI_", Sys.Date(),
                     ".jpeg"), width = 7, height = 6)

# lyr_data_mpspline_14c$est_1cm %>% 
#   tibble() %>% 
#   dplyr::filter(LD < 101) %>% 
#   dplyr::left_join(lyr_mpspline %>% 
#                      distinct(id,.keep_all = TRUE) %>% 
#                      dplyr::select(entry_name, id, lyr_obs_date_y), 
#                    by = "id") %>% 
#   mutate(sample_yr = cut(lyr_obs_date_y, dig.lab = 4,
#                         breaks = c(1899,1960,1984,1999,2022))) %>% 
#   group_by(UD, sample_yr) %>% 
#   mutate(median_pseudo = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$estimate,
#          lci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
#          uci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
#          n = n()) %>% 
#   ungroup() %>% 
#   ggplot(aes(y = UD)) +
#   geom_ribbon(aes(xmin = lci_median, xmax = uci_median, color = sample_yr,
#                   fill = sample_yr), alpha = 0.7) +
#   geom_line(aes(x = median_pseudo, color = sample_yr), 
#             size = 0.7, orientation = "y") +
#   geom_line(aes(x = n, color = sample_yr), linetype = "dashed", orientation = "y") +
#   theme_classic(base_size = 16) +
#   theme(axis.text = element_text(color = "black"),
#         panel.grid.major = element_line(color = "grey", linetype = "dotted",
#                                         size = 0.3),
#         panel.grid.minor = element_line(color = "grey", linetype = "dotted",
#                                         size = 0.2),
#         legend.position = c(0.2,0.8),
#         legend.background = element_blank()) +
#   scale_x_continuous("d14C", limits = c(-1000,375), expand = c(0,0),
#                      position = "top", breaks = seq(-1000,250,250)) +
#   scale_y_reverse("Depth [cm]", limits = c(105,0), expand = c(0,0),
#                   breaks = seq(0,200,50)) +
#   scale_color_viridis_d("Sampling year") +
#   scale_fill_viridis_d("Sampling year")
# ggsave(file = paste0("./Figure/ISRaD_14C_depth_mspline_sum_1_sampling_y_", Sys.Date(),
#                      ".jpeg"), width = 7, height = 6)
# 
# lyr_data_mpspline_14c$est_1cm %>% 
#   tibble() %>% 
#   dplyr::filter(LD != 201) %>% 
#   dplyr::left_join(lyr_mpspline %>% 
#                      distinct(id,.keep_all = TRUE) %>% 
#                      dplyr::select(entry_name, site_name, id, lyr_obs_date_y), 
#                    by = "id") %>% 
#   mutate(sample_yr = cut(lyr_obs_date_y, dig.lab = 4,
#                          breaks = c(1899,1960,1984,1999,2022))) %>% 
#   filter(sample_yr == "(1899,1960]") %>% 
#   count(entry_name, site_name, lyr_obs_date_y)

lyr_data_mpspline_14c$est_1cm %>% 
  tibble() %>% 
  dplyr::filter(LD < 101) %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order,
                                   pro_wrb_soil_order), 
                   by = "id") %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  mutate(#median_pseudo = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$estimate,
         #lci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         #uci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median = median(SPLINED_VALUE),
         mad = mad(SPLINED_VALUE),
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup() %>%
  filter(n_site > 4) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median), color = "red", size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = median-mad, xmax = median+mad) ,
                  color = "black", alpha = 0.2, fill = "red") +
  facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.9),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous("d14C", expand = c(0,0),
                     position = "top", breaks = seq(-1000,500,250)) +
  scale_y_reverse("Depth [cm]", limits = c(105,0), expand = c(0,0),
                  breaks = seq(0,200,50)) +
  coord_cartesian(xlim = c(-1000,155))
ggsave(file = paste0("./Figure/ISRaD_14C_depth_mspline_sum_soilt_1_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

lyr_data_mpspline_14c$est_1cm %>% 
  tibble() %>% 
  dplyr::filter(LD < 101) %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order), 
                   by = "id") %>% 
  drop_na(pro_usda_soil_order) %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  #group_by(pro_usda_soil_order, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup() %>%
  filter(n_site > 4) %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  summarise(n = n_distinct(id)) %>% view()

lyr_data_mpspline_14c$est_1cm %>% 
  tibble() %>% 
  dplyr::filter(LD != 101) %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order, 
                                   pro_country), 
                   by = "id") %>% 
  drop_na(pro_usda_soil_order) %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  mutate(n = n(),
         n_site = n_distinct(id)) %>% 
  filter(pro_usda_soil_order == "Oxisols") %>% 
  ungroup() %>%
  group_by(entry_name, pro_country) %>% 
  summarise(n = n_distinct(id))

lyr_data_mpspline_14c$est_1cm %>% 
  tibble() %>% 
  dplyr::filter(LD != 101) %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id, .keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, pro_KG_present_long), 
                   by = "id") %>% 
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(n = n()) %>% 
  count(ClimateZone) %>% 
  filter(UD == 1|
           UD == 101)

lyr_data_mpspline_14c$est_1cm %>% 
  tibble() %>% 
  dplyr::filter(LD != 101) %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id, .keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, pro_KG_present_long), 
                   by = "id") %>% 
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  filter(ClimateZone == "arid") %>% 
  group_by(entry_name) %>% 
  count(pro_KG_present_long)

lyr_data_mpspline_14c$est_1cm %>% 
  tibble() %>% 
  dplyr::filter(LD != 101) %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id, .keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, pro_KG_present_long), 
                   by = "id") %>% 
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(median_pseudo = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  ungroup() %>%
  filter(n > 4) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_pseudo, color = ClimateZone), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median, fill = ClimateZone, 
                  color = ClimateZone), alpha = 0.2) +
  geom_line(aes(x = n, color = ClimateZone), linetype = "dashed",  
            orientation = "y") +
  # facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous("d14C", expand = c(0,0),
                     position = "top", breaks = seq(-1000,500,250)) +
  scale_y_reverse("Depth [cm]", limits = c(100,0), expand = c(0,0),
                  breaks = seq(0,100,25)) +
  coord_cartesian(xlim = c(-1000,255))
ggsave(file = paste0("./Figure/ISRaD_14C_depth_mspline_sum_climate_1m_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

## Extract dominant clay type
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis

nc_clay <- nc_open("D:/Sophie/PhD/AfSIS_GlobalData/GlobalClay/Ito-Wagai_2016_Global_Clay/clay_fraction_5m_v1r1.nc4")

{
  sink("clay_fraction_5m_v1r1_metadata.txt")
  print(nc_clay)
  sink()
}

lon <- ncvar_get(nc_clay, "lon")
lat <- ncvar_get(nc_clay, "lat", verbose = F)
t <- ncvar_get(nc_clay, "time")

ndvi.array.top <- ncvar_get(nc_clay, "dominant clay mineral (topsoil)")
dim(ndvi.array.top) 
ndvi.array.bot <- ncvar_get(nc_clay, "dominant clay mineral (subsoil)")

nc_close(nc_clay) 

#values: 1 - kaolinite, 2 - smectite, 3 - vermiculite, 4 - illite/mica
r.top <- raster(t(ndvi.array.top), xmn=min(lon), xmx=max(lon), ymn=min(lat), 
            ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

plot(r.top)

clay_type_top <- raster::extract(r.top, cbind(lyr_mpspline$pro_long,
                                              lyr_mpspline$pro_lat))

#values: 1 - kaolinite, 2 - smectite, 3 - vermiculite, 4 - illite/mica
r.bot <- raster(t(ndvi.array.bot), xmn=min(lon), xmx=max(lon), ymn=min(lat), 
                ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

plot(r.bot)

clay_type_bot <- raster::extract(r.bot, cbind(lyr_mpspline$pro_long,
                                              lyr_mpspline$pro_lat))

lyr_mpspline_clay <- cbind(lyr_mpspline, clay_type_top, clay_type_bot) %>% 
  tibble() %>% 
  mutate(clay_type_top = case_when(
    clay_type_top == 1 ~ "low-activity clays",
    clay_type_top == 2 ~ "high-acivity clays",
    clay_type_top == 3 ~ "high-acivity clays",
    clay_type_top == 4 ~ "high-acivity clays"
  )) %>% 
  mutate(clay_type_bot = case_when(
    clay_type_bot == 1 ~ "low-activity clays",
    clay_type_bot == 2 ~ "high-acivity clays",
    clay_type_bot == 3 ~ "high-acivity clays",
    clay_type_bot == 4 ~ "high-acivity clays"
  )) %>%
  mutate(clay_type = case_when(
    lyr_bot <= 30 ~ clay_type_top,
    lyr_bot > 30 ~ clay_type_bot
  )) 

lyr_data_mpspline_14c$est_1cm %>% 
  tibble() %>% 
  dplyr::filter(LD != 201) %>% 
  dplyr::left_join(lyr_mpspline_clay %>% 
                     distinct(id, .keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, clay_type_bot, pro_KG_present_long), 
                   by = "id") %>%
  drop_na(clay_type_bot) %>% 
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  group_by(ClimateZone, clay_type_bot, UD) %>% 
  mutate(median_pseudo = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  ungroup() %>%
  filter(n > 4) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_pseudo, color = clay_type_bot), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median, fill = clay_type_bot, 
                  color = clay_type_bot), alpha = 0.2) +
  geom_line(aes(x = n, color = clay_type_bot), linetype = "dashed",  
            orientation = "y") +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.1,0.9),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous("d14C", expand = c(0,0),
                     position = "top", breaks = seq(-1000,500,250)) +
  scale_y_reverse("Depth [cm]", limits = c(205,0), expand = c(0,0),
                  breaks = seq(0,200,50)) +
  coord_cartesian(xlim = c(-1000,255))
ggsave(file = paste0("./Figure/ISRaD_14C_depth_mspline_sum_climate_clay_type_bot_1_", 
                     Sys.Date(), ".jpeg"), width = 11, height = 6)

lyr_data_mpspline_14c$est_1cm %>% 
  tibble() %>% 
  dplyr::filter(LD != 201) %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order, 
                                   pro_KG_present_long), 
                   by = "id") %>% 
  drop_na(pro_usda_soil_order) %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  #reclassify soil type Schuur_2001: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Schuur_2001",
                                       "Andisols")) %>%
  #reclassify soil type Guillet_1988: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Guillet_1988",
                                       "Andisols")) %>%
  #reclassify soil type Torn_1997: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Torn_1997",
                                       "Andisols")) %>%
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  mutate(clay_type = case_when(
    pro_usda_soil_order == "Oxisols" ~ "low-activity clays",
    pro_usda_soil_order == "Ultisols" ~ "low-activity clays",
    TRUE ~ "high-activity clays"
  )) %>% 
  group_by(ClimateZone, clay_type, UD) %>% 
  mutate(median_pseudo = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  ungroup() %>%
  filter(n > 4) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_pseudo, color = clay_type), size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median, fill = clay_type, 
                  color = clay_type), alpha = 0.2) +
  geom_line(aes(x = n, color = clay_type), linetype = "dashed",  
            orientation = "y") +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.1,0.9),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous("d14C", expand = c(0,0),
                     position = "top", breaks = seq(-1000,500,250)) +
  scale_y_reverse("Depth [cm]", limits = c(205,0), expand = c(0,0),
                  breaks = seq(0,200,50)) +
  coord_cartesian(xlim = c(-1000,255))
ggsave(file = paste0("./Figure/ISRaD_14C_depth_mspline_sum_climate_clay_type_1_", 
                     Sys.Date(), ".jpeg"), width = 11, height = 6)


## mspline CORG

lyr_data_mpspline_c <- lyr_mpspline %>% 
  dplyr::select(id, lyr_top, lyr_bot, CORG) %>% 
  mpspline_tidy(vlow = 0.01, vhigh = 60, lam = 1)

lyr_data_mpspline_c$tmse %>% 
  filter(ERROR_TYPE == "RMSE") %>% 
  summary()

lyr_data_mpspline_c$est_1cm %>% 
  filter(LD != 101) %>% 
  ggplot(aes(x = UD, y = SPLINED_VALUE)) +
  geom_path(aes(group = id), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(log(x)),
              fill = "lightblue") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]", expand = c(0,0), limits = c(0,105))

lyr_data_mpspline_c$est_1cm %>%
  tibble() %>% 
  dplyr::filter(LD < 101) %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order), 
                   by = "id") %>% 
  drop_na(pro_usda_soil_order) %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  mutate(median_pseudo = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_median = wilcox.test(SPLINED_VALUE, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup() %>%
  filter(n_site > 4) %>% 
  ggplot(aes(y = UD)) +
  geom_line(aes(x = median_pseudo), color = "red", size = 0.7, orientation = "y") +
  geom_ribbon(aes(xmin = lci_median, xmax = uci_median) ,
              color = "black", fill = "red", alpha = 0.2) +
  facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.9),
        legend.background = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
  scale_x_continuous("SOC", expand = c(0,0),
                     position = "top") +
  scale_y_reverse("Depth [cm]", limits = c(105,0), expand = c(0,0),
                  breaks = seq(0,200,50))
ggsave(file = paste0("./Figure/ISRaD_SOC_depth_mspline_sum_soilt_1_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

## Mspline 14C and CORG

mspline_14c_c <- lyr_data_mpspline_14c$est_1cm %>% 
  rename(lyr_14c = SPLINED_VALUE) %>% 
  full_join(lyr_data_mpspline_c$est_1cm %>% 
              rename(CORG = SPLINED_VALUE)) %>% 
  filter(LD < 101) %>% 
  tibble()


mspline_14c_c_soilt <- mspline_14c_c %>%
  tibble() %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order), 
                   by = "id") %>% 
  drop_na(pro_usda_soil_order) %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup() %>%
  filter(n_site > 4) 

p1 <- mspline_14c_c_soilt %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrange(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c), color = "#fee0d2") +
  facet_wrap(~pro_usda_soil_order) +
  theme_classic(base_size = 16) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]"))) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank())

p1 +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c), color = "#fee0d2") +
  geom_path()
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_sum_soilt_1m_", Sys.Date(),
                     ".jpeg"), width = 23, height = 15)

mspline_14c_c_soilt %>% 
  group_by(pro_usda_soil_order, UD) %>% 
  summarise(n = n_distinct(id)) %>% 
  filter(UD == 1| UD == 97) %>% view()

mspline_14c_c_climate <- mspline_14c_c %>%
  tibble() %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id, .keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_KG_present_long), 
                   by = "id") %>% 
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  group_by(ClimateZone, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup() %>%
  filter(n_site > 4)

p1 <- mspline_14c_c_climate %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  arrnge(UD) %>% 
  ggplot(aes(x = median_c, y = median_14c, color = ClimateZone)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c), alpha = 0.1) +
  # facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous("Delta14C") +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        legend.position = c(0.1,0.85))

p1 +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c), alpha = 0.1) +
  geom_path()
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_sum_climate_1m_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

mspline_14c_c_climate_clay <- mspline_14c_c %>%
  tibble() %>% 
  dplyr::left_join(lyr_mpspline %>% 
                     distinct(id,.keep_all = TRUE) %>% 
                     dplyr::select(entry_name, id, site_name, pro_usda_soil_order, 
                                   pro_KG_present_long), 
                   by = "id") %>% 
  drop_na(pro_usda_soil_order) %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  #reclassify soil type Schuur_2001: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Schuur_2001",
                                       "Andisols")) %>%
  #reclassify soil type Guillet_1988: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Guillet_1988",
                                       "Andisols")) %>%
  #reclassify soil type Torn_1997: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Torn_1997",
                                       "Andisols")) %>%
  mutate(ClimateZone = case_when(
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Polar") ~ "cold/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  mutate(clay_type = case_when(
    pro_usda_soil_order == "Oxisols" ~ "low-activity clays",
    pro_usda_soil_order == "Ultisols" ~ "low-activity clays",
    TRUE ~ "high-activity clays"
  )) %>% 
  group_by(ClimateZone, clay_type, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n(),
         n_site = n_distinct(site_name)) %>% 
  ungroup() %>%
  filter(n_site > 4)

p2 <- mspline_14c_c_climate_clay %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>%
  distinct(median_14c, .keep_all = TRUE) %>%
  ggplot(aes(x = median_c, y = median_14c, color = clay_type)) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c), alpha = 0.1) +
  facet_wrap(~ClimateZone) +
  theme_classic(base_size = 16) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_y_continuous("Delta14C") +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.background = element_blank(),
        legend.position = c(0.9,0.1))

p2 +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c), alpha = 0.1) +
  geom_path()
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_sum_climate_clay_1m_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

#### NOT IDEAL YET ####
mspline_14c_c %>% 
  group_by(UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$estimate,
         # lci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         # uci_c = wilcox.test(CORG, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  ungroup() %>% 
  dplyr::select(-c(id, lyr_14c, CORG)) %>% 
  distinct(median_14c, .keep_all = TRUE) %>% 
  ggplot(aes(x = median_c)) +
  geom_vline(xintercept = 0.511, linetype = "dashed", size = 0.5,
             color = "darkgrey") +
  geom_ribbon(aes(ymin = lci_14c, ymax = uci_14c), orientation = "x",
              fill = "#fee0d2", color = "black") +
  # geom_ribbon(aes(xmin = lci_c, xmax = uci_c),
  #             fill = "#fee0d2", color = "black") +
  geom_path(aes(y = median_14c), color = "red", size = 0.7) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dotted",
                                        size = 0.3),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted",
                                        size = 0.2),
        legend.position = c(0.2,0.1),
        legend.background = element_blank()) +
  scale_y_continuous("d14C", limits = c(-1000,125), expand = c(0,0),
                     breaks = seq(-1000,250,250)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") 
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_sum_1_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)

plotly::ggplotly(
  mspline_14c_c %>% 
    ggplot(aes(x = CORG, y = lyr_14c)) +
    geom_path(aes(group = id), alpha = 0.5) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_x_continuous(trans = "log10") 
)
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)


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

mspline_14c_c %>% 
  rename(lyr_14c_msp = lyr_14c,
         CORG_msp = CORG) %>% 
  left_join(lyr_mpspline, by = "id") %>% 
  drop_na(pro_usda_soil_order) %>% 
  filter(pro_usda_soil_order != "Aridisols") %>% 
  #reclassify soil type Schuur_2001: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Schuur_2001",
                                       "Andisols")) %>%
  #reclassify soil type Guillet_1988: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Guillet_1988",
                                       "Andisols")) %>%
  #reclassify soil type Torn_1997: all Andisols
  mutate(pro_usda_soil_order = replace(pro_usda_soil_order,
                                       entry_name == "Torn_1997",
                                       "Andisols")) %>%
  ggplot(aes(x = CORG_msp, y = lyr_14c_msp)) +
  geom_path(aes(group = id, color = pro_BIO12_mmyr_WC2.1), size = 1) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8,0.1)) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~pro_usda_soil_order) +
  scale_color_viridis_c("MAP [mm]", trans = "log10", direction = -1, 
                        limits = c(100,3000)) +
  guides(color = guide_colorbar(barwidth = 10, frame.colour = "black", 
                                ticks.linewidth = 2, direction = "horizontal",
                                title.position = "top"))
ggsave(file = paste0("./Figure/ISRaD_14C_SOC_mspline_soil_", Sys.Date(),
                     ".jpeg"), width = 11, height = 6)


# Cluster analysis
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


