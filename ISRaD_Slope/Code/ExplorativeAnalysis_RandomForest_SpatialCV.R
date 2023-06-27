# Data analysis: Profile data ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 08/03/2023 #

library(tidyverse)
library(ggpubr)
library(scales)
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(iml)
library(sf)
library(mlr3spatiotempcv) 
library(mlr3spatial)

##Spatial cross-validation ##
# Load filtered and splined data
lyr_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-03-09.csv")

# Check data
lyr_data %>% 
  dplyr::summarise(n_studies = n_distinct(entry_name),
                   n_sites = n_distinct(site_name),
                   n_profiles = n_distinct(id))

lyr_data %>% 
  dplyr::select(lyr_14c_msp, CORG_msp, UD, pro_MAT_mod, pro_AI, lyr_clay_mod) %>% 
  cor()

# Convert characters into factors
lyr_data$ClimateZone <- factor(lyr_data$ClimateZone,
                               levels = c("tundra/polar", "cold temperate", "arid",
                                          "warm temperate", "tropical"))

lyr_data$ClimateZoneAnd <- factor(lyr_data$ClimateZoneAnd,
                                  levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                             "warm temperate", "arid", "tropical"))

lyr_data %>% 
  group_by(ClimateZoneAnd) %>% 
  dplyr::summarise(n_studies = n_distinct(entry_name),
                   n_sites = n_distinct(site_name),
                   n_profiles = n_distinct(id))



## Random Forest Analysis ##
#Only use a subset of the data for now
set.seed(42)
split_id <- lyr_data %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>% 
  distinct(id) %>% 
  sample_frac(0.3)

split_data <- lyr_data %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  tibble::rowid_to_column() %>% 
  mutate(split = case_when(
    id %in% split_id$id ~ "test",
    !(id %in% split_id$id) ~ "train"
  )) 


lyr_data_sf <- sf::st_as_sf(split_data %>% filter(split == "test"), 
                            coords = c("pro_long", "pro_lat"), crs = 4326)

# Filter data (only keep variables of interest)
rf_data_14c <- lyr_data_sf %>% 
  arrange(id) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

## Set-up random forest
task_rf_14c <- as_task_regr_st(x = rf_data_14c %>% 
                              dplyr::select(-entry_name),
                            target = "lyr_14c_msp")

#Add id as group for CV (same id kept together)
task_rf_14c$set_col_roles("id", roles = "group")
task_rf_14c$groups
print(task_rf_14c)

lrn_rf_14c <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000, predict_type = "response")

#10 fold CV
resampling <- rsmp("cv", folds = 10)
set.seed(42)
resampling$instantiate(task_rf_14c)

autoplot(resampling, task_rf_14c, fold_id = c(1,2,5,10))

rr <- mlr3::resample(task = task_rf_14c, learner = lrn_rf_14c, resampling = resampling)

rr$aggregate(measures = msrs(c("regr.rmse", "regr.mse", "regr.rsq")))
rr_pred <- rr$prediction(predict_sets = "test")

rr_pred %>% 
  ggplot(aes(x = rr_pred$data$response, y = rr_pred$data$truth)) +
  geom_point() +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous(expression(paste("Observed ", Delta^14,"C [‰]")), 
                     limits = c(-1000,350), expand = c(0,0),
                     breaks = seq(-1000,250,250)) +
  scale_x_continuous(expression(paste("Predicted ", Delta^14,"C [‰]")), 
                     limits = c(-1000,350), expand = c(0,0),
                     breaks = seq(-1000,250,250))


