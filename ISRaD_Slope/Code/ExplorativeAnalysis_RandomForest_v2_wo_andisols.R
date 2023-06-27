# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 25/01/2023 #

#### RANDOM FOREST ####

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mlr3viz)
library(sf)
library(mlr3)
library(mlr3learners)
library(mlr3tuningspaces)
library(mlr3tuning)


### Load filtered and splined lyr data
mspline_14c_c_all <- read_csv("./Data/ISRaD_flat_splined_filled_2023-02-01.csv")

head(mspline_14c_c_all)
names(mspline_14c_c_all)

mspline_14c_c_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

mspline_14c_c_all %>% 
  group_by(ClimateZoneAnd) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  skimr::skim(UD, lyr_14c_msp, CORG_msp, pro_MAT_mod, 
              pro_GPP_Fluxcom_2001_2012_gC_m2d1, lyr_clay_mod, pro_AI)

mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>% 
  dplyr::select(UD, lyr_14c_msp, CORG_msp,pro_MAT_mod, 
                pro_GPP_Fluxcom_2001_2012_gC_m2d1, lyr_clay_mod, pro_AI) %>% 
  cor()

## Plot raw data
climate_all_14c_c_raw <- mspline_14c_c_all %>%
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>% 
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_14c_c_raw$ClimateZoneAnd <- factor(climate_all_14c_c_raw$ClimateZoneAnd,
                                               levels = c("cold temperate", 
                                                          "warm temperate", 
                                                          "arid", "tropical"))

depth_sum <- climate_all_14c_c_raw %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)


c1_14c_c_raw <- climate_all_14c_c_raw %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.81,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Observed SOC [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-500,125), breaks = seq(-500,100,100)) 
plot(c1_14c_c_raw)
ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_cwat_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

### Predict SOC and 14C

# Filter data (only keep variables of interest)
mspline_14c_c_and <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar")

mspline_14c_c_and %>% 
  summarise(n_samples = nrow(.),
            n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

# Create spatial object for cross-validation and resampling
mspline_14c_c_and_sf <- sf::st_as_sf(mspline_14c_c_and, 
                                     coords = c("pro_long", "pro_lat"), 
                                     crs = 4326)

rf_data_14c <- mspline_14c_c_and_sf %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 

rf_data_c <- mspline_14c_c_and_sf %>% 
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 

# Regression tree
# library(rpart)
# library(rpart.plot)
# 
# rt_14c <- rpart(lyr_14c_msp ~ UD + pro_MAT_mod + pro_GPP_Fluxcom_2001_2012_gC_m2d1 +
#           lyr_clay_mod + pro_AI, data = rf_data_14c, 
#           control = rpart.control(maxdepth = 5, minsplit = 500))
# 
# rt_c <- rpart(CORG_msp ~ UD + pro_MAT_mod + pro_GPP_Fluxcom_2001_2012_gC_m2d1 +
#                 lyr_clay_mod + pro_AI, data = rf_data_c, 
#               control = rpart.control(maxdepth = 5, minsplit = 500))
# 
# rpart.plot(rt_14c)
# rpart.plot(rt_c)

# Set-up random forest
# https://geocompr.robinlovelace.net/eco.html#eco
task_rf_14c <- mlr3spatiotempcv::as_task_regr_st(x = rf_data_14c %>% 
                                                   dplyr::select(-id, -entry_name), 
                                                 target = "lyr_14c_msp")

lrn_rf_14c <- lrn("regr.ranger", predict_type = "response",
                  importance = "permutation",
                  num.trees = 1000)

search_space_14c <- paradox::ps(
  mtry = paradox::p_int(lower = 1, upper = ncol(task_rf_14c$data()) - 1),
  sample.fraction = paradox::p_dbl(lower = 0.2, upper = 0.9),
  min.node.size = paradox::p_int(lower = 1, upper = 10)
)

# set.seed(42)
# resampling_sp = rsmp("repeated_spcv_coords", folds = 5)
# autoplot(resampling_sp, task_rf_14c, fold_id = c(1:4), size = 0.7)

set.seed(42)
autotuner_rf_14c <- mlr3tuning::AutoTuner$new(
  learner = lrn_rf_14c,
  resampling = mlr3::rsmp("repeated_spcv_coords", folds = 5), # spatial partitioning
  measure = mlr3::msr("regr.rmse"), # performance measure
  terminator = mlr3tuning::trm("evals", n_evals = 50), # specify 50 iterations
  search_space = search_space_14c, # predefined hyperparameter search space
  tuner = mlr3tuning::tnr("random_search") # specify random search
)


task_rf_c <- mlr3spatiotempcv::as_task_regr_st(x = rf_data_c %>% 
                                                 dplyr::select(-id, -entry_name), 
                                               target = "CORG_msp")

lrn_rf_c <- lrn("regr.ranger", predict_type = "response",
                importance = "permutation",
                num.trees = 1000)

search_space_c <- paradox::ps(
  mtry = paradox::p_int(lower = 1, upper = ncol(task_rf_c$data()) - 1),
  sample.fraction = paradox::p_dbl(lower = 0.2, upper = 0.9),
  min.node.size = paradox::p_int(lower = 1, upper = 10)
)

set.seed(42)
autotuner_rf_c <- mlr3tuning::AutoTuner$new(
  learner = lrn_rf_c,
  resampling = mlr3::rsmp("repeated_spcv_coords", folds = 5), # spatial partitioning
  measure = mlr3::msr("regr.rmse"), # performance measure
  terminator = mlr3tuning::trm("evals", n_evals = 50), # specify 50 iterations
  search_space = search_space_c, # predefined hyperparameter search space
  tuner = mlr3tuning::tnr("random_search") # specify random search
)

# Train model
set.seed(42)
# rf_model_14c <- lrn_rf_14c$train(task_rf_14c)
rf_model_14c_tuned <- autotuner_rf_14c$train(task_rf_14c)

# rf_model_14c$model
rf_model_14c_tuned$tuning_result
rf_model_14c_tuned$predict(task_rf_14c)

set.seed(42)
# rf_model_c <- lrn_rf_c$train(task_rf_c)
rf_model_c_tuned <- autotuner_rf_c$train(task_rf_c)

# rf_model_c$model
rf_model_c_tuned$tuning_result
rf_model_c_tuned$predict(task_rf_14c)

saveRDS(rf_model_14c_tuned, "./Data/rf_tuned_14c.rds")
saveRDS(rf_model_c_tuned, "./Data/rf_tuned_c.rds")

#############################

## Check model
autoplot(rf_model_14c_tuned$predict(task_rf_14c))
autoplot(rf_model_c_tuned$predict(task_rf_c))

rf_pred_14c_df <- data.frame(lyr_14c_pred = rf_model_14c$model$predictions)

rf_data_pred_14c <- cbind(mspline_14c_c_and, rf_pred_14c_df) %>% 
  tibble()

rf_pred_c_df <- data.frame(lyr_c_pred = rf_model_c$model$predictions)

rf_data_pred_c <- cbind(mspline_14c_c_and, rf_pred_c_df) %>% 
  tibble()

rf_data_pred_14c %>% 
  skimr::skim(UD, lyr_14c_msp, lyr_14c_pred, CORG_msp, pro_MAP_mod,
              pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
              lyr_clay_mod, pro_AI)

rf_data_pred_c %>% 
  skimr::skim(UD, lyr_14c_msp, lyr_c_pred, CORG_msp, pro_MAP_mod,
              pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
              lyr_clay_mod, pro_AI)

rf_data_pred_14c %>% 
  ggplot(aes(x = lyr_14c_pred, y = lyr_14c_msp)) +
  geom_point(shape = 21) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed 14C", limits = c(-1000,300), expand = c(0,0)) +
  scale_x_continuous("Predicted 14C", limits = c(-1000,300), expand = c(0,0)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_pred_obs_cwat_", Sys.Date(),
                     ".jpeg"), width = 6, height = 6)

rf_data_pred_c %>% 
  ggplot(aes(x = lyr_c_pred, y = CORG_msp)) +
  geom_point(shape = 21) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed C", limits = c(0,21), expand = c(0,0)) +
  scale_x_continuous("Predicted C", limits = c(0,21), expand = c(0,0)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_pred_obs_cwat_", Sys.Date(),
                     ".jpeg"), width = 6, height = 6)

data.frame(importance = rf_model_14c$importance()) %>% 
  rownames_to_column(var = "predictors") %>% 
  mutate(relImp = importance/sum(importance)*100) %>% 
  ggplot() +
  geom_bar(aes(x = reorder(predictors, -importance), y = relImp), 
           stat = "identity") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous("Relative importance (%; permutation)", expan = c(0,0),
                     limits = c(0,35))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_vrb_imp_cwat_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)

data.frame(importance = rf_model_c$importance()) %>% 
  rownames_to_column(var = "predictors") %>% 
  mutate(relImp = importance/sum(importance)*100) %>% 
  ggplot() +
  geom_bar(aes(x = reorder(predictors, -importance), y = relImp), 
           stat = "identity") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous("Relative importance (%; permutation)", expan = c(0,0),
                     limits = c(0,35))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_vrb_imp_cwat_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)

rf_data_pred_14c_c <- cbind(mspline_14c_c_and, rf_pred_14c_df, rf_pred_c_df) %>% 
  tibble()

climate_all_14c_c <- rf_data_pred_14c_c %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_14c_c$ClimateZoneAnd <- factor(climate_all_14c_c$ClimateZoneAnd,
                                           levels = c("cold temperate", 
                                                      "warm temperate", 
                                                      "arid", "tropical"))

depth_sum <- climate_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)


c1_14c_c <- climate_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.81,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Predicted SOC [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-500,125), breaks = seq(-500,100,100)) 
plot(c1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_climate_cwat_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

### Predictions with one fixed/changed variable

## PET/MAP
mspline_14c_c_and %>% 
  group_by(ClimateZoneAnd) %>% 
  mutate(pro_AI_mod_dry = 0.3 * pro_AI + pro_AI) %>% 
  mutate(pro_AI_mod_wet = pro_AI - (0.3 * pro_AI)) %>% 
  skimr::skim(pro_AI, pro_AI_mod_dry, pro_AI_mod_wet)

rf_data_14c_AI <- mspline_14c_c_and_sf %>% 
  # mutate(pro_AI = pro_AI - (0.30 * pro_AI)) %>% 
  mutate(pro_AI = 0.3 * pro_AI + pro_AI) %>% 
  # mutate(pro_AI = 1) %>%
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 

rf_data_c_AI <- mspline_14c_c_and_sf %>% 
  # mutate(pro_AI = pro_AI - (0.30 * pro_AI)) %>% 
  mutate(pro_AI = 0.3 * pro_AI + pro_AI) %>% 
  # mutate(pro_AI = 1) %>%  
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 

# rf_model_14c_AI <- lrn_rf_14c$predict_newdata(rf_data_14c_AI)
rf_model_14c_tuned_AI <- rf_model_14c_tuned$predict_newdata(rf_data_14c_AI)

# rf_model_c_AI <- lrn_rf_c$predict_newdata(rf_data_c_AI)
rf_model_c_tuned_AI <- rf_model_c_tuned$predict_newdata(rf_data_c_AI)

# rf_pred_14c_AI <- data.frame(lyr_14c_pred = rf_model_14c_AI$data$response)
# rf_pred_c_AI <- data.frame(lyr_c_pred = rf_model_c_AI$data$response)

rf_pred_14c_tuned_AI <- data.frame(lyr_14c_pred = rf_model_14c_tuned_AI$data$response)
rf_pred_c_tuned_AI <- data.frame(lyr_c_pred = rf_model_c_tuned_AI$data$response)

rf_data_pred_14c_c_AI <- cbind(mspline_14c_c_and, rf_pred_14c_tuned_AI, 
                               rf_pred_c_tuned_AI) %>% 
  tibble()

climate_all_14c_c_AI <- rf_data_pred_14c_c_AI %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_14c_c_AI$ClimateZoneAnd <- factor(climate_all_14c_c_AI$ClimateZoneAnd,
                                              levels = c("cold temperate", 
                                                         "warm temperate", 
                                                         "arid", "tropical"))

depth_sum <- climate_all_14c_c_AI %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)


c1_14c_c_AI <- climate_all_14c_c_AI %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.81,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Predicted SOC [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-500,125), breaks = seq(-500,100,100)) 
plot(c1_14c_c_AI)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_climate_AImod_1_cwat_", 
                     Sys.Date(), ".jpeg"), width = 6, height = 5.5)

# Plot deviation from original prediction
climate_all_14c_c_dev <- climate_all_14c_c_AI %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(UD, median_14c, median_c, ClimateZoneAnd) %>% 
  rename(median_14c_pred = median_14c,
         median_c_pred = median_c) %>% 
  full_join(climate_all_14c_c %>% 
              filter(ClimateZoneAnd != "volcanic soils") %>% 
              filter(n > 4 & n_rel > 33) %>% 
              dplyr::select(UD, median_14c, median_c, ClimateZoneAnd)) %>% 
  mutate(median_14c_dev = median_14c_pred - median_14c,
         median_c_dev = median_c_pred- median_c) 

depth_sum <- climate_all_14c_c_dev %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c_dev, median_14c_dev) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)

# climate_all_14c_c_dev %>% 
#   group_by(ClimateZoneAnd) %>% 
#   arrange(UD) %>% 
#   ggplot() +
#   geom_path(aes(x = median_c_dev, y = median_14c_dev, color = ClimateZoneAnd),
#             size = 2) +
#   geom_point(data = depth_sum, aes(x = median_c_dev, y = median_14c_dev), size = 2,
#              shape = 21, fill = "black", color = "white") +
#   theme_bw(base_size = 16) +
#   theme(axis.text = element_text(color = "black"),
#         legend.position = c(0.21,0.2),
#         legend.background = element_blank())

depth_sum %>% 
  group_by(ClimateZoneAnd) %>% 
  arrange(UD) %>% 
  ggplot(aes(x = median_c_dev, y = median_14c_dev)) +
  geom_path(aes(color = ClimateZoneAnd), size = 2) +
  geom_point(aes(fill = UD), size = 4, color = "white", shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.81,0.8),
        legend.background = element_blank()) +
  scale_fill_viridis_c("Depth [cm]", direction = -1, limits = c(0,100))

## GPP
mspline_14c_c_and %>% 
  mutate(pro_GPP_neg30 = pro_GPP_Fluxcom_2001_2012_gC_m2d1 - 
           (0.3 * pro_GPP_Fluxcom_2001_2012_gC_m2d1)) %>% 
  mutate(pro_GPP_pos30 = 0.3 * pro_GPP_Fluxcom_2001_2012_gC_m2d1 + 
           pro_GPP_Fluxcom_2001_2012_gC_m2d1) %>%
  group_by(ClimateZoneAnd) %>% 
  skimr::skim(pro_GPP_Fluxcom_2001_2012_gC_m2d1, pro_GPP_pos30, pro_GPP_neg30) 


rf_data_14c_GPP <- mspline_14c_c_and %>% 
  # mutate(pro_GPP_Fluxcom_2001_2012_gC_m2d1 = pro_GPP_Fluxcom_2001_2012_gC_m2d1 -
  #          (0.3 * pro_GPP_Fluxcom_2001_2012_gC_m2d1)) %>%
  # mutate(pro_GPP_Fluxcom_2001_2012_gC_m2d1 = 0.3 *
  #          pro_GPP_Fluxcom_2001_2012_gC_m2d1 + pro_GPP_Fluxcom_2001_2012_gC_m2d1) %>%
  mutate(pro_GPP_Fluxcom_2001_2012_gC_m2d1 = 10) %>%  
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 

rf_data_c_GPP <- mspline_14c_c_and %>% 
  # mutate(pro_GPP_Fluxcom_2001_2012_gC_m2d1 = pro_GPP_Fluxcom_2001_2012_gC_m2d1 -
  #          (0.3 * pro_GPP_Fluxcom_2001_2012_gC_m2d1)) %>%
  # mutate(pro_GPP_Fluxcom_2001_2012_gC_m2d1 = 0.3 *
  #          pro_GPP_Fluxcom_2001_2012_gC_m2d1 + pro_GPP_Fluxcom_2001_2012_gC_m2d1) %>%
  mutate(pro_GPP_Fluxcom_2001_2012_gC_m2d1 = 10) %>%  
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 

# rf_model_14c_GPP <- lrn_rf_14c$predict_newdata(rf_data_14c_GPP)
rf_model_14c_tuned_GPP <- rf_model_14c_tuned$predict_newdata(rf_data_14c_GPP)

# rf_model_c_GPP <- lrn_rf_c$predict_newdata(rf_data_c_GPP)
rf_model_c_tuned_GPP <- rf_model_c_tuned$predict_newdata(rf_data_c_GPP)

# rf_pred_14c_GPP <- data.frame(lyr_14c_pred = rf_model_14c_GPP$data$response)
# rf_pred_c_GPP <- data.frame(lyr_c_pred = rf_model_c_GPP$data$response)

rf_pred_14c_tuned_GPP <- data.frame(lyr_14c_pred = rf_model_14c_tuned_GPP$data$response)
rf_pred_c_tuned_GPP <- data.frame(lyr_c_pred = rf_model_c_tuned_GPP$data$response)

rf_data_pred_14c_c_GPP <- cbind(mspline_14c_c_and, rf_pred_14c_tuned_GPP, 
                                rf_pred_c_tuned_GPP) %>% 
  tibble()

climate_all_14c_c_GPP <- rf_data_pred_14c_c_GPP %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(lyr_c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_14c_c_GPP$ClimateZoneAnd <- factor(climate_all_14c_c_GPP$ClimateZoneAnd,
                                               levels = c("cold temperate", 
                                                          "warm temperate", 
                                                          "arid", "tropical"))

depth_sum <- climate_all_14c_c_GPP %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)


c1_14c_c_GPP <- climate_all_14c_c_GPP %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = ClimateZoneAnd), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Predicted SOC [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-500,125), breaks = seq(-500,100,100)) 
plot(c1_14c_c_GPP)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_climate_GPPmod_10_", 
                     Sys.Date(), ".jpeg"), width = 6, height = 5.5)

# mspline_14c_c_and %>% 
#   filter(pro_GPP_Fluxcom_2001_2012_gC_m2d1 == min(pro_GPP_Fluxcom_2001_2012_gC_m2d1)) %>% 
#   group_by(UD) %>% 
#   mutate(mean_c = mean(CORG_msp),
#          mean_14c = mean(lyr_14c_msp)) %>% 
#   ggplot() +
#   geom_path(aes(x = mean_c, y = mean_14c, color = ClimateZoneAnd), size = 2) +
#   scale_x_continuous(trans = "log10", limits = c(0.1,10), expand = c(0,0)) +
#   scale_y_continuous(limits = c(-1000,125), expand = c(0,0))

diff_14c_c <- rf_data_pred_14c_c_GPP %>% 
  mutate(CORG_diff = CORG_msp - lyr_c_pred,
         lyr_14c_diff = lyr_14c_msp - lyr_14c_pred) %>% 
  group_by(UD, ClimateZoneAnd) %>% 
  mutate(median_14c = median(lyr_14c_diff),
         median_c = median(CORG_diff))

diff_14c_c$ClimateZoneAnd <- factor(diff_14c_c$ClimateZoneAnd,
                                    levels = c("cold temperate", 
                                               "warm temperate", 
                                               "arid", "tropical"))

depth_sum <- diff_14c_c %>% 
  # filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(ClimateZoneAnd, UD, median_c, median_14c) %>% 
  dplyr::filter(UD == 1|
                  UD == 10|
                  UD == 20|
                  UD == 30|
                  UD == 40|
                  UD == 50|
                  UD == 60|
                  UD == 70|
                  UD == 80|
                  UD == 90|
                  UD == 99)


diff_14c_c_plot <- diff_14c_c %>% 
  arrange(UD) %>% 
  # filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), size = 2) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c, fill = UD, 
                                   color = ClimateZoneAnd), 
             size = 4, shape = 21) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("SOC raw - predicted", expand = c(0,0), limits = c(-2,2)) +
  scale_y_continuous("14C raw - predicted", expand = c(0,0), limits = c(-150,150)) +
  scale_fill_viridis_c(direction = -1, limits = c(0,100))
plot(diff_14c_c_plot)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_diff_climate_GPPmod_10_", 
                     Sys.Date(), ".jpeg"), width = 8, height = 5.5)
