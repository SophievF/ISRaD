# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 25/01/2023 #

#### RANDOM FOREST ####

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(iml)

### Load filtered and splined lyr data
mspline_14c_c_all <- read_csv("./Data/ISRaD_flat_splined_filled_2023-03-09.csv") %>% 
  arrange(UD, id)
  
head(mspline_14c_c_all)
names(mspline_14c_c_all)

mspline_14c_c_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

mspline_14c_c_all %>% 
  group_by(ClimateZoneAnd) %>% 
  # filter(ClimateZoneAnd != "volcanic soils") %>% 
  skimr::skim(UD, lyr_14c_msp, CORG_msp, pro_MAT_mod, 
              pro_GPP_Fluxcom_2001_2012_gC_m2d1, lyr_clay_mod, pro_AI)

mspline_14c_c_all %>% 
  # filter(ClimateZoneAnd != "volcanic soils") %>% 
  # filter(ClimateZoneAnd != "tundra/polar") %>% 
  dplyr::select(UD, lyr_14c_msp, CORG_msp, pro_MAT_mod, 
                lyr_clay_mod, pro_AI) %>% 
  cor()
  

### Predict SOC and 14C

# Filter data (only keep variables of interest)

rf_data_14c <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI) 

rf_data_c <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI) 

## Default set-up of RF 

# Set-up random forest
task_rf_14c <- as_task_regr(x = rf_data_14c %>% 
                              dplyr::select(-id, -entry_name), 
                            target = "lyr_14c_msp")

lrn_rf_14c <- lrn("regr.ranger", importance = "permutation",
              num.trees = 1000)

task_rf_c <- as_task_regr(x = rf_data_c %>% 
                            dplyr::select(-id, -entry_name), 
                          target = "CORG_msp")

lrn_rf_c <- lrn("regr.ranger", importance = "permutation",
                num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c <- lrn_rf_14c$train(task_rf_14c)

rf_model_14c$model

set.seed(42)
rf_model_c <- lrn_rf_c$train(task_rf_c)

rf_model_c$model

## Check model
rf_pred_14c_df <- rf_model_14c$predict(task_rf_14c)
lyr_14c_pred <- rf_pred_14c_df$data$response

rf_pred_c_df <- rf_model_c$predict(task_rf_c)
lyr_c_pred <- rf_pred_c_df$data$response

rf_data_pred_14c_c <- cbind(mspline_14c_c_all, lyr_14c_pred, lyr_c_pred) %>% 
  tibble()

rf_data_pred_14c_c %>% 
  skimr::skim(UD, lyr_14c_msp, lyr_14c_pred, CORG_msp, lyr_c_pred, pro_MAT_mod, 
              pro_GPP_Fluxcom_2001_2012_gC_m2d1, lyr_clay_mod, pro_AI)

rf_data_pred_14c_c %>% 
  ggplot(aes(x = lyr_14c_pred, y = lyr_14c_msp)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed 14C", limits = c(-1000,300), 
                     breaks = seq(-1000,250,250)) +
  scale_x_continuous("Predicted 14C", limits = c(-1000,300), 
                     breaks = seq(-1000,250,250)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_pred_obs_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

rf_data_pred_14c_c %>% 
  ggplot(aes(x = lyr_c_pred, y = CORG_msp)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed C", limits = c(0,55)) +
  scale_x_continuous("Predicted C", limits = c(0,55)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_pred_obs_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

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
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_vrb_imp_", Sys.Date(),
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
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_vrb_imp_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)



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
                                           levels = c("volcanic soils", 
                                                      "tundra/polar", 
                                                      "cold temperate", 
                                                      "warm temperate", 
                                                      "arid", "tropical"))

depth_sum <- climate_all_14c_c %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), linewidth = 2) +
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_pred_SOC_climate_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

mineral_all_14c_c <- rf_data_pred_14c_c %>%
  group_by(MineralType, UD) %>% 
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

depth_sum <- mineral_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(MineralType, UD, median_c, median_14c) %>% 
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


m1_14c_c <- mineral_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = MineralType), linewidth = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralType), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralType), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Predicted SOC [wt-%]", trans = "log10", 
                     limits = c(0.1,35), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(m1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_mineral_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

# Partial dependence plots
num_features <- c("UD", "pro_MAT_mod", "pro_GPP_Fluxcom_2001_2012_gC_m2d1", 
                  "lyr_clay_mod", "pro_AI")

model_14c <- Predictor$new(lrn_rf_14c, data = rf_data_14c %>% 
                             dplyr::select(-id, -entry_name))

effect_14c <- FeatureEffects$new(model_14c, method = "ale")

saveRDS(effect_14c, paste0(getwd(), "/Data/ISRaD_RF_all_default_ale_14c_", 
                           Sys.Date()))

pdp_14c <- plot(effect_14c, features = num_features)

plot(pdp_14c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_default_ale_", Sys.Date(),
                     ".jpeg"), width = 12, height = 5.5, plot = pdp_14c)

model_c <- Predictor$new(lrn_rf_c, data = rf_data_c %>% 
                           dplyr::select(-id, -entry_name))

effect_c <- FeatureEffects$new(model_c, method = "ale")

saveRDS(effect_c, paste0(getwd(), "/Data/ISRaD_RF_all_default_ale_c_", Sys.Date()))

pdp_c <- plot(effect_c, features = num_features)

plot(pdp_c) 
ggsave(file = paste0("./Figure/ISRaD_msp_predSOC_default_ale_", Sys.Date(),
                     ".jpeg"), width = 12, height = 5.5, plot = pdp_c)

## Train model with subset of profiles
# Random choice
set.seed(42)
test_id <- mspline_14c_c_all %>% 
  distinct(id) %>% 
  sample_frac(0.3)

split_data <- mspline_14c_c_all %>% 
  tibble::rowid_to_column() %>% 
  mutate(split = case_when(
    id %in% test_id$id ~ "test",
    !(id %in% test_id$id) ~ "train"
  )) 

split_data %>% 
  group_by(ClimateZoneAnd, split) %>% 
  summarise(n_profile = n_distinct(id))

splits <- split(split_data$rowid, split_data$split)

# Set-up random forest
task_rf_14c <- as_task_regr(x = rf_data_14c %>% 
                              dplyr::select(-id, -entry_name), 
                            target = "lyr_14c_msp")

lrn_rf_14c <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000)

task_rf_c <- as_task_regr(x = rf_data_c %>% 
                            dplyr::select(-id, -entry_name), 
                          target = "CORG_msp")

lrn_rf_c <- lrn("regr.ranger", importance = "permutation",
                num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c <- lrn_rf_14c$train(task_rf_14c, splits$train)

rf_model_14c$model

set.seed(42)
rf_model_c <- lrn_rf_c$train(task_rf_c, splits$train)

rf_model_c$model

## Check model
rf_pred_14c_df <-  lrn_rf_14c$predict(task_rf_14c, splits$test)
lyr_14c_pred <- rf_pred_14c_df$data$response

rf_pred_c_df <-  lrn_rf_c$predict(task_rf_c, splits$test)
lyr_c_pred <- rf_pred_c_df$data$response

rf_data_pred_test_14c_c <- cbind(split_data %>% 
                                   filter(split == "test"), 
                                 lyr_14c_pred, lyr_c_pred) %>% 
  tibble()

summary(lm(lyr_14c_msp ~ lyr_14c_pred, data = rf_data_pred_test_14c_c))
rf_pred_14c_df$score()
sqrt(rf_pred_14c_df$score())

rf_data_pred_test_14c_c %>% 
  ggplot(aes(x = lyr_14c_pred, y = lyr_14c_msp)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed 14C", limits = c(-1000,300), 
                     breaks = seq(-1000,250,250)) +
  scale_x_continuous("Predicted 14C", limits = c(-1000,300), 
                     breaks = seq(-1000,250,250)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_test_pred_obs_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

summary(lm(CORG_msp ~ lyr_c_pred, data = rf_data_pred_test_14c_c))
rf_pred_c_df$score()
sqrt(rf_pred_c_df$score())

rf_data_pred_test_14c_c %>% 
  ggplot(aes(x = lyr_c_pred, y = CORG_msp)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed C", limits = c(0,55)) +
  scale_x_continuous("Predicted C", limits = c(0,55)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_test_pred_obs_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

data.frame(importance = rf_model_14c$importance()) %>% 
  rownames_to_column(var = "predictors") %>% 
  mutate(relImp = importance/sum(importance)*100) %>% 
  ggplot() +
  geom_bar(aes(x = reorder(predictors, -importance), y = relImp), 
           stat = "identity") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous("Relative importance (%; permutation)", expand = c(0,0),
                     limits = c(0,45))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_test_vrb_imp_", Sys.Date(),
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
  scale_y_continuous("Relative importance (%; permutation)", expand = c(0,0),
                     limits = c(0,45))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_test_vrb_imp_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)

# Predict all data to reproduce figures (SOC ~ 14C)
rf_pred_14c_df <-  lrn_rf_14c$predict_newdata(rf_data_14c %>% 
                                                dplyr::select(-id, -entry_name))
lyr_14c_pred <- rf_pred_14c_df$data$response

rf_pred_c_df <- lrn_rf_c$predict_newdata(rf_data_c %>% 
                                           dplyr::select(-id, -entry_name))
lyr_c_pred <- rf_pred_c_df$data$response

rf_data_pred_all_14c_c <- cbind(mspline_14c_c_all, lyr_14c_pred, lyr_c_pred) %>% 
  tibble()

climate_all_14c_c <- rf_data_pred_all_14c_c %>%
  # filter(id %in% test_id$id) %>% 
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
                                           levels = c("volcanic soils", 
                                                      "tundra/polar", 
                                                      "cold temperate", 
                                                      "warm temperate", 
                                                      "arid", "tropical"))

depth_sum <- climate_all_14c_c %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), 
            linewidth = 2) +
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_test_climate_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

mineral_all_14c_c <- rf_data_pred_all_14c_c %>%
  group_by(MineralType, UD) %>% 
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

depth_sum <- mineral_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(MineralType, UD, median_c, median_14c) %>% 
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


m1_14c_c <- mineral_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = MineralType), size = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralType), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralType), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Predicted SOC [wt-%]", trans = "log10", 
                     limits = c(0.1,35), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(m1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_test_mineral_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

# Partial dependence plots
num_features <- c("UD", "pro_MAT_mod", "lyr_clay_mod", "pro_AI")

model_14c <- Predictor$new(rf_model_14c, rf_data_14c %>% 
                             dplyr::select(-id, -entry_name))

# model_14c <- Predictor$new(lrn_rf_14c, data = rf_train_data %>% 
#                              dplyr::select(-id, -entry_name, -CORG_msp))

# effect_14c <- FeatureEffects$new(model_14c, method = "ale")

# saveRDS(effect_14c_ice, paste0(getwd(), "/Data/ISRaD_RF_all_test_ice_14c_",
#                                Sys.Date()))

#pro_MAT
effect_14c_ice_MAT <- FeatureEffects$new(model_14c, method = "ice", 
                                         features = "pro_MAT_mod")

ice_mat <- mspline_14c_c_all %>% 
  dplyr::select(id, MineralType) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_14c_ice_MAT$effects$pro_MAT_mod$results, by = ".id",
             multiple = "all") %>%
  dplyr::group_by(MineralType, pro_MAT_mod) %>% 
  dplyr::mutate(pred_mean_14c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = pro_MAT_mod, y = pred_mean_14c, group = .id, color = MineralType)) +
  geom_path(linewidth = 1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("MAT [°C]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100),
                     expand = c(0,0))

#pro_AI
effect_14c_ice_AI <- FeatureEffects$new(model_14c, method = "ice", 
                                        features = "pro_AI")

ice_AI <- mspline_14c_c_all %>% 
  dplyr::select(id, MineralType) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_14c_ice_AI$effects$pro_AI$results, by = ".id",
             multiple = "all") %>%
  dplyr::group_by(MineralType, pro_AI) %>% 
  dplyr::mutate(pred_mean_14c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = pro_AI, y = pred_mean_14c, group = .id, color = MineralType)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("PET/MAP") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100),
                     expand = c(0,0))

#UD
effect_14c_ice_UD <- FeatureEffects$new(model_14c, method = "ice", 
                                        features = "UD")

ice_UD <- mspline_14c_c_all %>% 
  dplyr::select(id, MineralType) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_14c_ice_UD$effects$UD$results, by = ".id",
             multiple = "all") %>%
  dplyr::group_by(MineralType, UD) %>% 
  dplyr::mutate(pred_mean_14c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = UD, y = pred_mean_14c, group = .id, color = MineralType)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100),
                     expand = c(0,0))

#clay
effect_14c_ice_clay <- FeatureEffects$new(model_14c, method = "ice", 
                                          features = "lyr_clay_mod")

ice_clay <- mspline_14c_c_all %>% 
  dplyr::select(id, MineralType) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_14c_ice_clay$effects$lyr_clay_mod$results, by = ".id",
             multiple = "all") %>%
  dplyr::group_by(MineralType, lyr_clay_mod) %>% 
  dplyr::mutate(pred_mean_14c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = lyr_clay_mod, y = pred_mean_14c, group = .id, color = MineralType)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Clay content [%]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100), 
                     expand = c(0,0))

ggarrange(ice_UD, ice_mat, ice_AI, ice_clay, common.legend = TRUE)

ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_all_ice_mineral_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

# pdp_14c <- plot(effect_14c, features = num_features)
# 
# plot(pdp_14c)
# ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_test_ale_", Sys.Date(),
#                      ".jpeg"), width = 12, height = 5.5, plot = pdp_14c)
# 
model_c <- Predictor$new(rf_model_c, rf_data_c %>%
                           dplyr::select(-id, -entry_name))

# effect_c <- FeatureEffects$new(model_c, method = "ale")

# saveRDS(effect_c, paste0(getwd(), "/Data/ISRaD_RF_all_test_ale_c_", Sys.Date()))

#pro_MAT
effect_c_ice_MAT <- FeatureEffects$new(model_c, method = "ice", 
                                       features = "pro_MAT_mod")

ice_mat_c <- mspline_14c_c_all %>% 
  dplyr::select(id, MineralType) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_c_ice_MAT$effects$pro_MAT_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(MineralType, pro_MAT_mod) %>%
  dplyr::mutate(pred_mean_c = mean(.value)) %>%
  ungroup() %>%
  ggplot(aes(x = pro_MAT_mod, y = pred_mean_c, group = .id, color = MineralType)) +
  geom_path(linewidth = 1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("MAT [°C]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,25), expand = c(0,0))

#pro_AI
effect_c_ice_AI <- FeatureEffects$new(model_c, method = "ice", 
                                      features = "pro_AI")

ice_AI_c <- mspline_14c_c_all %>% 
  dplyr::select(id, MineralType) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_c_ice_AI$effects$pro_AI$results, by = ".id",
             multiple = "all") %>%
  dplyr::group_by(MineralType, pro_AI) %>% 
  dplyr::mutate(pred_mean_c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = pro_AI, y = pred_mean_c, group = .id, color = MineralType)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("PET/MAP") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,25), expand = c(0,0))

#UD
effect_c_ice_UD <- FeatureEffects$new(model_c, method = "ice", 
                                      features = "UD")

ice_UD_c <- mspline_14c_c_all %>% 
  dplyr::select(id, MineralType) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_c_ice_UD$effects$UD$results, by = ".id",
             multiple = "all") %>%
  dplyr::group_by(MineralType, UD) %>% 
  dplyr::mutate(pred_mean_c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = UD, y = pred_mean_c, group = .id, color = MineralType)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,25), expand = c(0,0))

#clay
effect_c_ice_clay <- FeatureEffects$new(model_c, method = "ice", 
                                        features = "lyr_clay_mod")

ice_clay_c <- mspline_14c_c_all %>% 
  dplyr::select(id, MineralType) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_c_ice_clay$effects$lyr_clay_mod$results, by = ".id",
             multiple = "all") %>%
  dplyr::group_by(MineralType, lyr_clay_mod) %>% 
  dplyr::mutate(pred_mean_c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = lyr_clay_mod, y = pred_mean_c, group = .id, color = MineralType)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Clay content [%]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,25), expand = c(0,0))

ggarrange(ice_UD_c, ice_mat_c, ice_AI_c, ice_clay_c, common.legend = TRUE)

ggsave(file = paste0("./Figure/ISRaD_msp_predSOC_all_ice_mineral_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

# pdp_c <- plot(effect_c, features = num_features)
# 
# plot(pdp_c) 
# ggsave(file = paste0("./Figure/ISRaD_msp_predSOC_test_ale_", Sys.Date(),
#                      ".jpeg"), width = 12, height = 5.5, plot = pdp_c)

## Random choice + w/o andisols
set.seed(42)
test_id <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  distinct(id) %>% 
  sample_frac(0.3)

split_data <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  tibble::rowid_to_column() %>% 
  mutate(split = case_when(
    id %in% test_id$id ~ "test",
    !(id %in% test_id$id) ~ "train"
  )) 

split_data %>% 
  group_by(ClimateZoneAnd, split) %>% 
  summarise(n_profile = n_distinct(id))

splits <- split(split_data$rowid, split_data$split)

rf_data_14c_wo_andi <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

rf_data_c_wo_andi <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

# Set-up random forest
task_rf_14c <- as_task_regr(x = rf_data_14c_wo_andi %>% 
                              dplyr::select(-id, -entry_name), 
                            target = "lyr_14c_msp")

lrn_rf_14c <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000)

task_rf_c <- as_task_regr(x = rf_data_c_wo_andi %>% 
                            dplyr::select(-id, -entry_name), 
                          target = "CORG_msp")

lrn_rf_c <- lrn("regr.ranger", importance = "permutation",
                num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c <- lrn_rf_14c$train(task_rf_14c, splits$train)

rf_model_14c$model

set.seed(42)
rf_model_c <- lrn_rf_c$train(task_rf_c, splits$train)

rf_model_c$model

## Check model
rf_pred_14c_df <-  lrn_rf_14c$predict(task_rf_14c, splits$test)
lyr_14c_pred <- rf_pred_14c_df$data$response

rf_pred_c_df <-  lrn_rf_c$predict(task_rf_c, splits$test)
lyr_c_pred <- rf_pred_c_df$data$response

rf_data_pred_test_14c_c <- cbind(split_data %>% 
                                   filter(split == "test"), 
                                 lyr_14c_pred, lyr_c_pred) %>% 
  tibble()

summary(lm(lyr_14c_msp ~ lyr_14c_pred, data = rf_data_pred_test_14c_c))
rf_pred_14c_df$score()
sqrt(rf_pred_14c_df$score())

rf_data_pred_test_14c_c %>% 
  ggplot(aes(x = lyr_14c_pred, y = lyr_14c_msp)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed 14C", limits = c(-1000,300), 
                     breaks = seq(-1000,250,250)) +
  scale_x_continuous("Predicted 14C", limits = c(-1000,300), 
                     breaks = seq(-1000,250,250)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_test_pred_obs_wo_Andi_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

summary(lm(CORG_msp ~ lyr_c_pred, data = rf_data_pred_test_14c_c))
rf_pred_c_df$score()
sqrt(rf_pred_c_df$score())

rf_data_pred_test_14c_c %>% 
  ggplot(aes(x = lyr_c_pred, y = CORG_msp)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed C", limits = c(0,55)) +
  scale_x_continuous("Predicted C", limits = c(0,55)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_test_pred_obs_wo_Andi_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

data.frame(importance = rf_model_14c$importance()) %>% 
  rownames_to_column(var = "predictors") %>% 
  mutate(relImp = importance/sum(importance)*100) %>% 
  ggplot() +
  geom_bar(aes(x = reorder(predictors, -importance), y = relImp), 
           stat = "identity") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous("Relative importance (%; permutation)", expand = c(0,0),
                     limits = c(0,35))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_test_vrb_imp_wo_Andi_", Sys.Date(),
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
  scale_y_continuous("Relative importance (%; permutation)", expand = c(0,0),
                     limits = c(0,36))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_test_vrb_imp_wo_Andi_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)

# Predict all data to reproduce figures (SOC ~ 14C)
rf_pred_14c_df <- lrn_rf_14c$predict_newdata(rf_data_14c_wo_andi %>% 
                                               dplyr::select(-id, -entry_name))
lyr_14c_pred <- rf_pred_14c_df$data$response

rf_pred_c_df <- lrn_rf_c$predict_newdata(rf_data_c_wo_andi %>% 
                                           dplyr::select(-id, -entry_name))
lyr_c_pred <- rf_pred_c_df$data$response

rf_data_pred_all_14c_c <- cbind(mspline_14c_c_all %>% 
                                  filter(ClimateZoneAnd != "volcanic soils"), 
                                lyr_14c_pred, lyr_c_pred) %>% 
  tibble()

climate_all_14c_c <- rf_data_pred_all_14c_c %>%
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
                                           levels = c("volcanic soils", 
                                                      "tundra/polar", 
                                                      "cold temperate", 
                                                      "warm temperate", 
                                                      "arid", "tropical"))

depth_sum <- climate_all_14c_c %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), linewidth = 2) +
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_test_climate_wo_AndiPolar_", 
                     Sys.Date(), ".jpeg"), width = 6, height = 5.5)

mineral_all_14c_c <- rf_data_pred_all_14c_c %>%
  group_by(MineralType, UD) %>% 
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

depth_sum <- mineral_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(MineralType, UD, median_c, median_14c) %>% 
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


m1_14c_c <- mineral_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = MineralType), linewidth = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralType), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralType), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Predicted SOC [wt-%]", trans = "log10", 
                     limits = c(0.1,35), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(m1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_test_mineral_wo_andi_", 
                     Sys.Date(), ".jpeg"), width = 6, height = 5.5)

# Partial dependence plots
num_features <- c("UD", "pro_MAT_mod", "lyr_clay_mod", "pro_AI")

model_14c <- Predictor$new(rf_model_14c, rf_data_14c_wo_andi %>%  
                             dplyr::select(-id, -entry_name))

# model_14c <- Predictor$new(lrn_rf_14c, data = rf_train_data %>% 
#                              dplyr::select(-id, -entry_name, -CORG_msp))

# effect_14c <- FeatureEffects$new(model_14c, method = "ice")

# saveRDS(effect_14c, paste0(getwd(), "/Data/ISRaD_RF_all_test_ale_14c_", 
#                            Sys.Date()))

# pdp_14c <- plot(effect_14c, features = num_features)
# 
# plot(pdp_14c)
# ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_test_pdp_wo_andi", Sys.Date(),
#                      ".jpeg"), width = 12, height = 5.5, plot = pdp_14c)

#pro_MAT
effect_14c_ice_MAT <- FeatureEffects$new(model_14c, method = "ice", 
                                         features = "pro_MAT_mod")

ice_mat <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  right_join(effect_14c_ice_MAT$effects$pro_MAT_mod$results, by = ".id") %>%
  dplyr::group_by(ClimateZoneAnd, pro_MAT_mod) %>% 
  dplyr::mutate(pred_mean_14c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = pro_MAT_mod, y = pred_mean_14c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("MAT [°C]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100),
                     expand = c(0,0))

#pro_AI
effect_14c_ice_AI <- FeatureEffects$new(model_14c, method = "ice", 
                                        features = "pro_AI")

ice_AI <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_14c_ice_AI$effects$pro_AI$results, by = ".id") %>%
  dplyr::group_by(ClimateZoneAnd, pro_AI) %>% 
  dplyr::mutate(pred_mean_14c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = pro_AI, y = pred_mean_14c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("PET/MAP") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100),
                     expand = c(0,0))

#UD
effect_14c_ice_UD <- FeatureEffects$new(model_14c, method = "ice", 
                                        features = "UD")

ice_UD <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_14c_ice_UD$effects$UD$results, by = ".id") %>%
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(pred_mean_14c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = UD, y = pred_mean_14c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100),
                     expand = c(0,0))

#clay
effect_14c_ice_clay <- FeatureEffects$new(model_14c, method = "ice", 
                                          features = "lyr_clay_mod")

ice_clay <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_14c_ice_clay$effects$lyr_clay_mod$results, by = ".id") %>%
  dplyr::group_by(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::mutate(pred_mean_14c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = lyr_clay_mod, y = pred_mean_14c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Clay content [%]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100), 
                     expand = c(0,0))

ggarrange(ice_UD, ice_mat, ice_AI, ice_clay, common.legend = TRUE)

ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_wo_Andi_ice_climate_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

model_c <- Predictor$new(rf_model_c, rf_data_c_wo_andi %>% 
                           dplyr::select(-id, -entry_name))

# model_c <- Predictor$new(lrn_rf_c, data = rf_train_data %>% 
#                            dplyr::select(-id, -entry_name, -lyr_14c_msp))

# effect_c <- FeatureEffects$new(model_c, method = "pdp")

# saveRDS(effect_c, paste0(getwd(), "/Data/ISRaD_RF_all_test_ale_c_", Sys.Date()))

# pdp_c <- plot(effect_c, features = num_features)
# 
# plot(pdp_c) 
# ggsave(file = paste0("./Figure/ISRaD_msp_predSOC_test_pdp_wo_andi_", Sys.Date(),
#                      ".jpeg"), width = 12, height = 5.5, plot = pdp_c)

#pro_MAT
effect_c_ice_MAT <- FeatureEffects$new(model_c, method = "ice", 
                                       features = "pro_MAT_mod")

ice_mat_c <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_c_ice_MAT$effects$pro_MAT_mod$results, by = ".id", 
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_MAT_mod) %>% 
  dplyr::mutate(pred_mean_c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = pro_MAT_mod, y = pred_mean_c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("MAT [°C]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,25), expand = c(0,0))

#pro_AI
effect_c_ice_AI <- FeatureEffects$new(model_c, method = "ice", 
                                      features = "pro_AI")

ice_AI_c <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_c_ice_AI$effects$pro_AI$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_AI) %>% 
  dplyr::mutate(pred_mean_c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = pro_AI, y = pred_mean_c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("PET/MAP") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,25), expand = c(0,0))

#UD
effect_c_ice_UD <- FeatureEffects$new(model_c, method = "ice", 
                                      features = "UD")

ice_UD_c <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_c_ice_UD$effects$UD$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(pred_mean_c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = UD, y = pred_mean_c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,25), expand = c(0,0))

#clay
effect_c_ice_clay <- FeatureEffects$new(model_c, method = "ice", 
                                        features = "lyr_clay_mod")

ice_clay_c <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_c_ice_clay$effects$lyr_clay_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::mutate(pred_mean_c = mean(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = lyr_clay_mod, y = pred_mean_c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Clay content [%]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,25), expand = c(0,0))

ggarrange(ice_UD_c, ice_mat_c, ice_AI_c, ice_clay_c, common.legend = TRUE)

ggsave(file = paste0("./Figure/ISRaD_msp_predSOC_Andi_ice_climate_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

## Random choice + w/o andisols, polar/tundra
set.seed(1234)
test_id <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>% 
  distinct(id) %>% 
  sample_frac(0.3)

split_data <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  tibble::rowid_to_column() %>% 
  mutate(split = case_when(
    id %in% test_id$id ~ "test",
    !(id %in% test_id$id) ~ "train"
  )) 

split_data %>% 
  group_by(ClimateZoneAnd, split) %>% 
  summarise(n_profile = n_distinct(id))

splits <- split(split_data$rowid, split_data$split)

rf_data_14c_wo_andi <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

rf_data_c_wo_andi <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, lyr_clay_mod, pro_AI)

# Set-up random forest
task_rf_14c <- as_task_regr(x = rf_data_14c_wo_andi %>% 
                              dplyr::select(-id, -entry_name), 
                            target = "lyr_14c_msp")

lrn_rf_14c <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000)

task_rf_c <- as_task_regr(x = rf_data_c_wo_andi %>% 
                            dplyr::select(-id, -entry_name), 
                          target = "CORG_msp")

lrn_rf_c <- lrn("regr.ranger", importance = "permutation",
                num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c <- lrn_rf_14c$train(task_rf_14c, splits$train)

rf_model_14c$model

set.seed(42)
rf_model_c <- lrn_rf_c$train(task_rf_c, splits$train)

rf_model_c$model

## Check model
rf_pred_14c_df <-  lrn_rf_14c$predict(task_rf_14c, splits$test)
lyr_14c_pred <- rf_pred_14c_df$data$response

rf_pred_c_df <-  lrn_rf_c$predict(task_rf_c, splits$test)
lyr_c_pred <- rf_pred_c_df$data$response

rf_data_pred_test_14c_c <- cbind(split_data %>% 
                                   filter(split == "test"), 
                                 lyr_14c_pred, lyr_c_pred) %>% 
  tibble()

summary(lm(lyr_14c_msp ~ lyr_14c_pred, data = rf_data_pred_test_14c_c))
rf_pred_14c_df$score()
sqrt(rf_pred_14c_df$score())

rf_data_pred_test_14c_c %>% 
  ggplot(aes(x = lyr_14c_pred, y = lyr_14c_msp)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed 14C", limits = c(-1000,300), 
                     breaks = seq(-1000,250,250)) +
  scale_x_continuous("Predicted 14C", limits = c(-1000,300), 
                     breaks = seq(-1000,250,250)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_test_pred_obs_wo_AndiPolar_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

summary(lm(CORG_msp ~ lyr_c_pred, data = rf_data_pred_test_14c_c))
rf_pred_c_df$score()
sqrt(rf_pred_c_df$score())

rf_data_pred_test_14c_c %>% 
  ggplot(aes(x = lyr_c_pred, y = CORG_msp)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed C", limits = c(0,45)) +
  scale_x_continuous("Predicted C", limits = c(0,45)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_test_pred_obs_wo_AndiPolar_", Sys.Date(),
                     ".jpeg"), width = 5, height = 5)

data.frame(importance = rf_model_14c$importance()) %>% 
  rownames_to_column(var = "predictors") %>% 
  mutate(relImp = importance/sum(importance)*100) %>% 
  ggplot() +
  geom_bar(aes(x = reorder(predictors, -importance), y = relImp), 
           stat = "identity") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous("Relative importance (%; permutation)", expand = c(0,0),
                     limits = c(0,40))
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_test_vrb_imp_wo_AndiPolar_", Sys.Date(),
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
  scale_y_continuous("Relative importance (%; permutation)", expand = c(0,0),
                     limits = c(0,40))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_test_vrb_imp_wo_AndiPolar_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)

# Predict all data to reproduce figures (SOC ~ 14C)
rf_pred_14c_df <- lrn_rf_14c$predict_newdata(rf_data_14c_wo_andi %>% 
                                                dplyr::select(-id, -entry_name))
lyr_14c_pred <- rf_pred_14c_df$data$response

rf_pred_c_df <- lrn_rf_c$predict_newdata(rf_data_c_wo_andi %>% 
                                           dplyr::select(-id, -entry_name))
lyr_c_pred <- rf_pred_c_df$data$response

rf_data_pred_all_14c_c <- cbind(mspline_14c_c_all %>% 
                                  filter(ClimateZoneAnd != "volcanic soils") %>% 
                                  filter(ClimateZoneAnd != "tundra/polar"), 
                                lyr_14c_pred, lyr_c_pred) %>% 
  tibble()

climate_all_14c_c <- rf_data_pred_all_14c_c %>%
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
                                           levels = c("volcanic soils", 
                                                      "tundra/polar", 
                                                      "cold temperate", 
                                                      "warm temperate", 
                                                      "arid", "tropical"))

depth_sum <- climate_all_14c_c %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = ClimateZoneAnd), linewidth = 2) +
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_test_climate_wo_AndiPolar_", 
                     Sys.Date(), ".jpeg"), width = 6, height = 5.5)

mineral_all_14c_c <- rf_data_pred_all_14c_c %>%
  group_by(MineralType, UD) %>% 
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

depth_sum <- mineral_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  dplyr::select(MineralType, UD, median_c, median_14c) %>% 
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


m1_14c_c <- mineral_all_14c_c %>% 
  filter(n > 4 & n_rel > 33) %>% 
  ggplot() + 
  geom_path(aes(x = median_c, y = median_14c, color = MineralType), linewidth = 2) +
  geom_errorbar(aes(ymin = lci_14c, ymax = uci_14c, x = median_c, 
                    color = MineralType), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = median_14c, 
                     color = MineralType), alpha = 0.3) +
  geom_point(data = depth_sum, aes(x = median_c, y = median_14c), size = 2,
             shape = 21, fill = "black", color = "white") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.21,0.2),
        legend.background = element_blank()) +
  scale_x_continuous("Predicted SOC [wt-%]", trans = "log10", 
                     limits = c(0.1,35), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(m1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_test_mineral_wo_andi_", 
                     Sys.Date(), ".jpeg"), width = 6, height = 5.5)

# Partial dependence plots
num_features <- c("UD", "pro_MAT_mod", "lyr_clay_mod", "pro_AI")

model_14c <- Predictor$new(rf_model_14c, rf_data_14c_wo_andi %>%  
                             dplyr::select(-id, -entry_name))

# model_14c <- Predictor$new(lrn_rf_14c, data = rf_train_data %>% 
#                              dplyr::select(-id, -entry_name, -CORG_msp))

# effect_14c <- FeatureEffects$new(model_14c, method = "ice")

# saveRDS(effect_14c, paste0(getwd(), "/Data/ISRaD_RF_all_test_ale_14c_", 
#                            Sys.Date()))

# pdp_14c <- plot(effect_14c, features = num_features)
# 
# plot(pdp_14c)
# ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_test_pdp_wo_andi", Sys.Date(),
#                      ".jpeg"), width = 12, height = 5.5, plot = pdp_14c)

#pro_MAT
# effect_14c_pdp_MAT <- FeatureEffects$new(model_14c, method = "pdp", 
#                                          features = "pro_MAT_mod")
# 
# pdp_mat <- effect_14c_pdp_MAT$effects$pro_MAT_mod$results %>% 
#   ggplot(aes(x = pro_MAT_mod, y = .value)) +
#   geom_path(linewidth = 1) +
#   theme_bw(base_size = 14) +
#   theme(axis.text = element_text(color = "black")) +
#   scale_x_continuous("MAT [°C]") +
#   scale_y_continuous("Mean predicted 14C", limits = c(-400,100),
#                      expand = c(0,0))

effect_14c_ice_MAT <- FeatureEffects$new(model_14c, method = "ice",
                                         features = "pro_MAT_mod")

ice_mat <- mspline_14c_c_all %>%
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>%
  rownames_to_column(var = ".id") %>%
  mutate(.id = as.integer(.id)) %>%
  full_join(effect_14c_ice_MAT$effects$pro_MAT_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_MAT_mod) %>%
  dplyr::mutate(pred_mean_14c = median(.value)) %>%
  ungroup() %>%
  ggplot(aes(x = pro_MAT_mod, y = pred_mean_14c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("MAT [°C]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100),
                     expand = c(0,0))

#pro_AI
# effect_14c_pdp_AI <- FeatureEffects$new(model_14c, method = "pdp", 
#                                          features = "pro_AI")
# 
# pdp_AI <- effect_14c_pdp_AI$effects$pro_AI$results %>% 
#   ggplot(aes(x = pro_AI, y = .value)) +
#   geom_path(linewidth = 1) +
#   theme_bw(base_size = 14) +
#   theme(axis.text = element_text(color = "black")) +
#   scale_x_continuous("PET/MAP") +
#   scale_y_continuous("Mean predicted 14C", limits = c(-400,100),
#                      expand = c(0,0))

effect_14c_ice_AI <- FeatureEffects$new(model_14c, method = "ice",
                                        features = "pro_AI")

ice_AI <- mspline_14c_c_all %>%
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>%
  rownames_to_column(var = ".id") %>%
  mutate(.id = as.integer(.id)) %>%
  full_join(effect_14c_ice_AI$effects$pro_AI$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_AI) %>%
  dplyr::mutate(pred_mean_14c = median(.value)) %>%
  ungroup() %>%
  ggplot(aes(x = pro_AI, y = pred_mean_14c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("PET/MAP") +
  scale_y_continuous("Mean predicted 14C", limits = c(-800,100),
                     expand = c(0,0))

#UD
effect_14c_ice_UD <- FeatureEffects$new(model_14c, method = "ice", 
                                        features = "UD")

ice_UD <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_14c_ice_UD$effects$UD$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(pred_mean_14c = median(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = UD, y = pred_mean_14c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-400,100),
                     expand = c(0,0))

#clay
effect_14c_ice_clay <- FeatureEffects$new(model_14c, method = "ice", 
                                          features = "lyr_clay_mod")

ice_clay <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_14c_ice_clay$effects$lyr_clay_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::mutate(pred_mean_14c = median(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = lyr_clay_mod, y = pred_mean_14c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Clay content [%]") +
  scale_y_continuous("Mean predicted 14C", limits = c(-400,100), 
                     expand = c(0,0))

ggarrange(ggarrange(ice_UD, ice_clay, common.legend = TRUE), 
          ggarrange(ice_mat, ice_AI, common.legend = TRUE), nrow = 2)

ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_wo_AndiPolar_ice_climate_depth20_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)


# model_c <- Predictor$new(lrn_rf_c, data = rf_train_data %>% 
#                            dplyr::select(-id, -entry_name, -lyr_14c_msp))

# effect_c <- FeatureEffects$new(model_c, method = "pdp")

# saveRDS(effect_c, paste0(getwd(), "/Data/ISRaD_RF_all_test_ale_c_", Sys.Date()))

# pdp_c <- plot(effect_c, features = num_features)
# 
# plot(pdp_c) 
# ggsave(file = paste0("./Figure/ISRaD_msp_predSOC_test_pdp_wo_andi_", Sys.Date(),
#                      ".jpeg"), width = 12, height = 5.5, plot = pdp_c)

model_c <- Predictor$new(rf_model_c, rf_data_c_wo_andi %>% 
                           dplyr::select(-id, -entry_name))

#pro_MAT
# effect_14c_pdp_MAT <- FeatureEffects$new(model_c, method = "pdp",
#                                          features = "pro_MAT_mod")
# 
# pdp_mat_c <- effect_14c_pdp_MAT$effects$pro_MAT_mod$results %>%
#   ggplot(aes(x = pro_MAT_mod, y = .value)) +
#   geom_path(linewidth = 1) +
#   theme_bw(base_size = 14) +
#   theme(axis.text = element_text(color = "black")) +
#   scale_x_continuous("MAT [°C]") +
#   scale_y_continuous("Mean predicted SOC", limits = c(0,6), expand = c(0,0))

effect_c_ice_MAT <- FeatureEffects$new(model_c, method = "ice",
                                       features = "pro_MAT_mod")

ice_mat_c <- mspline_14c_c_all %>%
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>%
  rownames_to_column(var = ".id") %>%
  mutate(.id = as.integer(.id)) %>%
  full_join(effect_c_ice_MAT$effects$pro_MAT_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_MAT_mod) %>%
  dplyr::mutate(pred_mean_c = median(.value)) %>%
  ungroup() %>%
  ggplot(aes(x = pro_MAT_mod, y = pred_mean_c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("MAT [°C]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,6), expand = c(0,0))

#pro_AI
# effect_14c_pdp_AI <- FeatureEffects$new(model_c, method = "pdp", 
#                                         features = "pro_AI")
# 
# pdp_AI_c <- effect_14c_pdp_AI$effects$pro_AI$results %>% 
#   ggplot(aes(x = pro_AI, y = .value)) +
#   geom_path(linewidth = 1) +
#   theme_bw(base_size = 14) +
#   theme(axis.text = element_text(color = "black")) +
#   scale_x_continuous("MAT [°C]") +
#   scale_y_continuous("Mean predicted SOC", limits = c(0,6), expand = c(0,0))

effect_c_ice_AI <- FeatureEffects$new(model_c, method = "ice",
                                      features = "pro_AI")

ice_AI_c <- mspline_14c_c_all %>%
  filter(ClimateZoneAnd != "volcanic soils") %>%
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>%
  rownames_to_column(var = ".id") %>%
  mutate(.id = as.integer(.id)) %>%
  full_join(effect_c_ice_AI$effects$pro_AI$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, pro_AI) %>%
  dplyr::mutate(pred_mean_c = median(.value)) %>%
  ungroup() %>%
  ggplot(aes(x = pro_AI, y = pred_mean_c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("PET/MAP") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,6), expand = c(0,0))

#UD
effect_c_ice_UD <- FeatureEffects$new(model_c, method = "ice", 
                                      features = "UD")

ice_UD_c <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_c_ice_UD$effects$UD$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(pred_mean_c = median(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = UD, y = pred_mean_c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Depth [cm]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,6), expand = c(0,0))

#clay
effect_c_ice_clay <- FeatureEffects$new(model_c, method = "ice", 
                                        features = "lyr_clay_mod")

ice_clay_c <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  full_join(effect_c_ice_clay$effects$lyr_clay_mod$results, by = ".id",
            multiple = "all") %>%
  dplyr::group_by(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::mutate(pred_mean_c = median(.value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = lyr_clay_mod, y = pred_mean_c, group = .id, color = ClimateZoneAnd)) +
  geom_path(linewidth = 1) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("Clay content [%]") +
  scale_y_continuous("Mean predicted SOC", limits = c(0,6), expand = c(0,0))

# ggarrange(ice_UD_c, pdp_mat_c, pdp_AI_c, ice_clay_c, common.legend = TRUE)

ggarrange(ggarrange(ice_UD_c, ice_clay_c, common.legend = TRUE), 
          ggarrange(ice_mat_c, ice_AI_c, common.legend = TRUE), nrow = 2)

ggsave(file = paste0("./Figure/ISRaD_msp_predSOC_AndiPolar_ice_climate_depth20_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

#Merge both PDP's
mat_ice <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(effect_14c_ice_MAT$effects$pro_MAT_mod$results %>% 
              dplyr::rename(pred_14c = .value), multiple = "all") %>% 
  left_join(effect_c_ice_MAT$effects$pro_MAT_mod$results %>% 
              dplyr::rename(pred_c = .value), multiple = "all") 

write_csv(mat_ice, file = paste0("./Data/ISRaD_RF_wo_andipol_ice_mat_14c_c_", 
                                  Sys.Date(), ".csv"))

ai_ice <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(effect_14c_ice_AI$effects$pro_AI$results %>% 
              dplyr::rename(pred_14c = .value), multiple = "all") %>% 
  left_join(effect_c_ice_AI$effects$pro_AI$results %>% 
              dplyr::rename(pred_c = .value), multiple = "all") 

write_csv(ai_ice, file = paste0("./Data/ISRaD_RF_wo_andipol_ice_ai_14c_c_", 
                                 Sys.Date(), ".csv"))

ud_ice <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(effect_14c_ice_UD$effects$UD$results %>% 
              dplyr::rename(pred_14c = .value), multiple = "all") %>% 
  left_join(effect_c_ice_UD$effects$UD$results %>% 
              dplyr::rename(pred_c = .value), multiple = "all") 

write_csv(ud_ice, file = paste0("./Data/ISRaD_RF_wo_andipol_ice_ud_14c_c_", 
                                Sys.Date(), ".csv"))

clay_ice <- mspline_14c_c_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  filter(ClimateZoneAnd != "tundra/polar") %>%
  dplyr::select(id, ClimateZoneAnd) %>% 
  rownames_to_column(var = ".id") %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(effect_14c_ice_clay$effects$lyr_clay_mod$results %>% 
              dplyr::rename(pred_14c = .value), multiple = "all") %>% 
  left_join(effect_c_ice_clay$effects$lyr_clay_mod$results %>% 
              dplyr::rename(pred_c = .value), multiple = "all") 

write_csv(clay_ice, file = paste0("./Data/ISRaD_RF_wo_andipol_ice_clay_14c_c_", 
                                  Sys.Date(), ".csv"))

clay_ice %>% 
  dplyr::group_by(ClimateZoneAnd, lyr_clay_mod) %>% 
  dplyr::mutate(pred_median_c = median(pred_c),
                pred_median_14c = median(pred_14c)) %>% 
  ungroup() %>% 
  arrange(lyr_clay_mod) %>% 
  ggplot(aes(y = pred_median_14c, x = pred_median_c, group = .id, color = ClimateZoneAnd)) +
  geom_point(aes(size = lyr_clay_mod)) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Mean predicted 14C", limits = c(-300,0), expand = c(0,0)) +
  scale_x_continuous("Mean predicted SOC", limits = c(0.3,1.6), position = "bottom",
                     expand = c(0,0)) +
  coord_trans(x = "log10", xlim = c(0.3,2)) +
  annotation_logticks(sides = "t", 
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm"))
ggsave(file = paste0("./Figure/ISRaD_msp_SOC_14c_AndiPolar_ice_climate_2_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)
