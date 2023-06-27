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

#Load filtered and splined lyr data
mspline_14c_c_all <- read_csv("./Data/ISRaD_flat_splined_filled_2023-01-25.csv")

head(mspline_14c_c_all)
names(mspline_14c_c_all)

mspline_14c_c_all %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

mspline_14c_c_all %>% 
  group_by(ClimateZoneAnd) %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
  skimr::skim(UD, lyr_14c_msp, CORG_msp, pro_MAP_mod,
              pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
              lyr_clay_mod, pro_AI)
  

#Filter data (only keep variables of interest)
rf_data <- mspline_14c_c_all %>% 
  
  dplyr::select(id, entry_name, UD, lyr_14c_msp, CORG_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 

rf_data %>% 
  skimr::skim()

## Set-up random forest
task_rf <- as_task_regr(x = rf_data %>% 
                          dplyr::select(-entry_name), 
                        target = "lyr_14c_msp")

lrn_rf <- lrn("regr.ranger", importance = "permutation",
              num.trees = 1000)

# Train model
set.seed(42)
# rf_model_split <- lrn_rf$train(task_rf, splits$train)

rf_model <- lrn_rf$train(task_rf)

rf_model$model

## Check model
rf_pred_df <- data.frame(lyr_14c_pred = rf_model$model$predictions)

rf_data_pred <- cbind(mspline_14c_c_all, rf_pred_df) %>% 
  tibble()

rf_data_pred %>% 
  skimr::skim(UD, lyr_14c_msp, lyr_14c_pred, CORG_msp, pro_MAP_mod,
              pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
              lyr_clay_mod, pro_AI)

rf_data_pred %>% 
  ggplot(aes(x = lyr_14c_pred, y = lyr_14c_msp)) +
  geom_point(shape = 21) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed 14C", limits = c(-1000,300), expand = c(0,0)) +
  scale_x_continuous("Predicted 14C", limits = c(-1000,300), expand = c(0,0)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_RF_pred_obs_", Sys.Date(),
       ".jpeg"), width = 6, height = 6)

# predictions <- lrn_rf$predict_newdata(rf_data)
# autoplot(predictions)
# 
# predictions <- lrn_rf$predict(task_rf, splits$test)
# autoplot(predictions)

data.frame(importance = rf_model$importance()) %>% 
  rownames_to_column(var = "predictors") %>% 
  mutate(relImp = importance/sum(importance)*100) %>% 
  ggplot() +
  geom_bar(aes(x = reorder(predictors, -importance), y = relImp), 
           stat = "identity") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous("Relative importance (%; permutation)", expan = c(0,0),
                     limits = c(0,25))
ggsave(file = paste0("./Figure/ISRaD_msp_RF_vrb_imp_", Sys.Date(),
                     ".jpeg"), width = 8, height = 6)

climate_all <- rf_data_pred %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all$ClimateZoneAnd <- factor(climate_all$ClimateZoneAnd,
                                     levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                "warm temperate", "arid", "tropical"))

depth_sum <- climate_all %>% 
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


c1 <- climate_all %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_depth_SOC_SOC_climate_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

mineral_all <- rf_data_pred %>%
  group_by(MineralType, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

depth_sum <- mineral_all %>% 
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

m1 <- mineral_all %>% 
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
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,35), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(m1)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_SOC_mineral_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

### Test influence of individual predictors on 14C ~ SOC

## Depth + SOC
rf_data_depth <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, CORG_msp) 

rf_data_depth %>% 
  skimr::skim()

# Set-up random forest

task_rf_depth <- as_task_regr(x = rf_data_depth %>% 
                                dplyr::select(-id, -entry_name), 
                              target = "lyr_14c_msp")

lrn_rf_depth <- lrn("regr.ranger", importance = "permutation",
                    num.trees = 1000)

# Train model
set.seed(42)
rf_model_depth <- lrn_rf_depth$train(task_rf_depth)

rf_model_depth$model

rf_depth_pred_df <- data.frame(lyr_14c_pred = rf_model_depth$model$predictions)

rf_depth_data_pred <- cbind(mspline_14c_c_all, rf_depth_pred_df) %>% 
  tibble()

climate_all_depth <- rf_depth_data_pred %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_depth$ClimateZoneAnd <- factor(climate_all_depth$ClimateZoneAnd,
                                           levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                      "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_depth %>% 
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


c1_depth <- climate_all_depth %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_depth)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_SOC_climate_depth_SOC_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

## MAT
rf_data_MAT <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, CORG_msp, pro_MAT_mod) 

rf_data_MAT %>% 
  skimr::skim()

# Set-up random forest

task_rf_MAT <- as_task_regr(x = rf_data_MAT %>% 
                              dplyr::select(-id, -entry_name), 
                            target = "lyr_14c_msp")

lrn_rf_MAT <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000)

# Train model
set.seed(42)
rf_model_MAT <- lrn_rf_MAT$train(task_rf_MAT)

rf_model_MAT$model

rf_MAT_pred_df <- data.frame(lyr_14c_pred = rf_model_MAT$model$predictions)

rf_MAT_data_pred <- cbind(mspline_14c_c_all, rf_MAT_pred_df) %>% 
  tibble()

climate_all_MAT <- rf_MAT_data_pred %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_MAT$ClimateZoneAnd <- factor(climate_all_MAT$ClimateZoneAnd,
                                         levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                    "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_MAT %>% 
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


c1_MAT <- climate_all_MAT %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_MAT)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_SOC_climate_MAT_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

## Clay content
rf_data_clay <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, CORG_msp, lyr_clay_mod) 

rf_data_clay %>% 
  skimr::skim()

# Set-up random forest
task_rf_clay <- as_task_regr(x = rf_data_clay %>% 
                               dplyr::select(-id, -entry_name), 
                             target = "lyr_14c_msp")

lrn_rf_clay <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000)

# Train model
set.seed(42)
rf_model_clay <- lrn_rf_clay$train(task_rf_clay)

rf_model_clay$model

rf_clay_pred_df <- data.frame(lyr_14c_pred = rf_model_clay$model$predictions)

rf_clay_data_pred <- cbind(mspline_14c_c_all, rf_clay_pred_df) %>% 
  tibble()

climate_all_clay <- rf_clay_data_pred %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_clay$ClimateZoneAnd <- factor(climate_all_clay$ClimateZoneAnd,
                                         levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                    "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_clay %>% 
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


c1_clay <- climate_all_clay %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_clay)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_SOC_climate_clay_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

## GPP
rf_data_GPP <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, CORG_msp, pro_GPP_Fluxcom_2001_2012_gC_m2d1) 

rf_data_GPP %>% 
  skimr::skim()

# Set-up random forest
task_rf_GPP <- as_task_regr(x = rf_data_GPP %>% 
                              dplyr::select(-id, -entry_name), 
                            target = "lyr_14c_msp")

lrn_rf_GPP <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000)

# Train model
set.seed(42)
rf_model_GPP <- lrn_rf_GPP$train(task_rf_GPP)

rf_model_GPP$model

rf_GPP_pred_df <- data.frame(lyr_14c_pred = rf_model_GPP$model$predictions)

rf_GPP_data_pred <- cbind(mspline_14c_c_all, rf_GPP_pred_df) %>% 
  tibble()

climate_all_GPP <- rf_GPP_data_pred %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_GPP$ClimateZoneAnd <- factor(climate_all_GPP$ClimateZoneAnd,
                                          levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                     "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_GPP %>% 
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


c1_GPP <- climate_all_GPP %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_GPP)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_SOC_climate_GPP_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

## PET/MAP
rf_data_AI <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, CORG_msp, pro_AI) 

rf_data_AI %>% 
  skimr::skim()

# Set-up random forest
task_rf_AI <- as_task_regr(x = rf_data_AI %>% 
                             dplyr::select(-id, -entry_name), 
                           target = "lyr_14c_msp")

lrn_rf_AI <- lrn("regr.ranger", importance = "permutation",
                 num.trees = 1000)

# Train model
set.seed(42)
rf_model_AI <- lrn_rf_AI$train(task_rf_AI)

rf_model_AI$model

rf_AI_pred_df <- data.frame(lyr_14c_pred = rf_model_AI$model$predictions)

rf_AI_data_pred <- cbind(mspline_14c_c_all, rf_AI_pred_df) %>% 
  tibble()

climate_all_AI <- rf_AI_data_pred %>%
  group_by(ClimateZoneAnd, UD) %>% 
  mutate(median_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_14c = wilcox.test(lyr_14c_pred, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
         lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
         uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
         n = n()) %>% 
  distinct(median_14c, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  mutate(n_rel = n * 100 / max(n))

climate_all_AI$ClimateZoneAnd <- factor(climate_all_AI$ClimateZoneAnd,
                                         levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                    "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_AI %>% 
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


c1_AI <- climate_all_AI %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
        axis.title = element_text(face = "bold"),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon [wt-%]", trans = "log10", 
                     limits = c(0.1,10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")), expand = c(0,0),
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_AI)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_SOC_climate_AI_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

### Predict SOC and 14C

# Filter data (only keep variables of interest)
rf_data_14c <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 

rf_data_c <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, CORG_msp,
                pro_MAT_mod, pro_GPP_Fluxcom_2001_2012_gC_m2d1,
                lyr_clay_mod, pro_AI) 


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
rf_pred_14c_df <- data.frame(lyr_14c_pred = rf_model_14c$model$predictions)

rf_data_pred_14c <- cbind(mspline_14c_c_all, rf_pred_14c_df) %>% 
  tibble()

rf_pred_c_df <- data.frame(lyr_c_pred = rf_model_c$model$predictions)

rf_data_pred_c <- cbind(mspline_14c_c_all, rf_pred_c_df) %>% 
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
ggsave(file = paste0("./Figure/ISRaD_msp_14CRF_pred_obs_", Sys.Date(),
                     ".jpeg"), width = 6, height = 6)

rf_data_pred_c %>% 
  ggplot(aes(x = lyr_c_pred, y = CORG_msp)) +
  geom_point(shape = 21) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous("Observed C", limits = c(0,55), expand = c(0,0)) +
  scale_x_continuous("Predicted C", limits = c(0,55), expand = c(0,0)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0("./Figure/ISRaD_msp_CORGRF_pred_obs_", Sys.Date(),
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

rf_data_pred_14c_c <- cbind(mspline_14c_c_all, rf_pred_14c_df, rf_pred_c_df) %>% 
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
                                           levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                      "warm temperate", "arid", "tropical"))

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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_pred_SOC_climate_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

## Look at individual predictors
# depth only
rf_data_14c_depth <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp) 

rf_data_c_depth <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, CORG_msp)

# Set-up random forest
task_rf_14c_depth <- as_task_regr(x = rf_data_14c_depth %>% 
                                    dplyr::select(-id, -entry_name), 
                                  target = "lyr_14c_msp")

lrn_rf_14c_depth <- lrn("regr.ranger", importance = "permutation",
                        num.trees = 1000)

task_rf_c_depth <- as_task_regr(x = rf_data_c_depth %>% 
                                  dplyr::select(-id, -entry_name), 
                                target = "CORG_msp")

lrn_rf_c_depth <- lrn("regr.ranger", importance = "permutation",
                      num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c_depth <- lrn_rf_14c_depth$train(task_rf_14c_depth)

rf_model_14c_depth$model

set.seed(42)
rf_model_c_depth <- lrn_rf_c_depth$train(task_rf_c_depth)

rf_model_c_depth$model

## Check model
rf_pred_14c_df_depth <- data.frame(lyr_14c_pred = rf_model_14c_depth$model$predictions)

rf_pred_c_df_depth <- data.frame(lyr_c_pred = rf_model_c_depth$model$predictions)

rf_data_pred_14c_c_depth <- cbind(mspline_14c_c_all, rf_pred_14c_df_depth, rf_pred_c_df_depth) %>% 
  tibble()

climate_all_14c_c_depth <- rf_data_pred_14c_c_depth %>%
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

climate_all_14c_c_depth$ClimateZoneAnd <- factor(climate_all_14c_c_depth$ClimateZoneAnd,
                                                 levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                            "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_14c_c_depth %>% 
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


c1_14c_c_depth <- climate_all_14c_c_depth %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c_depth)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_climate_depth_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

# MAT
rf_data_14c_MAT <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, pro_MAT_mod) 

rf_data_c_MAT <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, CORG_msp, pro_MAT_mod)

# Set-up random forest
task_rf_14c_MAT <- as_task_regr(x = rf_data_14c_MAT %>% 
                                  dplyr::select(-id, -entry_name), 
                                target = "lyr_14c_msp")

lrn_rf_14c_MAT <- lrn("regr.ranger", importance = "permutation",
                      num.trees = 1000)

task_rf_c_MAT <- as_task_regr(x = rf_data_c_MAT %>% 
                                dplyr::select(-id, -entry_name), 
                              target = "CORG_msp")

lrn_rf_c_MAT <- lrn("regr.ranger", importance = "permutation",
                    num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c_MAT <- lrn_rf_14c_MAT$train(task_rf_14c_MAT)

rf_model_14c_MAT$model

set.seed(42)
rf_model_c_MAT <- lrn_rf_c_MAT$train(task_rf_c_MAT)

rf_model_c_MAT$model

## Check model
rf_pred_14c_df_MAT <- data.frame(lyr_14c_pred = rf_model_14c_MAT$model$predictions)

rf_pred_c_df_MAT <- data.frame(lyr_c_pred = rf_model_c_MAT$model$predictions)

rf_data_pred_14c_c_MAT <- cbind(mspline_14c_c_all, rf_pred_14c_df_MAT, rf_pred_c_df_MAT) %>% 
  tibble()

climate_all_14c_c_MAT <- rf_data_pred_14c_c_MAT %>%
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

climate_all_14c_c_MAT$ClimateZoneAnd <- factor(climate_all_14c_c_MAT$ClimateZoneAnd,
                                               levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                          "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_14c_c_MAT %>% 
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


c1_14c_c_MAT <- climate_all_14c_c_MAT %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c_MAT)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_climate_MAT_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

# GPP
rf_data_14c_GPP <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, pro_GPP_Fluxcom_2001_2012_gC_m2d1) 

rf_data_c_GPP <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, CORG_msp, pro_GPP_Fluxcom_2001_2012_gC_m2d1)

# Set-up random forest
task_rf_14c_GPP <- as_task_regr(x = rf_data_14c_GPP %>% 
                                  dplyr::select(-id, -entry_name), 
                                target = "lyr_14c_msp")

lrn_rf_14c_GPP <- lrn("regr.ranger", importance = "permutation",
                      num.trees = 1000)

task_rf_c_GPP <- as_task_regr(x = rf_data_c_GPP %>% 
                                dplyr::select(-id, -entry_name), 
                              target = "CORG_msp")

lrn_rf_c_GPP <- lrn("regr.ranger", importance = "permutation",
                    num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c_GPP <- lrn_rf_14c_GPP$train(task_rf_14c_GPP)

rf_model_14c_GPP$model

set.seed(42)
rf_model_c_GPP <- lrn_rf_c_GPP$train(task_rf_c_GPP)

rf_model_c_GPP$model

## Check model
rf_pred_14c_df_GPP <- data.frame(lyr_14c_pred = rf_model_14c_GPP$model$predictions)

rf_pred_c_df_GPP <- data.frame(lyr_c_pred = rf_model_c_GPP$model$predictions)

rf_data_pred_14c_c_GPP <- cbind(mspline_14c_c_all, rf_pred_14c_df_GPP, rf_pred_c_df_GPP) %>% 
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
                                               levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                          "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_14c_c_GPP %>% 
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


c1_14c_c_GPP <- climate_all_14c_c_GPP %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c_GPP)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_climate_GPP_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)


# PET/MAP
rf_data_14c_AI <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, pro_AI) 

rf_data_c_AI <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, CORG_msp, pro_AI)

# Set-up random forest
task_rf_14c_AI <- as_task_regr(x = rf_data_14c_AI %>% 
                                 dplyr::select(-id, -entry_name), 
                               target = "lyr_14c_msp")

lrn_rf_14c_AI <- lrn("regr.ranger", importance = "permutation",
                      num.trees = 1000)

task_rf_c_AI <- as_task_regr(x = rf_data_c_AI %>% 
                                dplyr::select(-id, -entry_name), 
                              target = "CORG_msp")

lrn_rf_c_AI <- lrn("regr.ranger", importance = "permutation",
                    num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c_AI <- lrn_rf_14c_AI$train(task_rf_14c_AI)

rf_model_14c_AI$model

set.seed(42)
rf_model_c_AI <- lrn_rf_c_AI$train(task_rf_c_AI)

rf_model_c_AI$model

## Check model
rf_pred_14c_df_AI <- data.frame(lyr_14c_pred = rf_model_14c_AI$model$predictions)

rf_pred_c_df_AI <- data.frame(lyr_c_pred = rf_model_c_AI$model$predictions)

rf_data_pred_14c_c_AI <- cbind(mspline_14c_c_all, rf_pred_14c_df_AI, rf_pred_c_df_AI) %>% 
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
                                              levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                         "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_14c_c_AI %>% 
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


c1_14c_c_AI <- climate_all_14c_c_AI %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c_AI)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_climate_AI_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)

# clay
rf_data_14c_clay <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, lyr_14c_msp, lyr_clay_mod) 

rf_data_c_clay <- mspline_14c_c_all %>% 
  dplyr::select(id, entry_name, UD, CORG_msp, lyr_clay_mod)

# Set-up random forest
task_rf_14c_clay <- as_task_regr(x = rf_data_14c_clay %>% 
                                   dplyr::select(-id, -entry_name), 
                                 target = "lyr_14c_msp")

lrn_rf_14c_clay <- lrn("regr.ranger", importance = "permutation",
                       num.trees = 1000)

task_rf_c_clay <- as_task_regr(x = rf_data_c_clay %>% 
                                 dplyr::select(-id, -entry_name), 
                               target = "CORG_msp")

lrn_rf_c_clay <- lrn("regr.ranger", importance = "permutation",
                     num.trees = 1000)

# Train model
set.seed(42)
rf_model_14c_clay <- lrn_rf_14c_clay$train(task_rf_14c_clay)

rf_model_14c_clay$model

set.seed(42)
rf_model_c_clay <- lrn_rf_c_clay$train(task_rf_c_clay)

rf_model_c_clay$model

## Check model
rf_pred_14c_df_clay <- data.frame(lyr_14c_pred = rf_model_14c_clay$model$predictions)

rf_pred_c_df_clay <- data.frame(lyr_c_pred = rf_model_c_clay$model$predictions)

rf_data_pred_14c_c_clay <- cbind(mspline_14c_c_all, rf_pred_14c_df_clay, rf_pred_c_df_clay) %>% 
  tibble()

climate_all_14c_c_clay <- rf_data_pred_14c_c_clay %>%
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

climate_all_14c_c_clay$ClimateZoneAnd <- factor(climate_all_14c_c_clay$ClimateZoneAnd,
                                                levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                                           "warm temperate", "arid", "tropical"))

depth_sum <- climate_all_14c_c_clay %>% 
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


c1_14c_c_clay <- climate_all_14c_c_clay %>% 
  filter(ClimateZoneAnd != "volcanic soils") %>% 
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
                     limits = c(-1000,125), breaks = seq(-1000,0,250)) 
plot(c1_14c_c_clay)
ggsave(file = paste0("./Figure/ISRaD_msp_pred14C_predSOC_climate_clay_", Sys.Date(),
                     ".jpeg"), width = 6, height = 5.5)
