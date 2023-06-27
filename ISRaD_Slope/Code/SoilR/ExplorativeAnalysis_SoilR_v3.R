# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 10/02/2023 #

## Code and functions for modelling (SoilR) provided by Carlos Sierra

# devtools::install_github('MPIBGC-TEE/SoilR-exp/pkg')

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(SoilR)

### Load filtered and splined lyr data
lyr_msp_data <- read_csv("./Data/ISRaD_flat_splined_filled_2023-02-08.csv") %>%
  dplyr::select(entry_name, pro_name, pro_long, pro_lat, lyr_name, id,
                lyr_14c_msp, CORG_msp, pro_AI, pro_MAT_mod,
                pro_GPP_Fluxcom_2001_2012_gC_m2d1, lyr_clay_mod,
                ClimateZoneAnd, MineralType, UD, LD, lyr_obs_date_y)

skimr::skim_without_charts(lyr_msp_data)

lyr_msp_data %>%
  distinct(id, .keep_all = TRUE) %>%
  count(ClimateZoneAnd)

lyr_msp_data %>%
  distinct(id, .keep_all = TRUE) %>%
  count(MineralType)

lyr_msp_data %>%
  skimr::skim_without_charts(UD, lyr_14c_msp, CORG_msp, pro_MAT_mod, 
                             pro_GPP_Fluxcom_2001_2012_gC_m2d1, lyr_clay_mod, pro_AI)

### Set-up a depth-resolved model with SoilR (Carlos) ###
source("C:/Users/sfromm/Documents/GitHub/SoilR-exp/prototypes/VerticalTransfer/VTLM.R")
source("C:/Users/sfromm/Documents/GitHub/SoilR-exp/prototypes/VerticalTransfer/collect_by_layer.R")
source("C:/Users/sfromm/Documents/GitHub/SoilR-exp/prototypes/VerticalTransfer/collect_by_pool.R")
source("C:/Users/sfromm/Documents/GitHub/SoilR-exp/prototypes/VerticalTransfer/solveVTLM.R")

#### One-pool-model ####

### Different k-values

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
years <- seq(-48050,2019, by = 1)
# years <- seq(1941.5,2019, by = 0.5)
LitterInput <- 1  
# k1 <- 1/50
initialC <- 0 # You can put any number here. They will be ignored because steady-state will be assumed
initialF14 <- 0 # Same here. For the depth simulation the age of the steady-state carbon will be used to set initial F14C
AtmF14 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2021$NHZone2,
                         time.scale = "AD")[,1:2]
# AtmF14 <- Hua2021$NHZone2[,1:2]

k1_list <- list(k1.1 = 1/5,
                k1.2 = 1/10,
                k1.3 = 1/50,
                k1.4 = 1/100,
                k1.5 = 1/200,
                k1.6 = 1/500,
                k1.7 = 1/1000,
                k1.8 = 1/2000)

# Parameters for VTLM
depthLayers <- c(1:10)
rootInputs <- rep(0, times = length(depthLayers))
downwardTransferRate <- 0.005
upwardTransferRate <- 0

#if d is too small, e.g. 0.003:
#Error in solve.default(X) : 
  # system is computationally singular: reciprocal condition number = 9.1434e-19

## Step 2. Create SoilR pool model & Create Vertical Transport model
OPM_VTLM_k <- 
  map(k1_list, \(x) OnepModel14(
    t = years, 
    k = x, 
    C0 = initialC,
    F0_Delta14C = initialF14, 
    In = LitterInput,
    inputFc = AtmF14
  )) %>% 
  map(\(x) VTLM(Model = x, 
                lyrs = depthLayers, 
                latIn = rootInputs, 
                d = downwardTransferRate, 
                u = upwardTransferRate,
                vrm = 1,
                ivalF14 = "model"))

## Step 3. Solve the vertical transport model
Cdt_k <- OPM_VTLM_k %>% 
  map(\(x) getC(x))

C14dt_k <- OPM_VTLM_k %>% 
  map(\(x) getF14(x))

Cdt_2019_k <- Cdt_k %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "k_id") %>% 
  pivot_longer(!k_id, names_to = "layer", values_to = "c_stock") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2019_k <- C14dt_k %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "k_id") %>% 
  pivot_longer(!k_id, names_to = "layer", values_to = "lyr_14c") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2010_k <- C14dt_k %>%
  map(\(x) x[50061,]) %>%
  plyr::ldply(.id = "k_id") %>%
  pivot_longer(!k_id, names_to = "layer_v", values_to = "lyr_14c") %>%
  mutate(layer = rep(depthLayers, length(k1_list)))

k1 <- plyr::ldply(k1_list, .id = "k_id", data.frame) %>% 
  rename(k_value = X..i..)

OPM_output_k <- tibble(
  year = round(rep(tail(AtmF14$Year, 1), max(depthLayers)*length(k1_list)), 0),
  layer = rep(depthLayers, length(k1_list)),
  LitterInput = LitterInput,
  downwardTransferRate = downwardTransferRate,
  upwardTransferRate = upwardTransferRate,
  k_id = rep(k1$k_id, each = length(depthLayers)),
  k_value = round(rep(k1$k_value, each = length(depthLayers)), 4)
)

OPM_results_k <- OPM_output_k %>% 
  left_join(C14dt_2019_k) %>% 
  left_join(Cdt_2019_k) %>% 
  mutate(k_label = case_when(
    k_value == 1/5 ~ "1/5",
    k_value == 1/10 ~ "1/10",
    k_value == 1/50 ~ "1/50",
    k_value == 1/100 ~ "1/100",
    k_value == 1/200 ~ "1/200",
    k_value == 1/500 ~ "1/500",
    k_value == 1/1000 ~ "1/1000",
    k_value == 1/2000 ~ "1/2000"
  ))

OnePool_const_k <- OPM_results_k %>% 
  ggplot(aes(x = c_stock, y = lyr_14c)) +
  geom_path(aes(color = fct_reorder(k_label, -k_value)), linewidth = 2) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black", 
             color = "white") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.box = "horizontal") +
  scale_x_continuous("", limits = c(1e-14,5e+02), trans = "log10",
                     breaks = c(1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-08,
                                1e-07,1e-06,1e-05,1e-04,0.001,0.01,0.1,1,10,100), 
                     expand = c(0,0),
                     labels = c(1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-08,
                                1e-07,1e-06,1e-05,1e-04,0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")), limits = c(-250,300),
                     expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21),
                     labels = c("1", "2-10", "", "", "", "", "", "", "", "")) +
  scale_color_manual("k-value",
                     values = c("#f6e8c3", "#dfc27d", "#bf812d", "#8c510a", "#543005", 
                                "#7b3294", "#018571", "#80cdc1")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none")

# ggsave(file = paste0("./Figure/OnePoolModel_k_values_2019_", Sys.Date(), ".jpeg"), 
#        width = 12, height = 6)


### Different decomposition rates with depth

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
k1 <- 1/500
LitterInput <- 1
rootInputs <- rep(0, times = length(depthLayers))

# Parameters for VTLM
downwardTransferRate <- 0.005
upwardTransferRate <- 0

vkm_list <- list(vkm.0 = 1,
                 vkm.1 = seq(from  = 1, by = -1/length(depthLayers), 
                             length.out = length(depthLayers)),
                 vkm.2 = exp(log(1):(log(1)-(length(depthLayers)-1))/3),
                 vkm.3 = exp(log(1):(log(1)-(length(depthLayers)-1))))

## Step 2. Create SoilR pool model & Create Vertical Transport model
OPM <- OnepModel14(
  t = years, 
  k = k1, 
  C0 = initialC,
  F0_Delta14C = initialF14, 
  In = LitterInput,
  inputFc = AtmF14
)

OPM_VTLM_vrm <- vkm_list %>% 
  map(\(x) VTLM(Model = OPM, 
                lyrs = depthLayers, 
                latIn = rootInputs, 
                d = downwardTransferRate, 
                u = upwardTransferRate,
                vrm = x,
                ivalF14 = "model"))

## Step 3. Solve the vertical transport model
Cdt_vrm <- OPM_VTLM_vrm %>% 
  map(\(x) getC(x))

C14dt_vrm <- OPM_VTLM_vrm %>% 
  map(\(x) getF14(x))

Cdt_2019_vrm <- Cdt_vrm %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "vkm_id") %>% 
  pivot_longer(!vkm_id, names_to = "layer", values_to = "c_stock") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2019_vrm <- C14dt_vrm %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "vkm_id") %>% 
  pivot_longer(!vkm_id, names_to = "layer", values_to = "lyr_14c") %>% 
  mutate(layer = as.numeric(layer))

vrm <- plyr::ldply(vkm_list, .id = "vkm_id", data.frame) %>% 
  rename(vkm_value = X..i..) %>% 
  distinct(vkm_id, .keep_all = TRUE)

OPM_output_vrm <- tibble(
  year = round(rep(tail(AtmF14$Year, 1), max(depthLayers)*length(vkm_list)), 0),
  layer = rep(depthLayers, length(vkm_list)),
  LitterInput = LitterInput,
  vkm_id = rep(vrm$vkm_id, each = length(depthLayers)),
  vkm_value = rep(vrm$vkm_value, each = length(depthLayers)),
  upwardTransferRate = upwardTransferRate,
  downwardTransferRate = downwardTransferRate,
  k = k1
)

OPM_results_vrm <- OPM_output_vrm %>% 
  left_join(C14dt_2019_vrm, multiple = "all") %>% 
  left_join(Cdt_2019_vrm, multiple = "all")

OnePool_vkm <- OPM_results_vrm %>% 
  ggplot(aes(x = c_stock, y = lyr_14c)) +
  geom_path(aes(color = as.factor(vkm_id)), linewidth = 2) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.box = "horizontal") +
  scale_x_continuous("Soil organic carbon stocks; log-scale", limits = c(1,5e+02),
                     breaks = c(1,10,100), expand = c(0,0),
                     labels = c(1,10,100), trans = "log10",) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")), limits = c(-250,300),
                     expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21),
                     labels = c("1", "2-10", "", "", "", "", "", "", "", "")) +
  scale_color_manual("vertical k modifier", 
                     values = c("#7b3294", "#0570b0", "#74a9cf", "#bdc9e1")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none")

# ggsave(file = paste0("./Figure/OnePoolModel_vrm_values_2019_", Sys.Date(), ".jpeg"), 
#        width = 12, height = 6)

ggarrange(OnePool_const_k, OnePool_vkm, nrow = 2)
ggsave(file = paste0("./Figure/OnePoolModel_k_vkm_values_2019_", Sys.Date(), ".jpeg"),
       width = 10, height = 6)
  
### Different downward transfer

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
k1 <- 1/500
LitterInput <- 1

# Parameters for VTLM
d_list <- list(d.1 = 5,
               d.2 = 1,
               d.3 = 0.5,
               d.4 = 0.1,
               d.5 = 0.01,
               d.6 = 0.005,
               d.7 = 0.001,
               d.8 = 0.0005,
               d.9 = 0.0001)

## Step 2. Create SoilR pool model & Create Vertical Transport model
OPM <- OnepModel14(
  t = years, 
  k = k1, 
  C0 = initialC,
  F0_Delta14C = initialF14, 
  In = LitterInput,
  inputFc = AtmF14
)

OPM_VTLM_d <- d_list %>% 
  map(\(x) VTLM(Model = OPM, 
                lyrs = depthLayers, 
                latIn = rootInputs, 
                d = x, 
                u = upwardTransferRate,
                vrm = 1,
                ivalF14 = "model"))

## Step 3. Solve the vertical transport model
Cdt_d <- OPM_VTLM_d %>% 
  map(\(x) getC(x))

C14dt_d <- OPM_VTLM_d %>% 
  map(\(x) getF14(x))

Cdt_2019_d <- Cdt_d %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "dtr_id") %>% 
  pivot_longer(!dtr_id, names_to = "layer", values_to = "c_stock") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2019_d <- C14dt_d %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "dtr_id") %>% 
  pivot_longer(!dtr_id, names_to = "layer", values_to = "lyr_14c") %>% 
  mutate(layer = as.numeric(layer))

d <- plyr::ldply(d_list, .id = "dtr_id", data.frame) %>% 
  rename(dtr_value = X..i..)

OPM_output_d <- tibble(
  year = round(rep(tail(AtmF14$Year, 1), max(depthLayers)*length(d_list)), 0),
  layer = rep(depthLayers, length(d_list)),
  LitterInput = LitterInput,
  dtr_id = rep(d$dtr_id, each = length(depthLayers)),
  dtr_value = rep(d$dtr_value, each = length(depthLayers)),
  upwardTransferRate = upwardTransferRate,
  k = k1
)

OPM_results_d <- OPM_output_d %>% 
  left_join(C14dt_2019_d) %>% 
  left_join(Cdt_2019_d)

OPM_results_d %>% 
  ggplot(aes(x = c_stock, y = lyr_14c, label = dtr_value)) +
  geom_path(aes(color = as.factor(dtr_value)), linewidth = 2) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.box = "horizontal") +
  scale_x_continuous("Soil organic carbon stocks; log-scale", 
                     limits = c(1e-10,6e+02), trans = "log10",
                     breaks = c(1e-10,1e-9,1e-08,1e-07,1e-06,1e-05,1e-04,
                                0.001,0.01,0.1,1,10,100), 
                     expand = c(0,0),
                     labels = c(1e-10,1e-9,1e-08,1e-07,1e-06,1e-05,1e-04,
                                0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")), limits = c(-500,300),
                     expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_color_manual("vertical transport",
                     values = c("#c7eae5", "#5ab4ac", "#01665e", "#7b3294", 
                                "#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none")

ggsave(file = paste0("./Figure/OnePoolModel_d_values_2019_", Sys.Date(), ".jpeg"), 
       width = 12, height = 6)

### Different litter inputs

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
k1 <- 1/500

downwardTransferRate <- 0.005

LI_list <- list(LI.2 = 0.2,
                LI.3 = 0.5,
                LI.4 = 1,
                LI.5 = 2,
                LI.6 = 5)

## Step 2. Create SoilR pool model & Create Vertical Transport model
OPM_VTLM_LI <- 
  map(LI_list, \(x) OnepModel14(
    t = years, 
    k = k1, 
    C0 = initialC,
    F0_Delta14C = initialF14, 
    In = x,
    inputFc = AtmF14
  )) %>% 
  map(\(x) VTLM(Model = x, 
                lyrs = depthLayers, 
                latIn = rootInputs, 
                d = downwardTransferRate, 
                u = upwardTransferRate,
                vrm = 1,
                ivalF14 = "model"))

## Step 3. Solve the vertical transport model
Cdt_LI <- OPM_VTLM_LI %>% 
  map(\(x) getC(x))

C14dt_LI <- OPM_VTLM_LI %>% 
  map(\(x) getF14(x))

Cdt_2019_LI <- Cdt_LI %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "LI_id") %>% 
  pivot_longer(!LI_id, names_to = "layer", values_to = "c_stock") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2019_LI <- C14dt_LI %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "LI_id") %>% 
  pivot_longer(!LI_id, names_to = "layer", values_to = "lyr_14c") %>% 
  mutate(layer = as.numeric(layer))

LI <- plyr::ldply(LI_list, .id = "LI_id", data.frame) %>% 
  rename(LI_value = X..i..)

OPM_output_LI <- tibble(
  year = round(rep(tail(AtmF14$Year, 1), max(depthLayers)*length(LI_list)), 0),
  layer = rep(depthLayers, length(LI_list)),
  LitterInput = LitterInput,
  k = k1,
  downwardTransferRate = downwardTransferRate,
  upwardTransferRate = upwardTransferRate,
  LI_id = rep(LI$LI_id, each = length(depthLayers)),
  LI_value = rep(LI$LI_value, each = length(depthLayers))
)

OPM_results_LI <- OPM_output_LI %>% 
  left_join(C14dt_2019_LI) %>% 
  left_join(Cdt_2019_LI)

OnePool_Litter <- OPM_results_LI %>% 
  ggplot(aes(x = c_stock, y = lyr_14c)) +
  geom_path(aes(color = as.factor(LI_value)), linewidth = 2) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        legend.box = "vertical",
        plot.margin = margin(r = 15, t = 10, l = 10, b = 10)) +
  scale_x_continuous("Soil organic carbon stocks; log-scale", limits = c(1,1e+03),
                     breaks = c(1,10,100,1000), expand = c(0,0),
                     labels = c(1,10,100,1000), trans = "log10") +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")), limits = c(-250,100),
                     expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_color_manual("aboveground quantities", 
                     values = c("#dfc27d", "#a6611a", "#7b3294", 
                                "#018571", "#80cdc1")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none")
# ggsave(file = paste0("./Figure/OnePoolModel_LI_values_2019_", Sys.Date(), ".jpeg"), 
#        width = 12, height = 6)

### Different lateral C inputs

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for VTLM
k1 <- 1/500
LitterInput <- 1

#Root inputs all have the same sum ~ 1.58
RI_list <- list(RI.0 = rep(0, times = length(depthLayers)),
                RI.const = rep(0.158, times = length(depthLayers)),
                RI.lin = seq(from  = 0.288, by = -0.288/length(depthLayers), 
                                 length.out = length(depthLayers)),
                RI.exp.1 = exp((log(0.39):(log(0.39)-(length(depthLayers)-1)))/2),
                RI.exp.2 = exp(log(1):(log(1)-(length(depthLayers)-1))))

sum(rep(0.158, times = length(depthLayers)))
sum(seq(from  = 0.288, by = -0.288/length(depthLayers), 
        length.out = length(depthLayers)))
sum(exp((log(0.39):(log(0.39)-(length(depthLayers)-1)))/2))
sum(exp(log(1):(log(1)-(length(depthLayers)-1))))


plot(x = depthLayers, y = RI_list$RI.const, type = "l")
plot(x = depthLayers, y = RI_list$RI.lin, type = "l")
plot(x = depthLayers, y = RI_list$RI.exp.1, type = "l")
plot(x = depthLayers, y = RI_list$RI.exp.2, type = "l")


# RI_list <- list(RI.0 = rep(0, times = length(depthLayers)),
#                 RI.250_lind = seq(from = 250, to = 1, 
#                                   length.out = length(depthLayers)),
#                 RI_250_expd = exp(log(250):(log(250)-(length(depthLayers)-1))))

## Step 2. Create SoilR pool model & Create Vertical Transport model
OPM <- OnepModel14(
  t = years, 
  k = k1, 
  C0 = initialC,
  F0_Delta14C = initialF14, 
  In = LitterInput,
  inputFc = AtmF14
)

OPM_VTLM_RI <- RI_list %>% 
  map(\(x) VTLM(Model = OPM, 
                lyrs = depthLayers, 
                latIn = x, 
                d = downwardTransferRate, 
                u = upwardTransferRate,
                vrm = 1,
                ivalF14 = "model"))

## Step 3. Solve the vertical transport model
Cdt_RI <- OPM_VTLM_RI %>% 
  map(\(x) getC(x))

C14dt_RI <- OPM_VTLM_RI %>% 
  map(\(x) getF14(x))

Cdt_2019_RI <- Cdt_RI %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "RI_id") %>% 
  pivot_longer(!RI_id, names_to = "layer", values_to = "c_stock") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2019_RI <- C14dt_RI %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "RI_id") %>% 
  pivot_longer(!RI_id, names_to = "layer", values_to = "lyr_14c") %>% 
  mutate(layer = as.numeric(layer))

RI <- plyr::ldply(RI_list, .id = "RI_id", data.frame) %>% 
  rename(RI_value = X..i..) %>% 
  distinct(RI_id, .keep_all = TRUE)

OPM_output_RI <- tibble(
  year = round(rep(tail(AtmF14$Year, 1), max(depthLayers)*length(RI_list)), 0),
  layer = rep(depthLayers, length(RI_list)),
  RI_id = rep(RI$RI_id, each = length(depthLayers)),
  RI_value = rep(RI$RI_value, each = length(depthLayers)),
  downwardTransferRate = downwardTransferRate,
  upwardTransferRate = upwardTransferRate,
  k = k1,
  LitterInput = LitterInput
)

OPM_results_RI <- OPM_output_RI %>% 
  left_join(C14dt_2019_RI) %>% 
  left_join(Cdt_2019_RI)

OnePool_Roots <- OPM_results_RI %>% 
  ggplot(aes(x = c_stock, y = lyr_14c)) +
  geom_path(aes(color = as.factor(RI_id)), linewidth = 2) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        legend.box = "vertical",
        plot.margin = margin(r = 15, t = 10, l = 0, b = 10),
        axis.title.y = element_blank()) +
  scale_x_continuous("Soil organic carbon stocks; log-scale", limits = c(1,1e+03),
                     breaks = c(1,10,100,1000), expand = c(0,0),
                     labels = c(1,10,100,1000), trans = "log10") +
  scale_y_continuous("", limits = c(-250,100), expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_color_manual("belowground distributions", 
                     values = c("#7b3294", "#006837", "#31a354", 
                                "#78c679", "#c2e699")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none")
# ggsave(file = paste0("./Figure/OnePoolModel_RI_values_2019_", Sys.Date(), ".jpeg"), 
#        width = 12, height = 6)

ggarrange(OnePool_Litter, OnePool_Roots, nrow = 1, legend = "top")
ggsave(file = paste0("./Figure/OnePoolModel_litter_roots_2019_", Sys.Date(), ".jpeg"),
       width = 14, height = 6)

### Combine different model runs in one figure
## Inputs: different root inputs with same litter inputs 
years <- seq(-48050,2019, by = 1)
initialC <- 0 # You can put any number here. They will be ignored because steady-state will be assumed
initialF14 <- 0 # Same here. For the depth simulation the age of the steady-state carbon will be used to set initial F14C
AtmF14 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2021$NHZone2,
                         time.scale = "AD")[,1:2]
depthLayers <- c(1:10)
downwardTransferRate <- 0.005
upwardTransferRate <- 0
k1 <- 1/500
LitterInputs = 1

vrm.1 <- c(1.0, 0.8, 0.7, 0.6, 0.5, 0.45, 0.35, 0.25, 0.2, 0.1)

# Set up constant model values
RI_list <- list(RI.const = rep(0.158, times = length(depthLayers)),
                RI.lin = seq(from  = 0.288, by = -0.288/length(depthLayers),
                             length.out = length(depthLayers)),
                RI.exp = exp(log(1):(-(length(depthLayers)-1))))


# data.frame(RI_list, layer = depthLayers) %>% 
#   pivot_longer(!layer, names_to = "RootDist", values_to = "RootInputs") %>% 
#   dplyr::mutate(RootDist = factor(RootDist,
#                                   levels = c("RI.const", "RI.lin", "RI.exp"))) %>%
#   ggplot(aes(y = layer, x = RootInputs, color = RootDist)) +
#   geom_path(linewidth = 1.5) +
#   theme_bw(base_size = 16) +
#   theme(axis.text = element_text(color = "black")) +
#   scale_y_reverse("Depth layers", limits = c(10,1), expand = c(0,0),
#                   breaks = seq(1,10,1)) +
#   scale_x_continuous("Root inputs", limits = c(0,1), expand = c(0.01,0.01)) + 
#   scale_color_manual("Input distribution", 
#                      values = c("#fee391", "#fe9929", "#993404"))
# 
# ggsave(file = paste0("./Figure/OnePoolModel_RootInputs_Dist_", 
#                      Sys.Date(), ".jpeg"), width = 12, height = 6)

###Litter inputs = 1
## Step 2. Create SoilR pool model & Create Vertical Transport model
OPM_LI.1 <- OnepModel14(
  t = years, 
  k = k1, 
  C0 = initialC,
  F0_Delta14C = initialF14, 
  In = LitterInputs,
  inputFc = AtmF14
)

OPM_VTLM_RI_LI.1 <- RI_list %>% 
  map(\(x) VTLM(Model = OPM_LI.1, 
                lyrs = depthLayers, 
                latIn = x, 
                d = downwardTransferRate, 
                u = upwardTransferRate, 
                vrm = vrm.1,
                ivalF14 = "model"))

## Step 3. Solve the vertical transport model
Cdt_RI_LI.1 <- OPM_VTLM_RI_LI.1 %>% 
  map(\(x) getC(x))

C14dt_RI_LI.1 <- OPM_VTLM_RI_LI.1 %>% 
  map(\(x) getF14(x))

Cdt_2019_RI_LI.1 <- Cdt_RI_LI.1 %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "RI_id") %>% 
  pivot_longer(!RI_id, names_to = "layer", values_to = "c_stock") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2019_RI_LI.1 <- C14dt_RI_LI.1 %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "RI_id") %>% 
  pivot_longer(!RI_id, names_to = "layer", values_to = "lyr_14c") %>% 
  mutate(layer = as.numeric(layer))

OPM_results_RI_LI.1 <- C14dt_2019_RI_LI.1 %>% 
  left_join(Cdt_2019_RI_LI.1) %>% 
  mutate(LitterInput = 1)

OnePool_roots_sum <- OPM_results_RI_LI.1 %>% 
  ggplot(aes(x = c_stock, y = lyr_14c)) +
  geom_path(aes(color = as.factor(RI_id)), linewidth = 1.5) + 
  geom_path(aes(color = as.factor(RI_id)), linewidth = 1.5,
            arrow = arrow(type = "open"), show.legend = FALSE) + 
  geom_point(aes(shape = as.factor(layer)), size = 1.5, fill = "black",
             color = "white") +
  scale_color_manual("Root inputs", values = c("#fee391", "#fe9929", "#993404")) +
  scale_fill_discrete("Layer") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.15,0.7),
        plot.margin = margin(r = 15, t = 15, l = 15)) +
  scale_x_continuous("", expand = c(0,0),
                     breaks = c(50,100,500),
                     labels = c(50,100,500)) +
  coord_trans(x = "log10", xlim = c(30,501)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")),limits = c(-250,100),
                     expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  guides(shape = "none")

# ggsave(file = paste0("./Figure/OnePoolModel_Root_VRM_values_2019_", 
#                      Sys.Date(), ".jpeg"), width = 12, height = 6)

#Variable transport instead of variable root inputs
d_list <- list(d.1 = 0.01,
               d.2 = 0.005,
               d.3 = 0.003)

RootInputs <- exp(log(1):(-(length(depthLayers)-1)))
vrm.2 <- c(1.0, 0.9, 0.8, 0.7, 0.65, 0.6, 0.5, 0.4, 0.3, 0.2)

## Step 2. Create SoilR pool model & Create Vertical Transport model
OPM_VT <- OnepModel14(
  t = years, 
  k = k1, 
  C0 = initialC,
  F0_Delta14C = initialF14, 
  In = LitterInputs,
  inputFc = AtmF14
)

OPM_VTLM_VT <- d_list %>% 
  map(\(x) VTLM(Model = OPM_VT, 
                lyrs = depthLayers, 
                latIn = RootInputs, 
                d = x, 
                u = upwardTransferRate, 
                vrm = vrm.2, 
                ivalF14 = "model"))

## Step 3. Solve the vertical transport model
Cdt_VT <- OPM_VTLM_VT %>% 
  map(\(x) getC(x))

C14dt_VT <- OPM_VTLM_VT %>% 
  map(\(x) getF14(x))

Cdt_2019_VT <- Cdt_VT %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "d_id") %>% 
  pivot_longer(!d_id, names_to = "layer", values_to = "c_stock") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2019_VT <- C14dt_VT %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "d_id") %>% 
  pivot_longer(!d_id, names_to = "layer", values_to = "lyr_14c") %>% 
  mutate(layer = as.numeric(layer))

OPM_results_VT <- C14dt_2019_VT %>% 
  left_join(Cdt_2019_VT) 

OnePool_VT_sum <- OPM_results_VT %>% 
  ggplot(aes(x = c_stock, y = lyr_14c)) +
  geom_path(aes(color = as.factor(d_id)), linewidth = 1.5) + 
  geom_segment(aes(x = 170, xend = 390, y = 90, yend = 30), color = "#fee391", 
               linewidth = 1, arrow = arrow(type = "open")) +
  geom_segment(aes(x = 86, xend = 40, y = -105, yend = -220), color = "#993404", 
               linewidth = 1, arrow = arrow(type = "open")) +
  geom_point(aes(shape = as.factor(layer)), size = 1.5, fill = "black",
             color = "white") +
  scale_color_manual("Vertical transport", values = c("#2c7fb8", "#993404", "#7fcdbb")) +
  scale_fill_discrete("Layer") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.15,0.7),
        plot.margin = margin(r = 15, l = 15, b = 15, t = 5)) +
  scale_x_continuous("Soil organic carbon stocks; log-scale", expand = c(0,0),
                     breaks = c(50,100,500),
                     labels = c(50,100,500)) +
  coord_trans(x = "log10", xlim = c(30,501)) +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")),limits = c(-250,100),
                     expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  guides(shape = "none")

ggarrange(OnePool_roots_sum, OnePool_VT_sum, nrow = 2)
ggsave(file = paste0("./Figure/OnePoolModel_VT_values_2019_", 
                     Sys.Date(), ".jpeg"), width = 7, height = 7)

## Different approach
# OPM_LI_RI <- OPM_results_LI %>% 
#   rename(LI_lyr_14c = lyr_14c,
#          LI_c_stock = c_stock) %>% 
#   cbind(OPM_results_RI %>% 
#           rename(RI_lyr_14c = lyr_14c,
#                  RI_c_stock = c_stock) %>% 
#           dplyr::select(RI_id, RI_value, RI_lyr_14c, RI_c_stock)) %>% 
#   tibble()
# 
# library(ggnewscale)
# OPM_LI_RI %>% 
#   ggplot() +
#   geom_path(aes(x = LI_c_stock, y = LI_lyr_14c, color = as.factor(LI_value)),
#             linewidth = 1.5) +
#   scale_color_manual("Litter inputs", 
#                      values = c("#dfc27d", "#bf812d", "#762a83", "#8c510a", "#543005")) +
#   new_scale_color() +
#   geom_path(aes(x = RI_c_stock, y = RI_lyr_14c, color = as.factor(RI_id)),
#             linewidth = 1.5, linetype = "dashed") +
#   scale_color_manual("Root inputs", 
#                      values = c("#762a83", "#003c30", "#01665e", "#35978f", "#80cdc1")) +
#   theme_bw(base_size = 16) +
#   theme(axis.text = element_text(color = "black"),
#         panel.grid.minor.x = element_blank()) +
#   scale_x_continuous("C stocks; log-scale", expand = c(0,0),
#                      breaks = c(1,5,10,50,100,500)) +
#   coord_trans(x = "log10", xlim = c(1,800)) +
#   annotation_logticks(sides = "b", scaled = FALSE,
#                       short = unit(1.5,"mm"),
#                       mid = unit(3,"mm"),
#                       long = unit(4,"mm")) +
#   scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-250,75),
#                      expand = c(0,0))
# ggsave(file = paste0("./Figure/OnePoolModel_LitterRoot_values_2019_", 
#                      Sys.Date(), ".jpeg"), width = 15.36, height = 8.26, units = "in")
# 
# scale_x_continuous("Soil organic carbon [wt-%]; log-scaled", expand = c(0,0), 
#                    breaks = c(0.1,0.5,1,5),
#                    labels = c(0.1,0.5,1,5)) +
#   coord_trans(x = "log10", xlim = c(0.1,6)) +
#   annotation_logticks(sides = "b", scaled = FALSE,
#                       short = unit(1.5,"mm"),
#                       mid = unit(3,"mm"),
#                       long = unit(4,"mm")) +
# 
# OPM_results_LI %>% 
#   ggplot(aes(x = c_stock, y = lyr_14c, label = LI_value)) +
#   geom_path(aes(color = as.factor(LI_value)), linewidth = 1.5) +
#   geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
#   theme_bw(base_size = 14) +
#   theme(axis.text = element_text(color = "black")) +
#   scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0)) +
#   scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-500,300),
#                      expand = c(0,0)) +
#   scale_fill_discrete("Depth layer") +
#   scale_color_viridis_d("Litter inputs") +
#   geom_text(aes(x = 1e-9, y = -140), label = downwardTransferRate) +
#   geom_text(aes(x = 1e-10, y = -160), label = k1) +
#   annotate(geom = "text", x = 1e-11, y = -140, label = "downward transfer = ") +
#   annotate(geom = "text", x = 1e-11, y = -160, label = "k-value = ") +
#   annotate(geom = "text", x = 1e-11, y = -180, label = "no root inputs ") +
#   geom_text(data = OPM_results_LI %>% 
#               arrange(desc(layer)) %>% 
#               distinct(LI_value, .keep_all = TRUE), nudge_y = -15, nudge_x = -0.1) +
#   coord_cartesian(xlim = c(1e-15,1e+03))
# 
# OPM_results_RI %>% 
#   ggplot(aes(x = c_stock, y = lyr_14c, label = RI_id)) +
#   geom_path(aes(color = as.factor(RI_id)), linewidth = 1.5) +
#   geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
#   theme_bw(base_size = 14) +
#   theme(axis.text = element_text(color = "black")) +
#   scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0)) +
#   scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-500,300),
#                      expand = c(0,0)) +
#   scale_fill_discrete("Depth layer") +
#   scale_color_viridis_d("Root inputs") +
#   geom_text(aes(x = 1e-10, y = -180), label = k1) +
#   geom_text(aes(x = 1e-9, y = -160), label = downwardTransferRate) +
#   geom_text(aes(x = 1e-10, y = -140), label = LitterInput) +
#   annotate(geom = "text", x = 1e-11, y = -180, label = "k-value = ") +
#   annotate(geom = "text", x = 1e-11, y = -160, label = "downward transfer = ") +
#   annotate(geom = "text", x = 1e-11, y = -140, label = "C inputs = ") +
#   geom_text(data = OPM_results_RI %>% 
#               arrange(desc(layer)) %>% 
#               distinct(RI_id, .keep_all = TRUE), nudge_y = -10, nudge_x = -0.1) +  
#   # coord_cartesian(xlim = c(1e+00,1e+03), ylim = c(-200,100)) +
#   coord_cartesian(xlim = c(1e-15,1e+03))


### Two-pool-model ####

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
years <- seq(-48050,2019, by = 1)
# years <- seq(1941.5, 2019, by = 0.5)
LitterInputs <- 1  
k1 <- 1/50
k2 <- 1/2000
initialC <- c(0,0) # You can put any number here. They will be ignored because steady-state will be assumed
initialF14 <- c(0,0) # Same here. For the depth simulation the age of the steady-state carbon will be used to set initial F14C
AtmF14 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2021$NHZone2, 
                         time.scale = "AD")[,1:2]
# AtmF14 <- Hua2021$NHZone2[,1:2]

a21_list <- list(a21_5 = 0.05 * k1,
                 a_21_10 = 0.1 * k1,
                 a21_25 = 0.25 * k1,
                 a21_50 = 0.5 * k1,
                 a21_75 = 0.75 * k1)

# Parameters for VTLM
depthLayers <- c(1:10)
rootInputs <- exp((log(0.39):(log(0.39)-(length(depthLayers)-1)))/2)
# rootInputs <- exp(log(1):(log(1)-(length(depthLayers)-1)))
# rootInputs <- exp((log(0.5):(log(0.5)-(length(depthLayers)-1)))/3)

downwardTransferRate <- 0.001
upwardTransferRate <- 0

# vrm.1 <- c(1.00, 1.00, 0.8, 0.7, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)
vrm.1 <- c(1.0, 0.9, 0.8, 0.7, 0.65, 0.6, 0.5, 0.4, 0.3, 0.2)


a21_list <- list(a21_5 = 0.05 * k1,
                 a_21_10 = 0.1 * k1,
                 a21_25 = 0.25 * k1,
                 a21_50 = 0.5 * k1,
                 a21_75 = 0.75 * k1)


#Step 2. Create SoilR pool model
OPM_VTLM_a <- 
  map(a21_list, \(x) TwopSeriesModel14(
    t = years, 
    ks = c(k1, k2), 
    C0 = initialC,
    a21 = x,
    F0_Delta14C = initialF14, 
    In = LitterInputs,
    inputFc = AtmF14
  )) %>% 
  map(\(x) VTLM(Model = x, 
                lyrs = depthLayers, 
                latIn = rootInputs, 
                d = downwardTransferRate, 
                u = upwardTransferRate, 
                vrm = vrm.1,
                ivalF14 = "model"))

# Step 3. Solve the vertical transport model and aggregate by pool and layer
Cdt_a <- OPM_VTLM_a %>% 
  map(\(x) getC(x))

C14dt_a <- OPM_VTLM_a %>% 
  map(\(x) getF14(x))

L14t_a <- map2(Cdt_a, C14dt_a, 
               \(x, y) collect_by_layer(
                 x, 
                 y,
                 nlayer = max(depthLayers), 
                 npool = 2
               ))

L14t_2019_a <- L14t_a %>% 
  unlist(recursive = FALSE) %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "a_id_val") %>% 
  pivot_longer(!a_id_val, names_to = "layer", values_to = "values") %>% 
  mutate(layer = as.numeric(layer)) %>% 
  separate_wider_delim(a_id_val, ".", names = c("a_id", "measure")) %>% 
  pivot_wider(names_from = measure, values_from = values)

a21 <- plyr::ldply(a21_list, .id = "a_id", data.frame) %>% 
  rename(a_value = X..i..)

OPM_output_a <- tibble(
  year = round(rep(tail(AtmF14$Year, 1), max(depthLayers)*length(a21_list)), 0),
  layer = rep(depthLayers, length(a21_list)),
  LitterInput = LitterInputs,
  downwardTransferRate = downwardTransferRate,
  upwardTransferRate = upwardTransferRate,
  a_id = rep(a21$a_id, each = length(depthLayers)),
  a_value = rep(a21$a_value/k1, each = length(depthLayers))
)

OPM_results_a <- OPM_output_a %>% 
  left_join(L14t_2019_a) 

OPM_results_a %>% 
  ggplot(aes(x = C, y = F14C)) +
  geom_path(aes(color = as.factor(a_value*100)), linewidth = 1.5) +
  geom_point(aes(shape = as.factor(layer)), size = 1.5, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank()) +
  scale_x_continuous("Soil organic carbon stocks; log-scale", expand = c(0,0), 
                     limits = c(10,1000),
                     breaks = c(10,50,100,500,1000),
                     labels = c(10,50,100,500,1000)) +
  coord_trans(x = "log10") +
  annotation_logticks(sides = "b", scaled = FALSE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")), limits = c(-550,100),
                     expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_color_viridis_d("transfer to\nslow pool [%]") +
  guides(shape = "none")

ggsave(file = paste0("./Figure/TwoPoolModel_a_values_2019_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

# Solve using a new interface. Useful to see how the individual pools change with depth
targetYear <- which(years == 2019)

pool_results_a_2019 <- OPM_VTLM_a %>% 
  map(\(x) solveVTLM(VTLM = x,
                     npool = 2,
                     nlayer = length(depthLayers))) %>% 
  map_depth(2, \(x) x[,, targetYear]) %>% 
  unlist(recursive = FALSE) %>% 
  plyr::ldply(.id = "a_id_pool") %>% 
  rename(fastPool = "1",
         slowPool = "2") %>% 
  mutate(layer = rep(1:10, length(depthLayers))) %>% 
  separate_wider_delim(a_id_pool, ".", names = c("a_id", "measure")) %>% 
  pivot_longer(cols = c(fastPool, slowPool), names_to = "pool", values_to = "values") %>% 
  pivot_wider(names_from = "measure", values_from = "values")

pool_results_a_2019 %>% 
  left_join(OPM_output_a) %>% 
  ggplot(aes(x = C, y = C14)) +
  geom_path(aes(color = pool), linewidth = 1.5) +
  geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
  facet_wrap(~a_id) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0), 
                     limits = c(1e-02,1e+04)) +
  scale_y_continuous(expression(paste(Delta^14,"C")), limits = c(-550,150),
                     expand = c(0,0))
ggsave(file = paste0("./Figure/TwoPoolModel_a_values_2019_bothPools_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)





