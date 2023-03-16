# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 10/02/2023 #

## Code and functions for modelling (SoilR) provided by Carlos Sierra

devtools::install_github('MPIBGC-TEE/SoilR-exp/pkg')

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
                k1.3 = 1/20,
                k1.4 = 1/50,
                k1.5 = 1/100,
                k1.6 = 1/200,
                k1.7 = 1/500,
                k1.8 = 1/1000,
                k1.9 = 1/2000)

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
                vrm = 1))

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
    k_value == 1/20 ~ "1/20",
    k_value == 1/50 ~ "1/50",
    k_value == 1/100 ~ "1/100",
    k_value == 1/200 ~ "1/200",
    k_value == 1/350 ~ "1/350",
    k_value == 1/500 ~ "1/500",
    k_value == 1/1000 ~ "1/1000",
    k_value == 1/2000 ~ "1/2000"
  ))

OPM_results_k %>% 
  ggplot(aes(x = c_stock, y = lyr_14c, label = k_label)) +
  geom_path(aes(color = fct_reorder(k_label, -k_value)), linewidth = 1.5) +
  geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-500,300),
                     expand = c(0,0)) +
  scale_fill_discrete("Depth layer") +
  scale_color_viridis_d("k value", direction = -1) +
  geom_text(aes(x = 1e-9, y = -160), label = downwardTransferRate) +
  geom_text(aes(x = 1e-10, y = -140), label = LitterInput) +
  annotate(geom = "text", x = 1e-11, y = -160, label = "downward transfer = ") +
  annotate(geom = "text", x = 1e-11, y = -140, label = "C inputs = ") +
  annotate(geom = "text", x = 1e-11, y = -180, label = "no root inputs") +
  geom_text(data = OPM_results_k %>% 
              arrange(desc(layer)) %>% 
              distinct(k_label, .keep_all = TRUE), nudge_y = -10) +
  coord_cartesian(xlim = c(1e-15,1e+03))
ggsave(file = paste0("./Figure/OnePoolModel_k_values_2019_", Sys.Date(), ".jpeg"), 
       width = 12, height = 6)
  
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
                vrm = 1))

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
  geom_path(aes(color = as.factor(dtr_value)), linewidth = 1.5) +
  geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-500,300),
                     expand = c(0,0)) +
  scale_fill_discrete("Depth layer") +
  scale_color_viridis_d("downward transport") +
  geom_text(aes(x = 1e-10, y = -140), label = k1) +
  geom_text(aes(x = 1e-10, y = -160), label = LitterInput) +
  annotate(geom = "text", x = 1e-11, y = -140, label = "k-value = ") +
  annotate(geom = "text", x = 1e-11, y = -160, label = "C inputs = ") +
  annotate(geom = "text", x = 1e-11, y = -180, label = "no root inputs ") +
  geom_text(data = OPM_results_d %>% 
              arrange(desc(layer)) %>% 
              distinct(dtr_value, .keep_all = TRUE), nudge_y = -10, nudge_x = -0.2) +
  coord_cartesian(xlim = c(1e-15,1e+03))
ggsave(file = paste0("./Figure/OnePoolModel_d_values_2019_", Sys.Date(), ".jpeg"), 
       width = 12, height = 6)

# plot(AtmF14 %>% 
#        filter(Year.AD >= 1900), type = "l", lty = 1, ylim = c(-200,1000))
# matlines(OPM@times, C14dt_d$d.6, lty = 1, col = 2:11)

### Different decomposition rates

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
k1 <- 1/200
LitterInput <- 1
rootInputs <- rep(0, times = length(depthLayers))

# Parameters for VTLM
downwardTransferRate <- 0.005
upwardTransferRate <- 0

vrm_list <- list(vrm.0 = 1,
                 vrm.1 = seq(from  = 1, by = -1/length(depthLayers), 
                             length.out = length(depthLayers)),
                 vrm.2 = exp(log(1):(log(1)-(length(depthLayers)-1))/3),
                 vrm.3 = exp(log(1):(log(1)-(length(depthLayers)-1))))

## Step 2. Create SoilR pool model & Create Vertical Transport model
OPM <- OnepModel14(
  t = years, 
  k = k1, 
  C0 = initialC,
  F0_Delta14C = initialF14, 
  In = LitterInput,
  inputFc = AtmF14
)

OPM_VTLM_vrm <- vrm_list %>% 
  map(\(x) VTLM(Model = OPM, 
                lyrs = depthLayers, 
                latIn = rootInputs, 
                d = downwardTransferRate, 
                u = upwardTransferRate,
                vrm = x))

## Step 3. Solve the vertical transport model
Cdt_vrm <- OPM_VTLM_vrm %>% 
  map(\(x) getC(x))

C14dt_vrm <- OPM_VTLM_vrm %>% 
  map(\(x) getF14(x))

Cdt_2019_vrm <- Cdt_vrm %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "vrm_id") %>% 
  pivot_longer(!vrm_id, names_to = "layer", values_to = "c_stock") %>% 
  mutate(layer = as.numeric(layer))

C14dt_2019_vrm <- C14dt_vrm %>% 
  map(\(x) tail(x, 1)) %>% 
  plyr::ldply(.id = "vrm_id") %>% 
  pivot_longer(!vrm_id, names_to = "layer", values_to = "lyr_14c") %>% 
  mutate(layer = as.numeric(layer))

vrm <- plyr::ldply(vrm_list, .id = "vrm_id", data.frame) %>% 
  rename(vrm_value = X..i..) %>% 
  distinct(vrm_id, .keep_all = TRUE)

OPM_output_vrm <- tibble(
  year = round(rep(tail(AtmF14$Year, 1), max(depthLayers)*length(vrm_list)), 0),
  layer = rep(depthLayers, length(vrm_list)),
  LitterInput = LitterInput,
  vrm_id = rep(vrm$vrm_id, each = length(depthLayers)),
  vrm_value = rep(vrm$vrm_value, each = length(depthLayers)),
  upwardTransferRate = upwardTransferRate,
  downwardTransferRate = downwardTransferRate,
  k = k1
)

OPM_results_vrm <- OPM_output_vrm %>% 
  left_join(C14dt_2019_vrm) %>% 
  left_join(Cdt_2019_vrm)

OPM_results_vrm %>% 
  ggplot(aes(x = c_stock, y = lyr_14c, label = round(vrm_value, 2))) +
  geom_path(aes(color = as.factor(vrm_id)), linewidth = 1.5) +
  geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-500,300),
                     expand = c(0,0)) +
  scale_fill_discrete("Depth layer") +
  scale_color_viridis_d("vertical rate modifier") +
  geom_text(aes(x = 1e-10, y = -140), label = k1) +
  geom_text(aes(x = 1e-10, y = -160), label = LitterInput) +
  geom_text(aes(x = 1e-9, y = -200), label = downwardTransferRate) +
  annotate(geom = "text", x = 1e-11, y = -140, label = "k-value = ") +
  annotate(geom = "text", x = 1e-11, y = -160, label = "C inputs = ") +
  annotate(geom = "text", x = 1e-11, y = -180, label = "no root inputs ") +
  annotate(geom = "text", x = 1e-11, y = -200, label = "downward transport = ") +
  geom_text(data = OPM_results_vrm %>% 
              arrange(desc(layer)) %>% 
              distinct(vrm_value, .keep_all = TRUE), nudge_y = -10, nudge_x = -0.2) +
  coord_cartesian(xlim = c(1e-15,1e+03))
ggsave(file = paste0("./Figure/OnePoolModel_vrm_values_2019_", Sys.Date(), ".jpeg"), 
       width = 12, height = 6)

# plot(AtmF14 %>% 
#        filter(Year.AD >= 1900), type = "l", lty = 1, ylim = c(-200,1000))
# matlines(OPM@times, C14dt_d$d.6, lty = 1, col = 2:11)

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
                u = upwardTransferRate))

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
  downwardTransferRate = downwardTransferRate,
  upwardTransferRate = upwardTransferRate,
  LI_id = rep(LI$LI_id, each = length(depthLayers)),
  LI_value = rep(LI$LI_value, each = length(depthLayers))
)

OPM_results_LI <- OPM_output_LI %>% 
  left_join(C14dt_2019_LI) %>% 
  left_join(Cdt_2019_LI)

OPM_results_LI %>% 
  ggplot(aes(x = c_stock, y = lyr_14c, label = LI_value)) +
  geom_path(aes(color = as.factor(LI_value)), linewidth = 1.5) +
  geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-500,300),
                     expand = c(0,0)) +
  scale_fill_discrete("Depth layer") +
  scale_color_viridis_d("Litter inputs") +
  geom_text(aes(x = 1e-9, y = -140), label = downwardTransferRate) +
  geom_text(aes(x = 1e-10, y = -160), label = k1) +
  annotate(geom = "text", x = 1e-11, y = -140, label = "downward transfer = ") +
  annotate(geom = "text", x = 1e-11, y = -160, label = "k-value = ") +
  annotate(geom = "text", x = 1e-11, y = -180, label = "no root inputs ") +
  geom_text(data = OPM_results_LI %>% 
              arrange(desc(layer)) %>% 
              distinct(LI_value, .keep_all = TRUE), nudge_y = -15, nudge_x = -0.1) +
  coord_cartesian(xlim = c(1e-15,1e+03))
ggsave(file = paste0("./Figure/OnePoolModel_LI_values_2019_", Sys.Date(), ".jpeg"), 
       width = 12, height = 6)

### Different lateral C inputs

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for VTLM
k1 <- 1/500
LitterInput <- 1

#Root inputs all have the same sum ~ 1.58
RI_list <- list(RI.0 = rep(0, times = length(depthLayers)),
                RI.const = rep(0.158, times = length(depthLayers)),
                RI.lin = seq(from  = 0.29, by = -0.29/length(depthLayers), 
                                 length.out = length(depthLayers)),
                RI.exp.1 = exp((log(0.1):(log(0.1)-(length(depthLayers)-1)))/3),
                RI.exp.2 = exp(log(1):(log(1)-(length(depthLayers)-1))))

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
                u = upwardTransferRate))

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

OPM_results_RI %>% 
  ggplot(aes(x = c_stock, y = lyr_14c, label = RI_id)) +
  geom_path(aes(color = as.factor(RI_id)), linewidth = 1.5) +
  geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-500,300),
                     expand = c(0,0)) +
  scale_fill_discrete("Depth layer") +
  scale_color_viridis_d("Root inputs") +
  geom_text(aes(x = 1e-10, y = -180), label = k1) +
  geom_text(aes(x = 1e-9, y = -160), label = downwardTransferRate) +
  geom_text(aes(x = 1e-10, y = -140), label = LitterInput) +
  annotate(geom = "text", x = 1e-11, y = -180, label = "k-value = ") +
  annotate(geom = "text", x = 1e-11, y = -160, label = "downward transfer = ") +
  annotate(geom = "text", x = 1e-11, y = -140, label = "C inputs = ") +
  geom_text(data = OPM_results_RI %>% 
              arrange(desc(layer)) %>% 
              distinct(RI_id, .keep_all = TRUE), nudge_y = -10, nudge_x = -0.1) +  
  # coord_cartesian(xlim = c(1e+00,1e+03), ylim = c(-200,100)) +
  coord_cartesian(xlim = c(1e-15,1e+03))
ggsave(file = paste0("./Figure/OnePoolModel_RI_values_2019_", Sys.Date(), ".jpeg"), 
       width = 12, height = 6)

### Pseudo-realistic scenarios

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model

RI_list <- list(RI.arid = exp(log(0.6):(log(0.6)-(length(depthLayers)-1))/0.5),
                RI.temp = exp(log(1.1):(log(1.1)-(length(depthLayers)-1))),
                RI.trop = exp(log(1.1):(log(1.1)-(length(depthLayers)-1)))/2)

plot(exp(log(0.6):(log(0.6)-(length(depthLayers)-1))/0.5), type = "l")
plot(exp(log(1.1):(log(1.1)-(length(depthLayers)-1))), type = "l")
plot(exp(log(1.1):(log(1.1)-(length(depthLayers)-1)))/2)

LI_list <- list(LI.arid = 0.35,
                LI.temp = 1,
                LI.trop = 1.5)

d_list <- list(d.arid = 0.008,
               d.temp = 0.008,
               d.trop = 0.01)

k_list <- list(k.arid = 1/450,
               k.temp = 1/350,
               k.trop = 1/200)

vrm_list <- list(vrm.arid = seq(from  = 1, by = -0.09, 
                                length.out = length(depthLayers)),
                 vrm.temp = exp(log(1):(log(1)-(length(depthLayers)-1))/5),
                 vrm.trop = exp(log(1):(log(1)-(length(depthLayers)-1))/3))

## Arid
OPM_arid <- OnepModel14(
  t = years,
  k = k_list$k.arid,
  C0 = initialC,
  F0_Delta14C = initialF14,
  In = LI_list$LI.arid,
  inputFc = AtmF14
)

OPM_arid <- VTLM(
  Model = OPM_arid,
  lyrs = depthLayers,
  latIn = RI_list$RI.arid,
  d = d_list$d.arid,
  u = upwardTransferRate,
  vrm = vrm_list$vrm.arid
)

Cdt_arid <- getC(OPM_arid)
C14dt_arid <- getF14(OPM_arid)

Cdt_2019_arid <- tail(Cdt_arid, 1) %>% 
  plyr::ldply() %>% 
  rename("c_stock" = V1) %>% 
  mutate(layer = depthLayers,
         model = "arid")

C14dt_2019_arid <- tail(C14dt_arid, 1) %>% 
  plyr::ldply() %>% 
  rename("lyr_14c" = V1) %>% 
  mutate(layer = depthLayers,
         model = "arid")

OPM_results_arid <- Cdt_2019_arid %>% 
  left_join(C14dt_2019_arid)

## Temperate
OPM_temp <- OnepModel14(
  t = years,
  k = k_list$k.temp,
  C0 = initialC,
  F0_Delta14C = initialF14,
  In = LI_list$LI.temp,
  inputFc = AtmF14
)

OPM_temp <- VTLM(
  Model = OPM_temp,
  lyrs = depthLayers,
  latIn = RI_list$RI.temp,
  d = d_list$d.temp,
  u = upwardTransferRate,
  vrm = vrm_list$vrm.temp
)

Cdt_temp <- getC(OPM_temp)
C14dt_temp <- getF14(OPM_temp)

Cdt_2019_temp <- tail(Cdt_temp, 1) %>% 
  plyr::ldply() %>% 
  rename("c_stock" = V1) %>% 
  mutate(layer = depthLayers,
         model = "temperate")

C14dt_2019_temp <- tail(C14dt_temp, 1) %>% 
  plyr::ldply() %>% 
  rename("lyr_14c" = V1) %>% 
  mutate(layer = depthLayers,
         model = "temperate")

OPM_results_temp <- Cdt_2019_temp %>% 
  left_join(C14dt_2019_temp)

## Tropical
OPM_trop <- OnepModel14(
  t = years,
  k = k_list$k.trop,
  C0 = initialC,
  F0_Delta14C = initialF14,
  In = LI_list$LI.trop,
  inputFc = AtmF14
)

OPM_trop <- VTLM(
  Model = OPM_trop,
  lyrs = depthLayers,
  latIn = RI_list$RI.trop,
  d = d_list$d.trop,
  u = upwardTransferRate,
  vrm = vrm_list$vrm.trop
)

Cdt_trop <- getC(OPM_trop)
C14dt_trop <- getF14(OPM_trop)

Cdt_2019_trop <- tail(Cdt_trop, 1) %>% 
  plyr::ldply() %>% 
  rename("c_stock" = V1) %>% 
  mutate(layer = depthLayers,
         model = "tropical")

C14dt_2019_trop <- tail(C14dt_trop, 1) %>% 
  plyr::ldply() %>% 
  rename("lyr_14c" = V1) %>% 
  mutate(layer = depthLayers,
         model = "tropical")

OPM_results_trop <- Cdt_2019_trop %>% 
  left_join(C14dt_2019_trop)

OPM_results_all <- rbind(OPM_results_arid, OPM_results_temp, OPM_results_trop)

OPM_results_all %>% 
  ggplot(aes(x = c_stock, y = lyr_14c)) +
  geom_path(aes(color = as.factor(model)), linewidth = 1.5) +
  geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0),
                     labels = scales::comma) +
  scale_y_continuous(expression(paste(Delta^14,"C")),limits = c(-200,150),
                     expand = c(0,0)) +
  scale_fill_discrete("Depth layer") +
  scale_color_viridis_d("Model", option = "plasma", direction = -1) +
  coord_cartesian(xlim = c(1,300))
ggsave(file = paste0("./Figure/OnePoolModel_PseudoRealistic_values_2019_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)


### Two-pool-model ####

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
years <- seq(-48050,2019, by = 1)
# years <- seq(1941.5, 2019, by = 0.5)
LitterInput <- 1  
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
# rootInputs <- rep(0, times = length(depthLayers))
rootInputs <- exp((log(0.5):(log(0.5)-(length(depthLayers)-1)))/3)
downwardTransferRate <- 0.001
upwardTransferRate <- 0

plot(exp((log(0.5):(log(0.5)-(length(depthLayers)-1)))/3), type = "l")

# vrm.1 <- seq(from  = 1, by = -0.05, 
#              length.out = length(depthLayers))

vrm.1 <- c(1.00, 1.00, 0.8, 0.7, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)

# vrm.2 <- exp(log(1):(log(1)-(length(depthLayers)-1))/3)

#Step 2. Create SoilR pool model
OPM_VTLM_a <- 
  map(a21_list, \(x) TwopSeriesModel14(
    t = years, 
    ks = c(k1, k2), 
    C0 = initialC,
    a21 = x,
    F0_Delta14C = initialF14, 
    In = LitterInput,
    inputFc = AtmF14
  )) %>% 
  map(\(x) VTLM(Model = x, 
                lyrs = depthLayers, 
                latIn = rootInputs, 
                d = downwardTransferRate, 
                u = upwardTransferRate, 
                vrm = vrm.1,
                ivalF14 = "age"))

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
  LitterInput = LitterInput,
  downwardTransferRate = downwardTransferRate,
  upwardTransferRate = upwardTransferRate,
  a_id = rep(a21$a_id, each = length(depthLayers)),
  a_value = rep(a21$a_value/k1, each = length(depthLayers))
)

OPM_results_a <- OPM_output_a %>% 
  left_join(L14t_2019_a) 

OPM_results_a %>% 
  ggplot(aes(x = C, y = F14C, label = a_value*100)) +
  geom_path(aes(color = as.factor(a_value*100)), linewidth = 1.5) +
  geom_point(aes(fill = as.factor(layer)), shape = 21, size = 2.5) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_continuous("C stocks; log-scale", trans = "log10", expand = c(0,0)) +
  scale_y_continuous(expression(paste(Delta^14,"C")), limits = c(-500,125),
                     expand = c(0,0), breaks = c(0,-250,-500)) +
  scale_fill_discrete("Depth layer") +
  scale_color_viridis_d("transfer to\nslow pool [%]") +
  geom_text(aes(x = 1e+00, y = 80), label = downwardTransferRate) +
  geom_text(aes(x = 1e+00, y = 60), label = LitterInput) +
  geom_text(aes(x = 1e+00, y = 40), label = k1) +
  geom_text(aes(x = 1e+00, y = 20), label = k2) +
  geom_text(aes(x = 1e+00, y = 0), label = "exp(-x)") +
  annotate(geom = "text", x = 1e-01, y = 80, label = "downward transfer = ") +
  annotate(geom = "text", x = 1e-01, y = 60, label = "C inputs = ") +
  annotate(geom = "text", x = 1e-01, y = 40, label = "k1 = ") +
  annotate(geom = "text", x = 1e-01, y = 20, label = "k2 = ") +
  annotate(geom = "text", x = 1e-01, y = 0, label = "root inputs = ") +
  geom_text(data = OPM_results_a %>% 
              arrange(desc(layer)) %>% 
              distinct(a_value, .keep_all = TRUE), nudge_y = -10) +
  coord_cartesian(xlim = c(1e+01,1e+03))
ggsave(file = paste0("./Figure/TwoPoolModel_a_values_2019_roots_",
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
                     limits = c(1e+00,1e+04)) +
  scale_y_continuous(expression(paste(Delta^14,"C")), limits = c(-450,150),
                     expand = c(0,0))
ggsave(file = paste0("./Figure/TwoPoolModel_a_values_2019_bothPools_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

## Combine various model runs in one figure



