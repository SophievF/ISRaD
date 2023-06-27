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

### Set-up a depth-resolved model with SoilR (Carlos) ###
source("C:/Users/sfromm/Documents/GitHub/SoilR-exp/prototypes/VerticalTransfer/VTLM.R")
source("C:/Users/sfromm/Documents/GitHub/SoilR-exp/prototypes/VerticalTransfer/collect_by_layer.R")
source("C:/Users/sfromm/Documents/GitHub/SoilR-exp/prototypes/VerticalTransfer/collect_by_pool.R")
source("C:/Users/sfromm/Documents/GitHub/SoilR-exp/prototypes/VerticalTransfer/solveVTLM.R")

#### One-pool-model ####
p_empty <- ggplot() +
  theme_void()

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

k1_list <- list(k1.5 = 1/250,
                k1.6 = 1/500,
                k1.7 = 1/1000,
                k1.8 = 1/2000)

# Parameters for VTLM
depthLayers <- c(1:10)
rootInputs <- rep(0, times = length(depthLayers))
downwardTransferRate <- 0.005
upwardTransferRate <- 0

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
  left_join(C14dt_2019_k, multiple = "all") %>% 
  left_join(Cdt_2019_k, multiple = "all") %>% 
  mutate(k_label = case_when(
    k_value == 1/150 ~ "1/150",
    k_value == 1/250 ~ "1/250",
    k_value == 1/500 ~ "1/500",
    k_value == 1/1000 ~ "1/1000",
    k_value == 1/2000 ~ "1/2000"
  ))

OnePool_const_k <- OPM_results_k %>% 
  ggplot(aes(x = c_stock, y = lyr_14c)) +
  geom_path(aes(color = fct_reorder(k_label, -k_value), 
                linewidth = fct_reorder(k_label, -k_value))) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black", 
             color = "white") +
  theme_bw(base_size = 14) +
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.13,0.755),
        legend.background = element_blank(),
        plot.margin = margin(r = 20, t = 5, l = 0, b = 5)) +
  scale_x_continuous("", 
                     limits = c(0.5,1000), trans = "log10",
                     breaks = c(1,10,100,1000), 
                     expand = c(0,0),
                     labels = c(1,10,100,1000)) +
  scale_y_continuous("", limits = c(-350,125), expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_linewidth_manual("", values = c(2,3,2,2)) +
  scale_color_manual("k-value",
                     values = c("#a6611a", "#7b3294", 
                                "#018571", "#80cdc1")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none", linewidth = "none", 
         color = guide_legend(override.aes = list(linewidth = 2)))

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
  geom_path(aes(color = as.factor(vkm_id), linewidth = as.factor(vkm_id))) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.175,0.768),
        plot.margin = margin(r = 15, t = 5, l = 5, b = 5),
        legend.background = element_blank()) +
  scale_x_continuous("", 
                     limits = c(0.5,1000), trans = "log10",
                     breaks = c(1,10,100,1000), 
                     expand = c(0,0),
                     labels = c(1,10,100,1000)) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")),
                     limits = c(-350,125), expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_linewidth_manual("", values = c(3,2,2,2)) +
  scale_color_manual("vertical k modifier", 
                     values = c("#7b3294", "#0570b0", "#74a9cf", "#bdc9e1")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none", linewidth = "none", 
         color = guide_legend(override.aes = list(linewidth = 2)))
  

### Different litter inputs

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
k1 <- 1/500
rootInputs <- rep(0, times = length(depthLayers))

downwardTransferRate <- 0.005

LI_list <- list(LI.2 = 0.25,
                LI.3 = 0.5,
                LI.4 = 1,
                LI.5 = 2,
                LI.6 = 4)

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
  geom_path(aes(color = as.factor(LI_value),
                linewidth = as.factor(LI_value))) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.22,0.73),
        legend.background = element_blank(),
        plot.margin = margin(r = 20, t = 5, l = 0, b = 5)) +
  scale_x_continuous("", 
                     limits = c(0.5,1000), trans = "log10",
                     breaks = c(1,10,100,1000), 
                     expand = c(0,0),
                     labels = c(1,10,100,1000)) +
  scale_y_continuous("", limits = c(-350,125), expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_linewidth_manual("", values = c(2,2,3,2,2)) +
  scale_color_manual("aboveground C inputs", 
                     values = c("#dfc27d", "#a6611a", "#7b3294", 
                                "#018571", "#80cdc1")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none", linewidth = "none", 
         color = guide_legend(override.aes = list(linewidth = 2)))


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

# sum(rep(0.158, times = length(depthLayers)))
# sum(seq(from  = 0.288, by = -0.288/length(depthLayers), 
#         length.out = length(depthLayers)))
# sum(exp((log(0.39):(log(0.39)-(length(depthLayers)-1)))/2))
# sum(exp(log(1):(log(1)-(length(depthLayers)-1))))
# 
# 
# plot(x = depthLayers, y = RI_list$RI.const, type = "l")
# plot(x = depthLayers, y = RI_list$RI.lin, type = "l")
# plot(x = depthLayers, y = RI_list$RI.exp.1, type = "l")
# plot(x = depthLayers, y = RI_list$RI.exp.2, type = "l")


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
  geom_path(aes(color = as.factor(RI_id),
                linewidth = as.factor(RI_id))) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.22,0.73),
        legend.background = element_blank(),
        plot.margin = margin(r = 15, t = 5, l = 5, b = 5)) +
  scale_x_continuous("", 
                     limits = c(0.5,1000), trans = "log10",
                     breaks = c(1,10,100,1000), 
                     expand = c(0,0),
                     labels = c(1,10,100, 1000)) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")),
                     limits = c(-350,125), expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_linewidth_manual("", values = c(3,2,2,2,2)) +
  scale_color_manual("belowground C inputs", 
                     values = c("#7b3294", "#006837", "#31a354", 
                                "#78c679", "#c2e699")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none", linewidth = "none", 
         color = guide_legend(override.aes = list(linewidth = 2)))

### Different downward transfer

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
k1 <- 1/500
LitterInput <- 1
rootInputs <- rep(0, times = length(depthLayers))

# Parameters for VTLM
d_list <- list(d.3 = 0.02,
               d.5 = 0.01,
               d.6 = 0.005,
               d.7 = 0.0025)

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

OnePool_d <- OPM_results_d %>% 
  ggplot(aes(x = c_stock, y = lyr_14c, label = dtr_value)) +
  geom_path(aes(color = as.factor(dtr_value),
                linewidth = as.factor(dtr_value))) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.box = "horizontal",
        legend.position = c(0.18,0.758),
        legend.background = element_blank(),
        plot.margin = margin(r = 20, t = 5, l = 0, b = 5)) +
  scale_x_continuous("", 
                     limits = c(0.5,1000), trans = "log10",
                     breaks = c(1,10,100,1000), 
                     expand = c(0,0),
                     labels = c(1,10,100,1000)) +
  scale_y_continuous("", limits = c(-350,125), expand = c(0,0)) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_linewidth_manual("", values = c(2,2,3,2)) +
  scale_color_manual("vertical transport",
                     values = c("#80cdc1", "#018571", "#7b3294", "#a6611a")) +
  annotation_logticks(sides = "b", scaled = TRUE, short = unit(1.5,"mm"),
                      mid = unit(3,"mm"), long = unit(4,"mm")) +
  guides(shape = "none", linewidth = "none", 
         color = guide_legend(override.aes = list(linewidth = 2)))


### Two-pool-model ####

## Step 1. Prepare what you need for a SoilR pool model
# Parameters for one-pool model
years <- seq(-48050,2019, by = 1)
# years <- seq(1941.5, 2019, by = 0.5)
LitterInputs <- 1  
k1 <- 1/50
k2 <- 1/1250
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

downwardTransferRate <- 0.0025
upwardTransferRate <- 0

vrm.1 <- c(1.0, 0.9, 0.8, 0.7, 0.6, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)

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

TwoPool_a <- OPM_results_a %>% 
  ggplot(aes(x = C, y = F14C)) +
  geom_path(aes(color = as.factor(a_value*100)), linewidth = 2) +
  geom_point(aes(shape = as.factor(layer)), size = 2, fill = "black",
             color = "white") +
  theme_bw(base_size = 14) +
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.15,0.7),
        plot.margin = margin(r = 15, t = 0, l = 5, b = 10),
        legend.background = element_blank()) +
  scale_x_continuous("Soil organic carbon stocks; log-scale", 
                     limits = c(0.5,1000), trans = "log10",
                     breaks = c(1,10,100,1000), 
                     expand = c(0,0),
                     labels = c(1,10,100,1000)) +
  scale_y_continuous(expression(paste(Delta^14,"C [‰]")),
                     limits = c(-350,100), expand = c(0,0)) +
  annotation_logticks(sides = "b", scaled = TRUE,
                      short = unit(1.5,"mm"),
                      mid = unit(3,"mm"),
                      long = unit(4,"mm")) +
  scale_shape_manual("depth layer", values = c(25, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)) +
  scale_color_manual("transfer to\nslow pool [%]", 
                     values = c("#ffeda0", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")) +
  guides(shape = "none")

ggarrange(p_empty,  OnePool_const_k,  OnePool_vkm, OnePool_Litter, OnePool_Roots, 
          OnePool_d, TwoPool_a, nrow = 4, ncol = 2)
ggsave(file = paste0("./Figure/AllModels_2019_", Sys.Date(), ".png"),
       width = 12, height = 15)









