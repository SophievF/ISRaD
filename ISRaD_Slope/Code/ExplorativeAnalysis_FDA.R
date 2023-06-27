# Explore 14C profiles in ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 22/06/2022 #

## Functional data analysis (FDA) ##
#https://www.psych.mcgill.ca/misc/fda/ex-goods-a1.html

library(tidyverse)
library(fda)

#Load filtered lyr data
lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-06-28"))

lyr_fda <- lyr_data %>% 
  filter(lyr_bot <= 200) %>% 
  group_by(id) %>%
  #Filter for studies that have more than 2 depth layers
  filter(n() > 2) %>%
  arrange(depth, .by_group = TRUE) %>% 
  ungroup()

##Depth~14C; SOC~14C; built fda object with both options; compare to mspline
n_id <- length(unique(lyr_fda$id))
max_depth <- max(lyr_fda$depth)
n_points <- length(lyr_fda$id)

lyr_list <- lyr_fda %>%
  dplyr::select(id, depth, lyr_14c) %>% 
  group_by(id) %>% 
  summarise(across(everything(), list))
view(lyr_list)

ggplot(lyr_fda, aes(x = depth, y = lyr_14c, group = id, col = id)) +
  geom_path() +
  theme_classic() +
  theme(legend.position = "none")

knots <- c(seq(0, max_depth, 10))
n_knots <- length(knots)
n_order <- 4
n_basis <- length(knots) + n_order - 2
basis <- create.bspline.basis(rangeval = c(0, max_depth), n_basis)
plot(basis)

rgvals <- matrix(lyr_fda$depth, lyr_fda$lyr_14c, ncol = n_id)
str(rgvals)

## Functional PCA ##
#https://cran.r-project.org/web/packages/fdapace/index.html

library(fdapace)

lyr_list <- lyr_fda %>%
  dplyr::select(id, depth, lyr_14c) %>% 
  group_by(id) %>% 
  summarise(across(everything(), list))

FPCAsparse <- FPCA(lyr_list$lyr_14c, lyr_list$depth, list(plot = TRUE))

par(mfrow=c(1,2))
CreatePathPlot(FPCAsparse, subset = c(3,5,135), main = 'K = 11', pch = 4); grid()
CreatePathPlot(FPCAsparse, subset = c(3,5,135), K = 3, main = 'K = 3', pch = 4) ; grid()

par(mfrow=c(1,1))
CreateOutliersPlot(FPCAsparse, optns = list(K = 3, variant = 'KDE'))

CreateFuncBoxPlot(FPCAsparse, xlab = 'CORG', ylab = 'Delta14C', 
                  optns = list(K =3, variant='pointwise'))



