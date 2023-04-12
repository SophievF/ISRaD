library(dplyr)
library(ggplot2)
library(tidyr)

# setwd("C:/Users/benbu/OneDrive - soilbenchmark.com/Documents/Ben/Sophie Von Fromm")

original_data_example <- read.csv("./Data/original_data_example.csv")

pdp <- read.csv("./Data/PDP_14C_example.csv")

#Group by ClimateZoneAnd and then Summarise 95% of data
quantiles <- original_data_example %>% group_by(ClimateZoneAnd) %>% 
  dplyr::summarise(P5_MAT = quantile(MAT, 0.025),
                   P95_MAT = quantile(MAT, 0.975),
                   P5_PET_MAP = quantile(PET_MAP, 0.025),
                   P95_PET_MAP = quantile(PET_MAP, 0.975))

#Split the data based on the predictor so that we can create new vectors
#that can be interpolated onto
pdp_split <- split(pdp, f = pdp$predictor)

#create the higher resolution vectors to interpolate onto

#one for MAT
pred_value_high_res_mat <- seq(min(pdp_split$MAT$pred_value),
                               max(pdp_split$MAT$pred_value),
                               length.out = 100)

#one for PET_MAP
pred_value_high_res_pet <- seq(min(pdp_split$PET_MAP$pred_value),
                               max(pdp_split$PET_MAP$pred_value),
                               length.out = 100) #the higher this is the great the resolution

#Split the data again by ClimateZoneAnd
pdp_split_mat <- split(pdp_split$MAT, f = pdp_split$MAT$ClimateZoneAnd)
pdp_split_pet <- split(pdp_split$PET_MAP, f = pdp_split$PET_MAP$ClimateZoneAnd)

#Create a function that will compute a natural spline of the data based on the
#higher resolution vectors
spline_data <- function(data, new_x) {
  
 splined_data <- approx(x = data$pred_value,
                        y = data$pred_median_14c,
                        method = "linear",
                        xout = new_x) |> data.frame()
 
 names(splined_data) <- c("pred_value", "pred_median_14c")
 
 splined_data$ClimateZoneAnd <- unique(data$ClimateZoneAnd)
 
 splined_data$predictor <- unique(data$predictor)
 
 return(splined_data)
 
 }

#Apply the spline function to the MAT data
pdp_split_mat_splined <- lapply(pdp_split_mat, spline_data,
                                new_x = pred_value_high_res_mat)

#Apply the spline function to the PET_MAP data
pdp_split_pet_splined <- lapply(pdp_split_pet, spline_data,
                                new_x = pred_value_high_res_pet)

#Bind them all together
pdp_splined <- do.call(rbind, c(pdp_split_mat_splined,
                                pdp_split_pet_splined))

#Split them again by ClimateZoneAnd so that we can run the outlier loop
pdp_splined_split <- split(pdp_splined, f = pdp_splined$ClimateZoneAnd)

#Just a quick check that I'm dealing with data in the right order because
#I'm about to right a clunky for loop that uses data in the pdp list and in the
#quantiles data table.
identical(quantiles$ClimateZoneAnd, names(pdp_splined_split))

#Run the loop that assigns outliers
for (i in 1:length(pdp_splined_split )) {
  
  pdp_splined_split[[i]] <- pdp_splined_split[[i]] %>% mutate(outlier = case_when(
    predictor == "MAT" & (pred_value < quantiles$P5_MAT[i] | pred_value > quantiles$P95_MAT[i]) ~ TRUE,
    predictor == "PET_MAP" & (pred_value < quantiles$P5_PET_MAP[i] | pred_value > quantiles$P95_PET_MAP[i]) ~ TRUE,
    .default = FALSE
  ))
  
}

#Bind that data together for plotting
pdp_splined <- do.call(rbind, pdp_splined_split)

#Plot data: shows entire range for each climate zone PDP_14C_example %>%
ggplot(pdp_splined) +
#  geom_point(data = pdp_splined, aes(x = pred_value, y = pred_median_14c, color = ClimateZoneAnd),
#             size = 2) +
  geom_line(data = pdp_splined, aes(x = pred_value, y = pred_median_14c, color = ClimateZoneAnd),
            linewidth = 1.5, alpha = 0.1) +
  geom_line(data = pdp_splined[-which(pdp_splined$outlier == TRUE), ], aes(x = pred_value, y = pred_median_14c, color = ClimateZoneAnd),
            linewidth = 1.5) +
  geom_rug(data = original_data_example %>%
             pivot_longer(!ClimateZoneAnd, values_to = "pred_value", names_to = "predictor"),
           aes(x = pred_value, color = ClimateZoneAnd),
          sides = "b") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top") +
  facet_wrap(~predictor, scales = "free_x") +
  scale_x_continuous("", expand = c(0,0)) +
  scale_y_continuous(expression(paste("Median predicted  ", Delta^14, "C [‰]")),
                     expand = c(0,0), limits = c(-510,100))
