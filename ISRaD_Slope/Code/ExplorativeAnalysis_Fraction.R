# Compile ISRaD database #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 05/12/2022 #

## Fraction data ##

## NOT WORKING CURRENTLY ##

library(tidyverse)
library(ggpubr)
library(mpspline2)
library(ISRaD)

ISRaD_key <- readRDS("./Data/ISRaD_extra_2022-10-21")

names(ISRaD_key$layer)


frc_data_all <- ISRaD.flatten(ISRaD_key, "fraction")

frc_data <- frc_data_all %>% 
  drop_na(frc_14c) %>% 
  drop_na(frc_c_org) %>% 
  filter(frc_scheme == "density") %>% 
  tibble() %>% 
  unite("id", c(entry_name, site_name, pro_name), remove = FALSE) %>% 
  filter(lyr_obs_date_y.x > 1959) %>% 
  group_by(id) %>% 
  filter(n() > 2) %>% 
  ungroup() %>% 
  mutate(ClimateZoneAnd = case_when(
    entry_name == "Gentsch_2018" ~ "tundra/polar",
    pro_usda_soil_order == "Gelisols" ~ "tundra/polar",
    pro_usda_soil_order == "Andisols" ~ "volcanic soils",
    str_detect(pro_KG_present_long, "Tropical") ~ "tropical",
    str_detect(pro_KG_present_long, "Temperate") ~ "warm temperate",
    str_detect(pro_KG_present_long, "Cold") ~ "cold temperate",
    str_detect(pro_KG_present_long, "Polar") ~ "tundra/polar",
    str_detect(pro_KG_present_long, "Arid") ~ "arid",
  )) %>% 
  #remove entries that have duplicates/composite and not enough depth
  dplyr::filter(entry_name != "Czimczik_2010") %>% 
  #Remove duplicates
  distinct(id, frc_c_org, frc_14c, .keep_all = TRUE) %>%
  #Remove crytoturbated pockets in Gentsch_2018 and organic horizon
  dplyr::filter(is.na(lyr_hzn)|lyr_hzn != "Ajj_Ojj",
                is.na(lyr_hzn)|lyr_hzn != "O") %>% 
  #Remove incubation data from Vaughn_2018
  dplyr::filter(id != "Vaughn_2018_Barrow_C1C",
                id != "Vaughn_2018_Barrow_TC",
                id != "Vaughn_2018_Barrow_Z210C",
                id != "Vaughn_2018_Barrow_B3C",
                id != "Vaughn_2018_Barrow_A3C",
                id != "Vaughn_2018_Barrow_A4C",
                id != "Vaughn_2018_Barrow_H3C",
                id != "Vaughn_2018_Barrow_Z415C",
                id != "Vaughn_2018_Barrow_Z53C") 


#Generate frc_c_tot_filled column which includes all possible data for later calculation of C to N ratios
## USE FUNCTION FROM ISRaD INSTEAD?! THERE IS A SCRIPT FROM SHANE IN DEVSCRIPTS
frc <- dplyr::mutate(frc_data, frc_c_tot_filled = frc_c_tot, 
                     frc_c_tot_calc = frc_c_org) %>% 
  dplyr::mutate(., frc_c_tot_filled = replace(frc_c_tot_filled, is.na(frc_c_tot_filled), 
                                              frc_c_tot_calc[is.na(frc_c_tot_filled)])) %>% 
  dplyr::select(-frc_c_tot_calc)

#Ensure all possible c_to_n values included: (1) duplicate reported values ("filled"), (2) calculate values based on c and n data ("calc"), (3) merge reported and calculated
frc <- dplyr::mutate(frc, frc_c_to_n_filled = frc_c_to_n, frc_c_to_n_calc = frc_c_tot_filled/frc_n_tot) %>% 
  dplyr::mutate(., frc_c_to_n_filled = replace(frc_c_to_n_filled, is.na(frc_c_to_n_filled), 
                                               frc_c_to_n_calc[is.na(frc_c_to_n_filled)])) %>% 
  dplyr::select(-frc_c_to_n_calc)

#Fill in frc_c_perc, this reduces NA's 
frc_den2 <- frc %>% 
  dplyr::group_by(entry_name, site_name, pro_name, lyr_name, lyr_bot) %>% 
  dplyr::summarize(frc_c_perc_calc1 = sum(frc_c_tot_filled/100*frc_mass_perc/100)) %>% 
  ungroup() %>% 
  right_join(.,frc, by = c("entry_name", "site_name", "pro_name", "lyr_name", "lyr_bot")) %>% 
  dplyr::mutate(frc_c_perc_calc2 = ((frc_c_tot_filled/100*frc_mass_perc/100)/frc_c_perc_calc1)*100) %>% 
  dplyr::mutate(frc_c_perc_filled = ifelse(is.na(frc_c_perc), frc_c_perc_calc2, frc_c_perc))

names(frc_den2)

##Which fractions are multiples
frc_den_multiple <- frc_den2 %>%
  dplyr::group_by(entry_name, site_name, pro_name, lyr_name, lyr_bot, frc_property) %>%
  dplyr::summarize(count=length(frc_property)) %>%
  filter(count>1) 

#Weighted average for 13C, 15N, C:N, 14C, frc_c_tot_filled by proportion of C because its a good total org proxy, summed perc
group1 <- frc_den2 %>% 
  dplyr::group_by(entry_name, site_name, pro_name, lyr_name, lyr_bot, frc_property) %>%
  dplyr::summarize(frc_c_perc_sum = sum(frc_c_perc)) %>% 
  ungroup()

frc_data %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

frc_den_mult_comb <- left_join(frc_den2, group1, 
                               by = c("entry_name", "site_name", "pro_name", 
                                      "lyr_name", "lyr_bot", "frc_property")) %>%
  dplyr::mutate(c_weight = frc_c_perc_filled / frc_c_perc_sum) %>%
  dplyr::group_by(entry_name, site_name, pro_name, lyr_name, lyr_bot, frc_property) %>%
  dplyr::summarize(frc_c_perc_filled = sum(frc_c_perc_filled), 
                   count = length(frc_property), across(c(frc_13c, frc_15n, 
                                                          frc_14c, frc_c_to_n_filled, 
                                                          frc_c_tot_filled, frc_n_tot), 
                                                        ~weighted.mean(., w = c_weight))) %>% 
  ungroup(.)

# # #Any multiples left?
frc_den_multiple_test <- frc_den_mult_comb %>%
  dplyr::group_by(entry_name, site_name, pro_name, lyr_name, lyr_bot, frc_property) %>%
  dplyr::summarize(count=length(frc_property)) %>%
  filter(count>1) #none

frc_den_new <- frc_den2 %>%
  distinct(., across(c(entry_name, site_name, pro_name, lyr_name, lyr_bot, frc_property)), 
           .keep_all = TRUE) %>% #gets rid of duplicates, keeps the value from the first instance of a duplicate
  right_join(.,frc_den_mult_comb, 
             by = c("entry_name", "site_name", "pro_name", "lyr_name", "lyr_bot",
                    "frc_property")) %>% #join to combined variables (now denoted with .y)
  dplyr::mutate(frc_14c = ifelse(count > 1, frc_14c.y, frc_14c.x),
                frc_c_to_n_filled = ifelse(count > 1, frc_c_to_n_filled.y, frc_c_to_n_filled.x),
                frc_c_tot_filled = ifelse(count > 1, frc_c_tot_filled.y, frc_c_tot_filled.x),
                frc_n_tot = ifelse(count > 1, frc_n_tot.y, frc_n_tot.x),
                frc_c_perc_filled = ifelse(count > 1, frc_c_perc_filled.y, frc_c_perc_filled.x))%>%
  dplyr::select(-(c(contains(".x"), contains(".y"))))


### Apply mspline function

## mspline 14C
frc_data_mpspline_14c <- frc_den_new %>% 
  filter(frc_property == "heavy") %>% 
  group_by(entry_name, site_name, pro_name) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  dplyr::select(id, lyr_top, lyr_bot, frc_14c, frc_property) %>% 
  mpspline_tidy(vlow = -1000, lam = 0.5)

## mspline CORG
frc_data_mpspline_c <- frc_den_new %>% 
  filter(frc_property == "heavy") %>% 
  group_by(entry_name, site_name, pro_name) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  dplyr::select(id, lyr_top, lyr_bot, frc_c_org, frc_property) %>% 
  mpspline_tidy(vlow = -1000, lam = 0.5)

## 14C and SOC
frc_mspline_14c_c <- frc_data_mpspline_14c$est_1cm %>% 
  rename(frc_14c_msp = SPLINED_VALUE) %>% 
  full_join(frc_data_mpspline_c$est_1cm %>% 
              rename(CORG_msp = SPLINED_VALUE)) %>% 
  filter(LD < 101) %>% 
  tibble()

frc_mspline_14c_c_all <- frc_mspline_14c_c %>%
  dplyr::left_join(frc_den_new %>% 
                     distinct(id, .keep_all = TRUE), 
                   by = "id") %>% 
  group_by(id) %>% 
  arrange(UD) %>% 
  ungroup()

## Plotting
frc_mspline_14c_c_all %>% 
  group_by(ClimateZoneAnd, UD) %>% 
  summarize(median_14c = median(frc_14c_msp),
            mad_14c = mad(frc_14c_msp),
            median_CORG = median(CORG_msp),
            mad_CORG = mad(CORG_msp),
            n = n(),
            n_site = n_distinct(site_name),
            lower_CORG = median_CORG - mad_CORG,
            upper_CORG = median_CORG + mad_CORG) %>% 
  mutate(lower_CORG = if_else(lower_CORG <= 0.01, 0.01, lower_CORG)) %>% 
  ungroup() %>% 
  mutate(n_rel = n/max(n)) %>% 
  arrange(desc(median_14c)) %>% 
  filter(n_site > 4) %>%
  filter(n_rel > 1/9) %>%
  ggplot(aes(y = median_14c, x = median_CORG, color = ClimateZoneAnd)) +
  # geom_errorbarh(aes(xmin = lower_CORG, xmax = upper_CORG), color = "grey") +
  # geom_errorbar(aes(ymin = median_14c - mad_14c, ymax = median_14c + mad_14c), color = "grey") +
  geom_path(size = 2) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Delta 14C", expand = c(0,0), limits = c(-1000,200)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10", expand = c(0,0), limits = c(0.01, 8))




