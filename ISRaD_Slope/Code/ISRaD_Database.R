# Compile ISRaD database #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 01/06/2022 #

devtools::install_github('International-Soil-Radiocarbon-Database/ISRaD/Rpkg',
                         force = TRUE)

library(ISRaD)
library(tidyverse)
library(ggpubr)

ISRaD_dir <- "C:/Users/sfromm/Documents/GitHub/ISRaD_SvF/ISRaD_data_files"

#To compile ISRaD manually
ISRaD_comp <- compile(dataset_directory = ISRaD_dir, write_report = TRUE, 
                      write_out = TRUE, return = "list")

geo_dir <- "D:/Seafile/ISRaD_geospatial_data/ISRaD_extra_geodata"
ISRaD_extra <- ISRaD.extra(ISRaD_comp, geodata_directory = geo_dir)

ISRaD_key <- ISRaD.extra.geospatial.keys(ISRaD_extra, 
                                         geodata_keys = "D:/Seafile/ISRaD_geospatial_data/ISRaD_extra_keys")

names(ISRaD_key)

# To extract data from github
# ISRaD_extra <- ISRaD.getdata(directory = ISRaD_dir,
#                              dataset = "full", extra = TRUE,
#                              force_download = TRUE)

saveRDS(ISRaD_key, paste0(getwd(), "/Data/ISRaD_extra_", Sys.Date()))

# ISRaD_key <- readRDS("./Data/ISRaD_extra_2022-06-30")

lyr_data_all <- ISRaD.flatten(ISRaD_key, 'layer')

lyr_data_all %>% 
  count(entry_name)

names(lyr_data_all)

#Prepare and filter data
lyr_data <- lyr_data_all %>% 
  drop_na(lyr_14c) %>% 
  mutate(CORG = case_when(
    is.na(lyr_c_org) ~ lyr_c_tot,
    TRUE ~ lyr_c_org
  )) %>%
  drop_na(CORG) %>% 
  filter(lyr_top >= 0 &
           lyr_bot >= 0 &
           pro_land_cover != "wetland" &
           is.na(pro_peatland)) %>% 
  filter(lyr_top != "Inf",
         lyr_bot != "Inf") %>% 
  unite("id", c(entry_name, site_name, pro_name), remove = FALSE) %>% 
  filter(!grepl("peat|Peat", id)) %>% 
  #remove permafrost studies
  #filter(is.na(lyr_all_org_neg)) %>% 
  #remove Huang_1999: peat study not flagged
  filter(entry_name != "Huang_1999") %>% 
  #remove Huang_1996: peat study not flagged
  filter(entry_name != "Huang_1996") %>% 
  #remove Bol_1996: peat study not flagged
  filter(entry_name != "Bol_1996") %>% 
  #remove Baisden_2002: data is also in Baisden_2007
  filter(entry_name != "Baisden_2002") %>% 
  #remove Kogel-Knabner_2008: same data as in Eusterhues_2003
  filter(entry_name != "Kogel-Knabner_2008") %>% 
  #remove profile: Chalk River Laboratories (CRL):46,-77.4_563: wetland
  filter(pro_name != "Chalk River Laboratories (CRL):46,-77.4_563") %>% 
  #filter CORG >0
  filter(CORG > 0) %>% 
  #depth = 200
  #filter(lyr_bot <= 200) %>% 
  # group_by(id) %>%
  # #Filter for studies that have more than 2 depth layers
  # filter(n() > 2) %>%
  # ungroup() %>% 
  #calculate layer mid-depth
  mutate(depth = ((lyr_bot - lyr_top)/2) + lyr_top)

lyr_data %>% 
  count(id)

#Check for data that has same depth value for same id
lyr_data %>%
  dplyr::select(entry_name, id, lyr_name, depth, lyr_top, lyr_bot, lyr_14c, CORG) %>%
  group_by(id, depth) %>%
  filter(n() > 1) %>% view()

lyr_data_clean <- lyr_data %>% 
  #remove entries that have duplicates/composite and not enough depth
  dplyr::filter(entry_name != "Czimczik_2010") %>% 
  #Remove duplicates
  distinct(id, depth, CORG, lyr_14c, .keep_all = TRUE) %>%
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

lyr_data_clean %>%
  dplyr::select(entry_name, id, lyr_name, depth, lyr_top, lyr_bot, lyr_14c, CORG) %>%
  group_by(id, depth) %>%
  filter(n() > 1) 

lyr_data_clean %>% 
  count(entry_name)

lyr_data_clean %>% 
  ggplot(aes(x = depth, y = lyr_14c, group = entry_name)) + 
  geom_point(size = 3, shape = 21) +
  theme_bw(base_size = 16)

# Gap-fill missing global data with reported local data (or vice-versa)
lyr_data_clean %>% 
  count(pro_usda_soil_order)
lyr_data_clean %>% 
  count(pro_0.5_USDA_soilorder)
lyr_data_clean %>% 
  count(pro_KG_present_long)
lyr_data_clean %>% 
  count(pro_KG_present_short)

summary(lyr_data_clean$pro_BIO12_mmyr_WC2.1)

lyr_data_fill <- lyr_data_clean %>% 
  mutate(pro_BIO12_mmyr_WC2.1 = ifelse(is.na(pro_BIO12_mmyr_WC2.1), 
                                       pro_MAP, pro_BIO12_mmyr_WC2.1),
         pro_BIO1_C_WC2.1 = ifelse(is.na(pro_BIO1_C_WC2.1),
                                   pro_MAT, pro_BIO1_C_WC2.1)) %>% 
  #Use WRB to fill USDA: https://www.isric.org/sites/default/files/major_soils_of_the_world/annexes/index.pdf
  mutate(pro_usda_soil_order = ifelse(grepl("andosol|Andosol", pro_soil_taxon), "Andisols",
                                      pro_usda_soil_order),
         pro_usda_soil_order = ifelse(grepl("ferralsol|Ferralsol|Ferralols", pro_soil_taxon), "Oxisols",
                                      pro_usda_soil_order),
         pro_usda_soil_order = ifelse(grepl("podzol|Podzol", pro_soil_taxon), "Spodosols",
                                      pro_usda_soil_order),
         pro_usda_soil_order = ifelse(grepl("vertisol|Vertisol", pro_soil_taxon), "Vertisols",
                                      pro_usda_soil_order),
         pro_usda_soil_order = ifelse(grepl("Kastanozem|Chernozem|Phaeozem", pro_soil_taxon), "Mollisols",
                                      pro_usda_soil_order),
         pro_usda_soil_order = ifelse(grepl("luvisol|Luvisol", pro_soil_taxon), "Alfisols",
                                      pro_usda_soil_order)) %>% 
  mutate(pro_usda_soil_order = ifelse(is.na(pro_usda_soil_order),
                                    pro_0.5_USDA_soilorder, 
                                    pro_usda_soil_order)) %>% 
  #Fill missing soil type based on expert knowledge
  mutate(pro_usda_soil_order = ifelse((entry_name == "Lassey_1996" &
                                  is.na(pro_usda_soil_order)), "Inceptisols",
                               pro_usda_soil_order)) %>% 
  # Manually assign climate zone for Czimczik_Unpublished
  mutate(pro_KG_present_long = case_when(
    is.na(pro_KG_present_long) ~ "Polar, tundra",
    TRUE ~ pro_KG_present_long
  )) %>% 
  mutate(pro_KG_present_short = case_when(
    is.na(pro_KG_present_short) ~ "ET",
    TRUE ~ pro_KG_present_short
  ))

lyr_data_fill %>% 
  count(pro_usda_soil_order)

lyr_data_fill %>% 
  filter(is.na(pro_usda_soil_order)) %>% 
  count(entry_name)

lyr_data_fill %>% 
  filter(is.na(pro_KG_present_long)) %>% 
  count(entry_name)

summary(lyr_data_fill$pro_BIO12_mmyr_WC2.1)

saveRDS(lyr_data_fill, paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_", Sys.Date()))

