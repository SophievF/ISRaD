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

ISRaD_key <- readRDS("./Data/ISRaD_extra_2022-09-22")

lyr_data_all <- ISRaD.flatten(ISRaD_key, 'layer')

lyr_data_all %>% 
  count(entry_name)

names(lyr_data_all)

#Prepare and filter data
lyr_data <- lyr_data_all %>% 
  drop_na(lyr_14c) %>% 
  drop_na(lyr_c_org_filled) %>% 
  rename(CORG = lyr_c_org_filled) %>% 
  # mutate(CORG = case_when(
  #   is.na(lyr_c_org) ~ lyr_c_tot,
  #   TRUE ~ lyr_c_org
  # )) %>%
  # drop_na(CORG) %>% 
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
  #remove Diemont_1987: peat study
  filter(entry_name != "Diemont_1987") %>% 
  #remove OBrien_1986: profile from under a house
  filter(entry_name != "OBrien_1986") %>%
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
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))


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
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

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
  # Gap fill with global data product
  mutate(pro_usda_soil_order = ifelse(is.na(pro_usda_soil_order),
                                      pro_0.5_USDA_soilorder, 
                                      pro_usda_soil_order)) %>% 
  #Fill missing soil type based on expert knowledge
  mutate(pro_usda_soil_order = ifelse((entry_name == "Lassey_1996" &
                                  is.na(pro_usda_soil_order)), "Inceptisols",
                               pro_usda_soil_order)) %>% 
  #Fill missing soil type based on expert knowledge
  mutate(pro_usda_soil_order = ifelse((entry_name == "Dintwe_2022" &
                                         is.na(pro_usda_soil_order)), "Entisols",
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

## Add WRB soil classification
#https://data.isric.org/geonetwork/srv/ger/catalog.search#/metadata/5c301e97-9662-4f77-aa2d-48facd3c9e14

WRB_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/ISRIC_WRB_map.tif"
WRB_raster <- raster::raster(WRB_dir)
wrb_number <- raster::extract(WRB_raster, cbind(lyr_data_fill$pro_long,
                                                        lyr_data_fill$pro_lat))

wrb_legend <- read_csv("D:/Sophie/PhD/AfSIS_GlobalData/ISRIC_WRB_map_legend.csv") %>% 
  dplyr::select(Number, WRB_group) %>% 
  dplyr::rename(pro_250m_wrb_soil_order = WRB_group)

wrb_data <- wrb_number %>% 
  tibble() %>% 
  rename(Number = '.') %>% 
  left_join(wrb_legend) %>% 
  dplyr::select(-Number) 

lyr_data_fill_wrb <- cbind(lyr_data_fill, wrb_data) %>%
  tibble() %>% 
  #Manually fix missing values
  mutate(pro_250m_wrb_soil_order = ifelse(pro_name == "Turlock Lake_123",
                                          "Luvisols", pro_250m_wrb_soil_order)) %>% 
  mutate(pro_250m_wrb_soil_order = ifelse(pro_name == "Timmendorf Forest_Profile",
                                          "Luvisols", pro_250m_wrb_soil_order)) %>% 
  mutate(pro_250m_wrb_soil_order = ifelse(pro_name == "SC1",
                                          "Kastanozems", pro_250m_wrb_soil_order)) %>% 
  mutate(pro_250m_wrb_soil_order = ifelse(pro_name == "MES:46.9542,-99.2792_331",
                                          "Chernozems", pro_250m_wrb_soil_order)) %>% 
  mutate(pro_250m_wrb_soil_order = ifelse(pro_name == "Amapá T1-Savanna",
                                          "Ferralsols", pro_250m_wrb_soil_order)) %>% 
  mutate(pro_250m_wrb_soil_order = ifelse(pro_name == "Amapá T3-Forest",
                                          "Ferralsols", pro_250m_wrb_soil_order)) %>% 
  #gap fill existing data first and then use global products
  mutate(
    pro_wrb_soil_order = 
      ifelse(grepl("Nitisols|Nitisol|nitisol", pro_soil_taxon), 
             "Nitisols", 
             ifelse(grepl("Luvisols|Luvisol|luvisol", pro_soil_taxon), 
                    "Luvisols", 
                    ifelse(grepl("Acrisols|Acrisol|acrisol", pro_soil_taxon), 
                           "Acrisols",
                           ifelse(grepl("Cambisols|Cambisol|cambisol", pro_soil_taxon), 
                                  "Cambisols",
                                  ifelse(grepl("Vertisols|Vertisol|vertisol", pro_soil_taxon), 
                                         "Vertisols",
                                         ifelse(grepl("Leptosols|Leptosol|leptosol", pro_soil_taxon), 
                                                "Leptosols",
                                                ifelse(grepl("Ferralsols|Ferrasols|Ferralsol|ferralsol", pro_soil_taxon), 
                                                       "Ferralsols",
                                                       ifelse(grepl("Calcisols|Calcisol|calcisol", pro_soil_taxon), 
                                                              "Calcisols",
                                                              ifelse(grepl("Kastanozems|Kastanozem|kastanozem", pro_soil_taxon), 
                                                                     "Kastanozems",
                                                                     ifelse(grepl("Chernozems|Chernozem|chernozem", pro_soil_taxon), 
                                                                            "Chernozems",
                                                                            ifelse(grepl("Andosols|Andosol|andosol|ands", pro_soil_taxon), 
                                                                                   "Andosols",
                                                                                   ifelse(grepl("Lixisols|Lixisol|lixisol", pro_soil_taxon), 
                                                                                          "Lixisols",
                                                                                          ifelse(grepl("Podzols|Podzol|podzol", pro_soil_taxon), 
                                                                                                 "Podzols",
                                                                                                 ifelse(grepl("Phaeozems|Phaeozem|phaeozem", pro_soil_taxon), 
                                                                                                        "Phaeozems",
                                                                                                        ifelse(grepl("Fluvisols|Fluvisol|fluvisol", pro_soil_taxon), 
                                                                                                               "Fluvisols",
                                                                                                               ifelse(grepl("Stagnosols|Stagnosol|stagnosol", pro_soil_taxon), 
                                                                                                                      "Stagnosols",
                                                                                                                      ifelse(grepl("Alisols|Alisol|alisol", pro_soil_taxon), 
                                                                                                                             "Alisols",
                                                                                                                             pro_250m_wrb_soil_order)))))))))))))))))) %>%
  #Manually assignment based on expert knowledge
  mutate(pro_wrb_soil_order = case_when(
    entry_name == "Chiti_2010" ~ "Ferralsols",
    entry_name == "Chiti_2018b" ~ "Ferralsols",
    entry_name == "Crow_2015" & pro_usda_soil_order == "Andisols" ~ "Andosols",
    entry_name == "Czimczik_2005" ~ "Podzols",
    entry_name == "Czimczik_2007" ~ "Cryosols",
    entry_name == "Gentsch_2018" ~ "Cryosols",
    entry_name == "Grant_2022" ~ "Andosols",
    entry_name == "Horwath_2008" ~ "Cryosols",
    entry_name == "Horwath_Thesis_2007" ~ "Cryosols",
    entry_name == "Katsuno_2010" ~ "Andosols",
    entry_name == "Khomo_2017" & pro_usda_soil_order == "Andisols" ~ "Andosols",
    entry_name == "Kleber_2005_Lausche_Lausche_30" ~ "Andosols",
    entry_name == "Kondo_2010" ~ "Andosols",
    entry_name == "Kramer_2012" ~ "Andosols",
    entry_name == "Laskar_2012" & pro_usda_soil_order == "Spodosols" ~ "Podzols",
    entry_name == "O'Donnell_2011" ~ "Cryosols",
    entry_name == "Rasmussen_2018" & pro_usda_soil_order == "Andisols" ~ "Andosols",
    entry_name == "Torn_1997" & pro_usda_soil_order == "Andisols" ~ "Andosols",
    entry_name == "Trumbore_Harden_1997" & pro_usda_soil_order == "Gelisols" ~ "Cryosols",
    entry_name == "Trumbore_Harden_1997" & pro_usda_soil_order == "Spodosols" ~ "Podzols",
    TRUE ~ pro_wrb_soil_order
  )) 


lyr_data_fill_wrb %>% 
  filter(is.na(pro_250m_wrb_soil_order)) %>% 
  count(entry_name, pro_name)

lyr_data_fill_wrb %>% 
  count(entry_name, pro_name, pro_250m_wrb_soil_order,
        pro_wrb_soil_order, pro_soil_taxon) %>% view()

lyr_data_fill_wrb %>% 
  filter(is.na(pro_wrb_soil_order)) %>% 
  count(entry_name, pro_name)

lyr_data_fill_wrb %>% 
  summarise(n_studies = n_distinct(entry_name),
            n_sites = n_distinct(site_name),
            n_profiles = n_distinct(id))

saveRDS(lyr_data_fill_wrb, paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_", Sys.Date()))

