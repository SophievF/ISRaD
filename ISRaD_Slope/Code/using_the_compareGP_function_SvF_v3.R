###Update Gaussian Process Regression Script###
###2022-02-16###
###Based on Script from Xavier Raynaud###

#https://mycore.core-cloud.net/index.php/s/0oCu5YShl07Qd2I?path=%2F

source("GaussianProcess_functions_v2.R") # Update code. 

#This is to load your data table. Let's assume is creates an object called data
#devtools::install_github("International-Soil-Radiocarbon-Database/ISRaD/Rpkg", 
#                         ref = "master", force = TRUE)

library(ISRaD)
library(tidyverse)

ISRaD_extra <- ISRaD.getdata(directory = "IsRad", 
                             dataset = "full", extra = TRUE, 
                             force_download = F)

#Filter ISRaD data
lyr_data <- ISRaD.flatten(ISRaD_extra, 'layer') %>% 
  drop_na(lyr_14c_fill_extra) %>% 
  mutate(CORG = case_when(
    is.na(lyr_c_org) ~ lyr_c_tot,
    TRUE ~ lyr_c_org
  )) %>% 
  drop_na(CORG) %>%
  filter(lyr_top >= 0 &
           lyr_bot >= 0 &
           pro_land_cover != "wetland" &
           is.na(pro_peatland)) %>% 
  unite("id", c(entry_name, site_name, pro_name), remove = FALSE) %>% 
  filter(!grepl("peat|Peat", id)) %>% 
  mutate(depth = ((lyr_bot - lyr_top)/2)+lyr_top) %>%   # what does this line do ?
  rename(value = lyr_14c_fill_extra) %>%
  mutate(error = 0) %>%
  group_by(id) %>%
  mutate(norm_value = (value-min(value))/(max(value)-min(value))) %>%
  filter(n() > 2) %>%
  tibble() %>%
  group_by(id) %>%
  arrange(depth)  

#Visual exploration of the ISRaD data
alpha = 0.2
lyr_data %>% 
  ggplot() +
  #  geom_point(aes(x = value, y = depth, color = id), alpha = alpha) +
  geom_line(aes(x = norm_value, y = depth, group = id, color = id), orientation = "y", alpha = alpha) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous("Depth [cm]", trans = "reverse", limits = c(700,0),
                     expand = c(0,0)) +
  scale_x_continuous("Delta 14C")

lyr_data %>% 
  group_by(id) %>% 
  count()  

#Select all profiles that have alox measured
data <- lyr_data %>% 
  drop_na(lyr_al_ox) %>% 
  mutate(alox = case_when(
    lyr_al_ox == 0 ~ 0.001,
    TRUE ~ lyr_al_ox
  )) %>% 
  #remove data from Gentsch_2018 (causing problems: check why: depth issues - permafrost soil)
  filter(!grepl("Gentsch", id)) 

data %>% 
  group_by(id) %>% 
  count()%>%View(  )

#Data should be a dataframe or a tibble in long format, with at least one column
#to identify the profile ('id'), one column containing the depth ('depth'), one 
#containing the 14C data ('value') and one containing the measurement error ('error', can be empty). 

nsamples <- 1000

#ncores is the number of cpus you have on your computer. If you put a value 
#larger than the true number of cpus you have in your computer, it will be 
#reduced to the number of cpus minus 1. 

# The compareGP output has slightly changed: you can select whether you keep simulations or not (This make the d object waaaayyy smaller if you don't).
# The way to save the goodness of fit results has also changed. Instead of an nxn matrix, results are saved in the form of a long table with 3 columns, data, model and goodness of fit index. 
d <- compare_GP(data, 
                idcol = "id",
                xcol = "depth", 
                ycol = "norm_value", 
                errcol = "error", 
                nsamples = nsamples, 
                verbose = 0, 
                engine = "laGP_gaupro", 
                ncores = 20, keep_simulations = T)

#Save the file with a filename containing the date. For 600 profiles, the file is ~1Gb. 
#save(d, file = paste0("results_gsol_ISRaD_Alox", nsamples,"-", 
#                      format(Sys.time(), "%Y%m%d%H%M"), ".RData"))

#str(d)

# d$models$`Spielvogel_2008_National Park Bayrischer Wald_Leptic Cambisol`$plot()


alpha = 0.1
g = list()
profnames = c("Torn_1997_Kokee (4.1my)_Kokee (4.1my)_profile_1",
              "Spielvogel_2008_National Park Bayrischer Wald_Leptic Cambisol", 
              "Trumbore_1996_Trumbore Shaver_Trumbore Shaver 2 1959", 
              "Khomo_2017_GR-740-C_GR-740-C")

for (p in 1: length(profnames)) {
  g[[p]] = d$simuls[[profnames[p]]] %>% 
    data.frame() %>% 
    pivot_longer(cols = -1, names_to = "runs", values_to = "C14values") %>% 
    group_by(depth) %>% 
    mutate(mean = mean(C14values)) %>% 
    ungroup() %>% 
    ggplot() +
    geom_line(aes(y = depth, x = C14values, group = runs),alpha = alpha,
              color = "darkgrey", orientation = "y") +
    geom_line(aes(y = depth, x = mean),
              orientation = "y", color = "red", size = 0.7) +
    geom_point(data = lyr_data %>% 
                 filter(id == profnames[p]),
               aes(y = depth, x = lyr_14c)) +
    geom_errorbar(data = lyr_data %>% 
                    filter(id == profnames[p]),
                  aes(y = depth, xmax = lyr_14c*(1+error), xmin = lyr_14c*(1-error))) + 
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black")) +
    ggtitle(profnames[p])  +
    scale_y_reverse(limits = c(330,0), expand = c(0,0)) +
    scale_x_continuous(limits = c(-1300,550), expand = c(0,0))
}



print(ggpubr::ggarrange(g[[1]], g[[2]], g[[3]], g[[4]]) )
ggsave(paste0("GausianProcessRegression_TestRun_4Profiles_",nsamples,"runs_Alox_", format(Sys.time(), "%d%m%y"),".jpeg"), 
       width = 12, height = 9)

# Here, we're rebuilding the data x model goodness of fit matrix 
dists =  as_data_frame( d$results %>% dplyr::select(data, model, d) %>% pivot_wider(id_cols = c("data"), names_from=model, values_from = d )) %>% data.frame()
row.names(dists) = dists$data # Add rownames to the matrix
dists = dists[,-1] # remove first column (data)

hclust.method = "ward.D2"
dist.method = "euclidean"

clust = dists %>% dist(method = dist.method) %>% hclust(method = hclust.method)

as.dendrogram(clust) %>% dendextend::set("branches_k_color", k = 7) %>% tidygraph::as_tbl_graph() %>% ggraph::ggraph() + ggraph::geom_edge_elbow(aes())+ ggraph::geom_edge_elbow(aes(colour = col)) + ggraph::geom_node_text(aes(label=label,filter = leaf), angle = 90, hjust = 1, size = 3) + ggraph::theme_graph() + scale_y_continuous(limits = c(-8000,NA))

ggsave(paste0("dendrogram",nsamples,"runs_Alox_", format(Sys.time(), "%d%m%y"),".jpeg"), 
       width = 12, height = 9)

groups = cutree(clust, 7) #NbClust::NbClust(dists, distance = dist.method, method = hclust.method, index = "cindex" )$Best.partition

datagrp = data %>% left_join(cbind.data.frame(data.frame(id = names(groups), group = groups))) # We're joining the dataset with the group obtained from clustering

print(
  ggplot(datagrp, aes(y= depth, x = value)) + geom_point(aes(colour = id)) + geom_line(aes(colour = id)) + facet_grid(cols = vars(group))+ scale_y_reverse() +  theme(legend.position = "none")
)

ggsave(paste0("data_groups",nsamples,"runs_Alox_", format(Sys.time(), "%d%m%y"),".jpeg"), 
       width = 12, height = 9)
print(
ggplot(datagrp, aes(y= depth, x = value)) + geom_point(aes(colour = as.factor(group))) + geom_line(aes( colour = as.factor(group), group = as.factor(id))) +scale_y_reverse() + theme_classic() + scale_colour_manual(breaks = as.factor(1:7), values = colorspace::rainbow_hcl(7, c=90, l=50)[c(3,2,1,4,5,6,7)]) 
)
ggsave(paste0("data_groups_v2",nsamples,"runs_Alox_", format(Sys.time(), "%d%m%y"),".jpeg"), 
       width = 12, height = 9)
