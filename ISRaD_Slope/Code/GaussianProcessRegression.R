###Updated Gaussian Process Regression Script###
###2022-02-16###
###Based on Script from Xavier Raynaud###

source("./Code/GaussianProcess_functions_v2.R") # Update code. 

#This is to load your data table. Let's assume is creates an object called data
#devtools::install_github("International-Soil-Radiocarbon-Database/ISRaD/Rpkg", 
#                         ref = "master", force = TRUE)

library(ISRaD)
library(tidyverse)

lyr_data <- readRDS(paste0(getwd(), "/Data/ISRaD_lyr_data_filtered_2022-06-21"))

#Filter ISRaD data
lyr_data_gpr <- lyr_data %>% 
  filter(depth <= 200) %>% 
  mutate(error = 0.1) %>%
  group_by(id) %>%
  mutate(scale_value = scale(lyr_14c)) %>% 
  mutate(norm_value = (lyr_14c-min(lyr_14c))/(max(lyr_14c)-min(lyr_14c))) %>%
  filter(n() > 2) %>%
  ungroup() 

lyr_data_gpr %>% 
  count(id)

#Visual exploration of the ISRaD data
alpha = 0.2
lyr_data_gpr %>% 
  ggplot() +
  geom_line(aes(x = scale_value, y = depth, group = id, color = id), 
            orientation = "y", alpha = alpha) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous("Depth [cm]", trans = "reverse", expand = c(0,0)) +
  scale_x_continuous("Delta 14C")

#Select all profiles that have alox measured
data_alox <- lyr_data_gpr %>% 
  drop_na(lyr_al_ox) %>% 
  mutate(alox = case_when(
    lyr_al_ox == 0 ~ 0.001,
    TRUE ~ lyr_al_ox
  )) %>% 
  #remove data from Gentsch_2018 (causing problems: check why: depth issues - permafrost soil)
  filter(!grepl("Gentsch", id)) 

data_alox %>% 
  count(id)

#Data should be a dataframe or a tibble in long format, with at least one column
#to identify the profile ('id'), one column containing the depth ('depth'), one 
#containing the 14C data ('value') and one containing the measurement error ('error', can be empty). 

nsamples <- 100

#ncores is the number of cpus you have on your computer. If you put a value 
#larger than the true number of cpus you have in your computer, it will be 
#reduced to the number of cpus minus 1. 

# The compareGP output has slightly changed: you can select whether you keep simulations or not (This make the d object waaaayyy smaller if you don't).
# The way to save the goodness of fit results has also changed. Instead of an nxn matrix, results are saved in the form of a long table with 3 columns, data, model and goodness of fit index. 
d <- compare_GP(lyr_data_gpr, 
                idcol = "id",
                xcol = "CORG", 
                ycol = "lyr_14c", 
                errcol = "error", 
                nsamples = nsamples, 
                verbose = 0, 
                engine = "laGP_gaupro", 
                ncores = 20, keep_simulations = FALSE)

#Save the file with a filename containing the date. For 600 profiles, the file is ~1Gb.
save(d, file = paste0("./Data/GPR_TestRun_SOC_all_", nsamples,"_",
                      format(Sys.time(), "%Y%m%d"), ".RData"))

str(d)

d$models$`Spielvogel_2008_National Park Bayrischer Wald_Leptic Cambisol`$plot()

d$models$Lawrence_2015_EarlyEvansCreek1_CW04_1979$plot()

alpha <- 0.1
g <- list()
profnames <- c("Torn_1997_Kokee (4.1my)_Kokee (4.1my)_profile_1",
              "Spielvogel_2008_National Park Bayrischer Wald_Leptic Cambisol", 
              "Trumbore_1996_Trumbore Shaver_Trumbore Shaver 2 1959", 
              "Khomo_2017_GR-740-C_GR-740-C")

for (p in 1: length(profnames)) {
  g[[p]] = d$simuls[[profnames[p]]] %>% 
    data.frame() %>% 
    pivot_longer(cols = -1, names_to = "runs", values_to = "C14values") %>% 
    group_by(CORG) %>% 
    mutate(mean = mean(C14values),
           sd = sd(C14values),
           lci = t.test(C14values, conf.level = 0.95)$conf.int[1],
           uci = t.test(C14values, conf.level = 0.95)$conf.int[2]) %>% 
    ungroup() %>% 
    ggplot() +
    geom_line(aes(x = CORG, y = C14values, group = runs), alpha = alpha,
              color = "darkgrey", orientation = "y") +
    geom_line(aes(x = CORG, y = mean),
              orientation = "x", color = "red", size = 0.7) +
    geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd, x = CORG), alpha = 0.2,
                fill = "red", color = "black") +
    geom_point(data = data_alox %>% 
                 filter(id == profnames[p]),
               aes(x = CORG, y = lyr_14c)) +
    geom_errorbar(data = data_alox %>% 
                     filter(id == profnames[p]),
                   aes(x = CORG, 
                       ymax = lyr_14c*(1+error), 
                       ymin = lyr_14c*(1-error))) + 
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black")) +
    ggtitle(profnames[p])  +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(expand = c(0,0))
}

print(ggpubr::ggarrange(g[[1]], g[[2]], g[[3]], g[[4]]))

# ggsave(paste0("GausianProcessRegression_TestRun_4Profiles_", nsamples, "runs_Alox_", 
#               format(Sys.time(), "%d%m%y"),".jpeg"), 
#        width = 12, height = 9)

# Here, we're rebuilding the data x model goodness of fit matrix 
dists <- as_tibble(d$results %>% 
                     dplyr::select(data, model, d) %>% 
                     pivot_wider(id_cols = c("data"), 
                                 names_from = model, 
                                 values_from = d)) %>% 
  data.frame()

str(dists)
head(dists)

row.names(dists) <- dists$data # Add rownames to the matrix
dists <- dists[,-1] # remove first column (data)

hclust.method <- "ward.D2"
dist.method <- "euclidean"

clust <- dists %>% 
  dist(method = dist.method) %>% 
  hclust(method = hclust.method)

plot(clust)
summary(clust)

best_numbers <- NbClust::NbClust(dists, distance = dist.method, method = hclust.method, index = "cindex" )
best_numbers$All.index
best_numbers$Best.nc
best_numbers$Best.partition

as.dendrogram(clust) %>% 
  dendextend::set("branches_k_color", k = best_numbers$Best.nc[1]) %>% 
  tidygraph::as_tbl_graph() %>% 
  ggraph::ggraph() + 
  ggraph::geom_edge_elbow(aes()) + 
  ggraph::geom_edge_elbow(aes(colour = col)) + 
  ggraph::geom_node_text(aes(label = label, filter = leaf), angle = 90, 
                         hjust = 1, size = 3) + 
  ggraph::theme_graph() + 
  scale_y_continuous(limits = c(-900,NA))

# ggsave(paste0("dendrogram",nsamples,"runs_Alox_", 
#               format(Sys.time(), "%d%m%y"),".jpeg"), 
#        width = 12, height = 9)


groups <- cutree(clust, best_numbers$Best.nc[1]) 

data_grp <- data_alox %>% 
  left_join(cbind.data.frame(data.frame(id = names(groups), group = groups))) # We're joining the dataset with the group obtained from clustering

ggplot(data_grp, aes(x = CORG, y = lyr_14c)) + 
    geom_point(aes(colour = id)) + 
    geom_line(aes(colour = id), orientation = "y") + 
    facet_grid(cols = vars(group)) + 
    scale_x_continuous(trans = "log10") + 
    theme(legend.position = "none")

ggplot(data_grp, aes(x = CORG, y = lyr_14c)) + 
  geom_point(aes(colour = pro_KG_present_long)) + 
  geom_line(aes(group = id, colour = pro_KG_present_long), orientation = "y") + 
  facet_grid(cols = vars(group)) + 
  scale_x_continuous(trans = "log10")

ggplot(data_grp, aes(x = CORG, y = lyr_14c)) + 
  geom_point(aes(colour = factor(group))) + 
  geom_line(aes(group = id, colour = factor(group)), orientation = "y") + 
  facet_grid(cols = vars(pro_KG_present_long)) + 
  scale_x_continuous(trans = "log10")

ggplot(data_grp, aes(x = CORG, y = lyr_14c)) + 
  geom_point(aes(colour = factor(group))) + 
  geom_line(aes(group = id, colour = factor(group)), orientation = "y") + 
  facet_grid(cols = vars(pro_usda_soil_order)) + 
  scale_x_continuous(trans = "log10")

data_grp %>% 
  dplyr::select(id, CORG,lyr_14c, group,
                lyr_al_ox, pro_usda_soil_order, pro_BIO12_mmyr_WC2.1) %>% 
  view()

# ggsave(paste0("data_groups",nsamples,"runs_Alox_", 
#               format(Sys.time(), "%d%m%y"),".jpeg"), 
#        width = 12, height = 9)

ggplot(data_grp, aes(x = CORG, y = lyr_14c)) + 
  geom_point(aes(colour = as.factor(group))) + 
  geom_line(aes(colour = as.factor(group), group = as.factor(id)),
            orientation = "y") +
  scale_x_continuous(trans = "log10") + 
  theme_classic() + 
  scale_colour_manual(breaks = as.factor(1:3), 
                      values = colorspace::rainbow_hcl(3, c = 90, l = 50)[c(3,2,1,4,5,6,7)]) 

# ggsave(paste0("data_groups_v2",nsamples,"runs_Alox_", 
#               format(Sys.time(), "%d%m%y"),".jpeg"), 
#        width = 12, height = 9)
