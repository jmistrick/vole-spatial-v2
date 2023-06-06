
# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)

#clear environment
rm(list = ls())


#######----------------- VOLE CAPTURE DATA CLEANING -----------------------###############

### May 04 2023 version of vole capture and week recap are THE MOST UP-TO-DATE VERSIONS ###
                      ### May 10 2023 versions of FULLTRAP correspond ###

# #pull the data cleaning function
# source(here("01_FUNCTION_data_cleaning.R"))
#
# #load the FULL (2021 AND 2022) vole processing and WR data (pulled from Box)
#
# #clean 2021 data
# clean_data(processingdata = "vole_capture_data_05.04.23.csv",
#            WRdata= "week_recap_data_05.04.23.csv",
#            yr=2021,
#            fulltrap_output = "fulltrap21_05.10.23.rds")
#
# #clean 2022 data
# clean_data(processingdata = "vole_capture_data_05.04.23.csv",
#            WRdata= "week_recap_data_05.04.23.csv",
#            yr=2022,
#            fulltrap_output = "fulltrap22_05.10.23.rds")

#load the fulltrap 21 and 22 data (make sure it's the most recent version)
    ## NOTE ## all DP, DT, or S animals are still in the fulltrap datasets
ft21 <- readRDS(here("fulltrap21_05.10.23.rds")) %>% ###### edit on 5/11 - ONLY breeders, add this to cleaning if I like it
  filter(season_breeder=="breeder")
ft22 <- readRDS(here("fulltrap22_05.10.23.rds")) %>%
  filter(season_breeder=="breeder")


#######----------------- CONSTRUCT SPATIAL OVERLAP NETWORKS -----------------------###############

### CODE RUN ON 6.2.23 with full-year centroids ###

#call the functions
source(here("02_FUNCTIONS_construct_overlap_networks_seasontrt.R"))

# ##----- OPTION - LIFETIME CENTROIDS -----------
# #run for 2021 data
# centroid_calc(data = ft21,
#               centroids_file = "centroids21.rds")
#
# # centroids21 <- readRDS(here("centroids21.rds"))

generate_params(data = ft21,
                centroids_file = "centroids21.rds",
                params_file = "params21.rds")

# params21 <- readRDS(here("params21.rds"))
# centroids21 <- readRDS(here("centroids21.rds"))

create_overlap_networks(data = ft21,
                        centroids_file = "centroids21.rds",
                        params_file = "params21.rds",
                        networks_file = "overlapnets21.rds")

# overlapnets21 <- readRDS(here("overlapnets21.rds"))

calculate_network_metrics(data=ft21,
                          networks_file = "overlapnets21.rds",
                          netmets_file = "netmets21.rds")

# netmets21 <- readRDS(here("netmets21.rds"))

#run for 2022 data

# ##----- OPTION - LIFETIME CENTROIDS -----------
# centroid_calc(data = ft22,
#               centroids_file = "centroids22.rds")
#
# # centroids22 <- readRDS(here("centroids22.rds"))

generate_params(data = ft22,
                centroids_file = "centroids22.rds",
                params_file = "params22.rds")

# params22 <- readRDS(here("params22.rds"))
# centroids22 <- readRDS(here("centroids22.rds"))

create_overlap_networks(data = ft22,
                        params_file = "params22.rds",
                        centroids_file = "centroids22.rds",
                        networks_file = "overlapnets22.rds")

# overlapnets22 <- readRDS(here("overlapnets22.rds"))

calculate_network_metrics(data=ft22,
                          networks_file = "overlapnets22.rds",
                          netmets_file = "netmets22.rds")

# netmets22 <- readRDS(here("netmets22.rds"))

##-------------- PLOTTING (to keel you!) ----------------

# overlap_network_list <- readRDS(here("overlapnets21.rds"))
# overlap_network_list <- readRDS(here("overlapnets22.rds"))

######## PLOT MULTIPLE SITES ACROSS MONTHS #############

# for(i in 1:length(overlap_network_list)) {
#
#   png(filename = paste("NEW_spatial_overlap_", "breeder_wt0.1_", names(overlap_network_list)[[i]], "_2022", ".png", sep = ""),
#       width=10 , height=3, units="in", res=600)
#
#   par(mfrow = c(1,5))
#
#   for(j in 1:length(overlap_network_list[[i]])){
#
#     data <- overlap_network_list[[i]][[j]]
#
#     #plot adj matrix for network
#     g <- graph_from_adjacency_matrix(
#       data,
#       mode = c("undirected"),
#       weighted = TRUE,
#       diag = FALSE)
#
#     # #for thresholded edges option1
#     # #remove edges with weight less than threshold
#     # g2 <- delete.edges(g, which(E(g)$weight<=0.05))
#
#     #for thresholded edges option2
#     #remove edges with weight less than threshold
#     g2 <- delete.edges(g, which(E(g)$weight<=0.1))
#
#     #plot thresholded edges opt1 or opt2
#     plot(g2, vertex.size=5, vertex.label=NA,
#          main = paste(names(overlap_network_list[[i]])[j]))
#
#     # #for edges of varying thickness
#     # plot(g, vertex.size=5, vertex.label=NA,
#     #      # edge.width = ((E(g)$weight)*3),
#     #      main = paste(names(overlap_network_list[[i]])[j]))
#
#   }
#
#   dev.off()
#
# }
