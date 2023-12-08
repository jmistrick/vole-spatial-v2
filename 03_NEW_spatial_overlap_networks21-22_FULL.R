##### NEW VERSION 6/8/23 - includes BREEDERS and NONBREEDERS ##########


# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)

#clear environment
rm(list = ls())


#######----------------- VOLE CAPTURE DATA CLEANING -----------------------###############

### May 04 2023 version of vole capture and week recap data are THE MOST UP-TO-DATE VERSIONS ###
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
ft21 <- readRDS(here("fulltrap21_05.10.23.rds")) %>% ##breeders and nonbreeders are here
  unite(sb, sex, season_breeder, remove = FALSE)
ft22 <- readRDS(here("fulltrap22_05.10.23.rds")) %>% ##breeders and nonbreeders are here
  unite(sb, sex, season_breeder, remove = FALSE)

# #minimum number of animals in may=0, max=5
# ft21 %>% filter(month=="may") %>% group_by(site) %>% summarise(n = n_distinct(tag))
# ft22 %>% filter(month=="may") %>% group_by(site) %>% summarise(n = n_distinct(tag))


#######---------------- GENERATE H-R DISTRIBUTION PARAMETERS -----------------------###############
#######------------------ CREATE SPATIAL OVERLAP NETWORKS -----------------------###############
#######--------------------- CALCULATE NETWORK METRICS -----------------------###############

#call the functions
source(here("02_FUNCTIONS_construct_overlap_networks_SBT.R"))

#run for 2021 data

# ##----- OPTION - LIFETIME CENTROIDS -----------
# ## this function is only in the "...seasontrt.R" version of 02_FUNCTIONS - NOT what is pulled here
# #run for 2021 data
# centroid_calc(data = ft21,
#               centroids_file = "centroids21.rds")
#
# # centroids21 <- readRDS(here("centroids21.rds"))

generate_params(data = ft21,
                params_file = "params21_STSB.rds")

params21 <- readRDS(here("params21_STSB.rds"))

create_overlap_networks(data = ft21,
                        centroids_file = "centroids21_STSB.rds",
                        params_file = "params21_STSB.rds",
                        networks_file = "overlapnets21_STSB.rds")

# centroids21 <- readRDS(here("centroids21_STSB.rds"))
# overlapnets21 <- readRDS(here("overlapnets21_STSB.rds"))

calculate_network_metrics(data=ft21,
                          networks_file = "overlapnets21_STSB.rds",
                          netmets_file = "netmets21_bindeg.rds")

netmets21 <- readRDS(here("netmets21_bindeg.rds"))

test %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>%
  mutate(site = factor(site, levels=c("asema", "helmipollo", "hevonen",
                                      "ketunpesa", "kiirastuli", "mustikka",
                                      "kuoppa", "radio", "vaarinkorpi",
                                      "janakkala", "puro", "talo"))) %>%
  ggplot(aes(x=norm.wt.deg)) +
  geom_histogram() +
  facet_grid(site ~ month)


# netmets21 <- readRDS(here("netmets21_STSB.rds"))

#run for 2022 data

# ##----- OPTION - LIFETIME CENTROIDS -----------
# ## this function is only in the "...seasontrt.R" version of 02_FUNCTIONS - NOT what is pulled here
# centroid_calc(data = ft22,
#               centroids_file = "centroids22.rds")
#
# # centroids22 <- readRDS(here("centroids22.rds"))

generate_params(data = ft22,
                params_file = "params22_STSB.rds")

# params22 <- readRDS(here("params22_STSB.rds"))

create_overlap_networks(data = ft22,
                        params_file = "params22_STSB.rds",
                        centroids_file = "centroids22_STSB.rds",
                        networks_file = "overlapnets22_STSB.rds")

# centroids22 <- readRDS(here("centroids22_STSB.rds"))
# overlapnets22 <- readRDS(here("overlapnets22_STSB.rds"))

calculate_network_metrics(data=ft22,
                          networks_file = "overlapnets22_STSB.rds",
                          netmets_file = "netmets22_STSB.rds")

# netmets22 <- readRDS(here("netmets22_STSB.rds"))



#######----------------- NETWORKS PLOTTING (to keel you!) -----------------------###############

library(ggraph)
library(tidygraph)

overlap_network_list <- readRDS(here("overlapnets21_STSB.rds"))
# overlap_network_list <- readRDS(here("overlapnets22_STSB.rds"))

##requires ft files too
metadata <- readRDS(here("fulltrap21_05.10.23.rds")) %>% ##breeders and nonbreeders are here
  unite(sb, sex, season_breeder, remove = FALSE)
# data <- readRDS(here("fulltrap22_05.10.23.rds")) %>% ##breeders and nonbreeders are here
#   unite(sb, sex, season_breeder, remove = FALSE)

######## PLOT MULTIPLE SITES ACROSS MONTHS #############


for(i in 1:length(overlap_network_list)) {

  # set.seed(2111994)

  png(filename = paste("spatial_overlap_", "THEprettiest_", names(overlap_network_list)[[i]], "_2021", ".png", sep = ""),
      width=10, height=3, units="in", res=600)

  par(mfrow = c(1,5))

  for(j in 1:length(overlap_network_list[[i]])){

    data <- overlap_network_list[[i]][[j]]

    site.id <- names(overlap_network_list[i])
    month.id <- names(overlap_network_list[[i]][j])

    ## **PLOTTING FROM TIDYGRAPH OBJECT
    #metadata to get sex for node color
    #subset metadata
    netmeta <- metadata %>% #starts with fulltrap, need to get down to a 'traits' version
      filter(site==site.id & month==month.id) %>%
      group_by(tag) %>% slice(1) %>% #one entry per vole per month
      drop_na(sex) %>% #didn't build networks with animals with sex=NA
      drop_na(season_breeder) %>% #didn't build networks with animals with season_breeder=NA
    ###dropping both sex and season_breeder = NA should clean up sb so there are no NAs
      select(tag, sex, season_breeder, sb) #tag should be first column to match to adj mat

    ## **PLOTTING FROM TIDYGRAPH OBJECT
    #create tidygraph object from adj matrix
    tidyg <- as_tbl_graph(data, directed=FALSE) %>% left_join(netmeta, by=c('name'='tag'))
    tg <- tidyg %>% activate(edges) %>% arrange(weight) %>% filter(weight>0)

    # # **PLOTTING FROM iGRAPH OBJECT
    # #create graph from adj matrix (igraph verbs)
    # g <- graph_from_adjacency_matrix(
    #   data,
    #   mode = c("undirected"),
    #   weighted = TRUE,
    #   diag = FALSE)

    #--------------ignore this------------------
    # #for thresholded edges option1
    # #remove edges with weight less than threshold
    # g2 <- delete.edges(g, which(E(g)$weight<=0.05))

    # #for thresholded edges option2
    # #remove edges with weight less than threshold
    # g2 <- delete.edges(g, which(E(g)$weight<=0.1))

    # #plot thresholded edges opt1 or opt2
    # plot(g2, vertex.size=5, vertex.label=NA,
    #      main = paste(names(overlap_network_list[[i]])[j]))

    #keep the layout the same
    #https://www.kateto.net/wp-content/uploads/2016/01/NetSciX_2016_Workshop.pdf
    # l <- layout_with_fr(g2)

    # #for edges of varying thickness
    # plot(g2, vertex.size=5, vertex.label=NA,
    #      # layout=l,
    #      edge.width = ((E(g2)$weight)*5),
    #      edge.color = "#545454",
    #      main = paste(names(overlap_network_list[[i]])[j]))
    #-----------------------------------------------------------

   ## Matt M-S suggestion - color edges by weight, thickness by weight
    ##USING GGRAPH for plotting!

    # # **PLOTTING FROM iGRAPH OBJECT
    # #metadata to get sex for node color
    # #subset metadata
    # netmeta <- metadata %>% #starts with fulltrap, need to get down to a 'traits' version
    #   filter(site==site.id & month==month.id) %>%
    #   group_by(tag) %>% slice(1) %>% #one entry per vole per month
    #   drop_na(sex) %>% #didn't build networks with animals with sex=NA
    #   drop_na(season_breeder) #didn't build networks with animals with season_breeder=NA
    # ###dropping both sex and season_breeder = NA should clean up sb so there are no NAs
    # #set "sex" as a vertex attribute
    # g <- set.vertex.attribute(graph=g, name="sex", index=V(g), value=netmeta$sex)
    # #set "breeder" as a vertex attribute
    # g <- set.vertex.attribute(graph=g, name="breeder", index=V(g), value=netmeta$season_breeder)
    # #set "sb" as vertex attribute
    # g <- set.vertex.attribute(graph=g, name="sb", index=V(g), value=netmeta$sb)

    ggraph(tg, layout="fr") +
      geom_edge_link(aes(colour=weight, width=weight, group=I(1))) + # add edges to the plot (colored by weight)
      ## even with sorting the edges by weight, they won't necessarily plot with thickest on top
        ## unless you either 1) use 'geom_edge_link0() or 2) add 'group=I(1)' to the aes call
        ## Matt M-S is a lifesaver <3
      # geom_node_label(aes(label=name)) + # add nodes to the plot
      geom_node_point(aes(color=sb), size=7) +
      scale_edge_width(range=c(0,3), guide="none") + #scale edge width by weight
      scale_edge_colour_gradient(low="#F0F0F0", high="#000000") + # set the (gray)scale
      theme_void() +
      labs(title=paste(names(overlap_network_list[[i]])[j]))

    ## feed layout= a matrix with the node position - make sure it's in the same order as the nodes
    ## numeric matrix two column x, y for all the nodes
    ## ?? arc-ed edges or jitter to points to avoid overlap

  }

  dev.off()

}



######################

for(i in 1:length(overlap_network_list)) {

  set.seed(2111994)

  png(filename = paste("spatial_overlap_", "STSB_wtbyTHICK3_", names(overlap_network_list)[[i]], "_2021", ".png", sep = ""),
      width=10 , height=3, units="in", res=600)

  par(mfrow = c(1,5))

  for(j in 1:length(overlap_network_list[[i]])){

    data <- overlap_network_list[[i]][[j]]

    #plot adj matrix for network
    g <- graph_from_adjacency_matrix(
      data,
      mode = c("undirected"),
      weighted = TRUE,
      diag = FALSE)

    # #for thresholded edges option1
    # #remove edges with weight less than threshold
    # g2 <- delete.edges(g, which(E(g)$weight<=0.05))

    # #for thresholded edges option2
    # #remove edges with weight less than threshold
    # g2 <- delete.edges(g, which(E(g)$weight<=0.1))

    # #plot thresholded edges opt1 or opt2
    # plot(g2, vertex.size=5, vertex.label=NA,
    #      main = paste(names(overlap_network_list[[i]])[j]))

    #keep the layout the same
    #https://www.kateto.net/wp-content/uploads/2016/01/NetSciX_2016_Workshop.pdf
    l <- layout_with_fr(g)

    #for edges of varying thickness
    plot(g, vertex.size=5, vertex.label=NA,
         layout=l,
         edge.width = ((E(g)$weight)*3),
         edge.color = "#545454",
         main = paste(names(overlap_network_list[[i]])[j]))

  }

  dev.off()

}
