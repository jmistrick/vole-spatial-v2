# load packages
library(here)
library(tidyverse)
library(igraph)

#clear environment
rm(list = ls())


##---------------- LOAD THE DATA ----------------------

#load the fulltrap dataset (make sure it's the most recent version)
## NOTE ## all DP, DT, or S animals are still in the fulltrap dataset
fulltrap <- readRDS(file = "fulltrap_06.28.22.rds") %>%
  filter(month != "may") %>% #drop may data since not all sites had captures
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
  drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
  unite("sts", season, trt, sex, remove=FALSE) #add sts column to match params_summary

# load the network data
overlap_network_list <- readRDS(here("overlap_network_list_04.05.23.RDS"))


##--------------------------------------------

# Create list to store results
wt_net_mets_list <- list()

# Calculate network metrics
for(i in 1:length(overlap_network_list)){

  #length(overlap_network_list) == for each of the 12 sites

  print(names(overlap_network_list)[[i]])
  site <- list()

  for(j in 1:5){

    #(j in 1:5) == for the 5 trapping occasions

    data <- overlap_network_list[[i]][[j]]

    ids <- as.vector(row.names(data)) #tag ids for all the animals on the grid
    month <- rep(names(overlap_network_list[[i]])[j],length(ids)) #capture month

    #create WEIGHTED NETWORK from adjacency matrix
    inet <- graph_from_adjacency_matrix(data, mode="undirected", weighted = TRUE, diag = FALSE)

    #dataframe to hold results per month
    site[[j]] <- data.frame(ids, month)

    #network metrics to calculate
    site[[j]]$wt.deg <- igraph::strength(inet) #this is the sum of all degree weights for a node
    # site[[j]]$norm.wt.deg <- (igraph::strength(inet))/((igraph::gorder(inet))-1) #your strength/(total nodes-you)
    site[[j]]$avg.wt.deg <- rep(mean(site[[j]]$wt.deg), length(ids)) #calculate average wt degree for a site/occasion
    site[[j]]$n.components <- rep(igraph::count_components(inet), length(ids))
    site[[j]]$n.node <- rep(igraph::gorder(inet), length(ids))
    site[[j]]$wt.n.edge <- rep(sum(E(inet)$weight), length(ids)) ##### THIS IS WEIGHTED #####

    # site[[j]]$assort.wt.deg <- rep(igraph::assortativity_degree(inet, directed=FALSE), length(ids))

  }

  #write the list 'site' as a 1st-order item in wt_net_mets_list
  wt_net_mets_list[[i]] <- site

}

#name the 12 1st order elements of wt_nets_list as the sites
names(wt_net_mets_list) <- names(overlap_network_list)

#rename the sublist items (months) for each site
#all sites only have 5 months of data (May omitted)
for(i in 1:length(wt_net_mets_list)){
  names(wt_net_mets_list[[i]]) <- c("june", "july", "aug", "sept", "oct")
}

################################### ABOUT wt_net_mets_list #############################################
#the output of wt_net_mets_list is a list of 12 1st-order items, 1 per site
#under each site, there are 5 second-order items, 1 per trapping occasion (named by month)
#each of those 2nd-order items (months) is a df with tag ID, month, and all the network metrics...
#for ONLY THE ANIMALS captured on that grid, during that occasion
########################################################################################################


############## condense wt_net_mets_list to make it easier to use for analysis #######################

#collate results
#this collapses the 2nd order elements (network metrics for a single month) down to the 1st order element (site)
#so now wt_net_mets_list_summary is a list of 12 dfs, each is all the network metrics for a site across all the months

#make a list to store things
wt_net_mets_list_summary <- list()

#loop across all sites and collapse the dfs per occasion into one df for the site
for(i in 1:length(wt_net_mets_list)){

  #for all 12 sites
  summary <- do.call("rbind", wt_net_mets_list[[i]])
  wt_net_mets_list_summary[[i]] <- summary
}

#remove the row names of the dfs (its the month, but we already have that info)
for(i in 1:length(wt_net_mets_list_summary)){
  row.names(wt_net_mets_list_summary[[i]]) <- NULL
}

#name the 12 1st order elements as their sites
names(wt_net_mets_list_summary) <- names(overlap_network_list)

### HOT DAMN - wt_net_mets_list_summary is now a list of 12 elements where each element represents a site...
# and contains a df of all the individuals on the grid that month, their tag number,
#and all their network metrics

## make net_mets_list_summary into freiggein huge df
wt_net_mets_summary <- do.call(rbind.data.frame, wt_net_mets_list_summary)

#clean up the df
wt_net_mets_summary <- wt_net_mets_summary %>%
  rownames_to_column("name") %>% #row names are the sites, make that a column
  separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
  mutate(site = as.factor(site)) %>% #make site a factor
  rename(tag=ids)

#save it
saveRDS(wt_net_mets_summary, here("wt_net_mets_summary_04.05.23.rds"))
