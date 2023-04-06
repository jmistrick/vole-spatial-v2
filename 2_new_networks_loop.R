### UNLESS YOU NEED TO CHANGE SOMETHING - this code was run 4.5.23 and saved to RDS
    # load data using the "overlap_network_list..." R file

###LARGETRAP VERSION
#options here to run code using just the 'largetrap' animals
#essentially reproductive and or large (>15g or >17g) voles


# load packages
library(here)
library(tidyverse)
library(igraph)
# library(lubridate)
# library(janitor)

#clear environment
rm(list = ls())

##---------------- LOAD THE CAPTURE DATA ----------------------

#load the fulltrap dataset (make sure it's the most recent version)
  ## NOTE ## all DP, DT, or S animals are still in the fulltrap dataset
fulltrap <- readRDS(file = "fulltrap_06.28.22.rds") %>%
  filter(month != "may") %>% #drop may data since not all sites had captures
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
  drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
  unite("sts", season, trt, sex, remove=FALSE) #add sts column to match params_summary

# ###LARGETRAP VERSION
# fulltrap <- readRDS(file = "fulltrap_06.28.22.rds") %>%
#   filter(month != "may") %>% #drop may data since not all sites had captures
#   mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
#   drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
#   filter(mass >= 15 |
#            per=="1" | nip=="1" | preg=="1" |
#            test=="1") %>%
#   filter(caps_per_life > 1) %>%
#   unite("sts", season, trt, sex, remove=FALSE) #add sts column to match params_summary

#pull all the traps and their x, y coordinates, save as df with (trapID, x, y)
traps <- fulltrap %>% group_by(trap) %>% slice(1) %>%
  select(trap, x, y) %>%
  arrange(trap)

#save a vector of the site names (alphabetical order)
site_names <- unique(fulltrap$site) %>% sort()


##---------------- CREATE A LIST of CAPTURE DATA (nested by site, month) --------------------

# #split() makes a list consisting of individual data.frames based on a condition ('site' in this case)
# site_list <- split(fulltrap, fulltrap$site)
#
# #create new list to hold nested site, month capture
# sitemonth_list <- list()
#
# for(i in 1:length(site_list)){
#   # print(i)
#   temp_list <- split(site_list[[i]], site_list[[i]]$month)
#   sitemonth_list[[i]] <- temp_list
# }
#
# #name 1e list element as site (2e list elements are months May-Oct)
# names(sitemonth_list) <- site_names
#
# ### OUTPUT: SITEMONTH_LIST is a nested list of length 12
#   #1e level is all the sites (12)
#   #2e level is all the months per site (5) excluding May

##-------------- LOAD ADDITIONAL DATA ----------------------

# # Reading in required functions in wanelik_farine_functions.R
# source(here("wanelik_farine_functions.R"))
#
# #read in the params RDS
# params_summary <- readRDS(file = here("params_summary_04.04.23.rds"))
#
# # ###LARGETRAP VERSION
# # #read in the params RDS
# # params_summary <- readRDS(file = here("params_LARGETRAP_summary_04.05.23.rds"))

##--------------- (in a loop) GENERATE OVERLAP NETWORK ---------------------

# overlap_network_list <- list()
#
# for(i in 1:length(sitemonth_list)){
#
#   print(names(sitemonth_list[i])) #print the site name
#
#   overlap_network_list[[i]] <- list()
#
#   for(j in 1:5){
#     #1:5 for the five months (MAY DATA OMITTED)
#
#     print(names(sitemonth_list[[i]][j])) #print the month
#
#     data <- sitemonth_list[[i]][[j]] %>%
#       select(tag, sts, trap, x, y) %>%
#       group_by(tag, trap) %>%
#       mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
#       ungroup()
#
#     tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
#     tagsts <- data %>% group_by(tag) %>% slice(1) %>% select(tag, sts) #all tags with their STS (season/trt/sex)
#
#     # ###LARGETRAP VERSION (since Mustikka only has 1 animal in June)
#     # if(nrow(tags)=="1") {next}
#
#     traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
#     tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap
#
#     #creates VERY LARGE df of all tags for all traps
#     #new column of Det.obs - whether that animal was observed in that trap
#     #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
#     fulldata <- cbind(tags_rep, traps_rep) %>%
#       left_join(data, by=c("tag", "trap", "x", "y")) %>%
#       group_by(tag) %>% fill(sts, .direction="downup") %>%
#       mutate(Det.count = replace_na(Det.count, 0),
#              Det.obs = ifelse(Det.count==0, 0, 1)) %>%
#       arrange(tag, trap) %>%
#       select(tag, sts, x, y, Det.obs, Det.count) %>%
#       rename(x.trap = x,
#              y.trap = y,
#              Tag_ID = tag)
#
#     #----------------------- 3. Generating overlap network  ------------------------------#
#
#     # Recalculating (weighted) centroids
#     # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
#     matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
#     matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
#     matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
#     matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
#     matrix_dists_real2<-matrix_dists_real2[,-7] #drops the Det.count2 column
#     x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
#     y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
#     obs_centroids <- data.frame(Tag_ID=tagsts$tag, sts=tagsts$sts)
#     obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
#     obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]
#
#     # Creating a matrix of distances between each individual's observed centroid and each trap
#     matrix_dists_obs <- fulldata #only one entry per individual
#     matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #x coord of centroid of individual
#     matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #y coord of centroid of individual
#     matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
#     matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)
#
#     # # Removing unobserved individuals
#     # matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] #not needed, matrix is only of voles capped that month
#
#     # Calculating total number of detections per individual
#     ### TBH Janine isn't sure this is needed because it doesn't appear to be used anywhere, but a good check I guess
#     a <- tapply(matrix_dists_obs$Det.count, list(matrix_dists_obs$Tag_ID), sum)
#     matrix_dists_obs$Det.count.total <- a[match(matrix_dists_obs$Tag_ID, names(a))]
#
#     # Estimating home range profiles (negative sigmoidal curves) done separately by season/trt -> params_summary file
#
#     # Generating overlap network
#     as <- params_summary$a[match(obs_centroids$sts,params_summary$sts)]
#     bs <- params_summary$b[match(obs_centroids$sts,params_summary$sts)]
#
#     overlap_network_sep <- get_network_2D(obs_centroids$x[which(!is.na(obs_centroids$x))], obs_centroids$y[which(!is.na(obs_centroids$y))],as,bs)
#     rownames(overlap_network_sep)=colnames(overlap_network_sep)=na.omit(obs_centroids)$Tag_ID
#
#     overlap_network_list[[i]][[j]] <- overlap_network_sep
#
#   }
#
#
# }
#
#
# #name the 12 1st order elements of overlap_network_list as the sites
# names(overlap_network_list) <- names(sitemonth_list)
#
# #rename the sublist items (months) for each site
# for(i in 1:length(overlap_network_list)){
#   names(overlap_network_list[[i]]) <- c("june", "july", "aug", "sept", "oct")
# }
#
# # ###LARGETRAP VERSION (remove null 2e list item at Mustikka)
# # overlap_network_list$mustikka <- overlap_network_list$mustikka %>% discard(is.null)
#
# #save it since that code takes 5ever to run
# saveRDS(overlap_network_list, here("overlap_network_list_04.05.23.rds"))
#
# # ###LARGETRAP VERSION
# # #save it since that code takes 5ever to run
# # saveRDS(overlap_network_list, here("LARGETRAP_overlap_network_list_04.05.23.rds"))

# load the network data
overlap_network_list <- readRDS(here("overlap_network_list_04.05.23.RDS"))



##-------------- PLOTTING (to keel you!) ----------------


######## PLOT MULTIPLE SITES ACROSS MONTHS #############

for(i in 1:length(overlap_network_list)) {

  png(filename = paste("spatial_overlap_", "wt0.05_", names(overlap_network_list)[[i]], "_2021", ".png", sep = ""),
      width=10 , height=3, units="in", res=600)

  par(mfrow = c(1,5))

  for(j in 1:length(overlap_network_list[[i]])){ #adjust for LARGETRAP (since Mustikka has no June)

    data <- overlap_network_list[[i]][[j]]

    #plot adj matrix for network
    g <- graph_from_adjacency_matrix(
      data,
      mode = c("undirected"),
      weighted = TRUE,
      diag = FALSE)

    #remove edges with weight less than threshold
    g2 <- delete.edges(g, which(E(g)$weight<=0.05))

    # #remove edges with weight less than threshold
    # g2 <- delete.edges(g, which(E(g)$weight<=0.1))

    plot(g2, vertex.size=5, vertex.label=NA,
         main = paste(names(overlap_network_list[[i]])[j]))

  }

  dev.off()

}


# ############# ONE SITE ACROSS MULTIPLE MONTHS ###############
# ######## (intermediate version to get to the above) ##########
#
# #save each list element as its own list
# for(i in 1:length(overlap_network_list)){
#   data <- overlap_network_list[[i]] #pull the 1e list elements (sites)
#   assign(paste(names(overlap_network_list)[i]), data) #name it as the site name
# }
#
#
# #save multiple plots using par(mfrow)
# #https://stackoverflow.com/questions/18584722/save-multiple-plots-in-r-as-a-jpg-file-how
# png(filename = paste("spatial_overlap_", "helmipollo", "_2021", ".png", sep = ""),
#     width=10 , height=3, units="in", res=600)
#
# par(mfrow = c(1,5))
#
# for(i in 1:5){
#
#   data <- helmipollo[[i]]
#
#   #plot adj matrix for network
#   g <- graph_from_adjacency_matrix(
#     data,
#     mode = c("undirected"),
#     weighted = TRUE,
#     diag = FALSE)
#
#   g2 <- delete.edges(g, which(E(g)$weight<0.05))
#
#   plot(g2, vertex.size=5, vertex.label=NA,
#        main = paste(names(helmipollo)[i]))
#
# }
#
# dev.off()


##----------------------- END ---------------------------------
