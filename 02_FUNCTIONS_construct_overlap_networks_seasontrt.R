
## Function for generating the a and b parameters by season/trt/sex to define the distributions for vole HRs
## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
##          centroids_file = file name and extension in " " for the centroids_file (SEASONAL centroid)
##          params_file = file name and extension in " " for the params_file to be generated
## Output: dataframe of parameters (3 columns for seasontrtsex, a param, b param)

## also OPTION TO USE LIFETIME CENTROIDS generated in a separate function --> see commented out code in the fxn

## CENTROIDS here refers to SEASONAL CENTROIDS
generate_params <- function(data, centroids_file, params_file) {

  # ##---------------- LOAD THE LIFETIME CENTROIDS ----------------------
  #
  # obs_centroids <- readRDS(here(centroids_file))

  ##---------------- LOAD THE CAPTURE DATA ----------------------

  #load, clean the fulltrap dataset
  fulltrap <- data %>%
    filter(month != "may") %>% #drop may data since not all sites had captures
    mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
    drop_na(sex) #remove animals with sex=NA (since we can't assign then a HR)

  #pull all the traps and their x, y coordinates, save as df with (trapID, x, y)
  traps <- fulltrap %>% group_by(trap) %>% slice(1) %>%
    select(trap, x, y) %>%
    arrange(trap)

  #save a vector of the site names (alphabetical order)
  site_names <- unique(fulltrap$site) %>% sort()

  #save a vector of trts
  trt <- unique(fulltrap$trt) %>% sort()

  #save a vector of seasons
  season <- unique(fulltrap$season) %>% sort()

  ##---------------- CREATE A LIST of CAPTURE DATA (nested by season, trt) --------------------

  #split() makes a list consisting of individual data.frames based on a condition ('season' in this case)
  season_list <- split(fulltrap, fulltrap$season)

  #create new list to hold nested site, month capture
  seasontrt_list <- list()

  for(i in 1:length(season_list)){
    # print(i)
    temp_list <- split(season_list[[i]], season_list[[i]]$trt)
    seasontrt_list[[i]] <- temp_list
  }

  #name 1e list element as season (2e list elements are trts)
  names(seasontrt_list) <- season

  ### SEASONTRT_LIST is a nested list of length 2
  ### 1e level is season (summer/fall) (2)
  ### 2e level is all the trts per site (4)
  ### all grids per trt and months per season are combined together!!

  ##------------------------------------

# Following code from Wanelik & Farine 2022

  params_list <- list()

  #this is for SEASONAL CENTROIDS
  centroids_list <- list()

  for(i in 1:length(seasontrt_list)){

    print(i)

    params_list[[i]] <- list()

    #this is for SEASONAL CENTROIDS
    centroids_list[[i]] <- list()

    for(j in 1:4){
      #1:4 for the four treatments - not reproducible, but easy

      print(j)

      data <- seasontrt_list[[i]][[j]] %>%
        select(tag, sex, trap, x, y) %>%
        group_by(tag, trap) %>%
        mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
        ungroup()

      season_lab <- unique(seasontrt_list[[i]][[j]]$season)
      trt_lab <- unique(seasontrt_list[[i]][[j]]$trt)

      tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
      tagsex <- data %>% group_by(tag) %>% slice(1) %>% select(tag,sex) #all tags with their sex

      traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
      tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

      #VERY LARGE df of all tags for all traps
      #new column of Det.obs - whether that animal was observed in that trap
      #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
      fulldata <- cbind(tags_rep, traps_rep) %>%
        left_join(data, by=c("tag", "trap", "x", "y")) %>%
        group_by(tag) %>% fill(sex, .direction="downup") %>%
        mutate(Det.count = replace_na(Det.count, 0),
               Det.obs = ifelse(Det.count==0, 0, 1)) %>%
        arrange(tag, trap) %>%
        select(tag, sex, x, y, Det.obs, Det.count) %>%
        rename(x.trap = x,
               y.trap = y,
               Tag_ID = tag,
               Sex = sex)

      #----------------------- 3. Generating overlap network  ------------------------------#

      # OPTION TO REMOVE CENTROID CALCULATION HERE - USE LIFETIME CENTROID (comment out the following block of code)

      # Calculating (weighted) centroids (FOR each VOLE, in a SEASON)
      # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
      matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
      matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
      matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
      matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
      matrix_dists_real2<-matrix_dists_real2[,-7] #drops the Det.count2 column
      x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
      y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
      obs_centroids <- data.frame(Tag_ID=tagsex$tag, Sex=tagsex$sex)
      obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
      obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]

      #----END OPTION TO REMOVE CENTROID CALCULATION----

      # Creating a matrix of distances between each individual's observed centroid and each trap
      matrix_dists_obs <- fulldata #only one entry per individual
      matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #x coord of centroid of individual
      matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #y coord of centroid of individual
      matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
      matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)
      # # Removing unobserved individuals
      # matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] #not needed, matrix is only of voles capped that month

      # Calculating total number of detections per individual
      ### TBH Janine isn't sure this is needed because it doesn't appear to be used anywhere, but a good check I guess
      a <- tapply(matrix_dists_obs$Det.count, list(matrix_dists_obs$Tag_ID),sum)
      matrix_dists_obs$Det.count.total <- a[match(matrix_dists_obs$Tag_ID,names(a))]

      # Estimating home range profiles (negative sigmoidal curves) using this data
      # Example: Sex-specific curves (male home range profile = fit.males; female home range profile = fit.females)
      fit.males <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Sex=="M"),], family=binomial, control = list(maxit = 50))
      fit.females <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Sex=="F"),], family=binomial, control = list(maxit = 50))

      # Generating overlap network
      params <- data.frame(sex=c("M","F"), a=c(coef(fit.males)[1],coef(fit.females)[1]), b=c(coef(fit.males)[2],coef(fit.females)[2]))

      params <- params %>% mutate(season=season_lab,
                                  trt=trt_lab)

      #this is for SEASONAL CENTROIDS
      obs_centroids <- obs_centroids %>% mutate(season=season_lab,
                                                trt=trt_lab)

      params_list[[i]][[j]] <- params

      #this is for SEASONAL CENTROIDS
      centroids_list[[i]][[j]] <- obs_centroids
    }

  }

  #collate results
  #this collapses the 2nd order elements (params per trt) down to the 1st order element (season)

  #make a list to store things
  params_list_summary <- list()

  #this is for SEASONAL CENTROIDS
  centroids_list_summary <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(params_list)){

    #for both seasons
    summary <- do.call("rbind", params_list[[i]])
    params_list_summary[[i]] <- summary
  }

  #this is for SEASONAL CENTROIDS
  for(i in 1:length(centroids_list)){
    #for both seasons
    summary <- do.call("rbind", centroids_list[[i]])
    centroids_list_summary[[i]] <- summary
  }

  #collapse params_list_summary into a df
  params_summary <- do.call(rbind.data.frame, params_list_summary) %>%
    unite("sts", season, trt, sex)
  row.names(params_summary) <- NULL

  #this is for SEASONAL CENTROIDS
  centroids_summary <- do.call(rbind.data.frame, centroids_list_summary)  %>%
    unite("sts", season, trt, Sex)
  row.names(centroids_summary) <- NULL

  ##### OUTPUT: PARAMS_SUMMARY has the "a" and "b" parameters, calculated per season, per treatment, per sex
  ### params_summary has columns "sts" "a" and "b"

  #save to RDS file
  saveRDS(params_summary, file = here(params_file))
  #this is for SEASONAL CENTROIDS
  saveRDS(centroids_summary, file = here(centroids_file))

}

######------------------------------------------------------------------------

## Function for creating spatial overlap networks for individual voles for all the sites, months in a given year
## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
##          centroids_file = file name and extension in " " for the centroids_file (lifetime centroid)
##          params_file = file name and extension in "" for the params_file generated by generate_params function
##          networks_file = file name and extension in "" for the output file of adjacency matrices
## Output: list (saved to rds file) of all sites, months and an adjacency matrix for each of overlaps between voles

## LOAD A CENTROIDS FILE - LIKELY THE SEASONAL ONE - NOT THE LIFETIME ONE

create_overlap_networks <- function(data, centroids_file, params_file, networks_file){

  ##---------------- LOAD THE CAPTURE DATA ----------------------

  #load, clean  the fulltrap dataset
  fulltrap <- data %>%
    filter(month != "may") %>% #drop may data since not all sites had captures
    mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
    drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
    unite("sts", season, trt, sex, remove=FALSE) #add sts column to match params_summary

  #pull all the traps and their x, y coordinates, save as df with (trapID, x, y)
  traps <- fulltrap %>% group_by(trap) %>% slice(1) %>%
    select(trap, x, y) %>%
    arrange(trap)

  #save a vector of the site names (alphabetical order)
  site_names <- unique(fulltrap$site) %>% sort()


  ##---------------- CREATE A LIST of CAPTURE DATA (nested by site, month) --------------------

  #split() makes a list consisting of individual data.frames based on a condition ('site' in this case)
  site_list <- split(fulltrap, fulltrap$site)

  #create new list to hold nested site, month capture
  sitemonth_list <- list()

  for(i in 1:length(site_list)){
    # print(i)
    temp_list <- split(site_list[[i]], site_list[[i]]$month)
    sitemonth_list[[i]] <- temp_list
  }

  #name 1e list element as site (2e list elements are months May-Oct)
  names(sitemonth_list) <- site_names

  ### OUTPUT: SITEMONTH_LIST is a nested list of length 12
  #1e level is all the sites (12)
  #2e level is all the months per site (5) excluding May

  ##-------------- LOAD ADDITIONAL DATA ----------------------

  # Reading in required functions in wanelik_farine_functions.R
  source(here("wanelik_farine_functions.R"))

  # # OPTION: TO USE LIFETIME CENTROIDS
  # #read in the centroids RDS created in the centroid_calc() function

  # # OPTION: TO USE SEASONAL CENTROIDS
  # #read in the centroids RDS created in the generate_params() function
  obs_centroids_FULL <- readRDS(here(centroids_file))
  # #----END CENTROIDS----

  #read in the params RDS created in generate_params() function
  params_summary <- readRDS(file = here(params_file))

  ##--------------- (in a loop) GENERATE OVERLAP NETWORK ---------------------

  overlap_network_list <- list()

  for(i in 1:length(sitemonth_list)){

    print(names(sitemonth_list[i])) #print the site name

    overlap_network_list[[i]] <- list()

    for(j in 1:5){
      #1:5 for the five months (MAY DATA OMITTED)

      print(names(sitemonth_list[[i]][j])) #print the month

      data <- sitemonth_list[[i]][[j]] %>%
        select(tag, sts, trap, x, y) %>%
        group_by(tag, trap) %>%
        mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
        ungroup()

      #need these for either centroid option
      tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
      tagsts <- data %>% group_by(tag) %>% slice(1) %>% select(tag, sts) #all tags with their STS (season/trt/sex)

      # tagsex <- data %>% group_by(tag) %>% slice(1) %>% select(tag,sex) #all tags with their sex

      # skip creating network if there are 0 or 1 animals
      if(nrow(tags)=="0" | nrow(tags)=="1") {next}

      traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
      tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

      #the following creates VERY LARGE df of all tags for all traps
      #new column of Det.obs - whether that animal was observed in that trap
      #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
      fulldata <- cbind(tags_rep, traps_rep) %>%
        left_join(data, by=c("tag", "trap", "x", "y")) %>%
        group_by(tag) %>% fill(sts, .direction="downup") %>%
        mutate(Det.count = replace_na(Det.count, 0),
               Det.obs = ifelse(Det.count==0, 0, 1)) %>%
        arrange(tag, trap) %>%
        select(tag, sts, x, y, Det.obs, Det.count) %>%
        rename(x.trap = x,
               y.trap = y,
               Tag_ID = tag)

      #----------------------- 3. Generating overlap network  ------------------------------#

      # # OPTION FOR LIFETIME CENTROIDS
      # # #subset obs_centroids for only the animals on the grid that month (contains all animals in that season / year, depending on centroids file)
      # tag.vec <- as.vector(tags$tag)
      # obs_centroids <- obs_centroids_FULL %>% filter(Tag_ID %in% tag.vec) %>%
      #   left_join(tagsts, by=c("Tag_ID" = "tag")) #add the 'sts' column for later matchy-matchy

      # OPTION FOR SEASONAL CENTROIDS
      # subset obs_centroids_FULL to just those that month,
        #but make sure we're using the correct SEASON centroids for that vole (hence left join to tagsts)
      obs_centroids <- left_join(tagsts, obs_centroids_FULL, by=c("tag" = "Tag_ID", "sts"="sts")) %>%
        rename(Tag_ID = tag) #and change column name to match the rest of the code

      # Creating a matrix of distances between each individual's observed centroid and each trap
      matrix_dists_obs <- fulldata #only one entry per individual
      matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #x coord of centroid of individual
      matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #y coord of centroid of individual
      matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
      matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)

      # # Removing unobserved individuals
      # matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] #not needed, matrix is only of voles capped that month

      # Calculating total number of detections per individual
      ### TBH Janine isn't sure this is needed because it doesn't appear to be used anywhere, but a good check I guess
      a <- tapply(matrix_dists_obs$Det.count, list(matrix_dists_obs$Tag_ID), sum)
      matrix_dists_obs$Det.count.total <- a[match(matrix_dists_obs$Tag_ID, names(a))]

      # Estimating home range profiles (negative sigmoidal curves) done separately by season/trt -> params_summary file
      as <- params_summary$a[match(obs_centroids$sts,params_summary$sts)]
      bs <- params_summary$b[match(obs_centroids$sts,params_summary$sts)]

      overlap_network_sep <- get_network_2D(obs_centroids$x[which(!is.na(obs_centroids$x))], obs_centroids$y[which(!is.na(obs_centroids$y))],as,bs)
      rownames(overlap_network_sep)=colnames(overlap_network_sep)=na.omit(obs_centroids)$Tag_ID

      overlap_network_list[[i]][[j]] <- overlap_network_sep

    }


  }


  #name the 12 1st order elements of overlap_network_list as the sites
  names(overlap_network_list) <- names(sitemonth_list)

  #rename the sublist items (months) for each site
  for(i in 1:length(overlap_network_list)){
    names(overlap_network_list[[i]]) <- c("june", "july", "aug", "sept", "oct")
  }

  #in case there were NULL list elements (for a month with 0 or 1 animals)
  library(rlist)
  #a slick little function that removes NULL list elements (recursive=TRUE) to work through nested lists
  overlap_network_list <- list.clean(overlap_network_list, fun = is.null, recursive = TRUE)

  #save it since that code takes 5ever to run
  saveRDS(overlap_network_list, here(networks_file))


}

#####-------------------------------------------------------------------------------------

#### THIS CODE calculates the network metrics (mostly just weighted degree and its relatives)
# for all the voles at all the times
# output is "wt_net_mets_summary" file which can be used for downstream analysis

## Input: data = the FULL fulltrap df for that year
##        networks_file = networks file generated in create_overlap_networks
##        netmets_file = file to be generated, df of network metrics
## Output: netmets_file = dataframe of network metrics for every vole in every occasion it was captured


calculate_network_metrics <- function(data, networks_file, netmets_file){

  ##---------------- LOAD THE DATA ----------------------

  #load, clean the fulltrap dataset (make sure it's the most recent version)
  fulltrap <- data %>%
    filter(month != "may") %>% #drop may data since not all sites had captures
    mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
    drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
    unite("sts", season, trt, sex, remove=FALSE) #add sts column to match params_summary

  # load the network data
  overlap_network_list <- readRDS(here(networks_file))


  ##--------------------------------------------

  # Create list to store results
  wt_net_mets_list <- list()

  # Calculate network metrics
  for(i in 1:length(overlap_network_list)){

    #length(overlap_network_list) == for each of the 12 sites

    #for each site
    print(names(overlap_network_list[i]))
    site.id <- names(overlap_network_list[i])

    site <- list()

    for(j in 1:length(overlap_network_list[[i]])){

      #(j in 1:length(overlap_network_list[[i]])) == for the 4-5 trapping occasions
        #would be (j in 1:5) except that Kiirastuli 2022 only has 4 :P *thhhbbbbt*

      #for each month
      print(names(overlap_network_list[[i]][j]))
      month.id <- names(overlap_network_list[[i]][j])

      adjmat <- overlap_network_list[[i]][[j]]

      #create WEIGHTED NETWORK from adjacency matrix
      inet <- graph_from_adjacency_matrix(adjmat, mode="undirected", weighted = TRUE, diag = FALSE)

      ids <- get.vertex.attribute(inet, "name") #tag ids for all the animals on the grid
      month <- rep(names(overlap_network_list[[i]])[j],length(ids)) #capture month

      #dataframe to hold results per month
      site[[j]] <- data.frame(ids, month)

      #network metrics to calculate
      site[[j]]$wt.deg <- igraph::strength(inet) #this is the sum of all degree weights for a node
      # site[[j]]$norm.wt.deg <- (igraph::strength(inet))/((igraph::gorder(inet))-1) #your strength/(total nodes-you)
      site[[j]]$avg.wt.deg <- rep(mean(site[[j]]$wt.deg), length(ids)) #calculate average wt degree for a site/occasion

      ####### N COMPONENTS ISN'T WORKING RN
      # site[[j]]$n.component <- rep(igraph::count_components(inet), length(ids))

      site[[j]]$n.node <- rep(igraph::gorder(inet), length(ids))
      # site[[j]]$wt.n.edge <- rep(sum(E(inet)$weight), length(ids)) ##### THIS IS WEIGHTED #####


      ###### Calculate Male-degree, Female-degree #########

      #make a (new) network from the adj matrix
      g <- graph.adjacency(adjmat, mode="directed", weighted=TRUE, diag = FALSE)
          #technically, our adj matrix is symmetrical so the network is undirected,
          #BUT because I want count the number of ties going in/out from a given node,
          #So I need to have all the 'out' ties listed in one column of the edgelist (so effectively a directed edgelist)
          #So for these purposes, I built a directed network (but edge weights are symmetrical)

      # JK on second thought, don't do this because wt.deg isn't thresholded
      # #remove edges with weight less than threshold
      # g2 <- delete.edges(g, which(E(g)$weight<=0.05))

      #subset metadata
      metadata <- data %>% #starts with fulltrap, need to get down to a 'traits' version
        filter(site==site.id & month==month.id) %>%
        group_by(tag) %>% slice(1) %>% #one entry per vole per month
        drop_na(sex) #didn't build networks with animals with no sex

      #set "sex" as a vertex attribute
      g <- set.vertex.attribute(graph=g, name="sex", index=V(g), value=metadata$sex)
      #save vectors
      # focal_sex <- get.vertex.attribute(g, "sex") #female is 1, male is 2
      focal_id <- get.vertex.attribute(g, "name")

      ## code from Matt M-S for male strength/female strength
      sex_to <- get.vertex.attribute(g, "sex")[get.edgelist(g, names=FALSE)[,2]]
          #since vertices already have meaningful names, call names=FALSE to return vertex indices instead
      weight_to <- get.edge.attribute(g, "weight")[get.edgelist(g, names=FALSE)[,2]]
      degree_to_M <- strength(g, mode="out", weights=((sex_to == "M")*weight_to)) #this need to be mode="OUT" (the "to" individual)
      degree_to_F <- strength(g, mode="out", weights=((sex_to == "F")*weight_to)) #this need to be mode="OUT"


      #network metrics to calculate
      site[[j]]$focal_id <- focal_id
      # site[[j]]$focal_sex <- focal_sex
      site[[j]]$F.deg <- degree_to_F
      site[[j]]$M.deg <- degree_to_M

      # output <- as.data.frame(cbind(focal_id, focal_sex, degree_to_F, degree_to_M))

    }

    #write the list 'site' as a 1st-order item in wt_net_mets_list
    wt_net_mets_list[[i]] <- site

  }

  #name the 12 1st order elements of wt_nets_list as the sites
  names(wt_net_mets_list) <- names(overlap_network_list)

  #rename the sublist items (months) for each site
  #accounting for the fact that most sites have 5 months of data but some have 4
  for(i in 1:length(wt_net_mets_list)){
    ifelse( length(wt_net_mets_list[[i]]) == 5,
            names(wt_net_mets_list[[i]]) <- c("june", "july", "aug", "sept", "oct"),
            names(wt_net_mets_list[[i]]) <- c("july", "aug", "sept", "oct") )
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

  ## make net_mets_list_summary into freiggein huge df
  wt_net_mets_summary <- do.call(rbind.data.frame, wt_net_mets_list_summary)

  #clean up the df
  wt_net_mets_summary <- wt_net_mets_summary %>%
    rownames_to_column("name") %>% #row names are the sites, make that a column
    separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
    mutate(site = as.factor(site)) %>% #make site a factor
    rename(tag=ids)

  #save it
  saveRDS(wt_net_mets_summary, here(netmets_file))


}






















# ## Function for generating LIFETIME centroids
# ## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
# ##          centroids_file = file name and extension in " " for the centroids_file to be generated
# ## Output: dataframe of centroids (3 columns for tagid, sex, x, y)
#
# centroid_calc <- function(data, centroids_file) {
#
#   ##---------------- LOAD THE CAPTURE DATA ----------------------
#
#   #load, clean the fulltrap dataset
#   fulltrap <- data %>%
#     filter(month != "may") %>% #drop may data since not all sites had captures
#     mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
#     drop_na(sex) #remove animals with sex=NA (since we can't assign then a HR)
#
#   #pull all the traps and their x, y coordinates, save as df with (trapID, x, y)
#   traps <- fulltrap %>% group_by(trap) %>% slice(1) %>%
#     select(trap, x, y) %>%
#     arrange(trap)
#
#   ##---------------- CALCULATE CENTROIDS for all voles, across all captures --------------------
#
#   # Following code from Wanelik & Farine 2022
#
#   data <- fulltrap %>%
#     select(tag, sex, trap, x, y) %>%
#     group_by(tag, trap) %>%
#     mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
#     ungroup()
#
#   tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded
#   tagsex <- data %>% group_by(tag) %>% slice(1) %>% select(tag,sex) #all tags with their sex
#
#   traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
#   tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap
#
#   #VERY LARGE df of all tags for all traps
#   #new column of Det.obs - whether that animal was (ever) observed in that trap
#   #new column of Det.count - number of times animal was captured in that trap (LIFETIME)
#   fulldata <- cbind(tags_rep, traps_rep) %>%
#     left_join(data, by=c("tag", "trap", "x", "y")) %>%
#     group_by(tag) %>% fill(sex, .direction="downup") %>%
#     mutate(Det.count = replace_na(Det.count, 0),
#            Det.obs = ifelse(Det.count==0, 0, 1)) %>%
#     arrange(tag, trap) %>%
#     select(tag, sex, x, y, Det.obs, Det.count) %>%
#     rename(x.trap = x,
#            y.trap = y,
#            Tag_ID = tag,
#            Sex = sex)
#
#   # Calculating (weighted) centroids
#   # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
#   matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
#   matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
#   matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
#   matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
#   matrix_dists_real2<-matrix_dists_real2[,-7] #drops the Det.count2 column
#   x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
#   y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
#   obs_centroids <- data.frame(Tag_ID=tagsex$tag, Sex=tagsex$sex)
#   obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
#   obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]
#
#   #save it
#   saveRDS(obs_centroids, here(centroids_file))
#
# }



######------------------------------------------------------------------------
