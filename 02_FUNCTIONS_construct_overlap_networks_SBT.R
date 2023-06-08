######### THIS VERSION WAS USED FOR SEASON, SEX, BREEDING STATUS
## but basically the glms were REAL goofy for the nonbreeders and tbh who knows what they do anyway
## so lets leave them out, and just focus on the breeders and use the 'Season Trt Sex' version of this file
## but only with animals that were reproductively active
## 5/11/23 wishing I had a margarita in hand

## lol jk - the models were goofy because the juvs don't move / are only ever caught once...
## but I rewrote this code on 6/6 in a slightly different way because I'm an idiot and didn't realize I still had this


# data <- readRDS(here("fulltrap21_05.10.23.rds"))
# data <- readRDS(here("fulltrap22_05.10.23.rds"))
# params_file = "params21_stsb.rds"

## Function for generating the a and b parameters by season/trt/sex to define the distributions for vole HRs
## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
##          params_file = file name and extension in " " for the params_file to be generated
## Output: dataframe of parameters (3 columns for seasontrtsex, a param, b param)

generate_params <- function(data, params_file) {

  ##---------------- LOAD THE CAPTURE DATA ----------------------

  #load, clean the fulltrap dataset
  fulltrap <- data %>%
    filter(month != "may") %>% #drop may data since not all sites had captures
    mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
    drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
    drop_na(season_breeder) #remove animals with season_breeder=NA (since we can't assign then a HR)

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

  #save a vector of season_breeder
  season_breeder <- unique(fulltrap$season_breeder) %>% sort()

  ##---------------- CREATE A LIST of CAPTURE DATA (nested by season, trt) --------------------

  # #split() makes a list consisting of individual data.frames based on a condition ('season' in this case)
  # season_list <- split(fulltrap, fulltrap$season)
  #
  # #create new list to hold nested site, month capture
  # seasontrt_list <- list()
  #
  # for(i in 1:length(season_list)){
  #   # print(i)
  #   temp_list <- split(season_list[[i]], season_list[[i]]$trt)
  #   seasontrt_list[[i]] <- temp_list
  # }
  #
  # #name 1e list element as season (2e list elements are trts)
  # names(seasontrt_list) <- season
  #
  # ### SEASONTRT_LIST is a nested list of length 2
  # ### 1e level is season (summer/fall) (2)
  # ### 2e level is all the trts per site (4)
  # ### all grids per trt and months per season are combined together!!


  season_breeder_list <- lapply(split(fulltrap, fulltrap$season, drop = TRUE),
                                function(x) split(x, x[["season_breeder"]], drop = TRUE))

  #create new list to hold nested site, month capture
  sbt_list <- list()

  for(i in 1:length(season_breeder_list)){
    # for each season (summer, fall)
    season_list <- list()

    for(j in 1:length(season_breeder_list[[i]])){
      # for breeders and nonbreeders

      trt_list <- split(season_breeder_list[[i]][[j]], season_breeder_list[[i]][[j]]$trt)
      season_list[[j]] <- trt_list
    }

    names(season_list) <- season_breeder
    sbt_list[[i]] <- season_list

  }

  names(sbt_list) <- season

  ### SBT_LIST is a nested list of length 2
  ### 1e level is season (summer/fall) (2)
  ### 2e level is breeding status in that season (breeder / nonbreeder)
  ### 3e level is all the trts per site (4)
  ### all grids per trt and months per season are combined together!!

  ##------------------------------------

# Following code from Wanelik & Farine 2022

params_list <- list()

  for(i in 1:length(sbt_list)){

    print(i)

    params_list[[i]] <- list()

    for(k in 1:length(sbt_list[[i]])){

      print(k)

      params_list[[i]][[k]] <- list()

      for(j in 1:4){
        #1:4 for the four treatments - not reproducible, but easy

        print(j)

        data <- sbt_list[[i]][[k]][[j]] %>%
          select(tag, sex, trap, x, y) %>%
          group_by(tag, trap) %>%
          mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
          ungroup()

        season_lab <- unique(sbt_list[[i]][[k]][[j]]$season)
        breeder_lab <- unique(sbt_list[[i]][[k]][[j]]$season_breeder)
        trt_lab <- unique(sbt_list[[i]][[k]][[j]]$trt)

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

        ##SEASONAL CENTROIDS to calc SEASONAL PARAMS

        # Recalculating (weighted) centroids
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
        params <- data.frame(sex=c("M","F"),
                             a=c(coef(fit.males)[1],coef(fit.females)[1]),
                             b=c(coef(fit.males)[2],coef(fit.females)[2]))

        params <- params %>% mutate(season=season_lab, season_breeder=breeder_lab, trt=trt_lab)

        params_list[[i]][[k]][[j]] <- params

        }

    }

  }

  #collate results
  #this collapses the 2nd order elements (params per trt) down to the 1st order element (season)

  #make a list to store things
  params_list_summary1 <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(params_list)){

    #for both seasons
    summary <- do.call("rbind", params_list[[i]])
    params_list_summary1[[i]] <- summary
  }

  params_list_summary <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(params_list_summary1)){

    #for both seasons
    summary <- do.call("rbind", params_list_summary1[[i]])
    params_list_summary[[i]] <- summary
  }


  #collapse params_list_summary into a df
  params_summary <- do.call(rbind.data.frame, params_list_summary) %>%
    unite("stsb", season, trt, sex, season_breeder)
  row.names(params_summary) <- NULL

  ##### OUTPUT: PARAMS_SUMMARY has the "a" and "b" parameters, calculated per season, per treatment, per sex
  ### params_summary has columns "stsb" "a" and "b"

  #save to RDS file
  saveRDS(params_summary, file = here(params_file))

}

######------------------------------------------------------------------------


## Function for creating spatial overlap networks for individual voles for all the sites, months in a given year
## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
##          params_file = file name and extension in "" for the params_file generated by generate_params function
##          centroids_file = file name and extension in "" for the output file of monthly centroids for all voles
##          networks_file = file name and extension in "" for the output file of adjacency matrices
## Output: list (saved to rds file) of all sites, months and an adjacency matrix for each of overlaps between voles

# #clear environment
# rm(list = ls())
#
# ft21 <- readRDS(here("fulltrap21_05.10.23.rds"))
# data <- ft21
# params_file <- "params21_stsb.rds"
# networks_file <- "overlapnet21_stsb.rds"
# centroids_file <- "centroids21_monthly_bnb.rds"

create_overlap_networks <- function(data, params_file, centroids_file, networks_file){

  ##---------------- LOAD THE CAPTURE DATA ----------------------

  #load, clean  the fulltrap dataset
  fulltrap <- data %>%
    filter(month != "may") %>% #drop may data since not all sites had captures
    mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
    drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
    drop_na(season_breeder) %>% #remove animals with season_breeder
    unite("stsb", season, trt, sex, season_breeder, remove=FALSE) #add stsb column to match params_summary

  ##KEEPING BREEDERS AND NON-BREEDERS!

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

  #read in the params RDS created in generate_params() function
  params_summary <- readRDS(file = here(params_file))

  ##--------------- (in a loop) GENERATE OVERLAP NETWORK ---------------------

  overlap_network_list <- list()

  centroids_list <- list()

  for(i in 1:length(sitemonth_list)){

    print(names(sitemonth_list[i])) #print the site name

    overlap_network_list[[i]] <- list()

    centroids_list[[i]] <- list()

    for(j in 1:5){
      #1:5 for the five months (MAY DATA OMITTED)

      print(names(sitemonth_list[[i]][j])) #print the month

      data <- sitemonth_list[[i]][[j]] %>%
        select(tag, stsb, trap, x, y) %>%
        group_by(tag, trap) %>%
        mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
        ungroup()

      tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
      tagstsb <- data %>% group_by(tag) %>% slice(1) %>% select(tag, stsb) #all tags with their STSB (season/trt/sex/breeder)

      # skip creating network if there are 0 or 1 animals
      if(nrow(tags)=="0" | nrow(tags)=="1") {next}

      traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
      tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

      #creates VERY LARGE df of all tags for all traps
      #new column of Det.obs - whether that animal was observed in that trap
      #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
      fulldata <- cbind(tags_rep, traps_rep) %>%
        left_join(data, by=c("tag", "trap", "x", "y")) %>%
        group_by(tag) %>% fill(stsb, .direction="downup") %>%
        mutate(Det.count = replace_na(Det.count, 0),
               Det.obs = ifelse(Det.count==0, 0, 1)) %>%
        arrange(tag, trap) %>%
        select(tag, stsb, x, y, Det.obs, Det.count) %>%
        rename(x.trap = x,
               y.trap = y,
               Tag_ID = tag)

      #----------------------- 3. Generating overlap network  ------------------------------#

      ### MONTHLY CENTRIODS with SEASONAL A B PARAMS

      # Recalculating (weighted) centroids
      # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
      matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
      matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
      matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
      matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
      matrix_dists_real2<-matrix_dists_real2[,-7] #drops the Det.count2 column
      x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
      y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
      obs_centroids <- data.frame(Tag_ID=tagstsb$tag, stsb=tagstsb$stsb)
      obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
      obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]

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
      # a <- tapply(matrix_dists_obs$Det.count, list(matrix_dists_obs$Tag_ID), sum)
      # matrix_dists_obs$Det.count.total <- a[match(matrix_dists_obs$Tag_ID, names(a))]

      # Estimating home range profiles (negative sigmoidal curves) done separately by season/trt -> params_summary file

      # Generating overlap network
      as <- params_summary$a[match(obs_centroids$stsb,params_summary$stsb)]
      bs <- params_summary$b[match(obs_centroids$stsb,params_summary$stsb)]

      overlap_network_sep <- get_network_2D(obs_centroids$x[which(!is.na(obs_centroids$x))], obs_centroids$y[which(!is.na(obs_centroids$y))],as,bs)
      rownames(overlap_network_sep)=colnames(overlap_network_sep)=na.omit(obs_centroids)$Tag_ID

      overlap_network_list[[i]][[j]] <- overlap_network_sep

      #save monthly centroids (for plotting)
      centroids_list[[i]][[j]] <- obs_centroids %>% select(Tag_ID, x, y) %>% mutate(month=paste(names(sitemonth_list[[i]][j])),
                                                                                    site=paste(names(sitemonth_list[i])))

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



  #collate centroids results

  #make a list to store things
  centroids_summary1 <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(centroids_list)){

    #for all 12 sites
    summary <- do.call("rbind", centroids_list[[i]])
    centroids_summary1[[i]] <- summary
  }

  #name the 12 1st order elements as their sites
  names(centroids_summary1) <- names(sitemonth_list)

  ## make net_mets_list_summary into freiggein huge df
  centroids_summary <- do.call(rbind.data.frame, centroids_summary1)
  row.names(centroids_summary) <- NULL

  #saveRDS
  saveRDS(centroids_summary, here(centroids_file))

}

#####-------------------------------------------------------------------------------------

#### THIS CODE calculates the network metrics (mostly just weighted degree and its relatives)
# for all the voles at all the times
# output is "wt_net_mets_summary" file which can be used for downstream analysis

## Input: data = the FULL fulltrap df for that year
##        networks_file = networks file generated in create_overlap_networks
##        netmets_file = file to be generated, df of network metrics
## Output: netmets_file = dataframe of network metrics for every vole in every occasion it was captured

ft21 <- readRDS(here("fulltrap21_05.10.23.rds"))
data <- ft21
networks_file <- "overlapnet21_stsb.rds"
netmets_file <- "netmets21_stsb.rds"

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

  #create tag_sex df for assortativity by sex
  tag_sex <- fulltrap %>% group_by(tag) %>% slice(1) %>%
    select(tag, sex) %>%
    arrange(tag)


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
      # site[[j]]$avg.wt.deg <- rep(mean(site[[j]]$wt.deg), length(ids)) #calculate average wt degree for a site/occasion

      site[[j]]$n.node <- rep(igraph::gorder(inet), length(ids))
      # site[[j]]$wt.n.edge <- rep(sum(E(inet)$weight), length(ids)) ##### THIS IS WEIGHTED #####

      ####### N COMPONENTS ISN'T WORKING RN
      # site[[j]]$n.component <- rep(igraph::count_components(inet), length(ids))

      ### FOR ASSORTATIVITY - igraph doesn't do it with weighted degree - using assortnet in separate script


      ###### Calculate Male-degree, Female-degree #########

      ###### Calculate Breeder-degree, Nonbreeder-degree #########

      #make a (new) network from the adj matrix
      g <- graph_from_adjacency_matrix(adjmat, mode="directed", weighted=TRUE, diag = FALSE)
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
        drop_na(sex) %>% #didn't build networks with animals with no sex
        drop_na(season_breeder) #didn't build networks with animals with no breeder status

      #set "sex" as a vertex attribute
      g <- set.vertex.attribute(graph=g, name="sex", index=V(g), value=metadata$sex)
      #set "breeder" as a vertex attribute
      g <- set.vertex.attribute(graph=g, name="breeder", index=V(g), value=metadata$season_breeder)
      #save vectors
      # focal_sex <- get.vertex.attribute(g, "sex") #female is 1, male is 2
      # focal_breed <- get.vertex.attribute(g, "breeder") #breeder is 1, nonbreeder is 2
      focal_id <- get.vertex.attribute(g, "name")

      ## code from Matt M-S for male strength/female strength
      sex_to <- get.vertex.attribute(g, "sex")[get.edgelist(g, names=FALSE)[,2]]
          #since vertices already have meaningful names, call names=FALSE to return vertex indices instead
      weight_to <- get.edge.attribute(g, "weight")[get.edgelist(g, names=FALSE)[,2]]
      degree_to_M <- strength(g, mode="out", weights=((sex_to == "M")*weight_to)) #this need to be mode="OUT" (the "to" individual)
      degree_to_F <- strength(g, mode="out", weights=((sex_to == "F")*weight_to)) #this need to be mode="OUT"

      ## code from Matt M-S for breeder/nonbreeder strength
      breed_to <- get.vertex.attribute(g, "breeder")[get.edgelist(g, names=FALSE)[,2]]
      #since vertices already have meaningful names, call names=FALSE to return vertex indices instead
      weight_to <- get.edge.attribute(g, "weight")[get.edgelist(g, names=FALSE)[,2]]
      degree_to_b <- strength(g, mode="out", weights=((breed_to == "breeder")*weight_to)) #this need to be mode="OUT" (the "to" individual)
      degree_to_nb <- strength(g, mode="out", weights=((breed_to == "nonbreeder")*weight_to)) #this need to be mode="OUT"

      ##### TRYING SOMETHING 6/8 -- calculating weighted degree in the same way I calc M.deg/F.deg
      ## this at least does give me a degree measurement that = M.deg+F.deg
      ## but why isn't strength(inet) giving me the same result?
      nodestrength <- strength(g, mode="out", weights=(weight_to)) #this need to be mode="OUT" (the "to" individual)
      site[[j]]$strength <- nodestrength
      ##############################################

      #####-------------------------------------

      #network metrics to calculate
      site[[j]]$focal_id <- focal_id
      # site[[j]]$focal_sex <- focal_sex
      site[[j]]$F.deg <- degree_to_F
      site[[j]]$M.deg <- degree_to_M
      site[[j]]$b.deg <- degree_to_b
      site[[j]]$nb.deg <- degree_to_nb

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

  # #save it
  # saveRDS(wt_net_mets_summary, here(netmets_file))


}


