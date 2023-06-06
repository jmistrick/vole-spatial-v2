### GENERATE PARAMS for M / F breeders, non breeders in summer and fall


## Function for generating the a and b parameters by season/trt/sex to define the distributions for vole HRs
## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
##          centroids_file = file name and extension in " " for the centroids_file (SEASONAL centroid)
##          params_file = file name and extension in " " for the params_file to be generated
## Output: dataframe of parameters (3 columns for seasontrtsex, a param, b param)


# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)

#clear environment
rm(list = ls())



#load the fulltrap 21 and 22 data (make sure it's the most recent version)
## NOTE ## all DP, DT, or S animals are still in the fulltrap datasets
# ft21 <- readRDS(here("fulltrap21_05.10.23.rds")) %>% #includes breeders and nonbreeders
#   filter(sex=="M") #one sex at a time
#   # %>% filter(caps_per_season > 1)

ft22 <- readRDS(here("fulltrap22_05.10.23.rds")) %>% #includes breeders and nonbreeders
  filter(sex=="F") #one sex at a time
  # %>% filter(caps_per_season > 1)


# #visualize n captures per season
# ft21 %>% group_by(season, tag) %>% drop_na(season_breeder) %>% slice(1) %>%
#   # filter(trt=="fed_deworm") %>%
#   ggplot( aes(x=caps_per_season, fill=season_breeder)) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   scale_fill_manual(values=c("#69b3a2", "#FFBF00")) +
#   facet_wrap(~season) +
#   labs(fill="")
#
# #visualize n traps per season
# ft21 %>% group_by(season, tag) %>% drop_na(season_breeder) %>% slice(1) %>%
#   # filter(trt=="fed_deworm") %>%
#   ggplot( aes(x=traps_per_season, fill=season_breeder)) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   scale_fill_manual(values=c("#69b3a2", "#FFBF00")) +
#   facet_wrap(~season) +
#   labs(fill="")

# #n breeders/nonbreeders per treatment in season
# ft21 %>% group_by(season, tag) %>% drop_na(season_breeder) %>% slice(1) %>%
#   ungroup() %>% group_by(season, trt, season_breeder) %>%
#   summarise(n=length(tag))





#####------------------------------------------------------

# data = ft21
# centroids_file = "M_b-nb_centroids21.rds"
# params_file = "M_b-nb_params21.rds"

data = ft22
centroids_file = "F_b-nb_centroids22.rds"
params_file = "F_b-nb_params22.rds"


####---------------------------------


## CENTROIDS here refers to SEASONAL CENTROIDS
# generate_params <- function(data, centroids_file, params_file) {

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

      ### CHANGE FOR THIS VERSION: USING 'SEASON_BREEDER' instead of SEX to build the distributions

      data <- seasontrt_list[[i]][[j]] %>%
        select(tag, season_breeder, trap, x, y) %>%
        group_by(tag, trap) %>%
        mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
        ungroup()

      season_lab <- unique(seasontrt_list[[i]][[j]]$season)
      trt_lab <- unique(seasontrt_list[[i]][[j]]$trt)
      sex_lab <- unique(seasontrt_list[[i]][[j]]$sex)
      year_lab <- unique(seasontrt_list[[i]][[j]]$year)

      tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
      tagbreeder <- data %>% group_by(tag) %>% slice(1) %>% select(tag,season_breeder) #all tags with their breeding status

      traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
      tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

      #VERY LARGE df of all tags for all traps
      #new column of Det.obs - whether that animal was observed in that trap
      #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
      fulldata <- cbind(tags_rep, traps_rep) %>%
        left_join(data, by=c("tag", "trap", "x", "y")) %>%
        group_by(tag) %>% fill(season_breeder, .direction="downup") %>%
        mutate(Det.count = replace_na(Det.count, 0),
               Det.obs = ifelse(Det.count==0, 0, 1)) %>%
        arrange(tag, trap) %>%
        select(tag, season_breeder, x, y, Det.obs, Det.count) %>%
        rename(x.trap = x,
               y.trap = y,
               Tag_ID = tag,
               Breeder = season_breeder)

      #----------------------- 3. Generating overlap network  ------------------------------#

      # Calculating (weighted) centroids (FOR each VOLE, in a SEASON)
      # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
      matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
      matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
      matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
      matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
      matrix_dists_real2<-matrix_dists_real2[,-7] #drops the Det.count2 column
      x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
      y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
      obs_centroids <- data.frame(Tag_ID=tagbreeder$tag, Breeder=tagbreeder$season_breeder)
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
      fit.breeder <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Breeder=="breeder"),], family=binomial, control = list(maxit = 50))
      fit.nonbreeder <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Breeder=="nonbreeder"),], family=binomial, control = list(maxit = 50))

      # Generating overlap network
      params <- data.frame(breeder=c("breeder","nonbreeder"), a=c(coef(fit.breeder)[1],coef(fit.nonbreeder)[1]),
                           b=c(coef(fit.breeder)[2],coef(fit.nonbreeder)[2]))

      params <- params %>% mutate(season=season_lab,
                                  trt=trt_lab,
                                  sex=sex_lab,
                                  year=year_lab)

      #this is for SEASONAL CENTROIDS
      obs_centroids <- obs_centroids %>% mutate(season=season_lab,
                                                trt=trt_lab,
                                                sex=sex_lab,
                                                year=year_lab)

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
  params_summary <- do.call(rbind.data.frame, params_list_summary)
    # %>% unite("sts", season, trt, sex)
  # row.names(params_summary) <- NULL

  #this is for SEASONAL CENTROIDS
  centroids_summary <- do.call(rbind.data.frame, centroids_list_summary)
    # %>% unite("sts", season, trt, Sex)
  # row.names(centroids_summary) <- NULL

  ##### OUTPUT: PARAMS_SUMMARY has the "a" and "b" parameters, calculated per season, per treatment, per sex
  ### params_summary has columns "sts" "a" and "b"



####--------- SAVE OUTPUT ----------------------

#save to RDS file
saveRDS(params_summary, file = here(params_file))

#this is for SEASONAL CENTROIDS
saveRDS(centroids_summary, file = here(centroids_file))


#########################################################################################

Mcentroids21 <- readRDS(here("M_b-nb_centriods21.rds"))
Fcentroids21 <- readRDS(here("F_b-nb_centriods21.rds"))
Mcentroids22 <- readRDS(here("M_b-nb_centriods22.rds"))
Fcentroids22 <- readRDS(here("F_b-nb_centriods22.rds"))

bnb_centroids21.22 <- rbind(Mcentroids21, Mcentroids22, Fcentroids21, Fcentroids22)
saveRDS(bnb_centroids21.22, here("bnb_centroids21.22.rds"))

Mparams21  <- readRDS(here("M_b-nb_params21.rds"))
Mparams22  <- readRDS(here("M_b-nb_params22.rds"))
Fparams21  <- readRDS(here("F_b-nb_params21.rds"))
Fparams22  <- readRDS(here("F_b-nb_params22.rds"))

bnb_params21.22 <- rbind(Mparams21, Mparams22, Fparams21, Fparams22)
saveRDS(bnb_params21.22, here("bnb_params21.22.rds"))



