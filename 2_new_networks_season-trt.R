# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)

#clear environment
rm(list = ls())

##---------------- LOAD THE CAPTURE DATA ----------------------

#load the fulltrap dataset (make sure it's the most recent version)
  ## NOTE ## all DP, DT, or S animals are still in the fulltrap dataset
fulltrap <- readRDS(file = "fulltrap_06.28.22.rds") %>%
  drop_na(sex) #remove animals with sex=NA (since we can't assign then a HR)

########## LATER: will want to drop all the May data:
# filter(month != "may") %>% #drop may data since not all sites had captures
########## but that will require going back to 01_data_cleaning and changing the factor levels to remove "may"

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

# Reading in required functions in wanelik_farine_functions.R
source(here("wanelik_farine_functions.R"))

##------------------------------------

params_list <- list()

for(i in 1:length(seasontrt_list)){

  print(i)

  params_list[[i]] <- list()

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
    params <- data.frame(sex=c("M","F"), a=c(coef(fit.males)[1],coef(fit.females)[1]), b=c(coef(fit.males)[2],coef(fit.females)[2]))

    params <- params %>% mutate(season=season_lab,
                                trt=trt_lab)

    params_list[[i]][[j]] <- params
  }

}

#collate results
#this collapses the 2nd order elements (params per trt) down to the 1st order element (season)

#make a list to store things
params_list_summary <- list()

#loop across all sites and collapse the dfs per occasion into one df for the site
for(i in 1:length(params_list)){

  #for both seasons
  summary <- do.call("rbind", params_list[[i]])
  params_list_summary[[i]] <- summary
}

#collapse params_list_summary into a df
params_summary <- do.call(rbind.data.frame, params_list_summary) %>%
  unite("sts", season, trt, sex)
row.names(params_summary) <- NULL



##### OUTPUT: PARAMS_SUMMARY has the "a" and "b" parameters, calculated per season, per treatment, per sex
  ### params_summary has columns "sts" "a" and "b"

#save to RDS file
saveRDS(params_summary, file = here("params_summary_04.04.23.rds"))


