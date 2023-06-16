### THIS CODE tests the distance from centroid (HR size effectively?) for males/females and across treatments
  # SEPARATELY in summer and fall
  # in both seasons, fed treatments are different from unfed_control -- unfed_deworm (only diff in fall)
  # in SUMMER - males and females differ in HR size, NOT in fall (YAY! I FOUND THIS TOO! :D)

##sooooo basically I think it therefore makes sense to calculate the HR params separately per season/trt like I did


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
fulltrap <- readRDS(file = "fulltrap21_05.10.23.rds") %>%
  filter(month != "may") %>% #drop may data since not all sites had captures
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
  drop_na(sex) #remove animals with sex=NA (since we can't assign then a HR)


#pull all the traps and their x, y coordinates, save as df with (trapID, x, y)
traps <- fulltrap %>% group_by(trap) %>% slice(1) %>%
  select(trap, x, y) %>%
  arrange(trap)

# #save a vector of the site names (alphabetical order)
# site_names <- unique(fulltrap$site) %>% sort()
#
# #save a vector of trts
# trt <- unique(fulltrap$trt) %>% sort()
#
# #save a vector of seasons
# season <- unique(fulltrap$season) %>% sort()


##------------------------------------

# Reading in required functions in wanelik_farine_functions.R
source(here("wanelik_farine_functions.R"))

##------------ CALCULATE DISTANCES FOR EVERYONE, RUN GLM BY SEX, TRT ------------------------
##------------------------------------- SUMMER ----------------------------------------------

    data <- fulltrap %>%
      filter(season=="summer") %>%
      select(tag, trt, sex, trap, x, y) %>%
      group_by(tag, trap) %>%
      mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
      ungroup()

    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded
    tagdata <- data %>% group_by(tag) %>% slice(1) %>% select(tag, trt, sex) #all tags with their data

    traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
    tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

    #VERY LARGE df of all tags for all traps
    #new column of Det.obs - whether that animal was observed in that trap
    #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
    fulldata <- cbind(tags_rep, traps_rep) %>%
      left_join(data, by=c("tag", "trap", "x", "y")) %>%
      group_by(tag) %>%
      fill(c(sex, trt), .direction="downup") %>%
      ungroup() %>%
      mutate(Det.count = replace_na(Det.count, 0),
             Det.obs = ifelse(Det.count==0, 0, 1)) %>%
      arrange(tag, trap) %>%
      select(tag, sex, trt, x, y, Det.obs, Det.count) %>%
      rename(x.trap = x,
             y.trap = y,
             Tag_ID = tag)

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
    obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, trt=tagdata$trt)
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

    # Estimating home range profiles (negative sigmoidal curves) using this data
    fit.summer <- glm(Det.obs ~ Dist.log + sex + trt + Dist.log*sex + Dist.log*trt, data=matrix_dists_obs, family=binomial, control = list(maxit = 50))
    summary(fit.summer)

  ################# SHOULD THERE BE AN INTERACTION of SEX/TRT/ and DIST?

### SUMMER: YES - SIGNIFICANT INTERACTION Dist.log:sexM (compared to sexF)
      ### and Dist.log:trtfed_control and Dist.log:trtfed_deworm (compared to trtunfed_control)

##-----------------------------------------------------------------------------------------------------

    ##------------ CALCULATE DISTANCES FOR EVERYONE, RUN GLM BY SEX, TRT ------------------------
    ##------------------------------------- FALL ----------------------------------------------

    data <- fulltrap %>%
      filter(season=="fall") %>%
      select(tag, trt, sex, trap, x, y) %>%
      group_by(tag, trap) %>%
      mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
      ungroup()

    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded
    tagdata <- data %>% group_by(tag) %>% slice(1) %>% select(tag, trt, sex) #all tags with their data

    traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
    tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

    #VERY LARGE df of all tags for all traps
    #new column of Det.obs - whether that animal was observed in that trap
    #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
    fulldata <- cbind(tags_rep, traps_rep) %>%
      left_join(data, by=c("tag", "trap", "x", "y")) %>%
      group_by(tag) %>%
      fill(c(sex, trt), .direction="downup") %>%
      ungroup() %>%
      mutate(Det.count = replace_na(Det.count, 0),
             Det.obs = ifelse(Det.count==0, 0, 1)) %>%
      arrange(tag, trap) %>%
      select(tag, sex, trt, x, y, Det.obs, Det.count) %>%
      rename(x.trap = x,
             y.trap = y,
             Tag_ID = tag)

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
    obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, trt=tagdata$trt)
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

    # Estimating home range profiles (negative sigmoidal curves) using this data
    fit.fall <- glm(Det.obs ~ Dist.log + sex + trt + Dist.log*sex + Dist.log*trt, data=matrix_dists_obs, family=binomial, control = list(maxit = 50))
    summary(fit.fall)

    par(mfrow = c(2, 2))
    plot(fit.fall)

### FALL: ISH - NO INTERACTION Dist.log:sexM (compared to sexF) [matches what I found in OG analysis]
    ### SIGNIFICANT INTERACTION: (compared to trtunfed_control)
          ### Dist.log:unfed_deworm, Dist.log:trtfed_control, and Dist.log:trtfed_deworm
