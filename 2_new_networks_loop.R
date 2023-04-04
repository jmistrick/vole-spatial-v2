# load packages
library(here)
library(tidyverse)
# library(mosaic) #mosaic does stats stuff, useful but not necessary yet 11.29.21
library(CMRnet)
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

### SITEMONTH_LIST is a nested list of length 12
  ### 1e level is all the sites (12)
  ### 2e level is all the months per site (5) excluding May

##------------------------------------

# Reading in required functions in wanelik_farine_functions.R
source(here("wanelik_farine_functions.R"))

##------------------------------------

data <- fulltrap %>% filter(month=="sept" & site=="janakkala") %>%
  select(tag, sex, trap, x, y) %>%
  drop_na(sex) %>% #all animals need a sex since HR estimate is for males/females
  group_by(tag, trap) %>%
  mutate(Det.count = length(tag)) %>% slice(1) %>%
  ungroup()

tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag)
tagsex <- data %>% group_by(tag) %>% slice(1) %>% select(tag,sex)

traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE))
tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE))

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
matrix_dists_real2<-fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
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
matrix_dists_obs <- matrix_dists_real2
matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)]
matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)]
matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)
# Removing unobserved individuals
matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] #not needed, matrix is only of voles capped that month

# Calculating total number of detections per individual
a <- tapply(matrix_dists_obs$Det.count, list(matrix_dists_obs$Tag_ID),sum)
matrix_dists_obs$Det.count.total <- a[match(matrix_dists_obs$Tag_ID,names(a))]

# Estimating home range profiles (negative sigmoidal curves) using this data 
# Example: Sex-specific curves (male home range profile = fit.males; female home range profile = fit.females)
fit.males <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Sex=="M"),], family=binomial, control = list(maxit = 50))
fit.females <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Sex=="F"),], family=binomial, control = list(maxit = 50))

# Generating overlap network 
params <- data.frame(sex=c("M","F"), a=c(coef(fit.males)[1],coef(fit.females)[1]), b=c(coef(fit.males)[2],coef(fit.females)[2]))
as <- params$a[match(obs_centroids$Sex,params$sex)]
bs <- params$b[match(obs_centroids$Sex,params$sex)]

overlap_network_sep <- get_network_2D(obs_centroids$x[which(!is.na(obs_centroids$x))], obs_centroids$y[which(!is.na(obs_centroids$y))],as,bs)
rownames(overlap_network_sep)=colnames(overlap_network_sep)=na.omit(obs_centroids)$Tag_ID  


##--------------------------------------------------------------

#plot adj matrix for network
g <- graph_from_adjacency_matrix(
  overlap_network_sep,
  mode = c("undirected"),
  weighted = TRUE,
  diag = FALSE)

g2 <- delete.edges(g, which(E(g)$weight<0.05))

plot(g2, vertex.size=5, vertex.label=NA)

