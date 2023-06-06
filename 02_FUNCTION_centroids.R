
## Function for generating LIFETIME centroids
## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
##          centroids_file = file name and extension in " " for the centroids_file to be generated
## Output: dataframe of centroids (3 columns for tagid, sex, x, y)

centroid_calc <- function(data, centroids_file) {

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

  ##---------------- CALCULATE CENTROIDS for all voles, across all captures --------------------

# Following code from Wanelik & Farine 2022

    data <- fulltrap %>%
        select(tag, sex, trap, x, y) %>%
        group_by(tag, trap) %>%
        mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
        ungroup()

      tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded
      tagsex <- data %>% group_by(tag) %>% slice(1) %>% select(tag,sex) #all tags with their sex

      traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
      tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

      #VERY LARGE df of all tags for all traps
      #new column of Det.obs - whether that animal was (ever) observed in that trap
      #new column of Det.count - number of times animal was captured in that trap (LIFETIME)
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

      # Calculating (weighted) centroids
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

      #save it
      saveRDS(obs_centroids, here(centroids_file))

}
