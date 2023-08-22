### THIS IS THE FINAL VERSION MADE IN THE WEE HOURS OF JUNE 27th when all I wanted to do was go to bed but instead
    ## here I am eating blueberry bagels with peanut butter and coding at 2am

### THIS CODE tests the distance from centroid (HR size effectively?) for males/females, breeder/nonbreeder and across treatments
  # SEPARATELY in summer and fall
  # SEPARATELY for 2021 and 2022 (just change the fulltrap file that's being loaded in line 27 ish)

  #yes, differences by sex, breeder, trt in all the seasons

##sooooo basically I think it therefore makes sense to calculate the HR params separately per season/sex/breeder/trt like I did


###### jUST A NOTE 7/11/23 (fucking hell I'm running out of time) -- decided that fitting space use like this (sex, breeder, trt all in a season)
## is way better than whatevertf I was was doing before in pieces. SO I'm going to use this version for BOTH
## vole-spatial and vole-hanta -- I also added an interaction between food*helm because science
## THIS separately looks at sex and breeding (not 'functional group' as a 4 factor variable). I don't know if
## it matters, but this is easier to interpret so deal with it #aggressive

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
  drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
  drop_na(season_breeder)

# ## is small space use of non-reproductive voles in Summer 2021 real or an artifact of capture freq?
# fulltrap %>% group_by(season, season_breeder) %>%
#   summarise(mean = mean(traps_per_season),
#             sd = sd(traps_per_season),
#             min = min(traps_per_season),
#             max = max(traps_per_season))
#
# fulltrap %>% group_by(season, season_breeder) %>%
#   summarise(mean = mean(caps_per_season),
#             sd = sd(caps_per_season),
#             min = min(caps_per_season),
#             max = max(caps_per_season))


# fulltrap <- readRDS(file = "fulltrap22_05.10.23.rds") %>%
#   filter(month != "may") %>% #drop may data since not all sites had captures
#   mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
#   drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
#   drop_na(season_breeder)


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

##------------ CALCULATE DISTANCES FOR EVERYONE, RUN GLM BY SEX, BREEDER, TRT ------------------------
##------------------------------------- SUMMER ----------------------------------------------

    data <- fulltrap %>%
      filter(season=="summer") %>%
      # select(tag, trt, sex, season_breeder, trap, x, y) %>%
      select(tag, food_trt, helm_trt, sex, season_breeder, trap, x, y) %>%
      mutate(season_breeder = factor(season_breeder, levels=c("nonbreeder", "breeder"))) %>%
      group_by(tag, trap) %>%
      mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
      ungroup()

    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded
    tagdata <- data %>% group_by(tag) %>% slice(1) %>%
      # select(tag, trt, sex, season_breeder) #all tags with their data
      select(tag, food_trt, helm_trt, sex, season_breeder) #all tags with their data

    traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
    tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) %>% arrange(tag) #all tags replicated for each trap

    #VERY LARGE df of all tags for all traps
    #new column of Det.obs - whether that animal was observed in that trap
    #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
    fulldata <- cbind(tags_rep, traps_rep) %>%
      left_join(data, by=c("tag", "trap", "x", "y")) %>%
      group_by(tag) %>%
      # fill(c(sex, season_breeder, trt), .direction="downup") %>%
      fill(c(sex, season_breeder, food_trt, helm_trt), .direction="downup") %>%
      ungroup() %>%
      mutate(Det.count = replace_na(Det.count, 0),
             Det.obs = ifelse(Det.count==0, 0, 1)) %>%
      arrange(tag, trap) %>%
      # select(tag, sex, season_breeder, trt, x, y, Det.obs, Det.count) %>%
      select(tag, sex, season_breeder, food_trt, helm_trt, x, y, Det.obs, Det.count) %>%
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
    matrix_dists_real2<-matrix_dists_real2[,-9] #drops the Det.count2 column
    x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
    y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
    # obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, season_breeder=tagdata$season_breeder, trt=tagdata$trt)
    obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, season_breeder=tagdata$season_breeder, food_trt=tagdata$food_trt, helm_trt=tagdata$helm_trt)
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
    # fit.summer <- glm(Det.obs ~ Dist.log + sex + trt + season_breeder +
    #                     Dist.log*sex + Dist.log*trt + Dist.log*season_breeder, data=matrix_dists_obs, family=binomial, control = list(maxit = 50))
    fit.summer <- glm(Det.obs ~ Dist.log + sex + season_breeder +
                        food_trt + helm_trt + food_trt*helm_trt +
                        Dist.log*sex + Dist.log*season_breeder +
                        Dist.log*food_trt + Dist.log*helm_trt +
                        Dist.log*food_trt*helm_trt,
                      data=matrix_dists_obs, family=binomial, control = list(maxit = 50))

    summary(fit.summer)

  ################# SHOULD THERE BE AN INTERACTION of SEX/TRT/BREED and DIST?

    #yes


    library(visreg)
    library(cowplot)
    library(ggtext)


  sex <- visreg(fit.summer, "Dist.log", by="sex", scale='response', rug=FALSE,
                  gg=TRUE, overlay=TRUE) +
      scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
      scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
      scale_color_manual(values=c("#f282a7", "#00d0ff")) +
      scale_fill_manual(values=c("#f282a750", "#00d0ff50")) +
      theme_half_open() +
      theme(legend.position = "bottom",
            axis.title = element_text(size=18),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text = element_text(size=16),
            plot.margin = margin(t = 10, r = 20, b = 20, l = 10, unit = "pt")) +
      labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
      annotate(geom = "text", x=2.4, y=.85, size = 6,
               label = paste("p < 0.001" )) +
      annotate(geom = "text", x=2.4, y=.9, size = 6,
               label = paste("OR =",
                             round( exp(coef(summary(fit.summer))[8,1]), digits=3) )) +
      geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=sex),
                 size=3, alpha=0.1, shape=16)
    # geom_jitter(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=sex),
    #             width=0, height=0.025,
    #            size=3, alpha=0.2, shape=16) #jittered data points, kind of ugly

    #by repro
    repro <- visreg(fit.summer, "Dist.log", by="season_breeder", scale="response", rug=FALSE,
                    gg=TRUE, overlay=TRUE) +
      scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
      scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
      scale_color_manual(values=c("#8b64b9", "#e8ac65")) +
      scale_fill_manual(values=c("#8b64b950", "#e8ac6550")) +
      theme_half_open() +
      theme(legend.position = "bottom",
            axis.title = element_text(size=18),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text = element_text(size=16),
            plot.margin = margin(t = 10, r = 10, b = 20, l = 20, unit = "pt")) +
      labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
      annotate(geom = "text", x=2.4, y=.85, size = 6,
               label = paste("p < 0.001")) +
      annotate(geom = "text", x=2.4, y=.9, size = 6,
               label = paste("OR =",
                             round( exp(coef(summary(fit.summer))[9,1]), digits=3) )) +
      geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=season_breeder),
                 size=3, alpha=0.15, shape=16)
      # geom_jitter(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=season_breeder),
      #             width=0, height=0.025,
      #             size=3, alpha=0.2, shape=16)

    #by food
    food <- visreg(fit.summer, "Dist.log", by="food_trt", scale='response', rug=FALSE,
                   gg=TRUE, overlay=TRUE) +
      scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
      scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
      scale_color_manual(values=c("#794624", "#68b63e")) +
      scale_fill_manual(values=c("#79462450", "#68b63e50")) +
      theme_half_open() +
      theme(legend.position = "bottom",
            axis.title = element_text(size=18),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text = element_text(size=16),
            plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt")) +
      labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
      annotate(geom = "text", x=2.4, y=.85, size = 6,
               label = paste("p < 0.001")) +
      annotate(geom = "text", x=2.4, y=.9, size = 6,
               label = paste("OR =",
                             round( exp(coef(summary(fit.summer))[10,1]), digits=3) )) +
      geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=food_trt),
                 size=3, alpha=0.15, shape=16)
      # geom_jitter(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=food_trt),
      #             width=0, height=0.025,
      #             size=3, alpha=0.2, shape=16)

    #by worms
    worms <- visreg(fit.summer, "Dist.log", by="helm_trt", scale='response', rug=FALSE,
                    gg=TRUE, overlay=TRUE) +
      scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
      scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
      scale_color_manual(values=c("#808080", "#ffe048")) +
      scale_fill_manual(values=c("#80808070", "#ffe04850")) +
      theme_half_open() +
      theme(legend.position = "bottom",
            axis.title = element_text(size=18),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text = element_text(size=16),
            plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
      labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
      annotate(geom = "text", x=2.4, y=.85, size = 6,
               label = paste("p =",
                             round( coef(summary(fit.summer))[11,4], digits=3) )) +
      annotate(geom = "text", x=2.4, y=.9, size = 6,
               label = paste("OR =",
                             round( exp(coef(summary(fit.summer))[11,1]), digits=3) )) +
      geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=helm_trt),
                 size=3, alpha=0.1, shape=16)
      # geom_jitter(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=helm_trt),
      #             width=0, height=0.025,
      #             size=3, alpha=0.2, shape=16)


    png(filename="fitsummer21.png", height=12, width=16, units="in", res=600)
    plot_grid(sex, repro, food, worms,
              labels=NULL, nrow=2)
    dev.off()


##-----------------------------------------------------------------------------------------------------

    ##------------ CALCULATE DISTANCES FOR EVERYONE, RUN GLM BY SEX,BREEDER,TRT ------------------------
    ##------------------------------------- FALL ----------------------------------------------


    data <- fulltrap %>%
      filter(season=="fall") %>%
      # select(tag, trt, sex, season_breeder, trap, x, y) %>%
      select(tag, food_trt, helm_trt, sex, season_breeder, trap, x, y) %>%
      mutate(season_breeder = factor(season_breeder, levels=c("nonbreeder", "breeder"))) %>%
      group_by(tag, trap) %>%
      mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
      ungroup()

    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded
    tagdata <- data %>% group_by(tag) %>% slice(1) %>%
      # select(tag, trt, sex, season_breeder) #all tags with their data
      select(tag, food_trt, helm_trt, sex, season_breeder) #all tags with their data

    traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
    tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) %>% arrange(tag) #all tags replicated for each trap

    #VERY LARGE df of all tags for all traps
    #new column of Det.obs - whether that animal was observed in that trap
    #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
    fulldata <- cbind(tags_rep, traps_rep) %>%
      left_join(data, by=c("tag", "trap", "x", "y")) %>%
      group_by(tag) %>%
      # fill(c(sex, season_breeder, trt), .direction="downup") %>%
      fill(c(sex, season_breeder, food_trt, helm_trt), .direction="downup") %>%
      ungroup() %>%
      mutate(Det.count = replace_na(Det.count, 0),
             Det.obs = ifelse(Det.count==0, 0, 1)) %>%
      arrange(tag, trap) %>%
      # select(tag, sex, season_breeder, trt, x, y, Det.obs, Det.count) %>%
      select(tag, sex, season_breeder, food_trt, helm_trt, x, y, Det.obs, Det.count) %>%
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
    matrix_dists_real2<-matrix_dists_real2[,-9] #drops the Det.count2 column
    x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
    y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
    # obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, season_breeder=tagdata$season_breeder, trt=tagdata$trt)
    obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, season_breeder=tagdata$season_breeder, food_trt=tagdata$food_trt, helm_trt=tagdata$helm_trt)
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
    fit.fall <- glm(Det.obs ~ Dist.log + sex + season_breeder +
                      food_trt + helm_trt + food_trt*helm_trt +
                      Dist.log*sex + Dist.log*season_breeder +
                      Dist.log*food_trt + Dist.log*helm_trt +
                      Dist.log*food_trt*helm_trt,
                    data=matrix_dists_obs, family=binomial, control = list(maxit = 50))
    summary(fit.fall)

    ################# SHOULD THERE BE AN INTERACTION of SEX/TRT/BREED and DIST?


library(visreg) #https://pbreheny.github.io/visreg/articles/web/overlay.html
library(cowplot)
library(ggtext)
#by sex
# visreg(fit.fall, "Dist.log", by="sex", overlay=TRUE,
#        xlab="Log Distance from Center", ylab="Log Odds (Detection)")
sex <- visreg(fit.fall, "Dist.log", by="sex", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  scale_color_manual(values=c("#f282a7", "#00d0ff")) +
  scale_fill_manual(values=c("#f282a750", "#00d0ff50")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size=16),
        plot.margin = margin(t = 10, r = 20, b = 20, l = 10, unit = "pt")) +
  labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  annotate(geom = "text", x=2.4, y=.85, size = 6,
           label = paste("p =",
                         round( coef(summary(fit.fall))[8,4], digits=3) )) +
  annotate(geom = "text", x=2.4, y=.9, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(fit.fall))[8,1]), digits=3) )) +
  #https://github.com/pbreheny/visreg/issues/56 #color points by group separately
  geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=sex),
             size=3, alpha=0.15, shape=16)
  # geom_jitter(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=sex),
  #             width=0, height=0.025,
  #             size=3, alpha=0.2, shape=16)

#by repro
repro <- visreg(fit.fall, "Dist.log", by="season_breeder", scale="response", rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  scale_color_manual(values=c("#8b64b9", "#e8ac65")) +
  scale_fill_manual(values=c("#8b64b950", "#e8ac6550")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size=16),
        plot.margin = margin(t = 10, r = 10, b = 20, l = 20, unit = "pt")) +
  labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  annotate(geom = "text", x=2.4, y=.85, size = 6,
           label = paste("p =",
                         round( coef(summary(fit.fall))[9,4], digits=3) )) +
  annotate(geom = "text", x=2.4, y=.9, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(fit.fall))[9,1]), digits=3) )) +
  geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=season_breeder),
             size=3, alpha=0.15, shape=16)
  # geom_jitter(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=season_breeder),
  #             width=0, height=0.025,
  #             size=3, alpha=0.2, shape=16)

#by food
food <- visreg(fit.fall, "Dist.log", by="food_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  scale_color_manual(values=c("#794624", "#68b63e")) +
  scale_fill_manual(values=c("#79462450", "#68b63e50")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size=16),
        plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt")) +
  labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  annotate(geom = "text", x=2.4, y=.85, size = 6,
           label = paste("p =",
                         round( coef(summary(fit.fall))[10,4], digits=3) )) +
  annotate(geom = "text", x=2.4, y=.9, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(fit.fall))[10,1]), digits=3) )) +
  geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=food_trt),
             size=3, alpha=0.15, shape=16)
  # geom_jitter(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=food_trt),
  #             width=0, height=0.025,
  #             size=3, alpha=0.2, shape=16)

#by worms
worms <- visreg(fit.fall, "Dist.log", by="helm_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size=16),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  annotate(geom = "text", x=2.4, y=.85, size = 6,
           label = paste("p =",
                         round( coef(summary(fit.fall))[11,4], digits=3) )) +
  annotate(geom = "text", x=2.4, y=.9, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(fit.fall))[11,1]), digits=3) )) +
  geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=helm_trt),
             size=3, alpha=0.1, shape=16)
  # geom_jitter(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=helm_trt),
  #             width=0, height=0.025,
  #             size=3, alpha=0.2, shape=16)


png(filename="fitfall21.png", height=12, width=16, units="in", res=600)
plot_grid(sex, repro, food, worms,
          labels=NULL, nrow=2)
dev.off()





######## pretty OUTPUT MODEL SUMMARY #########

library(gtsummary) #https://www.danieldsjoberg.com/gtsummary/articles/tbl_regression.html

fit.summer %>% tbl_regression(exponentiate = TRUE,
                                pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
      bold_p(t = 0.10) %>%
      bold_labels() %>%
      italicize_levels() %>%
  gtsummary::as_tibble() %>%
  write.csv(here("fit_summer21.csv"))

fit.fall %>% tbl_regression(exponentiate = TRUE,
                              pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
      bold_p(t = 0.10) %>%
      bold_labels() %>%
      italicize_levels() %>% #stop here to get HTML output in RStudio
  gtsummary::as_tibble() %>%
  write.csv(here("fit_fall21.csv"))


# #https://rpubs.com/benhorvath/glm_diagnostics
# plot(density(resid(fit.summer, type='pearson')))
# plot(density(resid(fit.fall, type='pearson')))
#
# plot(density(rstandard(fit.summer, type='pearson')))
# plot(density(rstandard(fit.fall, type='pearson')))
#
# plot(density(resid(fit.summer, type='deviance')))
# plot(density(resid(fit.fall, type='deviance')))
#
# plot(density(rstandard(fit.summer, type='deviance')))
# plot(density(rstandard(fit.fall, type='deviance')))
#
# library(statmod)
# plot(density(qresid(fit.summer)))
# plot(density(qresid(fit.fall)))
#
# par(mfrow=c(1,2))
# plot(density(matrix_dists_obs$Det.obs), main='M_0 y_hat')
# lines(density(predict(fit.fall, type='response')), col='red')
#
#
# qqnorm(statmod::qresid(fit.summer)); qqline(statmod::qresid(fit.summer))
# qqnorm(statmod::qresid(fit.fall)); qqline(statmod::qresid(fit.fall))
