####### ------------ CriCkLES FOR BREEDERS / NON-BREEDERS -----------------------

# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)
library(ggforce) #for geom_circle in ggplot

#clear environment
rm(list = ls())

######## HOME RANGE DISTRIBUTION ###########
#negative sigmoidal curve calculated following Wanelik & Farine 2022: DOI 10.1007/s00265-022-03222-5
#where the declining probability (P) of an individual being detected at a distance (d)
#from the centroid of its home range is given by:

P(d) = 1 / (1 + e^(-a-bd))

#where a describes the overall size of the home range
#b describes the steepness of the edge of the home range
#and d is the logarithmic distance from the centroid

# p <- 0.01 #eg probability of detection 1%
# a <- params21[14,2]
# b <- params21[14,3]
#
# (log((1/p)-1) + a) / (-b)

######### ^^ this equation is what I'm using to calculate e.g. the distance from the centroid with 1% probability of finding the animal,
  ##just to plot some figure showing approx HR size to visualize the overlaps



#https://stackoverflow.com/questions/65089775/plotting-custom-functions-in-ggplot-with-variables-from-dataframe



####----------- LOAD DATA -----------------

#params
params21 <- readRDS(here("params21_STSB.rds"))
params22 <- readRDS(here("params22_STSB.rds"))



### another thing 6/30 for volespatial ms ###
## Figure for Vole Spatial manuscript ##

trt.labs <- as_labeller(c("unfed_control" = "Unfed-Control",
                          "unfed_deworm" = "Unfed-Deworm",
                          "fed_control" = "Fed-Control",
                          "fed_deworm" = "Fed-Deworm"))

season.labs <- as_labeller(c("summer" = "Summer",
                             "fall" = "Autumn"))

png(filename = here("spaceuse_sex_breed_label_0.01.png"), height=7, width = 16, units = "in", res=600)
readRDS(here("params21_STSB.rds")) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b),
         area = paste( round(pi*(rad_0.01^2)*10, digits = 2), "m\u00B2" )) %>%
separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
  unite(trt, food_trt, helm_trt) %>%
  mutate(x = case_when(sex=="F" ~ 4,
                       sex=="M" ~ 8.5),
         y = case_when(season_breeder=="breeder" ~ 11,
                       season_breeder=="nonbreeder" ~ 7),
         repro = case_when(season_breeder=="breeder" ~ "Reproductive",
                           season_breeder=="nonbreeder" ~ "Non-Reproductive"),
         y_lab = case_when(season_breeder=="breeder" ~ y,
                           season_breeder=="nonbreeder" & season=="summer" ~ y-1.75,
                           season_breeder=="nonbreeder" & season=="fall" ~ y-2.25)) %>%
  mutate(season = factor(season, levels=c("summer", "fall")),
         trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm")),
         repro = factor(repro, levels=c("Reproductive", "Non-Reproductive"))) %>%
  ggplot() +
  geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=sex, linetype=repro), linewidth=0.8) +
  # geom_point(aes(x=x, y=y), size=1) +
  geom_text(aes(x=x, y=y_lab, label=area), hjust=0.5, vjust=0.5, size=5) +
  # scale_color_manual(values=c("#f282a7", "#00d0ff")) +
  scale_fill_manual(values=c("#f282a780", "#00d0ff80")) +
  facet_grid(season~trt, labeller=labeller(trt=trt.labs, season=season.labs)) +
  coord_fixed() +
  labs(fill="Sex:", linetype="Reproductive Status:") +
  theme(legend.position="bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size=19),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) #adjusts order of legends to show sex first then repro
dev.off()
#also
ggsave(file = "spaceuse_sex_breed_label.eps")

########### end figure for Vole Spatial ms ######################




### just a lil thing 6/16 ###

#mean radius in 2021 by season, sex, breed
y21 <- readRDS(here("params21_STSB.rds")) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
  mutate(area = 2*pi*(rad_0.01^2)) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  group_by(season, breeder, helmtrt) %>%
  summarize(mean = mean(area),
            sd = sd(area))

#visualize M/F breeder/non size per trt from summer to fall
readRDS(here("params21_STSB.rds")) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  unite(trt, foodtrt, helmtrt) %>%
  ggplot(aes(x=season, y=rad_0.01, color=sex, shape=breeder)) +
           geom_point() +
  facet_wrap(~trt)


#mean radius in 2021 by season, sex, breed
y22 <- readRDS(here("params22_STSB.rds")) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
  mutate(area = 2*pi*(rad_0.01^2)) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  group_by(season, sex, breeder) %>%
  summarize(mean = mean(area))

#visualize M/F breeder/non size per trt from summer to fall
readRDS(here("params22_STSB.rds")) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  unite(trt, foodtrt, helmtrt) %>%
  ggplot(aes(x=season, y=rad_0.01, color=sex, shape=breeder)) +
  geom_point() +
  facet_wrap(~trt)


#mean across both years, all trts
table <- rbind(readRDS(here("params21_STSB.rds")), readRDS(here("params22_STSB.rds"))) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
  mutate(area = 2*pi*(rad_0.01^2)) %>%
  mutate(area = area*10) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  unite(trt, foodtrt, helmtrt) %>%
  group_by(season, breeder, sex) %>%
  summarize(mean = mean(area),
            sd = sd(area))


#############################################










#MONTHLY centroids
centroids21 <- readRDS(here("centroids21_STSB.rds")) %>% rename(tag = Tag_ID)
centroids22 <- readRDS(here("centroids22_STSB.rds")) %>% rename(tag = Tag_ID)

#load fulltrap data - pull tag, month, site
ft21 <- readRDS(here("fulltrap21_05.10.23.rds"))
trapdat21 <- ft21 %>% select(c(year, season, trt, site, month, tag, sex, season_breeder)) %>%
  filter(month != "may") %>%
  group_by(tag, month) %>% slice(1) %>%
  drop_na(sex) %>% drop_na(season_breeder)

ft22 <- readRDS(here("fulltrap22_05.10.23.rds"))
trapdat22 <- ft22 %>% select(c(year, season, trt, site, month, tag, sex, season_breeder)) %>%
  filter(month != "may") %>%
  group_by(tag, month) %>% slice(1) %>%
  drop_na(sex) %>% drop_na(season_breeder)



####----------- CREATE 'CIRCLES' df for PLOTTING -----------------

circles21 <- left_join(centroids21, trapdat21, by=c("tag", "month", "site")) %>%
  unite(stsb, season, trt, sex, season_breeder) %>% left_join(params21, by="stsb") %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #remove may from factor
  separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
  unite(trt, food_trt, helm_trt) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))

circles22 <- left_join(centroids22, trapdat22, by=c("tag", "month", "site")) %>%
  unite(stsb, season, trt, sex, season_breeder) %>% left_join(params22, by="stsb") %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #remove may from factor
  separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
  unite(trt, food_trt, helm_trt) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))

####### I SHOULD DEFINTELY CONFIRM THIS CALCULATION WITH SOMEONE BECAUSE I'M RULL DUMB AT MATH ################
#matt M-S said it was fine :)



####--------------real quick plots just to see----------------------

# library(ggforce) #for geom_circle in ggplot
#https://ggforce.data-imaginist.com/reference/geom_circle.html

# circles21 %>% filter(site=="vaarinkorpi" & month=="sept") %>%
#   ggplot() +
#   geom_point(aes(x=x, y=y, color=sex)) +
#   geom_circle( aes(x0=x, y0=y, r=rad_0.01, color=sex), alpha=0.5) +
#   geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
#             fill=NA, alpha = 0.4, color = "black", linetype=2) +
#   coord_fixed()


####-------------------------------------------


##########Janine 6.12 says: ########## THIS NEEDS TO BE UPDATED TO WHATER THE PUUV DATA IS THESE DAYS ####################

# ########### THE FOLLOWING IS PULLED FROM THE HANTA CLEANING FILE AND EDITED TO INCLUDE NONBREEDERS #############
# #########################################   LOAD & CLEAN PUUV IFA DATA   ########################################
#
# #load, clean, format PUUV IFA data
# #go up a level from current wd, then down to file for puuv data
# puuv_data <- read.csv(file="../volehantaR/puuv_ifa.csv") %>%
#   clean_names %>%
#   #populate column of FINAL PUUV status (result of second run if two runs were done, else result of first run)
#   mutate(FINAL_puuv = ifelse(is.na(puuv_confirm), as.character(puuv_initial), as.character(puuv_confirm))) %>%
#   mutate(samp_id = as.numeric(id),
#          date_initial = as_date(date_initial, format= "%m/%d/%Y"),
#          puuv_initial = as.factor(puuv_initial),
#          date_confirm = as_date(date_confirm, format= "%m/%d/%Y"),
#          puuv_confirm = as.factor(puuv_confirm),
#          FINAL_puuv = as.factor(FINAL_puuv)) %>%
#   drop_na(FINAL_puuv) %>%
#   dplyr::select(FINAL_puuv, samp_id) %>%
#   rename(puuv_ifa = FINAL_puuv) %>%
#   left_join(tagsampid21.22, by="samp_id") %>% select(!samp_id) %>%
#   drop_na(tag)
# #output is a df with year, month, tag, and PUUV status (0,1) (BOTH YEARS!)
#
# ############ REMOVE THE VOLES from puuv_data THAT SEROCONVERT POS TO NEG ##################
#
# #voles that seroconvert PUUV+ to PUUV-
# puuv_pos_neg <- puuv_data %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "1,\\s0")) #filter for animals that go from pos to neg
# #pull the PIT tags
# problemchildren <- puuv_pos_neg$tag
# #filter netmetsPUUV to remove 'problemchildren'
# puuv_data <- puuv_data %>%
#   filter(!tag %in% problemchildren)
#
#
# # add previous (lagged) degree (degree from previous month influences current PUUV status)
# # add 0,1 for serovert - animals that go 0-0 or 0-1
# # BUT! the previous month has to be in the same year (don't want 2021 fall to influence 2022 spring)
# puuv <- puuv_data %>% group_by(year, tag) %>%
#   arrange(month, .by_group = TRUE) %>%
#   mutate(prev_puuv = lag(puuv_ifa)) %>%
#   rename(curr_puuv = puuv_ifa) %>%
#   mutate(color_status = case_when(prev_puuv=="1" ~ "prev_pos",
#                                   curr_puuv=="1" ~ "new_pos",
#                                   curr_puuv=="0" ~ "neg")) %>%
#   select(year, month, tag, color_status)
#
# ## NOW PUUV has:
#   # breeders and nonbreeders
#   # is not tied to netmets data rn
#   # new color_status column of new/prev positives
#
# ###############################################################################



####------ ADD PUUV DATA FOR INFECTION STATUS ---------------

# #create the 'circles' dataframes with all the circle dimensions and infection colors
# circles21 <- circles21.22 %>% filter(year=="2021") %>%
#   left_join(puuv, by=c("year", "month", "tag")) %>%
#   mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) #remove may from factor
#
# circles22 <- circles21.22 %>% filter(year=="2022") %>%
#   left_join(puuv, by=c("year", "month", "tag")) %>%
#   mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) #remove may from factor



# # quick visual
# circles22 %>% filter(site=="luostari" & month=="july") %>%
#   ggplot() +
#   geom_point(aes(x=x, y=y, color=sex, size=2), show.legend = FALSE) +
#   geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=color_status, color=sex), alpha=0.5) +
#   scale_fill_manual(values=c("#999999", "#B8EC86", "#F19013"), na.value=NA) +
#   geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
#             fill=NA, alpha = 0.4, color = "black", linetype=2) +
#   coord_fixed()





######## PLOT MULTIPLE SITES ACROSS MONTHS #############

##---------------- CREATE A NESTED LIST of CIRCLES (nested by site, month) --------------------

#save a vector of the site names (alphabetical order)
site_names <- unique(circles21$site) %>% sort()

#split() makes a list consisting of individual data.frames based on a condition ('site' in this case)
circles21_site_list <- split(circles21, circles21$site)

#create new list to hold nested site, month capture
circles21_list <- list()

for(i in 1:length(circles21_site_list)){
  # print(i)
  temp_list <- split(circles21_site_list[[i]], circles21_site_list[[i]]$month)
  circles21_list[[i]] <- temp_list
}

#name 1e list element as site (2e list elements are months May-Oct)
names(circles21_list) <- site_names

### OUTPUT: CIRCLES##_LIST is a nested list of length 12
#1e level is all the sites (12)
#2e level is all the months per site (5) excluding May

####--------------------repeat for 2022 data----------------------

#save a vector of the site names (alphabetical order)
site_names <- unique(circles22$site) %>% sort()

#split() makes a list consisting of individual data.frames based on a condition ('site' in this case)
circles22_site_list <- split(circles22, circles22$site)

#create new list to hold nested site, month capture
circles22_list <- list()

for(i in 1:length(circles22_site_list)){
  # print(i)
  temp_list <- split(circles22_site_list[[i]], circles22_site_list[[i]]$month)
  circles22_list[[i]] <- temp_list
}

#name 1e list element as site (2e list elements are months May-Oct)
names(circles22_list) <- site_names

### OUTPUT: CIRCLES##_LIST is a nested list of length 12
#1e level is all the sites (12)
#2e level is all the months per site (5) excluding May


####----------------------------END----------------------------------


####------------ PLOT [JUST CIRCLES] IN A LOOP (per YEAR)---------------------------

library(gridExtra)

for(i in 1:length(circles21_list)) {

  png(filename = paste("ERRBODYcircles_monthlycentroids_", "rad0.01_", names(circles21_list)[[i]], "_2021", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles21_list[[i]])){

    data <- circles21_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=sex), show.legend=FALSE) +
      xlim(-2,14) + ylim(-2,13) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=sex), alpha=0.5) +
      scale_fill_manual(values=c("#f282a780", "#00d0ff80")) +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "black", linetype=2) +
      theme(legend.position = "bottom",
            axis.ticks = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank()) +
      labs(title=paste(names(circles21_list[[i]])[j]),
           fill="Sex") +
      coord_fixed()

  }

  do.call(grid.arrange, c(p, ncol=5, top=paste( c(names(circles21_list)[[i]]), "2021") )) #plot all the plots in list p

  dev.off()

}

####------------------repeat for 2022---------------------------

for(i in 1:length(circles22_list)) {

  png(filename = paste("ERRBODYcircles_", "rad0.01_", names(circles22_list)[[i]], "_2022", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles22_list[[i]])){

    data <- circles22_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=sex)) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=sex), alpha=0.5) +
      scale_fill_manual(values=c("#f282a780", "#00d0ff80")) +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "black", linetype=2) +
      theme(legend.position = "bottom") +
      labs(title=paste(names(circles22_list[[i]])[j])) +
      coord_fixed()

  }

  do.call(grid.arrange, c(p, ncol=5, top=paste( c(names(circles22_list)[[i]]), "2022") )) #plot all the plots in list p

  dev.off()

}

#####----------------------------------------------------------------------------


####------------ PLOT [INFECTION PRE - CURRENT] IN A LOOP (per YEAR)---------------------------

library(gridExtra)

for(i in 1:length(circles21_list)) {

  png(filename = paste("ERRBODYinfection_circles_", "rad0.01_", names(circles21_list)[[i]], "_2021", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles21_list[[i]])){

    data <- circles21_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=sex, size=2), show.legend = FALSE) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=color_status, color=sex), alpha=0.5) +
      scale_fill_manual(values=c("#CFCFCF", "#B8EC86", "#F19013"), na.value="white") +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "black", linetype=2) +
      theme(legend.position = "bottom") +
      labs(title=paste(names(circles21_list[[i]])[j])) +
      coord_fixed()

  }

  do.call(grid.arrange, c(p, ncol=5, top=paste( c(names(circles22_list)[[i]]), "2021") )) #plot all the plots in list p

  dev.off()

}

####------------------repeat for 2022---------------------------

for(i in 1:length(circles22_list)) {

  png(filename = paste("ERRBODYinfection_circles_", "rad0.01_", names(circles22_list)[[i]], "_2022", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles22_list[[i]])){

    data <- circles22_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=sex, size=2), show.legend = FALSE) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=color_status, color=sex), alpha=0.5) +
      scale_fill_manual(values=c("#CFCFCF", "#B8EC86", "#F19013"), na.value="white") +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "black", linetype=2) +
      theme(legend.position = "bottom") +
      labs(title=paste(names(circles22_list[[i]])[j])) +
      coord_fixed()

  }

  do.call(grid.arrange, c(p, ncol=5, top=paste( c(names(circles22_list)[[i]]), "2022") )) #plot all the plots in list p

  dev.off()

}

#####----------------------------------------------------------------------------


