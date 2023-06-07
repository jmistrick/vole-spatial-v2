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



# p <- 0.01
# a <- params21[14,2]
# b <- params21[14,3]
#
# (log((1/0.01)-1) + a) / (-b)

#https://stackoverflow.com/questions/65089775/plotting-custom-functions-in-ggplot-with-variables-from-dataframe



####----------- LOAD DATA -----------------

#params
params21.22 <- readRDS(here("bnb_params21.22.rds"))

#centroids
centroids21.22 <- readRDS(here("bnb_centroids21.22.rds")) %>%
  rename(breeder = Breeder,
         tag = Tag_ID)

#load fulltrap data - pull tag, month, site
ft21 <- readRDS(here("fulltrap21_05.10.23.rds"))
tagmonth21 <- ft21 %>% select(c(year, season, site, month, tag)) %>%
  filter(month != "may") %>%
  group_by(tag, month) %>% slice(1)

ft22 <- readRDS(here("fulltrap22_05.10.23.rds"))
tagmonth22 <- ft22 %>% select(c(year, season, site, month, tag)) %>%
  filter(month != "may") %>%
  group_by(tag, month) %>% slice(1)

ft21.22 <- rbind(ft21, ft22)
tagsampid21.22 <- ft21.22 %>% select(year, month, tag, samp_id) %>%
  distinct(tag, samp_id, .keep_all = TRUE) #remove duplicate rows
tagmonth21.22 <- rbind(tagmonth21, tagmonth22)

circleparts <- left_join(centroids21.22, params21.22, by=c("year", "season", "trt", "sex", "breeder"))

circles21.22 <- left_join(tagmonth21.22, circleparts, by=c("year", "season", "tag")) %>%
  drop_na(x) %>% drop_na(a) %>% #drop rows without circle parameters (may animals, animals w/o sex, breeder, trap etc)
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))

####### I SHOULD DEFINTELY CONFIRM THIS CALCULATION WITH SOMEONE BECAUSE I'M RULL DUMB AT MATH ################






####--------------real quick plots just to see-----------------------------------------

# library(ggforce) #for geom_circle in ggplot
#https://ggforce.data-imaginist.com/reference/geom_circle.html


# circles21.22 %>% filter(year=="2021" & site=="vaarinkorpi" & month=="sept") %>%
#   ggplot() +
#   geom_point(aes(x=x, y=y, color=sex)) +
#   geom_circle( aes(x0=x, y0=y, r=rad_0.01, color=sex), alpha=0.5) +
#   geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
#             fill=NA, alpha = 0.4, color = "black", linetype=2) +
#   coord_fixed()


####-------------------------------------------



########### THE FOLLOWING IS PULLED FROM THE HANTA CLEANING FILE AND EDITED TO INCLUDE NONBREEDERS #############
#########################################   LOAD & CLEAN PUUV IFA DATA   ########################################

#load, clean, format PUUV IFA data
#go up a level from current wd, then down to file for puuv data
puuv_data <- read.csv(file="../volehantaR/puuv_ifa.csv") %>%
  clean_names %>%
  #populate column of FINAL PUUV status (result of second run if two runs were done, else result of first run)
  mutate(FINAL_puuv = ifelse(is.na(puuv_confirm), as.character(puuv_initial), as.character(puuv_confirm))) %>%
  mutate(samp_id = as.numeric(id),
         date_initial = as_date(date_initial, format= "%m/%d/%Y"),
         puuv_initial = as.factor(puuv_initial),
         date_confirm = as_date(date_confirm, format= "%m/%d/%Y"),
         puuv_confirm = as.factor(puuv_confirm),
         FINAL_puuv = as.factor(FINAL_puuv)) %>%
  drop_na(FINAL_puuv) %>%
  dplyr::select(FINAL_puuv, samp_id) %>%
  rename(puuv_ifa = FINAL_puuv) %>%
  left_join(tagsampid21.22, by="samp_id") %>% select(!samp_id) %>%
  drop_na(tag)
#output is a df with year, month, tag, and PUUV status (0,1) (BOTH YEARS!)

############ REMOVE THE VOLES from puuv_data THAT SEROCONVERT POS TO NEG ##################

#voles that seroconvert PUUV+ to PUUV-
puuv_pos_neg <- puuv_data %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
  dplyr::select(year, month, tag, puuv_ifa) %>%
  summarise(status_time = toString(puuv_ifa)) %>%
  filter(str_detect(status_time, "1,\\s0")) #filter for animals that go from pos to neg
#pull the PIT tags
problemchildren <- puuv_pos_neg$tag
#filter netmetsPUUV to remove 'problemchildren'
puuv_data <- puuv_data %>%
  filter(!tag %in% problemchildren)


# add previous (lagged) degree (degree from previous month influences current PUUV status)
# add 0,1 for serovert - animals that go 0-0 or 0-1
# BUT! the previous month has to be in the same year (don't want 2021 fall to influence 2022 spring)
puuv <- puuv_data %>% group_by(year, tag) %>%
  arrange(month, .by_group = TRUE) %>%
  mutate(prev_puuv = lag(puuv_ifa)) %>%
  rename(curr_puuv = puuv_ifa) %>%
  mutate(color_status = case_when(prev_puuv=="1" ~ "prev_pos",
                                  curr_puuv=="1" ~ "new_pos",
                                  curr_puuv=="0" ~ "neg")) %>%
  select(year, month, tag, color_status)

## NOW PUUV has:
  # breeders and nonbreeders
  # is not tied to netmets data rn
  # new color_status column of new/prev positives

###############################################################################



#create the 'circles' dataframes with all the circle dimensions and infection colors
circles21 <- circles21.22 %>% filter(year=="2021") %>%
  left_join(puuv, by=c("year", "month", "tag")) %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) #remove may from factor

circles22 <- circles21.22 %>% filter(year=="2022") %>%
  left_join(puuv, by=c("year", "month", "tag")) %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) #remove may from factor



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

  png(filename = paste("ERRBODYcircles_", "rad0.01_", names(circles21_list)[[i]], "_2021", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles21_list[[i]])){

    data <- circles21_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=sex)) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=sex), alpha=0.5) +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "black", linetype=2) +
      theme(legend.position = "bottom") +
      labs(title=paste(names(circles21_list[[i]])[j])) +
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

