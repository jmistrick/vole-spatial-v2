# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)

#clear environment
rm(list = ls())

params21 <- readRDS(here("params21.rds"))
params22 <- readRDS(here("params22.rds"))

# p <- 0.01
# a <- params21[14,2]
# b <- params21[14,3]
#
# (log((1/0.01)-1) + a) / (-b)

#https://stackoverflow.com/questions/65089775/plotting-custom-functions-in-ggplot-with-variables-from-dataframe







# params21 <- params21 %>% separate_wider_delim(sts, delim="_", names=c("season", "f_trt", "h_trt", "sex")) %>%
#   unite(trt, f_trt, h_trt) %>%
#   mutate(rad_0.001 = (log((1/0.001)-1) + a) / (-b))
#
# params22 <- params22 %>% separate_wider_delim(sts, delim="_", names=c("season", "f_trt", "h_trt", "sex")) %>%
#   unite(trt, f_trt, h_trt) %>%
#   mutate(rad_0.001 = (log((1/0.001)-1) + a) / (-b))

params21 <- params21 %>% separate_wider_delim(sts, delim="_", names=c("season", "f_trt", "h_trt", "sex")) %>%
  unite(trt, f_trt, h_trt) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))

params22 <- params22 %>% separate_wider_delim(sts, delim="_", names=c("season", "f_trt", "h_trt", "sex")) %>%
  unite(trt, f_trt, h_trt) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))

# params21 <- params21 %>% separate_wider_delim(sts, delim="_", names=c("season", "f_trt", "h_trt", "sex")) %>%
#   unite(trt, f_trt, h_trt) %>%
#   mutate(rad_0.05 = (log((1/0.05)-1) + a) / (-b))
#
# params22 <- params22 %>% separate_wider_delim(sts, delim="_", names=c("season", "f_trt", "h_trt", "sex")) %>%
#   unite(trt, f_trt, h_trt) %>%
#   mutate(rad_0.05 = (log((1/0.05)-1) + a) / (-b))



#### OKAY... BUT NOW I FUCKED THINGS UP bc these lifetimes centroids make things dumb in the glmer so I went back to seasonal centroids

centroids21 <- readRDS(here("centroids21.rds"))
centroids22 <- readRDS(here("centroids22.rds"))


netmets21 <- readRDS(here("netmets21.rds"))
netmets22 <- readRDS(here("netmets22.rds"))

circles21 <- netmets21 %>% select(c(site, month, tag)) %>%
  left_join(centroids21, by=c("tag" = "Tag_ID")) %>%
  left_join(read.csv(here("grid_trts.csv")), by="site") %>% unite(trt, food_trt, helm_trt) %>%
  rename(sex = Sex) %>%
  mutate(season = case_when(month == "june" ~ "summer",
                            month == "july" ~ "summer",
                            month == "aug" ~ "summer",
                            month == "sept" ~ "fall",
                            month == "oct" ~ "fall")) %>%
  left_join(params21, by=c("season", "trt", "sex")) %>%
  mutate(month=factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>%
  select(-c(a, b))

library(ggforce)
#https://ggforce.data-imaginist.com/reference/geom_circle.html

# circles21 %>% filter(site=="hevonen" & month=="july") %>%
#   ggplot() +
#   geom_point(aes(x=x, y=y)) +
#   geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
#             fill=NA, alpha = 0.4, color = "black", linetype=2) +
#   geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=sex), alpha=0.5) +
#   coord_fixed()


circles22 <- netmets22 %>% select(c(site, month, tag)) %>%
  left_join(centroids22, by=c("tag" = "Tag_ID")) %>%
  left_join(read.csv(here("grid_trts.csv")), by="site") %>% unite(trt, food_trt, helm_trt) %>%
  rename(sex = Sex) %>%
  mutate(season = case_when(month == "june" ~ "summer",
                            month == "july" ~ "summer",
                            month == "aug" ~ "summer",
                            month == "sept" ~ "fall",
                            month == "oct" ~ "fall")) %>%
  left_join(params22, by=c("season", "trt", "sex")) %>%
  mutate(month=factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>%
  select(-c(a, b))

# circles22 %>% filter(site=="hevonen" & month=="july") %>%
#   ggplot() +
#   geom_point(aes(x=x, y=y)) +
#   geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=sex), alpha=0.5) +
#   geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
#             fill=NA, alpha = 0.4, color = "black", linetype=2) +
#   coord_fixed()


####-------------------------------------------


#go up a level from current wd, then down to file
netmets_puuv <- readRDS(file="../volehantaR/netmets_puuv_06.05.23.rds")

nm_puuv21 <- netmets_puuv %>% filter(year=="2021") %>%
  select(c(year, site, month, tag, prev_curr_puuv)) %>%
  separate_wider_delim(prev_curr_puuv, names=c("prev_puuv", "curr_puuv"), delim="-") %>%
  mutate(color_status = ifelse(prev_puuv=="1", "prev_pos",
                               ifelse(curr_puuv=="1", "new_pos", "neg"))) %>%
  ungroup() %>%
  select(!c(year, prev_puuv, curr_puuv))

nm_puuv22 <- netmets_puuv %>% filter(year=="2022") %>%
  select(c(year, site, month, tag, prev_curr_puuv)) %>%
  separate_wider_delim(prev_curr_puuv, names=c("prev_puuv", "curr_puuv"), delim="-") %>%
  mutate(color_status = ifelse(prev_puuv=="1", "prev_pos",
                               ifelse(curr_puuv=="1", "new_pos", "neg")))  %>%
  ungroup() %>%
  select(!c(year, prev_puuv, curr_puuv))


circles21 <- circles21 %>% left_join(nm_puuv21, by=c("site", "month", "tag"))
circles22 <- circles22 %>% left_join(nm_puuv22, by=c("site", "month", "tag"))




# circles22 %>% filter(site=="hevonen" & month=="july") %>%
#   ggplot() +
#   geom_point(aes(x=x, y=y, color=sex, size=2), show.legend = FALSE) +
#   geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=color_status, color=sex), alpha=0.5) +
#   scale_fill_manual(values=c("#CFCFCF", "#B8EC86", "#F19013"), na.value="white") +
#   geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
#             fill=NA, alpha = 0.4, color = "black", linetype=2) +
#   coord_fixed()





############## TO DO NEXT #####################
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

  png(filename = paste("circles_", "rad0.01_", names(circles21_list)[[i]], "_2021", ".png", sep = ""),
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

  png(filename = paste("circles_", "rad0.01_", names(circles22_list)[[i]], "_2022", ".png", sep = ""),
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

  png(filename = paste("infection_circles_", "rad0.01_", names(circles21_list)[[i]], "_2021", ".png", sep = ""),
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

  png(filename = paste("infection_circles_", "rad0.01_", names(circles22_list)[[i]], "_2022", ".png", sep = ""),
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


