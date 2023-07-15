# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)
library(ggridges)

#clear environment
rm(list = ls())


ft21 <- readRDS(here("fulltrap21_05.10.23.rds")) %>% ##breeders and nonbreeders are here
  filter(month != "may") %>%
  drop_na(sex) %>% drop_na(season_breeder) %>%
  unite(sb, sex, season_breeder, remove = FALSE)

# # summary stats on number of voles, multicaps, mean caps per life etc
# n_distinct(ft21$tag) #742 individuals
# onepertag <- ft21 %>% group_by(tag) %>% slice(1) %>% ungroup()
# onepertag %>% summarise(mean = mean(caps_per_life),
#             sd = sd(caps_per_life),
#             min = min(caps_per_life),
#             max = max(caps_per_life))
# multicap <- ft21 %>% filter(caps_per_life >= 2)
# n_distinct(multicap$tag) #448 individuals (59.5%) captured at least 2x

tagsex <- ft21 %>% group_by(tag) %>% slice(1) %>% select(tag, sex)

netmets21 <- readRDS(here("netmets21_bindeg.rds")) %>% left_join(read.csv(here('grid_trts.csv')), by="site") %>%
  unite(trt, food_trt, helm_trt) %>%
  left_join(tagsex, by="tag")

# #mean degree by treatment, month
# netmets21 %>%
#   group_by(month, trt) %>%
#   summarise(mean = mean(norm.wt.deg),
#             median = median(norm.wt.deg))



month.labs <- as_labeller(c("june" = "June", "july" = "July", "aug" = "Aug", "sept" = "Sept", "oct" = "Oct"))

# #ridgeline plot of weighted degree (density), each site separate
# netmets21 %>%
#   mutate(site = fct_relevel(site, "janakkala", "puro", "talo",
#                             "kuoppa", "radio", "vaarinkorpi",
#                             "ketunpesa", "kiirastuli", "mustikka",
#                             "asema", "helmipollo", "hevonen")) %>%
#   mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
#   mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))) %>%
#   ggplot( aes(y=site, x=wt.deg,  fill=trt)) +
#   # stat_density_ridges(quantile_lines=TRUE,
#   #                     quantile_fun=function(x,...)median(x),
#   #                     alpha=0.6) +
#   geom_density_ridges(alpha=0.75, stat="binline", scale=1) + #stat_binline plots histogram
#   # geom_density_ridges(alpha=0.6) + #density plot
#   facet_wrap(~month, labeller=month.labs, ncol=5) +
#   xlim(-1, 7) +
#   scale_fill_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A"),
#                     name = "Treatment",
#                     labels = c("Unfed - Control", "Unfed - Deworm", "Fed - Control", "Fed - Deworm")) +
#   theme(axis.ticks.y=element_blank(),
#         axis.text.y=element_blank(),
#         legend.position="none",
#         legend.text=element_text(size=9),
#         axis.title=element_text(size=15),
#         panel.spacing = unit(0.2, "lines"),
#         strip.text.x = element_text(size = 10)) +
#   labs(x="Spatial Overlap Strength", y="Replicate Population", fill="Treatment")


#overlay density plots
# netmets21 %>%
#   mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
#   mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))) %>%
#   ggplot(aes(x=wt.deg, color=trt)) +
#   geom_density(alpha=.25, size=1) +
#   facet_wrap(~month, ncol=5)


#ridgeline plot (density) of weighted degree, sites combined by trt
wtdeg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  ggplot( aes(y=trt, x=wt.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  # geom_density_ridges(alpha=0.75, stat="binline", bins=55, scale=0.9) + #histogram
  # geom_density_ridges(alpha=0.6, scale=0.9) + #density plot
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  # facet_grid(sex~month) + #by sex
  # xlim(0,4) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  # theme_half_open() + panel_border() + background_grid(major="y") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Weighted Degree",
       fill="Treatment") #y="Density", title="Histogram of Weighted Degree Distribution"


#binary degree (threshold wt > 0.01)
bindeg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  ggplot( aes(y=trt, x=bin.01.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  # geom_density_ridges(alpha=0.75, stat="binline", bins=55, scale=0.9) + #histogram
  # geom_density_ridges(alpha=0.6, scale=0.9, rel_min_height=0.01) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  # facet_grid(sex~month) + #by sex
  # xlim(0, 40) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  # theme_half_open() + panel_border() + background_grid(major="y") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Unweighted Degree", y="Count",
       fill="Treatment") #title="Histogram of Binary Degree Distribution"

#'normalized' binary degree (# voles you overlap with/# you could overlap with)
normbindeg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  mutate(norm.bin.deg = bin.01.deg/(n.node-1)) %>%
  ggplot( aes(y=trt, x=norm.bin.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  # geom_density_ridges(alpha=0.75, stat="binline", bins=55, scale=0.9) + #histogram
  # geom_density_ridges(alpha=0.6, scale=0.9, rel_min_height = 0.01) + #density plot
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  # facet_grid(sex~month) + #by sex
  # xlim(0, 1) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  # theme_half_open() + panel_border() + background_grid(major="y") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.title = element_text(size=14),
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Normalized Unweighted Degree", y="Count",
       fill="Treatment") #title="Histogram of NORMALIZED Binary Degree Distribution"


library(cowplot)
png(filename = 'degree_distributions-flip.png', width = 12, height = 12, units = "in", res=600)
plot_grid(wtdeg, bindeg, normbindeg,
          ncol = 1)
dev.off()



#mean and sd norm bin deg
netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))) %>%
  mutate(season = case_when(month=="june" ~ "summer",
                            month=="july" ~ "summer",
                            month=="aug" ~ "summer",
                            month=="sept" ~ "fall",
                            month=="oct" ~ "fall")) %>%
  mutate(norm.bin.deg = bin.01.deg/(n.node-1)) %>%
  group_by(season, trt) %>%
  summarise(mean = mean(norm.bin.deg),
            sd = sd(norm.bin.deg))


netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))) %>%
  ggplot( aes(y=trt, x=norm.wt.deg,  fill=trt)) +
  # stat_density_ridges(quantile_lines=TRUE,
  #                     quantile_fun=function(x,...)median(x),
  #                     alpha=0.6) +
  geom_density_ridges(alpha=0.75, stat="binline", bins=55, scale=0.9) +
  #geom_density_ridges(alpha=0.6, bandwidth=1) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  # xlim(-1, 5) +
  scale_fill_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A"),
                    name = "Treatment",
                    labels = c("Unfed - Control", "Unfed - Deworm", "Fed - Control", "Fed - Deworm")) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=9),
        axis.title=element_text(size=15),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 10)) +
  labs(x="Normalized Wt Degree", y="Count",
       fill="Treatment", title="Histogram of Normalized Wt Degree Distribution")









###############################################################################
## code pulled for volenetworks R project, 02.1_trapping_networks_viz.R file ##
###############################################################################

library(kableExtra)

#to do the network metrics, we need a table that only has one entry per site per month
#(otherwise its a weighted average which isn't right)
netmets_onepersite <- netmets21 %>%
  group_by(month, trt, site) %>%
  slice(1)

#Network Size
sumstats_netsize <- netmets_onepersite %>%
  group_by(month, trt) %>%
  summarise(Mean = round(mean(n.node),2),
            SD = round(sd(n.node),2),
            Min = min(n.node),
            Max = max(n.node)) %>%
  mutate(trt=case_when(trt == "unfed_control" ~ "Unfed-Control",
                       trt == "unfed_deworm" ~ "Unfed-Deworm",
                       trt == "fed_control" ~ "Fed-Control",
                       trt == "fed_deworm" ~ "Fed-Deworm")) %>%
  mutate(month = case_when(month=="june" ~ "June",
                           month=="july" ~ "July",
                           month=="aug" ~ "Aug",
                           month=="sept" ~ "Sept",
                           month=="oct" ~ "Oct")) %>%
  mutate(month = as.factor(month),
         month = factor(month, levels=c("June", "July", "Aug", "Sept", "Oct"))) %>%
  mutate(trt = as.factor(trt),
         trt = factor(trt, levels=c("Unfed-Control", "Unfed-Deworm", "Fed-Control", "Fed-Deworm"))) %>%
           rename(c(Treatment=trt,  Month=month)) %>%
  arrange(Month, Treatment)

kbl(sumstats_netsize) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = F)



#Weighted Degree
sumstats_wtdeg <- netmets21 %>%
  group_by(month, trt) %>%
  summarise(N = length(tag),
            Mean = round(mean(wt.deg),2),
            SD = round(sd(wt.deg),2),
            Min = round(min(wt.deg),2),
            Max = round(max(wt.deg),2)) %>%
  mutate(trt=case_when(trt == "unfed_control" ~ "Unfed-Control",
                       trt == "unfed_deworm" ~ "Unfed-Deworm",
                       trt == "fed_control" ~ "Fed-Control",
                       trt == "fed_deworm" ~ "Fed-Deworm")) %>%
  mutate(month = case_when(month=="june" ~ "June",
                           month=="july" ~ "July",
                           month=="aug" ~ "Aug",
                           month=="sept" ~ "Sept",
                           month=="oct" ~ "Oct")) %>%
  mutate(month = as.factor(month),
         month = factor(month, levels=c("June", "July", "Aug", "Sept", "Oct"))) %>%
  mutate(trt = as.factor(trt),
         trt = factor(trt, levels=c("Unfed-Control", "Unfed-Deworm", "Fed-Control", "Fed-Deworm"))) %>%
  rename(c(Treatment=trt,  Month=month)) %>%
  arrange(Month, Treatment)

kbl(sumstats_wtdeg) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = F)


#Unweighted Degree
sumstats_deg <- netmets21 %>%
  group_by(month, trt) %>%
  summarise(N = length(tag),
            Mean = round(mean(bin.01.deg),2),
            SD = round(sd(bin.01.deg),2),
            Min = min(bin.01.deg),
            Max = max(bin.01.deg)) %>%
  mutate(trt=case_when(trt == "unfed_control" ~ "Unfed-Control",
                       trt == "unfed_deworm" ~ "Unfed-Deworm",
                       trt == "fed_control" ~ "Fed-Control",
                       trt == "fed_deworm" ~ "Fed-Deworm")) %>%
  mutate(month = case_when(month=="june" ~ "June",
                           month=="july" ~ "July",
                           month=="aug" ~ "Aug",
                           month=="sept" ~ "Sept",
                           month=="oct" ~ "Oct")) %>%
  mutate(month = as.factor(month),
         month = factor(month, levels=c("June", "July", "Aug", "Sept", "Oct"))) %>%
  mutate(trt = as.factor(trt),
         trt = factor(trt, levels=c("Unfed-Control", "Unfed-Deworm", "Fed-Control", "Fed-Deworm"))) %>%
  rename(c(Treatment=trt,  Month=month)) %>%
  arrange(Month, Treatment)

kbl(sumstats_deg) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = F)


#NORMALIZED Unweighted Degree
sumstats_normdeg <- netmets21 %>%
  mutate(norm.bin.deg = bin.01.deg/(n.node-1)) %>%
  group_by(month, trt) %>%
  summarise(N = length(tag),
            Mean = round(mean(norm.bin.deg),2),
            SD = round(sd(norm.bin.deg),2),
            Min = round(min(norm.bin.deg),2),
            Max = round(max(norm.bin.deg),2)) %>%
  mutate(trt=case_when(trt == "unfed_control" ~ "Unfed-Control",
                       trt == "unfed_deworm" ~ "Unfed-Deworm",
                       trt == "fed_control" ~ "Fed-Control",
                       trt == "fed_deworm" ~ "Fed-Deworm")) %>%
  mutate(month = case_when(month=="june" ~ "June",
                           month=="july" ~ "July",
                           month=="aug" ~ "Aug",
                           month=="sept" ~ "Sept",
                           month=="oct" ~ "Oct")) %>%
  mutate(month = as.factor(month),
         month = factor(month, levels=c("June", "July", "Aug", "Sept", "Oct"))) %>%
  mutate(trt = as.factor(trt),
         trt = factor(trt, levels=c("Unfed-Control", "Unfed-Deworm", "Fed-Control", "Fed-Deworm"))) %>%
  rename(c(Treatment=trt,  Month=month)) %>%
  arrange(Month, Treatment)

kbl(sumstats_normdeg) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = F)



######################### 12.2.22 brief foray into model fitting for network size #######################################

library(lme4)

#network size per site per month
net_size <- netmets21 %>%
  group_by(month, site) %>%
  slice(1) %>%
  dplyr::select(site,trt,month,n.node) %>%
  separate(trt, c("food", "worms")) %>%
  rename(Site = site,
         network_size = n.node) %>%
  #make month numeric so it's easier to see change through time
  mutate(Month = case_when(month=="june" ~ 0,
                           month=="july" ~ 1,
                           month=="aug" ~ 2,
                           month=="sept" ~ 3,
                           month=="oct" ~ 4)) %>%
  #make food and worms numeric so the reference group is what we want
  mutate(Food_Addition = case_when(food=="fed" ~ 1,
                                   food=="unfed" ~ 0),
         Helminth_Removal = case_when(worms=="control" ~ 0,
                                      worms=="deworm" ~ 1))

# #run a linear model
# netsize_mod <- lm(n.node ~ food + worms + month , data=net_size)
# summary(netsize_mod)
# plot(netsize_mod) #ehh things might be a little funky
# #formal test!
# shapiro.test(residuals(netsize_mod)) #yeah, residuals aren't normal
#
# #just to check: look at the distribution of the network sizes
# hist(net_size$n.node) #kind of skewed
#
# #so what if we log-transform n.node?
# #run a linear model
# netsize_mod <- lm(log(n.node) ~ food + worms + month , data=net_size)
# summary(netsize_mod) #looks better (estimates are scaled 0-1)
# plot(netsize_mod) #looks better
# # but confirm with a formal test!
# shapiro.test(residuals(netsize_mod)) #yaas the pvalue is large, we gucci
#
# #and confirm
# hist(log(net_size$n.node)) #yep, more normal shaped


#next step, see if we can make it run in a lmer
netsize_mod <- lmer(log(network_size) ~ Food_Addition + Helminth_Removal + Month + (1|Site), data=net_size)
summary(netsize_mod)

#summary table of mixed model output
#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
library(sjPlot)
tab_model(netsize_mod)

## MODEL DIAGNOSTICS - CODE FROM JASMINE (this does the same as below but baseR)
plot(netsize_mod)
plot(resid(netsize_mod))
qqnorm(resid(netsize_mod))
qqline(resid(netsize_mod))
hist(resid(netsize_mod))
#construct CI
confint(netsize_mod)
#chi square for pvalues (Ben BOLKER SAID TO DO THIS)
drop1(netsize_mod, test="Chisq")
#easy code for figures
library(visreg)
visreg(netsize_mod) #plotting partial residuals, change one parameter, control the others
#shows the figures for each of the predictors
#shows the predicted values based on the model
##use this after running the model (look at the model visually in addition to the model summary)


## MODEL DIAGNOSTICS - https://dfzljdn9uc3pi.cloudfront.net/2020/9522/1/MixedModelDiagnostics.html <- #MattSilkfortheWin
#STEP 1: plot standardized residuals vs fitted values
#resid vs fitted check for homoscedasticity
plot(netsize_mod, resid(., scaled=TRUE) ~ fitted(.),
     abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals") #looks like a nice cloud

#STEP 1.5: plot residuals vs explanatory values, tests for linearity
#In our case we plot standardised residuals against food
plot(netsize_mod, resid(., scaled=TRUE) ~ Food_Addition,
     abline = 0, pch=16, xlab="Food",ylab="Standardised residuals")
#and plot standardised residuals against worms
plot(netsize_mod, resid(., scaled=TRUE) ~ Helminth_Removal,
     abline = 0, pch=16, xlab="Worms",ylab="Standardised residuals")

#STEP 2: Residuals split by random effect groupings (check for homosce... in RE group too)
#We can use a default plot of the standardised residuals versus fitted
# #But add | <random effect> to see it per level of the random effect
# plot(netsize_mod, resid(., scaled=TRUE) ~ fitted(.)| site,
#      abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals")
#   #hard to say since only 5 points per group
#You can also compare the difference in residuals between random effect levels using default boxplots
plot(netsize_mod, as.factor(Site) ~ resid(., scaled=TRUE),
     abline=0,pch=16,xlab="Standardised residuals",ylab="Sites")

#STEP 3: QQ plot to make sure residuals are normally distributed
#qq plot - check that residuals are normally distributed
qqnorm(resid(netsize_mod),pch=16)
qqline(resid(netsize_mod)) #meh, decent?

#STEP 4: Leverage and Cook's Distance - check for INFLUENTIAL DATA POINTS
#Calculate leverage
lev<-hat(model.matrix(netsize_mod))
#Plot leverage against standardised residuals
plot(resid(netsize_mod,type="pearson")~lev,las=1,ylab="Standardised residuals",xlab="Leverage")
#Calculate Cook's Distance
cd<-cooks.distance(netsize_mod)
max(cd) #If Cook's distance is greater than 1 this highlights problematic datapoints
#Plot leverage and Cook's distance together
par(mfrow=c(1,1))
plot(lev,pch=16,col="red",ylim=c(0,1.2),las=1,ylab="Leverage/Cook's distance value")
points(cd,pch=17,col="blue")
points(x=150,y=1.1,pch=16,col="red")
points(x=150,y=0.9,pch=17,col="blue")
text(x=155,y=1.1,"Leverage",adj=c(0,0.5))
text(x=155,y=0.9,"Cook's distance",adj=c(0,0.5))

#STEP 5: Inspect the RE distribution (should be Gaussian distributed)
#It is possible to do this using a histogram but this can be unclear if there are few random effect levels as in our example
#A dotplot is a more effective way of doing this.
lattice::dotplot(ranef(netsize_mod, condVar=TRUE))

#################################################end###################################################################
