#Effect of deworm treatment on helminth presence and helminth infection intensity
#Code developed by Jasmine Veitch
#Testing on FULL 2021 FEC dataset (1031 entries)

#This code is associated with two .csv files from Jasmine (sent 11.09.23)
  # "FEC_master_data_02.01.23.csv"
  # "Helminth Data.csv"

#load packages
library(here)
library(tidyverse)
library(lme4)
# library(nlme)
library(lmerTest) #like lme4 but with pvalues

#clear environment
rm(list = ls())

## CODE DEVELOPED BY JASMINE S.M. VEITCH ##

#meta data for FEC data
moremetadata <- readRDS(here("fulltrap21_05.10.23.rds")) %>% ##breeders and nonbreeders are here
  group_by(samp_id) %>% slice(1) %>%
  mutate(tag = as.numeric(tag)) %>%
  dplyr::select(c("site", "trt", "food_trt", "helm_trt",
                  "month", "occasion",
                  "tag", "samp_id", "head", "sex", "current_breeder", "season_breeder"))

#FEC data for 2021
all2021FEC <- read_csv(here("Helminth Data.csv")) %>%
  select("samp_id", "nematode.y.n", "nematode.epg") %>%
  # mutate(site = as.factor(site),
  #        food_trt = as.factor(food_trt),
  #        helm_trt = as.factor(helm_trt),
  #        occasion = as.numeric(case_when((month == "may") ~ 1,
  #                                        (month == "june") ~ 2,
  #                                        (month == "july") ~ 3,
  #                                        (month == "aug") ~ 4,
  #                                        (month == "sept") ~ 5,
  #                                        (month == "oct") ~ 6,
  #                                        TRUE ~ NA))) %>%
  #remove three duplicate entries (details to follow)
  filter(samp_id != 274) %>% #vole 226170 sampled twice in July, remove second sample ID (eggs detected in both samples)
  filter(samp_id != 340) %>% #vole 219980 sampled twice in Aug, remove first sample ID (no eggs in first sample, eggs in second)
  distinct(samp_id, .keep_all = TRUE) %>% #remove duplicate entry for sample id 508
  left_join(moremetadata, by=c("samp_id")) %>%
  group_by(tag) %>%
  arrange(occasion, .by_group = TRUE) %>%
  mutate(capt_nbr = row_number()) %>% #add column for the capture number (ie the first, second, third capture of given vole)
  mutate(pre_post = as.factor(case_when(capt_nbr == 1 ~ "pre",
                                         capt_nbr > 1 ~ "post"))) %>%
  ungroup()



#meta data for FEC data
newmeta <- readRDS(here("fulltrap21_05.10.23.rds")) %>% ##breeders and nonbreeders are here
  group_by(samp_id) %>% slice(1) %>%
  dplyr::select(c("site", "trt", "food_trt", "helm_trt",
                  "season", "month", "occasion",
                  "tag", "samp_id", "sex", "current_breeder", "season_breeder",
                  "head"))

#better FEC data? from Sarah Dec 19 2023
FECdata <- read_csv(here("FEC data for Janine 12-19-23.csv"))

newFEC <- FECdata %>%
  filter("exclude from FEC analysis" != 1) %>% #remove samples with ?weird tube mass?
  rename("samp_id" = "Sample ID",
         "nematode.epg" = "nematode epg") %>%
  select(c("samp_id", "nematode.epg", "nematode.y.n")) %>%
  inner_join(newmeta, by="samp_id") %>%
  drop_na(nematode.epg) %>% drop_na(nematode.y.n) %>%
  mutate(nematode.epg = as.numeric(nematode.epg),
       nematode.y.n = as.factor(nematode.y.n)) %>%
  #remove three duplicate entries (details to follow)
  filter(samp_id != 274) %>% #vole 226170 sampled twice in July, remove second sample ID (eggs detected in both samples)
  filter(samp_id != 340) %>% #vole 219980 sampled twice in Aug, remove first sample ID (no eggs in first sample, eggs in second)
  filter(samp_id != 712) %>% #vole 219809 was sampled twice in Sept, remove first sample ID (eggs detected in both samples)
  filter(samp_id != 64) %>% #remove negative FEC value?
  group_by(tag) %>% arrange(occasion, .by_group = TRUE) %>%
  mutate(capt_nbr = row_number()) %>% #add column for the capture number (ie the first, second, third capture of given vole)
  mutate(pre_post = as.factor(case_when(capt_nbr == 1 ~ "pre",
                                        capt_nbr > 1 ~ "post"))) %>%
  drop_na(sex) %>% drop_na(season_breeder)

  # mutate(diff = (nematode.epg - helig.epg)) %>% #yes, helig and nematode are different
  # filter(diff > 0)


jasmine <- all2021FEC %>% select(tag, samp_id, month, nematode.epg, nematode.y.n)
sarah <-  newFEC %>% select(tag, samp_id, month, nematode.epg, nematode.y.n) %>%
  mutate(tag = as.numeric(tag))

check <- full_join(sarah, jasmine, by=c("tag","samp_id","month")) %>% arrange(tag, month, samp_id) %>%
  mutate(diff = nematode.epg.x - nematode.epg.y) %>%
  filter(diff >0)

############################################################################################



## what we're doing:
# P - marginal
# S - not significant
# find the best P model (wrt sex, head, occasion, etc.) and make I model match
# (deworm:prepost is the focal predictor)
#report the presence:
#p = presence
#adjust the optimizer: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
p.prepost <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
                   data = all2021FEC, family = binomial,
                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) #n=1020 (would be 1031 but missing sex, repro data for some)

NEWp.prepost <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
                         data = newFEC, family = binomial,
                         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(p.prepost)
summary(NEWp.prepost)
#and the intensity (only count >0)
## run it on only counts with eggs (intensity)
i.prepost <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
                  data = subset(all2021FEC, nematode.epg > 0)) #n=433 (would be 437 but missing sex, repro data for some)

NEWi.prepost <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
                            data = subset(newFEC, nematode.epg > 0))

summary(i.prepost)
summary(NEWi.prepost)


#summary tables

#summary table of mixed model output
#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
library(sjPlot)
#save it https://stackoverflow.com/questions/67280933/how-to-save-the-output-of-tab-model
tab_model(p.prepost, file="p.prepost_summary.doc")
tab_model(i.prepost, file="i.prepost_summary.doc")



## YEAh, whatever - ignore this mostly, just pick a good model for P based on biology and make the I model the same
# ## HOWEVER, based on AIC - the model with all the predictors is not the best fit
# ## a rather non-scientific (ie no model selection, not a full AIC comparison)
# dat <- all2021FEC %>% drop_na(sex) %>% drop_na(season_breeder) %>% drop_na(head)
# im1 <- lmer(log(nematode.epg+1) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
#                   data = subset(dat, nematode.epg > 0)) #n=437
# im2 <- lmer(log(nematode.epg+1) ~ helm_trt + pre_post + helm_trt:pre_post + sex + occasion + (1|tag),
#                   data = subset(dat, nematode.epg > 0)) #n=437
# im3 <- lmer(log(nematode.epg+1) ~ helm_trt + pre_post + helm_trt:pre_post + sex + occasion + season_breeder + (1|tag),
#                   data = subset(dat, nematode.epg > 0)) #n=437
# im4 <- lmer(log(nematode.epg+1) ~ helm_trt + pre_post + helm_trt:pre_post + sex + occasion + scale(head) + (1|tag),
#             data = subset(dat, nematode.epg > 0)) #n=437
# AIC(im1, im2, im3, im4) #m1 and m2 are best and equivalent
#
# pm1 <- glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
#                    data = dat, family = binomial) #n=1031
# pm2 <- glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + (1|tag),
#                    data = dat, family = binomial) #n=1031
# pm3 <- glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
#                    data = dat, family = binomial) #n=1031
# pm4 <- glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + scale(head) + (1|tag),
#              data = dat, family = binomial) #n=1031
# pm5 <- glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + scale(head) + (1|tag),
#              data = dat, family = binomial) #n=1031
# AIC(pm1, pm2, pm3, pm4, pm5) #m2 and m3 are best and equivalent
#
# ## so based on this, the best fit models for presence vs intensity differ,
#   #best for P=occasion + sex + scale(head)
#   #best for I=[base model]   OR   I=sex + occasion
#
# ## Essentially - the best pvalue for P is sex+occasion+season_breeder but
#                 #the best fit model is sex+occasion+scale(head)
# ## while the best fit for I is the base model or sex+occasion and
#         #the best pvalue is base model (worst pvalue is sex+occasion+scale(head))



## more new stuff Jan 9 2024 - visualize the interaction effects

library(visreg) #https://pbreheny.github.io/visreg/articles/web/overlay.html
library(cowplot)
library(ggtext)

# #plot looks more like boxplots if pre_post is categorical
# visreg(p.prepost, "pre_post", by="helm_trt", overlay=TRUE,
#        xlab="Treatment Stage", ylab="Log Odds (Infected)")

# #change pre_post to numeric to get a regression line
newFEC.num <- newFEC %>% mutate(pre_post = as.numeric(pre_post))
p.prepost.num <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
                           data = newFEC.num, family = binomial,
                           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
# visreg(p.prepost.num, "pre_post", by="helm_trt", overlay=TRUE,
#        xlab="Treatment Stage", ylab="Log Odds (Infected)")

#the same (continuous pre_post), but prettier
visreg(p.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
                gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="Likelihood of Helminth Infection") +
  annotate(geom = "text", x=2, y=.81, size = 6,
           label = paste("p =",
                         round( coef(summary(p.prepost.num))[7,4], digits=3) )) +
  annotate(geom = "text", x=2, y=.87, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(p.prepost.num))[7,1]), digits=3) )) +
  geom_jitter(data=newFEC.num, aes(x=pre_post, y=nematode.y.n, color=helm_trt),
              width=0.1, height=0.05,
              size=3, alpha=0.2, shape=16)


# #again, looks like a boxplot if pre_post is categorical
# visreg(i.prepost.lmer, "pre_post", by="helm_trt", overlay=TRUE,
#        xlab="Treatment Stage", ylab="Infection Intensity ( ln(EPG) )")

#pre_post as numeric to look like a regression line
#visreg doesn't like it when the model has a subset command in it #https://github.com/pbreheny/visreg/issues/99
witheggs2021FEC.num <- subset(newFEC.num, nematode.epg > 0)
i.prepost.num <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
                  data = witheggs2021FEC.num)
# visreg(i.prepost.num, "pre_post", by="helm_trt", overlay=TRUE,
#        xlab="Treatment Stage", ylab="Infection Intensity ( ln(EPG) )", gg=TRUE, scale="response")

#the same (continuous pre_post), but prettier
visreg(i.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="Infection Intensity ( ln(EPG) )") +
  annotate(geom = "text", x=1.9, y=4.3, size = 6,
           label = paste("p =",
                         round( coef(summary(i.prepost.num))[7,5], digits=3) )) +
  annotate(geom = "text", x=1.9, y=4.4, size = 6,
           label = paste("Est. =",
                         round( coef(summary(i.prepost.num))[7,1], digits=3) ))
#the scale of the actual datapoints goes up to 5000 epg (or 2000 if you remove one point) and messes with the plot
  # geom_jitter(data=witheggs2021FEC.num, aes(x=pre_post, y=nematode.epg, color=helm_trt),
  #             width=0.1, height=0.05,
  #             size=3, alpha=0.2, shape=16)









### NEW! 12 Dec 2023 (updated Dec 15 with Sarah Budischak)
################################################################################################
############## here we're running two models, prevalence and intensity #########################
########### models have both pre- and post- measurements, interaction term #####################
######## as requested by JAE AE for volespatial v3 *finalplzgod* resubmission ##################
################################################################################################

#this code runs on the 'all2021FEC' dataset (models combine pre- and post- data with an interaction by helm_trt) ...
    #using 'helm_trt' and 'pre_post' (both FACTORS) as predictors

############ HELMINTH PRESENCE (p=presence) ##################

# #predictors: only treatment and prepost
# p.prepost <- glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
#                    data = all2021FEC, family = binomial) #n=1031
# summary(p.prepost) #b=-0.5153, p=0.0648.
# #add occasion and sex (better)
# p.prepost <- glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + (1|tag),
#                    data = all2021FEC, family = binomial) #n=1031
# summary(p.prepost) #b=-0.53158, p=0.05950.
# #add scale(head, center=TRUE, scale=TRUE) (WORSE)
# #add current_breeder (WORSE)

#add occasion, sex, and season_breeder (best)
p.prepost <- glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
                   data = all2021FEC, family = binomial,
                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) #n=1031
summary(p.prepost) #b=-0.55645, p=0.05048.

############# presence GLMM model diagnostics ############
#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = p.prepost)
plot(simulationOutput) #qq plot and residual vs fitted
testDispersion(simulationOutput) #formal test for overdispersion
testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)
plotResiduals(simulationOutput, all2021FEC$deworm)
plotResiduals(simulationOutput, all2021FEC$pre_post)
#### OVERALL: DHARMa looks good #####



########## HELMINTH INTENSITY (i=intensity) #############

library(lmerTest) #package is basically lme4 but provides p-values in the model summary

## okay, so going back to helminth intensity (take out the count=0 animals)
## run model on only counts with eggs (HELMINTH INTENSITY)
i.prepost <- lmer(log(nematode.epg+1) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
                  data = subset(all2021FEC, nematode.epg > 0)) #n=437
summary(i.prepost)
#b=-0.30372 p=0.2379
#as you add sex, occasion, season_breeder, pvalue for deworm*prepost gets larger :/
#deworm*prepost p=0.198 with only deworm, prepost params


    # ## @future_Janine FYSI: if you use lme4::lmer(), you don't get a pvalue
    # ## (this is why using lmerTest::lmer() is better)
    # #i = intensity
    # i.prepost <- lme4::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + (1|tag),
    #                         data = subset(all2021FEC, nematode.epg > 0)) #n=437
    # # i.prepost <- lme4::lmer(log(nematode.epg) ~ helm_trt*pre_post + (1|tag), data = subset(all2021FEC, nematode.epg > 0)) #equivalent model
    # summary(i.prepost) #b=-0.32790 (but no pvalue provided)
    #
    # #to get a pvalue (says Ben Bolker:
    # ##https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode)
    # library(nlme)
    # i.prepost2 <- lme(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex, random=~1|tag,
    #                   data=subset(all2021FEC, nematode.epg > 0))
    # # i.prepost2 <- lme(log(nematode.epg) ~ helm_trt*pre_post, random=~1|tag, data=subset(all2021FEC, nematode.epg > 0))
    # anova(i.prepost2) #p=0.2225
    # #tbh this all is very cumbersome and annoying so use lmerTest and it's much easier

#model diagnostics
plot(i.prepost)
shapiro.test(resid(i.prepost)) #this is *more* reasonable

######## intensity LMM model diagnostics ############
#LMM model diagnostics > https://goodekat.github.io/redres/
# devtools::install_github("goodekat/redres")
library(redres)

# creates a plot of the conditional studentized residuals versus the fitted values
plot_redres(i.prepost, type="std_cond")
plot_redres(i.prepost, type="pearson_cond") #Pearson conditional residuals
# creates a residual quantile plot for the error term
plot_resqq(i.prepost) #dece
# creates normal quantile plots for each random effect
plot_ranef(i.prepost) #also dece
###### OVERALL: meh, plots seem alright


# ## trying something else... (developed with S. Budischak)
# ## run model on all counts (including 0's) - effectively, HELMINTH ABUNDANCE
# i.prepost <- lmer(log(nematode.epg+1) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
#                   data = all2021FEC) #n=1031
# summary(i.prepost)
# #model diagnostics
# plot(i.prepost)
# shapiro.test(resid(i.prepost)) #p-value is VERY small - violating assumption of normality
# #more diagnostics
# library(redres)
# # creates a residual quantile plot for the error term
# plot_resqq(i.prepost) #YIKES
# # creates normal quantile plots for each random effect
# plot_ranef(i.prepost) #YIKES
# ##RESULT: nope, helminth abundance is not the solution


############################################# END the good stuff #################################################





###########################################################################################
###### here is another thing that didn't really work, developed with S Budischak ##########
###########################################################################################
#### tried running models with worm abundance (combo of presence and intensity)
#### P * I = A
#### and with Poisson and NegBinomial error families (to deal with non-normality of data)

## Result: Significant effect of deworm*prepost in the Poisson model
        ## no effect of deworm*prepost in negative binomial models
        ## and the neg bin models have much better fit than the Poisson
## so not a great option because it doesn't prove what we want :(

library(glmmTMB)

#create a parameter for worm abundance, needs to be integer, thus the *10
#glmmTMB also requires the predictor (ie tag) to be a factor
all2021FEC <- newFEC %>%
  mutate(worms = round( (log(nematode.epg+1))*10, 0 )) %>%
  mutate(tagfac = as.factor(tag))

hist(all2021FEC$worms) #fairly Poisson

#worms is an ABUNDANCE variable now - "at the population scale, how infected is each population?"
#Prev*Intensity = Abundance
  #FEWER infected animals AND infected animals have LOWER intensity
  #and this happens AFTER deworming
mpois <- glmmTMB(worms ~ helm_trt + pre_post + occasion + sex + season_breeder + helm_trt:pre_post + (1|tag),
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                         data=all2021FEC, ziformula=~1, family=poisson)

mnbin1 <- glmmTMB(worms ~ helm_trt + pre_post + occasion + sex + season_breeder + helm_trt:pre_post + (1|tag),
                  control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                         data=all2021FEC, ziformula=~1, family=nbinom1)

mnbin2 <- glmmTMB(worms ~ helm_trt + pre_post + occasion + sex + season_breeder + helm_trt:pre_post + (1|tag),
                         data=all2021FEC, ziformula=~1, family=nbinom2)


summary(mpois) #signif effect of occasion and deworm*prepost
summary(mnbin1) #no signif effects
summary(mnbin2) #no signif effects
#strange that there's no effect of deworming or occasion in these models.. l
  #makes Sarah wonder if they're actually fitting properly

#save it
tab_model(mpois, file="abundance_model_poisson.doc")


#model diagnostics
plot(resid(mpois, type="pearson") ~ fitted(mpois))
plot(resid(mnbin1, type="pearson") ~ fitted(mnbin1)) #looks a little better than mpois
plot(resid(mnbin2, type="pearson") ~ fitted(mnbin2)) #looks a little better than mpois

qqnorm(resid(mpois, type="pearson"))
qqnorm(resid(mnbin1, type="pearson")) #looks a little better than mpois
qqnorm(resid(mnbin2, type="pearson")) #looks a little better than mpois

library(lmtest) #basically AIC for non normal distribution models
lrtest(mpois, mnbin1, mnbin2) #still shows that negbin models are better than poisson




############# presence GLMM model diagnostics ############
#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = mpois)
simulationOutput <- simulateResiduals(fittedModel = mnbin1)
simulationOutput <- simulateResiduals(fittedModel = mnbin2)
plot(simulationOutput) #qq plot and residual vs fitted
testDispersion(simulationOutput) #formal test for overdispersion
testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)
plotResiduals(simulationOutput, all2021FEC$helm_trt)
plotResiduals(simulationOutput, all2021FEC$pre_post)
plotResiduals(simulationOutput, all2021FEC$sex)
plotResiduals(simulationOutput, all2021FEC$season_breeder)
plotResiduals(simulationOutput, all2021FEC$occasion)
#### OVERALL: DHARMa looks good #####






###########################################################################################





# #############################################################################################
# ########### everything here down runs four separate models using all2021FEC #################
# ############## (this is what I used for the volespatial v2 resubmission) ####################
# #############################################################################################
#
# #Subset all2021FEC into first and subsequent captures (for modeling)
# firstcap <- all2021FEC %>% arrange(tag, occasion) %>%
#   distinct(tag, .keep_all = TRUE) #697 entries
# subseqcap <- all2021FEC %>% arrange(tag, occasion) %>%
#   filter(capt_nbr > 1) #334 entries
#
# ### This is essentially the full dataset of FEC data from 2021 - I elected to use this data for the deworm efficacy analysis
# ### since it's the most complete. The results are the same if I subset the FEC data for only the voles in the network analysis (code follows)
# ### but that honestly didn't make a ton of sense since we're just interested in whether the deworm analysis worked, not necessarily
# ### how effective it was for the particular animals that have network data.
#
# #First capture models (pre-treatment)
# #p = presence, i = intensity
# #intensity models are subsetted to only captures infected with worms
# p.first <- glm(nematode.y.n ~ helm_trt, data = firstcap, family = binomial) #n=697
# summary(p.first) #b=-0.2869, p=0.0643.
# #marginal significance, weak effect of treatment on helminth presence at first capture, deworm has slightly lower presence
#
# i.first <- lm(log(nematode.epg) ~ helm_trt, data = subset(firstcap, nematode.epg > 0)) #n=287
# summary(i.first) #b=-0.06278, p=0.697
# #nonsignificant effect of treatment on helminth infection intensity at first capture
#
# #Subsequent capture models (post-treatment)
# #p = presence, i = intensity
# p.subseqcap <- glmer(nematode.y.n ~ helm_trt + (1|tag), data = subseqcap, family = binomial) #n=334
# summary(p.subseqcap) #b=-0.8105, p=0.000621***
# #strong, significant effect of treatment on helminth presence at subsequent captures, deworm has lower presence
#
# i.subseqcap <- lmer(log(nematode.epg) ~ helm_trt + (1|tag), data = subset(subseqcap, nematode.epg > 0)) #n=150
# summary(i.subseqcap) #b=-0.4153
# #to get a pvalue (says Ben Bolker:
#   ##https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode)
# # library(nlme)
# i.subseqcap2 <- lme(log(nematode.epg) ~ helm_trt, random=~1|tag, data=subset(subseqcap, nematode.epg > 0))
# anova(i.subseqcap2) #p=0.0453*
# #moderate, significant effect of treatment on helminth infection intensity at subsequent captures, deworm has lower helminth epg


###############################still doing the four models thing########################################
# # Also ran the models on just the animals that have both FEC AND netmets data
#     # it's a smaller dataset but shows the same results
#     # this was not included anywhere in the v2 submission, just for Janine's own knowledge
#
# #FEC data from Sarah
# helmdat <- read_csv(here("DEWORM_ANALYSIS", "Helminth Data.csv")) %>%
#   #remove three duplicate entries (details to follow)
#   filter(samp_id != 274) %>% #vole 226170 sampled twice in july, remove second sample ID (eggs detected in both samples)
#   filter(samp_id != 340) %>% #vole 219980 sampled twice in aug, remove first sample ID (no eggs in first sample, eggs in second)
#   distinct(samp_id, .keep_all = TRUE)
#
# #Netmets 2021 data - one entry per vole, per month
# # subset of fulltrap: these are only the animals that were included in the networks
# netmets21 <- readRDS(here("netmets21_STSB.rds")) %>%
#   mutate(tag = as.numeric(tag))
#
# #join FEC data and netmets, keep only entries that have both network and FEC data
# nmFEC <- left_join(netmets21, helmdat, by=(c("site", "month", "tag"))) %>%
#   drop_na(nematode.epg) %>% #drop entries that have netmets data but no FEC
#   mutate(occasion = as.numeric(case_when((month == "may") ~ 1,
#                                   (month == "june") ~ 2,
#                                   (month == "july") ~ 3,
#                                   (month == "aug") ~ 4,
#                                   (month == "sept") ~ 5,
#                                   (month == "oct") ~ 6,
#                                   TRUE ~ NA))) %>%
#   group_by(tag) %>%
#   mutate(capt_nbr = row_number()) %>% #add column for the capture number
#   ungroup()
#
# #Subset into first and subsequent captures
# firstcap <- nmFEC %>% arrange(tag, occasion) %>%
#   distinct(tag, .keep_all = TRUE) #681 entries
# subseqcap <- nmFEC %>% arrange(tag, occasion) %>%
#   filter(capt_nbr > 1) #318 entries
#
# #First capture models (pre-treatment)
# #p = presence, i = intensity
# #intensity models are subsetted to only captures infected with worms
# p.first <- glm(nematode.y.n ~ helm_trt, data = firstcap, family = binomial)
# summary(p.first) #b=-0.2624, p=0.0957
#
# i.first <- lm(log(nematode.epg) ~ helm_trt, data = subset(firstcap, nematode.epg > 0))
# summary(i.first) #b=-0.09286, p=0.578
#
# #Subsequent capture models (post-treatment)
# p.subseqcap <- glmer(nematode.y.n ~ helm_trt + (1|tag), data = subseqcap, family = binomial)
# summary(p.subseqcap) #b=-0.8690, p=0.000288***
#
# i.subseqcap <- lmer(log(nematode.epg) ~ helm_trt + (1|tag), data = subset(subseqcap, nematode.epg > 0))
# summary(i.subseqcap) #b=-0.4966, p=0.0168*
#
# #to get a pvalue (says Ben Bolker:
# ##https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode)
# library(nlme)
# i.subseqcap2 <- lme(log(nematode.epg) ~ helm_trt, random=~1|tag, data=subset(subseqcap, nematode.epg > 0))
# anova(i.subseqcap2)
#########################################################################################

