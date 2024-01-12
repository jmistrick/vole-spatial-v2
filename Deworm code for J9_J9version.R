#Effect of deworm treatment on helminth presence and helminth infection intensity
## Code developed with input from JASMINE S.M. VEITCH & SARAH A. BUDISCHAK ##
#Testing on FULL 2021 FEC dataset (1025 entries)

#This code runs with FEC data from Sarah Budischak "FEC data for Janine 12-19-23.csv"

#load packages
library(here)
library(tidyverse)
library(lme4)
library(lmerTest) #like lme4 but with pvalues

# #clear environment
# rm(list = ls())


######################################### LOAD DATA ###################################################

#metadata for FEC data
metadata <- readRDS(here("fulltrap21_05.10.23.rds")) %>% ##breeders and nonbreeders are here
  group_by(samp_id) %>% slice(1) %>% #remove duplicate entries per occasion (ie week recaps)
  dplyr::select(c("site", "trt", "food_trt", "helm_trt",
                  "season", "month", "occasion",
                  "tag", "samp_id", "sex", "current_breeder", "season_breeder",
                  "head"))

#full (2021-23) FEC data from Sarah v.Dec 19 2023
FECdata <- read_csv(here("FEC data for Janine 12-19-23.csv"))

#pull 2021 FEC data, combine with metadata, keep only complete entries
#helm_trt, pre_post, and nematode.y.n are all factors
all2021FEC <- FECdata %>%
  filter("exclude from FEC analysis" != 1) %>% #remove samples with feces wt less than 0.005
  rename("samp_id" = "Sample ID",
         "nematode.epg" = "nematode epg") %>%
  select(c("samp_id", "nematode.epg", "nematode.y.n")) %>%
  inner_join(metadata, by="samp_id") %>%
  drop_na(nematode.epg) %>% drop_na(nematode.y.n) %>% #remove entries missing FEC data
  mutate(nematode.epg = as.numeric(nematode.epg),
       nematode.y.n = as.factor(nematode.y.n)) %>%
  #remove three duplicate entries (details to follow)
  filter(samp_id != 274) %>% #vole 226170 sampled twice in July, remove second sample ID (eggs detected in both samples)
  filter(samp_id != 340) %>% #vole 219980 sampled twice in Aug, remove first sample ID (no eggs in first sample, eggs in second)
  filter(samp_id != 712) %>% #vole 219809 was sampled twice in Sept, remove first sample ID (eggs detected in both samples)
  filter(samp_id != 64) %>% #remove negative FEC value
  group_by(tag) %>% arrange(occasion, .by_group = TRUE) %>%
  mutate(capt_nbr = row_number()) %>% #add column for the capture number (ie the first, second, third capture of given vole)
  mutate(pre_post = as.factor(case_when(capt_nbr == 1 ~ "pre",
                                        capt_nbr > 1 ~ "post"))) %>%
  # drop_na(sex) %>% drop_na(season_breeder) %>% #model has 10ish voles with sex or repro=NA but that's not a model parameter so it's fine
  ungroup()

################################################################################################
############################ RUN PREVALENCE & INTENSITY MODELS #################################
################################################################################################
############## here we're running two models, prevalence and intensity #########################
########### models have both pre- and post- measurements, interaction term #####################
######## as requested by JAE AE for volespatial v3 *finalplzgod* resubmission ##################
################################################################################################

#this code runs on the 'all2021FEC' dataset (models combine pre- and post- data with an interaction by helm_trt) ...
  #using 'helm_trt' and 'pre_post' (both FACTORS) as predictors

#p = prevalence model
#if needed: adjust the optimizer: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
p.prepost <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                   data = all2021FEC, family = binomial) #n=1035
summary(p.prepost) #trt*prepost b=-0.671 (log odds) p=0.0215*
exp(coef(summary(p.prepost))[,1]) #trt*prepost b=0.511 (odds ratio)

# ###### GLMM model diagnostics ######
# #https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# library(DHARMa)
# #calculate residuals (then run diagnostics on these)
# simulationOutput <- simulateResiduals(fittedModel = p.prepost)
# #qq plot and residual vs fitted
# plot(simulationOutput) #qq look good, scaled residuals look good
# #formal test for overdispersion
# testDispersion(simulationOutput) #looks good, pvalue is large
# #formal test for zero inflation (common type of overdispersion)
# testZeroInflation(simulationOutput) #looks good, pvalue is large
# plotResiduals(simulationOutput, all2021FEC$helm_trt) #looks good
# plotResiduals(simulationOutput, all2021FEC$pre_post) #looks good
# #### OVERALL: DHARMa looks good #####

#------------------------------------------------

#i = intensity (only samples with nematode.epg >0)
i.prepost <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                  data = subset(all2021FEC, nematode.epg > 0)) #n=440
summary(i.prepost) #trt*prepost b=-0.740 p=0.0053**

    # ## @future_Janine FYSI: if you use lme4::lmer(), you don't get a pvalue
    # ## (this is why using lmerTest::lmer() is better)
    # #i = intensity
    # i.prepost <- lme4::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
    #                         data = subset(all2021FEC, nematode.epg > 0))
    # summary(i.prepost) #log-odds=... (but no pvalue provided)
    #
    # #to get a pvalue (says Ben Bolker:
    # ##https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode)
    # library(nlme)
    # i.prepost2 <- lme(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder, random=~1|tag,
    #                   data=subset(all2021FEC, nematode.epg > 0))
    # anova(i.prepost2) #p-value=...
    # #tbh this all is very cumbersome and annoying so use lmerTest and it's much easier

# #model diagnostics
# plot(i.prepost) #okay... but some patterning
# shapiro.test(resid(i.prepost)) #p=0.00985 pretty small, some issues of non-normality
#
# ######## intensity LMM model diagnostics ############
# #LMM model diagnostics > https://goodekat.github.io/redres/
# # devtools::install_github("goodekat/redres")
# library(redres)
# # creates a plot of the conditional studentized residuals versus the fitted values
# plot_redres(i.prepost, type="std_cond")
# plot_redres(i.prepost, type="pearson_cond") #Pearson conditional residuals
# # creates a residual quantile plot for the error term
# plot_resqq(i.prepost) #looks good
# # creates normal quantile plots for each random effect
# plot_ranef(i.prepost) #also looks good
# ###### OVERALL: model diagnostics look good

#------------------------------------------------

#summary table of mixed model output
#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
library(sjPlot)
#save it https://stackoverflow.com/questions/67280933/how-to-save-the-output-of-tab-model
tab_model(p.prepost, file="p.prepost_summary.doc")
tab_model(i.prepost, file="i.prepost_summary.doc")


#################################### VISUALIZE MODEL OUTPUT ############################################

## visualize the interaction effects

library(visreg) #https://pbreheny.github.io/visreg/articles/web/overlay.html
library(cowplot)
library(ggtext)

### PREVALENCE MODEL ###

#change pre_post to numeric to get a regression line in visreg output
all2021FEC.num <- all2021FEC %>%
  mutate(pre_post = as.numeric(pre_post)) %>%
  mutate(nematode.y.n_num = case_when(nematode.y.n == "0" ~ 0,
                                      nematode.y.n == "1" ~ 1)) #nematode.y.n from factor to numeric to plot points on continous y axis
p.prepost.num <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                           data = all2021FEC.num, family = binomial)

#visualize infection probability, pretty
png(filename="prev_model_vis.png", height=4, width=6, units="in", res=600)
visreg(p.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
                gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(breaks=c(1,2),
                     expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=11),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="Likelihood of Helminth Infection") +
  annotate(geom = "text", x=1.35, y=1.05, size = 4,
           label = paste("Trt Grp*Trt Stage")) +
  annotate(geom = "text", x=1.35, y=0.97, size = 4,
           label = paste("OR =",
                         round( exp(coef(summary(p.prepost.num))[4,1]), digits=3) )) +
  annotate(geom = "text", x=1.35, y=.90, size = 4,
           label = paste("p =",
                         round( coef(summary(p.prepost.num))[4,4], digits=3) )) +
  geom_jitter(aes(x=pre_post, y=nematode.y.n_num, color=helm_trt),
              data=all2021FEC.num,
              width=0.1, height=0.05,
              size=3, alpha=0.2, shape=16)
dev.off()

# #ALTERNATIVE, visreg plot looks more like boxplots if pre_post is categorical
# visreg(p.prepost, "pre_post", by="helm_trt", overlay=TRUE,
#        xlab="Treatment Stage", ylab="Log Odds (Infected)")


### INTENSITY MODEL ###

#pre_post as numeric to look like a regression line
#visreg doesn't like it when the model has a subset command in it #https://github.com/pbreheny/visreg/issues/99
witheggs2021FEC.num <- subset(all2021FEC.num, nematode.epg > 0)
i.prepost.num <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                  data = witheggs2021FEC.num)

#visualize infection intensity, pretty
png(filename="int_model_vis.png", height=4, width=6, units="in", res=600)
visreg(i.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(breaks=c(1,2),
                     expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="ln(eggs per gram of feces)") +
  annotate(geom = "text", x=1.35, y=8.6, size = 4,
           label = paste("Trt Grp*Trt Stage")) +
  annotate(geom = "text", x=1.35, y=8.0, size = 4,
           label = paste("\u03B2 =",
                         round( coef(summary(i.prepost.num))[4,1], digits=3) )) +
  annotate(geom = "text", x=1.35, y=7.5, size = 4,
           label = paste("p =",
                         round( coef(summary(i.prepost.num))[4,5], digits=3) )) +
  geom_jitter(data=witheggs2021FEC.num, aes(x=pre_post, y=log(nematode.epg), color=helm_trt),
              width=0.05, size=3, alpha=0.2, shape=16)
dev.off()

# #ALTERNATIVE, visreg plot looks like a boxplot if pre_post is categorical
# witheggs2021FEC <- subset(all2021FEC, nematode.epg > 0)
# i.prepost.visreg <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
#                                 data = witheggs2021FEC)
# visreg(i.prepost.visreg, "pre_post", by="helm_trt", overlay=TRUE,
#        xlab="Treatment Stage", ylab="Infection Intensity ( ln(EPG) )")




###########################################################################################
########################### SINGLE MODEL, HELMINTH ABUNDANCE ##############################
###########################################################################################

### NOT using this model for the v3.R1 manuscript submission, just for shiggles

# ## run single model on all counts (including 0's) - ie, HELMINTH ABUNDANCE (developed with S. Budischak)
# a.prepost <- lmer(log(nematode.epg+1) ~ helm_trt + pre_post + helm_trt:pre_post + occasion + sex + season_breeder + (1|tag),
#                   data = all2021FEC) #n=1035
# summary(a.prepost)
#
# #model diagnostics
# plot(a.prepost) #cloud looks good (zeros in line)
# shapiro.test(resid(a.prepost)) #p-value is VERY small - violating assumption of normality
# #more diagnostics
# # library(redres)
# # creates a residual quantile plot for the error term
# plot_resqq(a.prepost) #YIKES
# # creates normal quantile plots for each random effect
# plot_ranef(a.prepost) #YIKES
# ##RESULT: helminth abundance can't be run as an lmer


#### models with worm abundance (combo of presence and intensity, P * I = A)
#### with Poisson and NegBinomial error families (to deal with non-normality of data)

## Result: NegBin model is better fit than Poisson, significant interaction term (trt*prepost)

library(glmmTMB) #https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf

#create a parameter for worm abundance, needs to be integer: thus the *10, round to 0 decimal
#glmmTMB also requires the predictor (ie tag) to be a factor
all2021FEC.abund <- all2021FEC %>%
  mutate(worms = round( (log(nematode.epg+1))*10, 0 )) %>%
  mutate(tagfac = as.factor(tag))

# hist(all2021FEC.abund$worms) #definitely not normal, zero-inflated Poisson/NegBin for sure

#worms is an ABUNDANCE variable now - "at the population scale, how infected is each population?"
#Prevalence * Intensity = Abundance
  #FEWER infected animals AND infected animals have LOWER intensity
  #and this happens AFTER deworming
a.pois <- glmmTMB(worms ~ helm_trt + pre_post + occasion + sex + season_breeder + helm_trt:pre_post + (1|tag),
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                         data=all2021FEC.abund, ziformula=~1, family=poisson)

a.nbin1 <- glmmTMB(worms ~ helm_trt + pre_post + occasion + sex + season_breeder + helm_trt:pre_post + (1|tag),
                  control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                         data=all2021FEC.abund, ziformula=~1, family=nbinom1)

a.nbin2 <- glmmTMB(worms ~ helm_trt + pre_post + occasion + sex + season_breeder + helm_trt:pre_post + (1|tag),
                  control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                         data=all2021FEC.abund, ziformula=~1, family=nbinom2)


summary(a.pois) #signif effect of occasion and deworm*prepost
summary(a.nbin1) #signif effect of occasion and deworm*prepost
summary(a.nbin2) #signif effect of occasion and deworm*prepost

#model diagnostics
plot(resid(a.pois, type="pearson") ~ fitted(a.pois)) #clear pattern to residuals
plot(resid(a.nbin1, type="pearson") ~ fitted(a.nbin1)) #two clouds, slight pattern
plot(resid(a.nbin2, type="pearson") ~ fitted(a.nbin2)) #two clouds, but no pattern

qqnorm(resid(a.pois, type="pearson")) #pretty curvy
qqnorm(resid(a.nbin1, type="pearson")) #looks a little better than a.pois
qqnorm(resid(a.nbin2, type="pearson")) #looks a little better than a.pois

library(lmtest) #model comparison for non normal distribution models
lrtest(a.pois, a.nbin1, a.nbin2) #NegBin models are better fit than Poisson, negbin models are identical


############# abundance model diagnostics ############
#DHARMa works for glmmTMB too > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
# simulationOutput <- simulateResiduals(fittedModel = a.pois)
# simulationOutput <- simulateResiduals(fittedModel = a.nbin1)
simulationOutput <- simulateResiduals(fittedModel = a.nbin2)

#qq plot and residual vs fitted
plot(simulationOutput)
#formal test for overdispersion
testDispersion(simulationOutput)
#formal test for zero inflation (common type of overdispersion)
testZeroInflation(simulationOutput)
#plots of residuals for individual predictors
plotResiduals(simulationOutput, all2021FEC$helm_trt)
plotResiduals(simulationOutput, all2021FEC$pre_post)
plotResiduals(simulationOutput, all2021FEC$sex)
plotResiduals(simulationOutput, all2021FEC$season_breeder)
plotResiduals(simulationOutput, all2021FEC$occasion)
#### OVERALL: DHARMa looks ...meh for NegBin #####


#save model summary ouput
tab_model(a.pois, file="abundance_model_poisson.doc")
tab_model(a.nbin1, file="abundance_model_negbin1.doc")
tab_model(a.nbin2, file="abundance_model_negbin2.doc")

############################################## END #############################################
