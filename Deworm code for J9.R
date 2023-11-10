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

#clear environment
rm(list = ls())

## CODE DEVELOPED BY JASMINE S.M. VEITCH ##

#FEC data for 2021
all2021FEC <- read_csv(here("DEWORM_ANALYSIS", "Helminth Data.csv")) %>%
  mutate(site = as.factor(site),
         food_trt = as.factor(food_trt),
         helm_trt = as.factor(helm_trt),
         occasion = as.numeric(case_when((month == "may") ~ 1,
                                         (month == "june") ~ 2,
                                         (month == "july") ~ 3,
                                         (month == "aug") ~ 4,
                                         (month == "sept") ~ 5,
                                         (month == "oct") ~ 6,
                                         TRUE ~ NA))) %>%
  #remove three duplicate entries (details to follow)
  filter(samp_id != 274) %>% #vole 226170 sampled twice in July, remove second sample ID (eggs detected in both samples)
  filter(samp_id != 340) %>% #vole 219980 sampled twice in Aug, remove first sample ID (no eggs in first sample, eggs in second)
  distinct(samp_id, .keep_all = TRUE) %>% #remove duplicate entry for sample id 508
  group_by(tag) %>%
  mutate(capt_nbr = row_number()) %>% #add column for the capture number (ie the first, second, third capture of given vole)
  ungroup()

#Subset all2021FEC into first and subsequent captures (for modeling)
firstcap <- all2021FEC %>% arrange(tag, occasion) %>%
  distinct(tag, .keep_all = TRUE) #697 entries
subseqcap <- all2021FEC %>% arrange(tag, occasion) %>%
  filter(capt_nbr > 1) #334 entries

### This is essentially the full dataset of FEC data from 2021 - I elected to use this data for the deworm efficacy analysis
### since it's the most complete. The results are the same if I subset the FEC data for only the voles in the network analysis (code follows)
### but that honestly didn't make a ton of sense since we're just interested in whether the deworm analysis worked, not necessarily
### how effective it was for the particular animals that have network data.

#First capture models (pre-treatment)
#p = presence, i = intensity
#intensity models are subsetted to only captures infected with worms
p.first <- glm(nematode.y.n ~ helm_trt, data = firstcap, family = binomial)
summary(p.first) #b=-0.2869, p=0.0643.
#marginal significance, weak effect of treatment on helminth presence at first capture, deworm has slightly lower presence

i.first <- lm(log(nematode.epg) ~ helm_trt, data = subset(firstcap, nematode.epg > 0))
summary(i.first) #b=-0.06278, p=0.697
#nonsignificant effect of treatment on helminth infection intensity at first capture

#Subsequent capture models (post-treatment)
#p = presence, i = intensity
p.subseqcap <- glmer(nematode.y.n ~ helm_trt + (1|tag), data = subseqcap, family = binomial)
summary(p.subseqcap) #b=-0.8105, p=0.000621***
#strong, significant effect of treatment on helminth presence at subsequent captures, deworm has lower presence

i.subseqcap <- lmer(log(nematode.epg) ~ helm_trt + (1|tag), data = subset(subseqcap, nematode.epg > 0))
summary(i.subseqcap) #b=-0.4153
#to get a pvalue (says Ben Bolker:
  ##https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode)
library(nlme)
i.subseqcap2 <- lme(log(nematode.epg) ~ helm_trt, random=~1|tag, data=subset(subseqcap, nematode.epg > 0))
anova(i.subseqcap2) #p=0.0453*
#moderate, significant effect of treatment on helminth infection intensity at subsequent captures, deworm has lower helminth epg



#########################################################################################
# # Also ran the models on just the animals that have both FEC AND netmets data
#     # it's a smaller dataset but shows the same results
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

