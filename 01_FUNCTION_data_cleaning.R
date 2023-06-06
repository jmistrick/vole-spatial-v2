
########### TO RUN THIS CODE IN THIS FILE ################
# #clean 2021 data
# processingdata = "vole_capture_data_05.04.23.csv"
# WRdata= "week_recap_data_05.04.23.csv"
# yr=2021
# fulltrap_output = "fulltrap21_05.10.23.rds"
#
# #clean 2022 data
# processingdata = "vole_capture_data_05.04.23.csv"
# WRdata= "week_recap_data_05.04.23.csv"
# yr=2022
# fulltrap_output = "fulltrap22_05.10.23.rds"
#########################################################


## This function cleans and combines the vole processing and the week recap data (input as .csv files)
## Inputs: processing data = "filename.csv" - (fed into here()) THIS SHOULD BE THE FULL 2021, 2022 CLEANED DATA
##          WRdata = "filename.csv" - same as above
##          year = number (2021 or 2022) to get the year you want
##          fulltrap_output = "filename.rds" - where to save the final 'fulltrap' file
## Output: fulltrap file with all the capture/WR data combined, plus new stuff calculated for desired year

clean_data <- function(processingdata, WRdata, yr, fulltrap_output) {

  ###############################   ENTRY & CLEANING VOLE CAPTURE DATA   ##################################

  #load data
  voledata <- read.csv(here(processingdata))
  #june 3 2022 version - has corrected tag IDs and ambiguous sexes resolved (2021 season)
  #oct 3 2022 version - has weirdly large head for 226016 removed (2021 season)
  #feb 28 2023 version - has cleaned/updated PIT tags, vole sexes (2022 season)

  #clean names
  voledata <- voledata %>%
    clean_names %>%
    #rename column
    rename(session = sess, samp_id = id) %>%
    #mutate columns to their proper formats
    mutate(year = as.numeric(year),
           occasion = as.numeric(occasion), #as.factor(occasion),
           site = as.character(site),
           food_trt = as.factor(food_trt),
           helm_trt = as.factor(helm_trt),
           date = as_date(date, format= "%m/%d/%Y"),
           session = as.numeric(session),
           trap = as.character(trap),
           tag = as.character(tag),
           samp_id = as.numeric(samp_id), #as.factor(samp_id),
           sex = as.factor(sex),
           ow = as.factor(ow),
           per = as.factor(per),
           nip = as.factor(nip),
           preg = as.factor(preg),
           test = as.factor(test),
           head = as.numeric(head),
           mass = as.numeric(mass),
           ticks = as.numeric(ticks),
           ear_ed = as.factor(ear_ed),
           saliva_sr = as.factor(saliva_sr),
           smear_bs = as.factor(smear_bs),
           bloodspin_bc = as.factor(bloodspin_bc),
           bloodrna_br = as.factor(bloodrna_br),
           fecalegg_fe = as.factor(fecalegg_fe),
           fecalrna_fr = as.factor(fecalrna_fr),
           deworm = as.factor(deworm),
           new = as.factor(new),
           fate = as.factor(fate),
           handler = as.factor(handler),
           notes = as.character(notes)) %>%
    #filter for just one year
    filter(year==yr) %>%
    #remove uusi 2022 data
    filter(site!="uusi") %>%
    #remove any animals without a trap location
    filter(!is.na(trap)) %>%
    #remove any animals without a tag id
    filter(!is.na(tag)) %>%
    #remove animals found dead when setting/supplementing
    filter(!session == "0") %>%
    filter(!is.na(session))

  ########## THE ABOVE CODE removes all animals that have missing data or were found dead NOT during a trapping occasion ##################
  # voledata %>% filter(is.na(trap))
  # #219825 - euthanized 9.1.21 during field course (caught off the grid)
  # #219917 - euthanized 9.1.21 during field course (caught off the grid)
  # voledata %>% filter(is.na(tag)) #these are fine to remove, nothing important here
  # voledata %>% filter(session=="0")
  # #these are animals found dead when setting - but I don't know when they went into the trap (assumed right after we left)
  # #they aren't technically session 4 captures, though so I don't know what I'd do with them if I left them in
  # voledata %>% filter(is.na(session)) #these animals were found DT when supplementing

  ########## RIGHT NOW, all of the DP, DT, and S animals (except those that were part of ^ above) are still in the dataset ################
  # #remove animals euthanized for terminal sampling
  # filter(!fate == "S")
  # #remove animals DP or DT
  # filter(!fate == "DP") %>%
  # filter(!fate == "DT")

  #check for spelling errors, extra groups, weird data
  # df$column[df$column == "old_value"] <- new_value #find and replace
  # unique(voledata$site)
  # unique(voledata$year)
  # unique(voledata$occasion)
  # unique(voledata$session)
  # unique(voledata$food_trt)
  # unique(voledata$helm_trt)
  # unique(voledata$sex)
  # unique(voledata$fate)
  # unique(voledata$ow)
  # unique(voledata$per)
  # unique(voledata$nip)
  # unique(voledata$preg)
  # unique(voledata$test)
  # range(voledata$head, na.rm = TRUE)
  # range(voledata$mass, na.rm = TRUE)
  # range(voledata$ticks, na.rm = TRUE)


  ################################ create a time column to replace session (for CMRnet) ##################################
  voledata <-
    voledata %>%
    mutate(time = case_when(
      session == "1" ~ "06:00:00",
      session == "2" ~ "18:00:00",
      session == "3" ~ "06:00:00",
      session == "4" ~ "18:00:00",)) %>%
    #so we didn't actually check traps 12hrs apart (more like 6am, 4pm) but this is cleaner for the distance moved analysis
    relocate(time, .after = session)

  #turn the time into a lubridate time
  voledata$time <- hms(voledata$time)
  #combine date and time columns into a lubridate time
  voledata$date_time <- ymd_hms(paste(voledata$date, voledata$time))
  #move date and time around
  voledata <-
    voledata %>%
    relocate(date, .after = year) %>%
    relocate(session, .after = occasion) %>%
    relocate(date_time, .after = date) %>%
    dplyr::select(-time)
  ######################################################### END ############################################################

  ################################ convert trap number to a grid coordinate ####################################
  voledata<-
    voledata %>%
    #separate trap ID (letter/number) into column of the letter (x) and number (y)
    separate(trap, into = c("x", "y"), sep = "(?<=[A-Z])(?=[0-9])", remove=FALSE) %>%
    #recode each letter as a number
    #this makes a new column with the letter turned into its corresponding number
    mutate(new_x = case_when(
      x == "A" ~ 1,
      x == "B" ~ 2,
      x == "C" ~ 3,
      x == "D" ~ 4,
      x == "E" ~ 5,
      x == "F" ~ 6,
      x == "G" ~ 7,
      x == "H" ~ 8,
      x == "I" ~ 9,
      x == "J" ~ 10,
      x == "K" ~ 11)) %>%
    #move that new column right after the original letter column
    relocate(new_x, .after = x) %>%
    #remove that original letter column
    dplyr::select(-x) %>%
    #rename the new number x column
    rename(x = new_x) %>%
    #make sure y and x are both numeric, not characters
    mutate(x = as.numeric(x), y = as.numeric(y))


  voledata<-
    voledata %>%
    #adjust the trap number (y) so it corresponds to a grid number
    #if remainder of x/2 is 0 (if x is even), multiply y by 2 - if else multiply by 2 then subtract 1
    mutate(y_new = ifelse((x %% 2) == 0, ((voledata$y)*2), (((voledata$y)*2)-1))) %>%
    relocate(y_new, .after = y) %>%
    #remove the old y column and rename the new one
    dplyr::select(-y) %>%
    rename(y = y_new)
  ######################################################### END ##################################################################

  #remove un-needed columns from datatable - make a smaller version for easier analysis
  trap <- voledata %>%
    dplyr::select(-id_as_number, -date, -ow, -ticks, -ear_ed, -saliva_sr, -smear_bs,
                  -bloodspin_bc, -bloodrna_br, -fecalegg_fe, -fecalrna_fr, -deworm, -new, -handler, -notes)

  ##### SAMPLE_ID and YEAR are back in the data (May 2023)


  #############################  LOAD AND CLEAN WEEK RECAP DATA   ###############################################

  #load the data
  wr_data <- read.csv(here(WRdata))
  #june 3 2022 version is the most up-to-date, occ 5 and 6 have been checked and sexes corrected (2021 season)
  #oct 3 2022 version only had updates to the processing data, wr data remained unchanged from june 28 version (2021 season)
  #feb 28 2023 version has cleaned/updated PIT tags/vole sexes (2022 season)

  #several columns (sex, per, nip, preg, test, fate, handler) have "" or "not noted" --> change these to NA
  wr_data$sex[wr_data$sex == "not noted"] <- NA
  wr_data$sex[wr_data$sex == ""] <- NA
  wr_data$per[wr_data$per == "not noted"] <- NA
  wr_data$per[wr_data$per == ""] <- NA
  wr_data$nip[wr_data$nip == "not noted"] <- NA
  wr_data$nip[wr_data$nip == ""] <- NA
  wr_data$preg[wr_data$preg == "not noted"] <- NA
  wr_data$preg[wr_data$preg == ""] <- NA
  wr_data$test[wr_data$test == "not noted"] <- NA
  wr_data$test[wr_data$test == ""] <- NA
  wr_data$fate[wr_data$fate == ""] <- NA
  wr_data$handler[wr_data$handler == ""] <- NA

  #clean names
  wr_data <- wr_data %>%
    remove_empty(which="rows") %>%
    clean_names %>%
    #rename column
    rename(session = sess) %>%
    #mutate columns to their proper formats
    mutate(year = as.numeric(year),
           occasion = as.numeric(occasion), #as.factor(occasion),
           site = as.character(site),
           food_trt = as.factor(food_trt),
           helm_trt = as.factor(helm_trt),
           date = as_date(date, format= "%m/%d/%Y"),
           session = as.numeric(session),
           trap = as.character(trap),
           tag = as.character(tag),
           sex = as.factor(sex),
           per = as.factor(per),
           nip = as.factor(nip),
           preg = as.factor(preg),
           test = as.factor(test),
           fate = as.factor(fate),
           handler = as.factor(handler),
           notes = as.character(notes)) %>%
    #filter for desired year
    filter(year==yr) %>%
    #remove uusi data
    filter(site!="uusi") %>%
    filter(!is.na(trap)) %>% #remove NA trap
    filter(!is.na(tag)) #remove NA tag

  # #check for missing data
  # wr_data %>% filter(is.na(session)) #no missing sessions
  # wr_data %>% filter(is.na(date)) #no missing dates

  #check for spelling errors, extra groups, weird data
  # unique(wr_data$year)
  # unique(wr_data$occasion)
  # unique(wr_data$site)
  # unique(wr_data$food_trt)
  # unique(wr_data$helm_trt)
  # unique(wr_data$session)
  # unique(wr_data$sex)
  # unique(wr_data$per)
  # unique(wr_data$nip)
  # unique(wr_data$preg)
  # unique(wr_data$test)
  # unique(wr_data$fate)

  ############ DEAD ANIMALS ARE STILL IN THE DATASET ###################
  # #remove animals euthanized for terminal sampling
  # filter(!fate == "S") %>%
  # #remove animals DP or DT
  # filter(!fate == "DP") %>%
  # filter(!fate == "DT")

  ################################ create a time column to replace session (for CMRnet) ##################################
  wr_data <-
    wr_data %>%
    mutate(time = case_when(
      session == "1" ~ "06:00:00",
      session == "2" ~ "18:00:00",
      session == "3" ~ "06:00:00",
      session == "4" ~ "18:00:00",)) %>%
    relocate(time, .after = session)

  #turn the time into a lubridate time
  wr_data$time <- hms(wr_data$time)
  #combine date and time columns into a lubridate time
  wr_data$date_time <- ymd_hms(paste(wr_data$date, wr_data$time))
  #move date and time around
  wr_data <-
    wr_data %>%
    relocate(date, .after = year) %>%
    relocate(session, .after = occasion) %>%
    relocate(date_time, .after = date) %>%
    dplyr::select(-time)
  ######################################################### END ############################################################

  ################################ convert trap number to a grid coordinate ####################################
  wr_data<-
    wr_data %>%
    separate(trap, into = c("x", "y"), sep = "(?<=[A-Z])(?=[0-9])", remove=FALSE) %>%
    #recode each letter as a number
    #this makes a new column with the letter turned into its corresponding number
    mutate(new_x = case_when(
      x == "A" ~ 1,
      x == "B" ~ 2,
      x == "C" ~ 3,
      x == "D" ~ 4,
      x == "E" ~ 5,
      x == "F" ~ 6,
      x == "G" ~ 7,
      x == "H" ~ 8,
      x == "I" ~ 9,
      x == "J" ~ 10,
      x == "K" ~ 11)) %>%
    #move that new column right after the original letter column
    relocate(new_x, .after = x) %>%
    #remove that original letter column
    dplyr::select(-x) %>%
    #rename the new number x column
    rename(x = new_x) %>%
    #make sure y and x are both numeric, not characters
    mutate(x = as.numeric(x), y = as.numeric(y))


  wr_data<-
    wr_data %>%
    #adjust the trap number (y) so it corresponds to a grid number
    #if remainder of x/2 is 0 (if x is even), multiply y by 2 - if else multiply by 2 then subtract 1
    #not sure why, but this code doesn't work if you string it together with the above code... the values for y are all wrong
    mutate(y_new = ifelse((x %% 2) == 0, ((wr_data$y)*2), (((wr_data$y)*2)-1))) %>%
    relocate(y_new, .after = y) %>%
    #remove the old y column and rename the new one
    dplyr::select(-y) %>%
    rename(y = y_new)
  ######################################################### END ############################################################

  #remove un-needed columns from data table - make a smaller version for easier analysis
  recap <- wr_data %>%
    dplyr::select(-date, -handler, -notes)


  #############################  COMBINE AND CLEAN FULLTRAP DATA   ###############################################

  #compare columns and their classes
  # compare_df_cols(trap, recap)

  #combining trap and recap dataframes
  fulltrap <- bind_rows(trap, recap)

  #change food_trt to "fed" "unfed"
  fulltrap$food_trt <- fct_recode(fulltrap$food_trt, "fed"="supplement", "unfed"="control")
  fulltrap$food_trt <-fct_relevel(fulltrap$food_trt, c("unfed", "fed"))
  #create trt column
  fulltrap <- fulltrap %>%
    unite(trt, food_trt, helm_trt, sep = "_", remove = FALSE)
  #and then relevel it
  fulltrap$trt <-fct_relevel(fulltrap$trt, c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))

  #add month column to replace occasion - keep occasion, session since later code might need numbers
  fulltrap <- fulltrap %>% mutate(month = case_when(occasion == "1" ~ "may",
                                                    occasion == "2" ~ "june",
                                                    occasion == "3" ~ "july",
                                                    occasion == "4" ~ "aug",
                                                    occasion == "5" ~ "sept",
                                                    occasion == "6" ~ "oct")) %>%
    mutate(month = factor(month, levels=c("may", "june", "july", "aug", "sept", "oct"))) %>%
    relocate(month, .after=date_time)

  #FIRSTCAP column: give a 1 if the capture is the first occurrence of that tag, else 0
  fulltrap <- fulltrap %>%
    unite(occ.sess, occasion, session, sep = ".", remove = FALSE) %>% #make a new occ.sess column so I can do things in order
    group_by(tag) %>%
    mutate(firstcap = ifelse(occ.sess == min(occ.sess), 1, 0)) %>%
    mutate(firstcap = factor(firstcap)) %>%
    ungroup()
  #this will record all the WRs as 0 as well

  ################ END FIRSTCAP/NEW ###########################

  # #in case igraph is already running, it masks "%--%" and the lubridate code won't run
  # # library(needs)
  # # prioritize(lubridate)
  # #create a measure of TIME KNOWN ALIVE (also date of first cap, last capture)
  # fulltrap <- fulltrap %>%
  #   group_by(tag) %>%
  #   arrange(date_time) %>%
  #   mutate(first_seen = min(date_time),
  #          last_seen = max(date_time)) %>%
  #   mutate(days_known_alive = round( as.duration(first_seen %--% last_seen) / ddays(1) , digits=2) ) %>%
  #   ungroup()

  #create caps_per_life column - number of captures of that individual
  fulltrap <- fulltrap %>%
    group_by(tag) %>%
    mutate(caps_per_life = length(tag)) %>%
    ungroup()
  #count of number of unique traps per animal in lifetime
  fulltrap <- fulltrap %>%
    group_by(tag) %>%
    mutate(traps_per_life = length(unique(trap))) %>%
    relocate(traps_per_life, .after=caps_per_life) %>%
    ungroup()
  #create a RECAPPED column (binary - were you recapped or not)
  fulltrap <- fulltrap %>%
    group_by(tag) %>%
    mutate(recapped = case_when(
      caps_per_life == "1" ~ 0,
      caps_per_life > "1" ~ 1)) %>%
    mutate(recapped = factor(recapped)) %>%
    relocate(recapped, .before=caps_per_life) %>%
    ungroup()

  #season column
  fulltrap <- fulltrap %>%
    mutate(season = ifelse(month=="sept" | month=="oct", "fall", "summer")) %>%
    mutate(season = factor(season, levels=c("summer", "fall"))) %>%
    relocate(season, .after=month)
  #caps_per_season
  fulltrap <- fulltrap %>%
    group_by(tag, season) %>%
    mutate(caps_per_season = length(tag)) %>%
    mutate(res = ifelse(caps_per_season >= 5, "resident", "nonresident")) %>% #add resident status
    mutate(res = factor(res)) %>%
    relocate(res, .after = traps_per_life) %>%
    ungroup()
  ################ NOTE: RESIDENT = 5 caps per SEASON -- NOT 5 caps per LIFE #############
  #traps_per_season
  fulltrap <- fulltrap %>%
    group_by(tag, season) %>%
    mutate(traps_per_season = length(unique(trap))) %>%
    ungroup()

  #caps_per_occ column - number of captures per occasion
  fulltrap <- fulltrap %>%
    group_by(occasion, tag) %>%
    mutate(caps_per_occ = length(tag)) %>%
    ungroup()
  #count of number of unique traps per animal per occasion
  ##n_distinct() is a dplyr wrapper for length(unique())
  fulltrap <- fulltrap %>%
    group_by(tag, occasion) %>%
    mutate(traps_per_occ  = length(unique(trap))) %>%
    ungroup()

  #make a column (current_breeder) of "breeder"/"nonbreeder" based on 1 for per/nip/preg or test IN THAT CAPTURE
  #NAs persist
  #then make a new column (ever_breeder) based on current_breeder column
  #where: if ever the animal was breeder during any capture, it gets "breeder"
  #use NA.rm in this to make sure NAs are ignored and get replaced with the right word
  #HOWEVER - if animal NEVER had breeding condition recorded in any capture, ever_breeder will be NA
  fulltrap <- fulltrap %>%
    #current_breeder in this occasion
    mutate(current_breeder = case_when(per=="1" | nip=="1" | preg == "1" ~ "breeder",
                                       test == "1" ~ "breeder",
                                       per=="0" & nip=="0" & preg == "0" ~ "nonbreeder",
                                       test == "0" ~ "nonbreeder")) %>%
    relocate(current_breeder, .after=sex) %>%
    #season_breeder if ever in the season (June-Aug, Sept-Oct)
    group_by(season, tag) %>%
    mutate(season_breeder_init = ifelse(all(is.na(current_breeder)), NA,
                                 ifelse(any(current_breeder == "breeder", na.rm = TRUE), "breeder", "nonbreeder"))) %>%
    relocate(season_breeder_init, .after = current_breeder) %>%
    #ever_breeder if ever in June-Oct
    group_by(tag) %>%
    mutate(ever_breeder = ifelse(all(is.na(current_breeder)), NA,
                                 ifelse(any(current_breeder == "breeder", na.rm = TRUE), "breeder", "nonbreeder"))) %>%
    relocate(ever_breeder, .after = season_breeder_init) %>%
    #adjust fall non-breeders to make sure they're Fall breeders if they were in the Summer
    mutate(season_breeder = ifelse(season=="fall" & season_breeder_init=="nonbreeder" & ever_breeder=="breeder",
                                    "breeder", season_breeder_init)) %>%
    mutate(season_breeder = factor(season_breeder, levels=c("breeder", "nonbreeder"))) %>%
    relocate(season_breeder, .after=season_breeder_init)
    # select(!c(current_breeder, ever_breeder, season_breeder_init))

  #reorganize this unholy mess
  fulltrap <- fulltrap %>%
    relocate(c(tag,samp_id,firstcap,trap,x,y,caps_per_occ,traps_per_occ), .after=helm_trt)

  #aight, I only want NAs if we DO NOT have that information at all in any way, so...
  #fill within an occasion for data not collected on WRs but collected during processing
  #ie: sex, breeding condition (per,nip,preg,test), current breeder, head, mass
  fulltrap <- fulltrap %>%
    group_by(tag, occasion) %>%
    arrange(session, .by_group = TRUE) %>%
    fill(c(samp_id, sex, per, nip, preg, test,
           season_breeder, mass, head), .direction="downup") %>% #fill missing data within an occasion
    ungroup()

  #### 2021 VOLES WITH SEX=NA (DUE TO FUCKUPS and other things and it's okay, I'm not mad about it)
  # 226128 # 226211 # 226769 # 219682
  ## some 2022 voles (at least 2) have sex=NA because we couldn't determine the correct sex


  #############################################################
  ## Body Condition Score Dropped (see 01_2021_data_cleaning)##
  #############################################################

  #based on some research on things, I've decided that 17.5g is the juv-subadult/adult cutoff
  #not bad, September breeders are v tiny regardless, but changing the cutoff to 15g makes it better other months

  # fulltrap <- fulltrap %>%
  #   mutate(adult = ifelse(mass>17.5, "adult", "juv-subA")) %>%
  #   relocate(adult, .after=sex)
  #
  # check <- fulltrap %>%
  #   group_by(month, current_breeder, adult) %>%
  #   summarise(n = length(tag))
  #
  # check %>%
  #   drop_na( c(adult,current_breeder) ) %>%
  #   ggplot(aes(x=current_breeder, y=n, fill=adult)) +
  #   geom_bar(stat="identity", position="dodge") +
  #   facet_wrap(~month)

  ###############################################################################################


  #save fulltrap so I can pull it for other scripts

  # Save fulltrap to a rdata file
  saveRDS(fulltrap, file = here(fulltrap_output))

  # Restore fulltrap from the rdata file
  # fulltrap <- readRDS(file = OUTPUT FILE NAME GOES HERE)

  #################### END ########################


}













########################################################################
############# Extra stuff, just for reference and posterity ############
########################################################################

######### A few checks on the full data set to make sure everything is happy #########

########## FOR 2022: TO FIND THE DISCREPANCIES BETWEEN FIRSTCAP and "NEW" column -- FIXED 5.4.23 ##############
# check <- fulltrap22 %>% select(occ.sess, tag, new, firstcap) %>%
#   drop_na(new) %>%
#   mutate(new = case_when(new == "1" ~ 1,
#                          new == "0" ~ 0)) %>%
#   mutate(firstcap = case_when(firstcap == "1" ~ 1,
#                          firstcap == "0" ~ 0)) %>%
#   mutate(confirm = new + firstcap) %>%
#   filter(confirm == 1)
#
# check2 <- fulltrap22 %>% filter(tag %in% check$tag) %>%
#   group_by(tag) %>%
#   arrange(occ.sess, .by_group = TRUE)

# fulltrap22 <- fulltrap22 %>% select(!new)

################################################

#there are also a number of animals that were on the vole processing sheet twice because they were
#captured, weighed, and released due to time issues (mostly occasion 6)
#I need to find these animals and make sure they only have 0 or 1 for the first occasion they were seen,
#not both if they have two processing entries

#this is why I created the 'firstcap' column like above

#confirm that firstcap is doing what we want...
# first <- fulltrap %>% filter(firstcap == "1")
# notfirst <- fulltrap %>% filter(firstcap == "0")
# nrow(first) + nrow(notfirst)
#06.01.22 -- yes, this is good, first + notfirst = fulltrap

#NOTE: a few animals were processed TWICE in a single occasion - make sure things are working properly for them:

#############################################

#make sure there isn't any weird fuckery with animals on multiple grids or with multiple entries per occ.sess

######## check to see how many unique grids animals are on ########
# fulltrap %>%
#   group_by(tag) %>%
#   summarise(site_unique = n_distinct(site)) %>%
#   filter(site_unique > 1) %>%
#   ungroup()
#if this is 0, you're good

######## check to make sure animals don't have more entries than unique occ.sess combos ########
# check <- fulltrap %>% group_by(tag) %>% summarise(count = length(occ.sess), unique = n_distinct(occ.sess))
# check <- check %>% mutate(diff=unique-count) %>% filter(diff != 0)
# tags <- check$tag
# check <- fulltrap %>% filter(tag %in% tags) %>% arrange(tag, occ.sess)
#if this is 0, you're good


######## check to make sure everyone only has one sex #########
######## these issues were checked and then resolved in the raw data #######
# fulltrap$sex <- as.character(fulltrap$sex)
# check <- fulltrap %>% group_by(tag) %>%
#   summarise(n = unique(sex)) %>%
#   filter(!is.na(n))
#save the tag IDs of animals with multiple sexes
# duplist <- as.vector(check$tag[duplicated(check$tag)])
#how many?
# length(duplist) #0 voles!
# sexswaps <- fulltrap %>%
#   filter(tag %in% duplist) %>%
#   arrange(tag, occ.sess)
#write as csv to share
# write.csv(sexswaps, here("confirm_sex.csv"), row.names=FALSE)

############################################### END CHECKING for FUCKERY #######################################################

