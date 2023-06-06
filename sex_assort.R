#clear environment
rm(list = ls())

ft21 <- readRDS(here("fulltrap21_05.10.23.rds")) %>% ###### edit on 5/11 - ONLY breeders, add this to cleaning if I like it
  filter(season_breeder=="breeder")
data = ft21
networks_file = "overlapnets21.rds"


sex_assort <- function(data, networks_file, netmets_file){

  ##---------------- LOAD THE DATA ----------------------

  #running using assortnet (D.Farine) because it can work with weighted networks
  library(assortnet)

  #load, clean the fulltrap dataset (make sure it's the most recent version)
  fulltrap <- data %>%
    filter(month != "may") %>% #drop may data since not all sites had captures
    mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
    drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
    unite("sts", season, trt, sex, remove=FALSE) #add sts column to match params_summary

  # load the network data
  overlap_network_list <- readRDS(here(networks_file))

  #create tag_sex df for assortativity by sex
  tag_sex <- fulltrap %>% group_by(tag) %>% slice(1) %>%
    select(tag, sex) %>%
    arrange(tag)


  ##--------------------------------------------

  # Create list to store results
  sex_assort_list <- list()

  # Calculate network metrics
  for(i in 1:length(overlap_network_list)){

    #length(overlap_network_list) == for each of the 12 sites

    #for each site
    print(names(overlap_network_list[i]))
    site.id <- names(overlap_network_list[i])

    site <- list()

    for(j in 1:length(overlap_network_list[[i]])){

      #(j in 1:length(overlap_network_list[[i]])) == for the 4-5 trapping occasions
      #would be (j in 1:5) except that Kiirastuli 2022 only has 4 :P *thhhbbbbt*

      #for each month
      print(names(overlap_network_list[[i]][j]))
      month <- names(overlap_network_list[[i]][j])

      #dataframe to hold results per month
      site[[j]] <- data.frame(month)

      adjmat <- overlap_network_list[[i]][[j]]
      diag(adjmat) <- 0 #matrix diagonal is NA - assortnet needs it to be 0

      #create WEIGHTED NETWORK from adjacency matrix
      inet <- graph_from_adjacency_matrix(adjmat, mode="undirected", weighted = TRUE, diag = FALSE)

      ids <- get.vertex.attribute(inet, "name") #tag ids for all the animals on the grid

      ### FOR ASSORTATIVITY
      #filter tag_sex for only ids caught this site/occ
      ids_sex <- tag_sex %>% filter(tag %in% ids)
      #calculate assortativity
      out <- assortment.discrete(adjmat, as.vector(ids_sex$sex), weighted=TRUE, SE=FALSE, M=1, na.rm=FALSE)
        #out is a list, $r has the assort coef across all individuals, $mixing_matrix has ppn of edges by sex

      check <- n_distinct(ids_sex$sex) #how many sexes are represented? (to catch site/month when only M or F are present)

      mat <- out$mixing_matrix #save just the mixing matrix
      #pull the percent of M/F overlap, F/F, and M/M
        #check is to make sure NA is input if there was only one sex of breeders in that month
      site[[j]]$fm <- ifelse(check==1, NA, mat["M","F"]*2) #double the fm overlaps since network is undirected
      site[[j]]$ff <- ifelse(check==1, NA, mat["F","F"])
      site[[j]]$mm <- ifelse(check==1, NA, mat["M","M"])

    }

    #write the list 'site' as a 1st-order item in wt_net_mets_list
    sex_assort_list[[i]] <- site

  }

  #name the 12 1st order elements of wt_nets_list as the sites
  names(sex_assort_list) <- names(overlap_network_list)

  #rename the sublist items (months) for each site
  #accounting for the fact that most sites have 5 months of data but some have 4
  for(i in 1:length(sex_assort_list)){
    ifelse( length(sex_assort_list[[i]]) == 5,
            names(sex_assort_list[[i]]) <- c("june", "july", "aug", "sept", "oct"),
            names(sex_assort_list[[i]]) <- c("july", "aug", "sept", "oct") )
  }


  ############## condense sex_assort_list to make it easier to use for analysis #######################

  #collate results
  #this collapses the 2nd order elements (network metrics for a single month) down to the 1st order element (site)
  #so now wt_net_mets_list_summary is a list of 12 dfs, each is all the network metrics for a site across all the months

  #make a list to store things
  sex_assort_summary <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(sex_assort_list)){

    #for all 12 sites
    summary <- do.call("rbind", sex_assort_list[[i]])
    sex_assort_summary[[i]] <- summary
  }

  #remove the row names of the dfs (its the month, but we already have that info)
  for(i in 1:length(sex_assort_summary)){
    row.names(sex_assort_summary[[i]]) <- NULL
  }

  #name the 12 1st order elements as their sites
  names(sex_assort_summary) <- names(overlap_network_list)

  ## make sex_assort_summary into freiggein huge df
  sex_assort_summary <- do.call(rbind.data.frame, sex_assort_summary)

  #clean up the df
  sex_assort_summary <- sex_assort_summary %>%
    rownames_to_column("name") %>% #row names are the sites, make that a column
    separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
    mutate(site = as.factor(site)) #make site a factor

  #save it
  # saveRDS()


}



#visualize
sex_assort_summary %>%
  left_join(read.csv(here("grid_trts.csv")), by="site") %>% unite(trt, food_trt, helm_trt) %>%
  pivot_longer(-c(site, trt, month), names_to = "group", values_to = "pct") %>%
  mutate(month = factor(month, levels=c('june', "july", "aug", "sept", "oct"))) %>%
  ggplot(aes(x=month, y=pct, color=group)) +
  geom_point() +
  facet_wrap(~trt)
