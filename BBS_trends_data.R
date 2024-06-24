#### Download and prep NA BBS data ####
# Required libraries
library(bbsBayes2) # devtools::install_github("bbsBayes/bbsBayes2")
#library(bbsBayes)
library(tidyr)
library(dplyr)
library(dtplyr)
library(ade4)
library(ape)
library(rotl)
library(sf)
library(spdep)
library(MRFtools) # devtools::install_github("eric-pedersen/MRFtools")

# Download NA BBS data up to 2022 using
# bbsBayes2 package; this only needs to be done once as the data will be stored
# and used for any future calls from bbsBayes functions
if(!bbsBayes2::have_bbs_data(quiet = TRUE)){
 bbsBayes2::fetch_bbs_data() # only downloads if doesn't already exist locally
}

# Load the full bbs dataset
all_dat <- bbsBayes2::load_bbs_data()

counts <- all_dat$birds %>%  # positive counts of each species during every BBS survey
  dplyr::select(route_data_id, aou, species_total) #dropping all but the critical columns

sampling_events <- all_dat$routes %>%  # date, time, starting location, for every BBS survey since 1966
  dplyr::select(country, state, st_abrev, route_name, bcr,
         country_num, state_num,route,
         latitude, longitude,
         route_data_id, #this is the critical unique id for a sampling event
         year, month, day,obs_n)

species_list <- all_dat$species %>% # species list
  dplyr::filter(unid_combined == FALSE,
                !(grepl("^(unid.)", english) | grepl("^(hybrid)", english) )) %>%  #this is specific to the BBS,
  dplyr::select(aou, english, french, genus, species) %>%  # the TRUE option combines 13 taxonomic units that have been split or
  dplyr::mutate(latin = paste(genus, species)) # lumped over the history of the BBS, this approach (== FALSE) retains the taxonomic
# units, exactly as they are stored in the BBS database

luni <- function(x){
  y <- length(unique(x))
}

# strata loop to limit species lists by strata ----------------------------
#setting limits on the strata that are included
more_than_n_routes_stratum <- 10 # More than this number of BBS routes in the stratum
more_than_n_surveys_stratum <- 200 # More than this number of surveys (routes * years) in the stratum

strata_list <- sampling_events %>% 
  dplyr::group_by(bcr) %>% 
  dplyr::summarise(.,
            n_surveys = dplyr::n(),
            n_routes = luni(route_name)) %>% 
  dplyr::filter(n_routes > more_than_n_routes_stratum, #filter on the total number of routes (sampling locations)
         n_surveys > more_than_n_surveys_stratum) #filter on the total number of sampling events
#seems high, but it's only ~4 surveys annually (low for such large regions)
# above drops BCR1 and BCR3
                        
data_all <- NULL # empty object to facilitate bind_rows() at the end of the next loop

# this strata-specific approach removes the species by strata
# combinations in the data that are exclusively 0-values. 
# Species are only included in a given strata, if they meet the 
# following inclusion criteria:
more_than_n_counts_species <- 100 
# More than this number of surveys (routes * years) on which the species has been observed
# this is probably too low... only ~2 counts/year 

for(reg in strata_list$bcr){
  samp_ev_select <- sampling_events %>% 
    dplyr::filter(bcr == reg)
  
  counts_select <- counts %>% 
    dplyr::filter(route_data_id %in% samp_ev_select$route_data_id)
  
  # counting number of observations of each species
  species_inc <- counts_select %>% 
    dplyr::group_by(aou) %>% 
    dplyr::summarise(.,
              n_obs = dplyr::n()) #n_obs represents the number of non-zero counts across routes and year
  
  #identifying which species to keep (n_obs > more_than_n_counts_species)
  species_keep <- species_inc %>% 
    dplyr::filter(n_obs > more_than_n_counts_species) %>%  
    dplyr::select(aou)
  
  # complete set of species by surveys
  samp_ev_select <- samp_ev_select %>% 
    expand_grid(., species_keep)
  # above creates a complete set of species by survey combinations
  
  #zero-filling
  data_sel <- samp_ev_select %>% 
    dplyr::full_join(.,counts_select,
              by = c("route_data_id",
                     "aou")) %>% 
    dplyr::mutate(species_total = ifelse(is.na(species_total),
                                  0,
                                  species_total))
  # above creates the full set of observations with appropriate zero values for
  # the species in species_keep at all surveys conducted in the region
  
  
  data_all <- dplyr::bind_rows(data_all,
                        data_sel)
  rm(data_sel)
}  

data_all <- data_all %>% 
  dplyr::left_join(.,species_list,
            by = "aou")
n_species = length(unique(data_all$english))

dir.create('BBS_data', recursive = TRUE, showWarnings = FALSE)
saveRDS(data_all, "BBS_data/All_BBS_data_by_BCR.rds")



data_all <- readRDS("BBS_data/All_BBS_data_by_BCR.rds")
all_data <- data_all %>% 
  rename(route_id = route_name,
         sp_latin = latin,
         count = species_total) %>% 
  dplyr::mutate(strata_name = paste0("BCR", bcr)) %>%  #to match the strata names in the bcr map below
  dplyr::filter(year > 1991)
## 12 Million rows in above

head(all_data)
length(unique(all_data$english))
# 476 species


# Aggregate route-level counts to sums within strata ----------------------


# Calculate some useful measures for approximating sampling effort
# which we can use for filtering the data
all_data %>%
  # Group by species and stratum to take the sum of counts
  # per year; use dtplyr to convert to data.table code for faster
  # grouping operations
  lazy_dt() %>%
  dplyr::group_by(strata_name, year) %>%
  # Calculate number of route observations per region, per year
  # to help form an offset of sampling effort
  dplyr::mutate(n_records = dplyr::n_distinct(route_id)) %>%
  dplyr::ungroup() %>%
  # Aggregate species' counts by region and by year
  dplyr::group_by(sp_latin, strata_name, year) %>%
  dplyr::mutate(count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::select(count, year, sp_latin, strata_name, n_records) %>%
  dplyr::distinct() %>%
  # Calculate the number of records > 1 for each species, which
  # we may use later to define a threshold for limiting the size
  # of the data
  dplyr::group_by(sp_latin) %>%
  dplyr::mutate(n_nonzero = length(which(count > 0))) %>%
  dplyr::ungroup() %>%
  as_tibble() -> all_data_sub
length(unique(all_data_sub$sp_latin))

# For now, keep only the 100 most commonly-observed species for starters
# so we can get the pipeline working and tested
all_data_sub %>%
  dplyr::select(sp_latin, n_nonzero) %>%
  dplyr::distinct() %>%
  dplyr::arrange(dplyr::desc(n_nonzero)) %>%
  dplyr::slice_head(n = 100) %>%
  dplyr::pull(sp_latin) -> sp_keep

all_data_sub %>%
  dplyr::filter(sp_latin %in% sp_keep) -> all_data_sub
length(unique(all_data_sub$sp_latin))
min(all_data_sub$n_nonzero)

# Complete missing combinations with zeros (safe to ignore these warnings)
# # warnings relate to region*year combinations with no survey-data
all_data_sub %>%
  tidyr::complete(tidyr::nesting(strata_name),
                  sp_latin, year,
                  fill = list(count = 0)) %>%
  dplyr::group_by(strata_name, year) %>%
  dplyr::mutate(n_records = ifelse(is.na(n_records), 
                                   max(n_records, na.rm = TRUE),
                                   n_records)) %>%
  dplyr::mutate(n_records = ifelse(is.finite(n_records),
                                   n_records,
                                   0)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(count = ifelse(n_records == 0,
                               NA,
                               count)) -> all_data_sub

tmp <- all_data_sub %>% 
  filter(is.na(count))
# na values for count reflect region*year combinations with no surveys
# this is now empty because of alternate data preparation

# Add functional trait information using the EltonTraits bird
# dataset of Wilman et al 2014 (10.1890/13-1917.1)
temp <- tempfile()
download.file('https://ndownloader.figshare.com/files/5631081',
              temp)
traits <- read.table(temp, header = TRUE,
                     fill = TRUE,
                     quote = '"',
                     stringsAsFactors = FALSE,
                     sep = "\t")
unlink(temp)


# Join the data together and filter out any species that don't have
# associated trait information
all_data_sub %>%
  dplyr::left_join(traits, by = c('sp_latin' = 'Scientific')) %>%
  dplyr::filter(!is.na(PassNonPass)) -> all_data_sub
length(unique(all_data_sub$sp_latin))

# Add phylogenetic information using the Open Tree of Life
tx_search <- rotl::tnrs_match_names(names = unique(all_data_sub$sp_latin),
                              context_name = "All life")
ott_in_tree <- rotl::ott_id(tx_search)[
  rotl::is_in_tree(rotl::ott_id(tx_search))]
ott_remove <- setdiff(names(ott_in_tree), unique(all_data_sub$sp_latin))
ott_in_tree <- ott_in_tree[!names(ott_in_tree) %in% ott_remove]
tr <- rotl::tol_induced_subtree(ott_ids = ott_in_tree,
                          label_format = 'name')
tr <- ape::compute.brlen(tr)

# Save the full phylogenetic tree and trait database, which may be useful
# for later plotting
dir.create('data', recursive = TRUE, showWarnings = FALSE)
saveRDS(tr, "./data/tree.rds")
saveRDS(traits, "./data/traits.rds")

tr <- readRDS("./data/tree.rds")
traits <- readRDS("./data/traits.rds")

# Filter the BBS data to only keep those species with phylogenetic information
all_data_sub %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sp_latin = sub(' ', '_', sp_latin)) %>%
  dplyr::filter(sp_latin %in% tr$tip.label) -> mod_data

# Set levels of species names to match tree order
mod_data$sp_latin <- factor(mod_data$sp_latin, 
                            levels = tr$tip.label)
setdiff(mod_data$sp_latin, tr$tip.label)
setdiff(tr$tip.label, mod_data$sp_latin)
length(tr$tip.label)

# Create the phylogenetic penalty matrix
phylo_penalty <- MRFtools::mrf_penalty(tr)
setdiff(mod_data$sp_latin, dimnames(phylo_penalty)[[1]])

# Now create the functional dendrogram and its associated
# penalty matrix. First gather the desired trait variables 
# that will be used for creating ecological dendrograms. 
# Here, we choose variables to represent species' proportional 
# use of seven different foraging habitat categories and 
# proportional use of ten diet categories, as these categories 
# together should give a reasonable representation of species' 
# local ecological niches
hab_dat <- mod_data %>%
  dplyr::select(matches(c('sp_latin', colnames(traits)[24:30]))) %>%
  dplyr::distinct() %>%
  as.data.frame()
hab_dat <- hab_dat[order(match(hab_dat$sp_latin, tr$tip.label)),]
rownames(hab_dat) <- hab_dat[,1]; hab_dat <- hab_dat[,-1]

diet_dat <- mod_data %>%
  dplyr::select(matches(c('sp_latin', colnames(traits)[10:19]))) %>%
  dplyr::distinct() %>%
  as.data.frame()
diet_dat <- diet_dat[order(match(diet_dat$sp_latin, tr$tip.label)),]
rownames(diet_dat) <- diet_dat[,1]; diet_dat <- diet_dat[,-1]

source('Functions/utilities.R')
func_penalty <- preptrait_penalty(trait_dfs = list(hab_dat, diet_dat), 
                                  prep_types = c('prep.fuzzy', 'prep.fuzzy'),
                                  ordering = tr$tip.label)
setdiff(mod_data$sp_latin, dimnames(func_penalty)[[1]])

# Load the BCR strata shapefile
bcr_shp <- bbsBayes2::load_map("bcr") %>% 
  filter(strata_name %in% unique(mod_data$strata_name))
# bcr_shp <- sf::st_read("maps/BBS_BCR_strata.shp", 
#                         layer = "BBS_BCR_strata")

# Ensure there are no differences in the stratum names
setdiff(mod_data$strata_name, bcr_shp$strata_name)

# Ensure strata_name is also a factor in the data
mod_data$strata_name <- factor(mod_data$strata_name, 
                         levels = bcr_shp$strata_name)

# Create the binary neighborhood adjacency matrix that defines
# which polygons are the neighbours of each focal polygon
strat_penalty <- spdep::poly2nb(bcr_shp, row.names = as.character(bcr_shp$strata_name))
names(strat_penalty) <- as.character(bcr_shp$strata_name)#attr(strat_penalty, "region.id")

# Create a non-trait random intercept penalty that can allow us to use
# reduced-rank MRFs for random effects based purely on the species name
sp_penalty <- MRFtools::mrf_penalty(mod_data$sp_latin,
                                    type = 'individual')

# Final tidying of data for modelling; remove un-needed trait columns
# and use a log(x+1) transformation for the number of records per 
# space-time combination, which will be helpful to include as an offset.
# We also add the log of the empirical mean count for each species*region
# combination, which will be useful for informing species' spatial average counts
# and help us to hopefully avoid the need for an extremely large random effect
# to capture this nuisance variation
logmean = function(x){
  if(all(is.na(x))){
    x <- 0
  }
  
  x <- x[!is.na(x)]
  
  if(all(x == 0)){
    out <- 0
  } else if(mean(x) == 0){
    out <- 0
  } else {
    out <- log(mean(x))
  }
  return(out)
}

mod_data %>%
  dplyr::mutate(sp_latin_phy = sp_latin,
                sp_latin_func = sp_latin) %>%
  dplyr::select(-c(n_nonzero:Record.Comment)) %>%
  
  # Number of observations per region x time combo
  dplyr::mutate(n_records = log(n_records + 1)) %>%
  
  # Standardizing counts to zero-centred, unit variance; this is 
  # probably a necessary evil to ensure we can validate out of sample
  # predictions of trend trajectories
  dplyr::group_by(sp_latin, strata_name) %>%
  dplyr::mutate(mean_est = mean(count, na.rm = TRUE),
                sd_est = sd(count, na.rm = TRUE)) %>%
  dplyr::mutate(count_sc = (count - mean_est) / sd_est) %>%
  dplyr::ungroup() -> mod_data

# Drop species * region combinations that were all zeros 
# to remove the stable-at-zero trajectories
sp_by_bcr <- mod_data %>% 
  group_by(sp_latin,strata_name) %>% 
  summarise(n_obs = n(),
            n_zero = length(which(count == 0))) %>%
  dplyr::mutate(keep = n_zero < n_obs) %>%
  dplyr::ungroup()

mod_data %>%
  dplyr::left_join(sp_by_bcr) %>%
  dplyr::filter(keep == TRUE) %>%
  # Add a factor to capture any remaining variation in species' 
  # spatially-varying intercepts
  dplyr::mutate(strat_sp  = interaction(strata_name, sp_latin)) %>%
  droplevels() -> mod_data
rm(sp_by_bcr)

# Save all data objects
save(mod_data, phylo_penalty, 
     func_penalty, strat_penalty,
     sp_penalty,
     file = "./data/model_objects.rda")

