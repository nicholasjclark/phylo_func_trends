#### Download and prep NA BBS data ####
# Required libraries
library(bbsBayes)
library(tidyr)
library(dplyr)
library(dtplyr)
library(ade4)
library(ape)
library(rotl)
library(sf)
library(spdep)
library(MRFtools) # devtools::install_github("eric-pedersen/MRFtools")

# Download NA BBS data up to 2021 using Adam Smith's 
# bbsBayes package; this only needs to be done once as the data will be stored
# and used for any future calls from bbsBayes functions
bbsBayes::fetch_bbs_data()

# Stratify counts by the Bird Conservation Regions, which are large but
# useful polygons that can help with dimension reduction
strat <- bbsBayes::stratify(by = "bcr")

# Clean and prepare data
strat$bird_strat %>%
  dplyr::left_join(strat$route_strat %>%
                     dplyr::mutate(route_id = 
                                     factor(paste0(statenum, '_', Route)))) %>%
  dplyr::left_join(strat$species_strat %>%
                     dplyr::mutate(AOU = as.integer(aou))) %>%
  dplyr::mutate(sp_latin = paste(genus, species),
                count = StopTotal) %>%
  dplyr::select(count, route_id, Year, Month, 
                Day, ObsN, strat_name,
                sp_latin) %>%
  dplyr::filter(is.finite(Year)) %>%
  dplyr::filter(is.finite(Month)) %>%
  dplyr::filter(is.finite(Day)) %>%
  dplyr::filter(is.finite(ObsN)) %>%
  dplyr::filter(is.finite(count)) %>%
  
  # Keep data from 1992 onward; for no real reason other than 
  # to limit the size of the data
  dplyr::filter(Year > 1991) %>%
  janitor::clean_names() %>%
  dplyr::mutate(ST_12 = strat_name) %>%
  dplyr::select(-strat_name) -> all_data
head(all_data)
NROW(all_data)
length(unique(all_data$sp_latin))

# Calculate some useful measures for approximating sampling effort
# which we can use for filtering the data
all_data %>%
  # Group by species and stratum to take the sum of counts
  # per year; use dtplyr to convert to data.table code for faster
  # grouping operations
  lazy_dt() %>%
  dplyr::group_by(ST_12, year) %>%
  # Calculate number of route observations per region, per year
  # to help form an offset of sampling effort
  dplyr::mutate(n_records = dplyr::n_distinct(route_id)) %>%
  dplyr::ungroup() %>%
  # Aggregate species' counts by region and by year
  dplyr::group_by(sp_latin, ST_12, year) %>%
  dplyr::mutate(count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::select(count, year, sp_latin, ST_12, n_records) %>%
  dplyr::distinct() %>%
  # Calculate the number of records > 1 for each species, which
  # we may use later to define a threshold for limiting the size
  # of the data
  dplyr::group_by(sp_latin) %>%
  dplyr::mutate(n_nonzero = length(which(count > 1))) %>%
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
all_data_sub %>%
  tidyr::complete(tidyr::nesting(ST_12),
                  sp_latin, year,
                  fill = list(count = 0)) %>%
  dplyr::group_by(ST_12, year) %>%
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
dir.create('data')
saveRDS(tr, "./data/tree.rds")
saveRDS(traits, "./data/traits.rds")

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
bcr_shp <- sf::st_read("maps/BBS_BCR_strata.shp", 
                        layer = "BBS_BCR_strata")

# Ensure there are no differences in the stratum names
setdiff(mod_data$ST_12, bcr_shp$ST_12)

# Ensure ST_12 is also a factor in the data
mod_data$ST_12 <- factor(mod_data$ST_12, 
                         levels = bcr_shp$ST_12)

# Create the binary neighborhood adjacency matrix that defines
# which polygons are the neighbours of each focal polygon
strat_penalty <- spdep::poly2nb(bcr_shp, row.names = bcr_shp$ST_12)
names(strat_penalty) <- attr(strat_penalty, "region.id")

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
  
  # Empirical mean (logged) of counts per species x region combo
  dplyr::group_by(sp_latin, ST_12) %>%
  dplyr::mutate(sp_st_mean = logmean(count)) %>%
  dplyr::ungroup() %>%
  
  # Use both of these to form the offset so we can avoid 
  # specifying large random effects just to estimate spatial
  # variation in species' average counts
  dplyr::mutate(offset = sp_st_mean + n_records) -> mod_data

# Save all data objects
save(mod_data, phylo_penalty, 
     func_penalty, strat_penalty,
     sp_penalty,
     file = "./data/model_objects.rda")
