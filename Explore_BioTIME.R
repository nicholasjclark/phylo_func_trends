#########
# Explore BIOTIME
#
# Data can be found here: https://zenodo.org/records/5026943 
#########


# specify dir -------------------------------------------------------------

dir <- '~/Work/Research/Projects/phylo_func_trends/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(janitor)


# read in data ------------------------------------------------------------

biotime <- read.csv(paste0(dir, 'Data/L0/BioTIME/BioTIMEQuery_24_06_2021.csv')) %>%
  janitor::clean_names()

biotime_meta <- read.csv(paste0(dir, 'Data/L0/BioTIME/BioTIMEMetadata_24_06_2021.csv')) %>% 
  janitor::clean_names()


# process data ------------------------------------------------------------

#modified from Johnson et al: https://github.com/GitTFJ/correlated_effect_model/blob/main/code/data_compile.Rmd
biotime2 <- biotime %>% dplyr::select(study_id, latitude, longitude, year, plot,
                              genus_species, sum_allrawdata_abundance) %>%
  dplyr::left_join(dplyr::select(biotime_meta, study_id, abundance_type),
            dplyr::join_by(study_id)) %>% 
  dplyr::filter(abundance_type != "Presence/Absence" & !is.na(abundance_type)) %>% 
  dplyr::mutate(dataset_id = "BioTIME",
         site = paste0(study_id, '_', plot),
         country = NA,
         unit = dplyr::case_match(
           abundance_type,
           "Count" ~ "n_individuals",
           "MeanCount" ~ "mean_n_individuals",
           "Density" ~ "density"
         )) %>% 
  dplyr::rename(date = year, species = genus_species, abundance = sum_allrawdata_abundance) %>% 
  dplyr::select(dataset_id, site, country, latitude, longitude, date, species, abundance, unit)


# explore data ------------------------------------------------------------

str(biotime_meta)
hist(biotime_meta$number_of_species, breaks = 1000)
dplyr::filter(biotime_meta, number_of_species > 4000)
