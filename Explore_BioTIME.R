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
  dplyr::select(dataset_id, study_id, site, country, latitude, longitude, date, species, abundance, unit)


# explore data ------------------------------------------------------------

str(biotime_meta)

#how many studies of each taxonomic group
ggplot(biotime_meta, aes(taxa)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
#log(number of locations) - many are single location
hist(log(biotime_meta$number_lat_long), breaks = 50)

#log(number of species)
hist(log(biotime_meta$number_of_species), breaks = 50)

#filter based on arbitrary thresholds
biotime_meta_f <- dplyr::filter(biotime_meta, 
              number_lat_long >= 10,
              number_of_species >= 10,
              #time points
              data_points >= 25,
              taxa %in% c('Birds',
                          'Mammals',
                          'Amphibians',
                          'Terrestrial plants',
                          'Fish'))

#how many studies of each taxonomic group for filtered data
ggplot(biotime_meta_f, aes(taxa)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

str(biotime_meta_f)
biotime_meta_f$title
biotime_meta_f$methods

