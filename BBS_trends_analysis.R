#### Post-processing ####
library(mgcv)
library(dplyr)
library(ggplot2)
library(marginaleffects)

# Source a few utility functions
source('Functions/utilities.R')

# Source the fitted models and data objects
load("./data/model_objects.rda")
mod <- readRDS("./models/mod.rds")

# Don't run full mgcv summaries; they are extremely computationally
# expensive for these models
summary(mod, re.test = FALSE)
AIC(mod)

mod2 <- readRDS("./models/mod2.rds")
summary(mod2, re.test = FALSE)
AIC(mod2)

# Calculate phylogenetic contributions to species' nonlinear trend
# estimates
phylo_derivs <- compute_derivs(mod_data,
                               model = mod,
                               type = 'phylogenetic',
                               n_cores = 8)
plot_trait_conts(trait_derivs = phylo_derivs,
                 type = 'phylogenetic')

# Proportion of trends in which phylogenetic information contributes at
# least 50% to the estimated wiggliness of the function compared to 
# other information (i.e. spatial and non-phylogenetic variation)
phylo_derivs %>%
  # Filter out those with very small wiggliness estimates, as 
  # the contributions don't add any meaningful information here
  # because the estimated trend is flat
  dplyr::filter(wiggliness >= quantile(wiggliness, probs = 0.15)) %>%
  dplyr::summarise(prop = 100 * length(which(trait_conts >= 0.5)) / 
                     dplyr::n())

# Repeat for the functional relatedness model
func_derivs <- compute_derivs(mod_data,
                               model = mod2,
                               type = 'functional',
                               n_cores = 8)
plot_trait_conts(trait_derivs = func_derivs,
                 type = 'functional')
func_derivs %>%
  dplyr::filter(wiggliness >= quantile(wiggliness, probs = 0.15)) %>%
  dplyr::summarise(prop = 100 * length(which(trait_conts >= 0.5)) / 
                     dplyr::n())

# Plot some predictions against truth (note, these predictions will not include
# contributions from the residual AR process)
plot_sp_trends(model = mod,
               data = mod_data,
               species = 'Hirundo_rustica',
               type = 'response',
               median_records = FALSE)
plot_sp_trends(model = mod,
               data = mod_data,
               species = levels(mod_data$sp_latin)[8],
               type = 'response',
               median_records = FALSE)
plot_sp_trends(model = mod,
               data = mod_data,
               species = levels(mod_data$sp_latin)[9],
               type = 'response',
               median_records = FALSE)
plot_sp_trends(model = mod,
               data = mod_data,
               species = levels(mod_data$sp_latin)[30],
               type = 'response',
               median_records = FALSE)
plot_sp_trends(model = mod,
               data = mod_data,
               species = levels(mod_data$sp_latin)[65],
               type = 'response',
               median_records = FALSE)

# Plot some expected trends by fixing the offsets and asking
# how the model would expect species' abundances to change over time
plot_sp_trends(model = mod,
               data = mod_data,
               species = 'Hirundo_rustica',
               type = 'expected',
               median_records = TRUE)
plot_sp_trends(model = mod,
               data = mod_data,
               species = levels(mod_data$sp_latin)[8],
               type = 'expected',
               median_records = TRUE)
plot_sp_trends(model = mod,
               data = mod_data,
               species = levels(mod_data$sp_latin)[9],
               type = 'expected',
               median_records = TRUE)
plot_sp_trends(model = mod,
               data = mod_data,
               species = levels(mod_data$sp_latin)[30],
               type = 'expected',
               median_records = TRUE)
plot_sp_trends(model = mod,
               data = mod_data,
               species = levels(mod_data$sp_latin)[2],
               type = 'expected',
               median_records = TRUE)

# Some attempt at an 'average' temporal trend for the last
# 25 years of data. Use only 65% of species
plot_av_trend(model = mod, mod_data = mod_data, 
              prop_species = 0.65)

# Will need to:
# 1. fit simpler models that ignore all relationship information, as 
#    well as phylogenetic/functional slopes (linear) models for comparisons
# 2. leave 10% of combos (strata x species) out and generate
#    predictions; compute CRPS from each model; run several CV folds
# 3. assess variation in prediction accuracy across space and 
#    across the trees
# 4. visualise the proportional change in variance explained 
#    (as calculated above) on the tree to see if there are clusters
#    of species that depend more heavily on phylogenetic relationships
# 5. assess support for nonlinearity of trends; determine when trends
#    were accelerating / decelerating most rapidly



