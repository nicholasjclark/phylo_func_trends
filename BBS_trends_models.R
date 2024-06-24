#### Fit candidate phylogenetic and functional smoothing models
# to the full dataset ####
library(mgcv)
library(dplyr)

# Load data objects
load("./data/model_objects.rda")

# Source a few utility functions
source('Functions/utilities.R')

# Randomly drop 10% of species * region combinations for using as a 
# validation set
all_sp_combos <- mod_data %>%
  dplyr::select(sp_latin, strata_name) %>%
  dplyr::distinct()

validation_combos <- all_sp_combos %>%
  dplyr::slice_sample(prop = 0.1) %>%
  dplyr::mutate(training = 'No')

mod_data <- validation_combos %>%
  dplyr::right_join(mod_data) %>%
  dplyr::mutate(training = ifelse(is.na(training), 'Yes', 'No'))

training_data <- mod_data %>%
  dplyr::filter(training == 'Yes')

validation_data <- mod_data %>%
  dplyr::filter(training == 'No')
rm(validation_combos, all_sp_combos)

# Save the training splits
save(training_data,
        validation_data,
        file = "./data/training_split.rds")

# Some useful details on MRF smooths:
# https://stats.stackexchange.com/questions/638522/gam-model-with-spatial-account-via-mrf

# Fit a phylogenetic trend model uses a Gaussian observation
# process. The model is a decomposition that includes terms capturing 
# marginal spatial and spatiotemporal fields, as well as contributions
# from species' phylogenetic relationships. Higher-order interaction terms
# allow each species' spatiotemporal field to be informed by phylogeny
ptm <- proc.time()
mod <- bam(count_sc ~ 
             
             ## First order effects ##
             # No need for a global intercept
             0 + 
             
             # Account for variation in number of routes per region per year
             n_records + 
             
             # Primary smooth of year
             s(year, bs = 'cr', k = 10) +
             
             ## Second order effects ##
             
             # Non-phylogenetic average spatiotemporal pattern
             ti(year, strata_name, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(nb = strat_penalty)),
                k = c(10, 25)) +
             
             # Phylogenetically-informed average temporal trends
             ti(year, sp_latin_phy, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(penalty = phylo_penalty)),
                k = c(10, 25)) +
             
             # Non-phylogenetic average temporal trends
             ti(year, sp_latin, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(penalty = sp_penalty)),
                k = c(10, 25)) +
             
             ## Third order effects ##
             
             # Phylogenetically-informed spatiotemporal effects
             ti(year, strata_name, sp_latin_phy,
                bs = c('cr', 'mrf', 'mrf'),
                xt = list(list(penalty = NULL),
                          list(nb = strat_penalty),
                          list(penalty = phylo_penalty)),
                k = c(10, 25, 25)) +

             # Non-phylogenetic spatiotemporal effects
             ti(year, strata_name, sp_latin,
                bs = c('cr', 'mrf', 'mrf'),
                xt = list(list(penalty = NULL),
                          list(nb = strat_penalty),
                          list(penalty = sp_penalty)),
                k = c(10, 25, 25)),
           family = gaussian(),
           data = training_data,
           method = 'fREML',
           select = TRUE,
           drop.unused.levels = FALSE,
           discrete = TRUE,
           
           # Adjust appropriately; I'm using a beefy i9 processor
           # and this model takes ~90 - 120 minutes to complete
           nthreads = 12)
runtime <- proc.time() - ptm
gc()
mod$runtime <- runtime

# Reduce model size and add draws of posterior coefficients
mod <- post_process(mod)

# Save the model object
dir.create('models')
saveRDS(mod, "./models/mod.rds")

# Now for two benchmark variants of model 1. Benchmark 1 ignores the phylogenetic 
# components but still allows for nonlinear trends
ptm <- proc.time()
mod_bench1 <- bam(count_sc ~ 
                    0 + 
                    n_records + 
                    s(year, bs = 'cr', k = 10) +
                    ti(year, strata_name, 
                       bs = c('cr', 'mrf'),
                       xt = list(list(penalty = NULL), 
                                 list(nb = strat_penalty)),
                       k = c(10, 25)) +
                    ti(year, sp_latin, 
                       bs = c('cr', 'mrf'),
                       xt = list(list(penalty = NULL), 
                                 list(penalty = sp_penalty)),
                       k = c(10, 25)) +
                    ti(year, strata_name, sp_latin,
                       bs = c('cr', 'mrf', 'mrf'),
                       xt = list(list(penalty = NULL),
                                 list(nb = strat_penalty),
                                 list(penalty = sp_penalty)),
                       k = c(10, 25, 25)),
                  family = gaussian(),
                  data = training_data,
                  method = 'fREML',
                  select = TRUE,
                  drop.unused.levels = FALSE,
                  discrete = TRUE,
                  nthreads = 12)
runtime <- proc.time() - ptm
gc()
mod_bench1$runtime <- runtime
mod_bench1 <- post_process(mod_bench1)
saveRDS(mod_bench1, "./models/mod_bench1.rds")

# The second benchmark allows for phylogenetic and non-phylogenetic 
# slopes but assumes the trend is linear
ptm <- proc.time()
mod_bench2 <- bam(count_sc ~ 
                    0 + 
                    n_records + 
                    # Feeding in large smoothing penalties will regularise
                    # smooths back to linear functions
                    s(year, bs = 'cr', k = 3, sp = .Machine$double.xmax) +
                    ti(year, strata_name, 
                       bs = c('cr', 'mrf'),
                       xt = list(list(penalty = NULL), 
                                 list(nb = strat_penalty)),
                       k = c(3, 25),
                       sp = c(.Machine$double.xmax, -1)) +
                    ti(year, sp_latin_phy, 
                       bs = c('cr', 'mrf'),
                       xt = list(list(penalty = NULL), 
                                 list(penalty = phylo_penalty)),
                       k = c(13, 25),
                       sp = c(.Machine$double.xmax, -1)) +
                    ti(year, sp_latin, 
                       bs = c('cr', 'mrf'),
                       xt = list(list(penalty = NULL), 
                                 list(penalty = sp_penalty)),
                       k = c(3, 25),
                       sp = c(.Machine$double.xmax, -1)) +
                    ti(year, strata_name, sp_latin_phy,
                       bs = c('cr', 'mrf', 'mrf'),
                       xt = list(list(penalty = NULL),
                                 list(nb = strat_penalty),
                                 list(penalty = phylo_penalty)),
                       k = c(3, 25, 25),
                       sp = c(.Machine$double.xmax, -1, -1)) +
                    ti(year, strata_name, sp_latin,
                       bs = c('cr', 'mrf', 'mrf'),
                       xt = list(list(penalty = NULL),
                                 list(nb = strat_penalty),
                                 list(penalty = sp_penalty)),
                       k = c(3, 25, 25),
                       sp = c(.Machine$double.xmax, -1, -1)),
                  family = gaussian(),
                  data = training_data,
                  method = 'fREML',
                  drop.unused.levels = FALSE,
                  discrete = TRUE,
                  nthreads = 12)
runtime <- proc.time() - ptm
gc()
mod_bench2$runtime <- runtime
mod_bench2 <- post_process(mod_bench2)
saveRDS(mod_bench2, "./models/mod_bench2.rds")

# Now a second model that uses functional relationships in place
# of phylogenetic relationships
ptm <- proc.time()
mod2 <- bam(count_sc ~ 
             0 + 
             n_records + 
             s(year, bs = 'cr', k = 10) +
             ti(year, strata_name, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(nb = strat_penalty)),
                k = c(10, 25)) +
             ti(year, sp_latin_func, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(penalty = func_penalty)),
                k = c(10, 25)) +
             ti(year, sp_latin, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(penalty = sp_penalty)),
                k = c(10, 25)) +
             ti(year, strata_name, sp_latin_func, 
                bs = c('cr', 'mrf', 'mrf'),
                xt = list(list(penalty = NULL),
                          list(nb = strat_penalty), 
                          list(penalty = func_penalty)),
                k = c(10, 25, 25)) +
             ti(year, strata_name, sp_latin, 
                bs = c('cr', 'mrf', 'mrf'),
                xt = list(list(penalty = NULL),
                          list(nb = strat_penalty), 
                          list(penalty = sp_penalty)),
                k = c(12, 25, 25)),
           family = gaussian(),
           data = training_data,
           method = 'fREML',
           select = TRUE,
           drop.unused.levels = FALSE,
           discrete = TRUE,
           nthreads = 12)
runtime <- proc.time() - ptm
gc()
mod2$runtime <- runtime
mod2 <- post_process(mod2)
saveRDS(mod2, "./models/mod2.rds")
