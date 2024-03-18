#### Fit candidate phylogenetic and functional smoothing models
# to the full dataset ####
library(mgcv)
library(MRFtools)
library(dplyr)
library(ggplot2)
library(marginaleffects)
library(patchwork)

# Load data objects
load("./data/model_objects.rda")

# Some useful details on MRF smooths:
# https://stats.stackexchange.com/questions/638522/gam-model-with-spatial-account-via-mrf

# Fit a phylogenetic trend model using a Poisson observation
# process. The model is a decomposition that includes terms capturing 
# marginal spatial and spatiotemporal fields, as well as contributions
# from species' phylogenetic relationships. Higher-order interaction terms
# allow each species' spatiotemporal field to be informed by phylogeny
mod <- bam(count ~ 
             ## First order effects ##
             
             # Offset to account for variation in number of routes
             # per region per year AND species' empirical means over
             # regions
             offset(offset) + 
             
             # Primary smooth of year
             s(year, bs = 'cr', k = 10) +
             
             # Phylogenetically-informed intercepts
             s(sp_latin_phy, 
               bs = 'mrf', 
               xt = list(penalty = phylo_penalty),
               k = 20) +
             
             # Primary average spatial field
             s(ST_12, 
               bs = 'mrf', 
               xt = list(nb = strat_penalty),
               k = 20) +
             
             ## Second order effects ##
             
             # Average spatiotemporal trend
             ti(year, ST_12, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(nb = strat_penalty)),
                k = c(10, 20)) +
             
             # Phylogenetically-informed spatial fields
             ti(ST_12, sp_latin_phy, 
                bs = c('mrf', 'mrf'),
                xt = list(list(nb = strat_penalty), 
                          list(penalty = phylo_penalty)),
                k = c(20, 20)) +
             
             # Phylogenetically-informed average temporal trends
             ti(year, sp_latin_phy, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(penalty = phylo_penalty)),
                k = c(10, 20)) +
             
             ## Third order effects ##
             
             # Phylogenetically-informed spatiotemporal effects
             ti(year, ST_12, sp_latin_phy, 
                bs = c('cr', 'mrf', 'mrf'),
                xt = list(list(penalty = NULL),
                          list(nb = strat_penalty), 
                          list(penalty = phylo_penalty)),
                k = c(10, 20, 20)),
           
           # Residual AR1 parameter to apply to working residuals;
           # we expect some autocorrelation but unfortunately cannot estimate
           # this parameter. 0.3 seems like a reasonable compromise
           rho = 0.3,
           
           # Logical variable telling bam() where each time series
           # begins, ensuring the AR1 process is applied appropriately
           # (only works if data are arranged properly)
           AR.start = mod_data$year == min(mod_data$year),
           
           # Poisson observations because we don't expect each species
           # to share the same level of overdispersion, as would be assumed
           # in a Negative Binomial or Tweedie model
           family = poisson(),
           data = mod_data,
           method = 'fREML',
           drop.unused.levels = FALSE,
           discrete = TRUE,
           
           # Adjust appropriately; I'm using a beefy i9 processor
           # and this model takes ~15 - 20 minutes to complete
           nthreads = 15)

# Add draws of posterior coefficients so this only needs to be done once, and
# then save the model object
dir.create('models')
mod$coef_posterior <- rmvn(500,
                           mu = coef(mod),
                           V = mod$Vp)
saveRDS(mod, "./models/mod.rds")

# Now a second model that uses functional relationships in place
# of phylogenetic relationships
mod2 <- bam(count ~ 
             offset(offset) + 
             s(year, bs = 'cr', k = 10) +
             s(sp_latin_func, 
               bs = 'mrf', 
               xt = list(penalty = func_penalty),
               k = 20) +
             s(ST_12, 
               bs = 'mrf', 
               xt = list(nb = strat_penalty),
               k = 20) +
             ti(year, ST_12, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(nb = strat_penalty)),
                k = c(10, 20)) +
             ti(ST_12, sp_latin_func, 
                bs = c('mrf', 'mrf'),
                xt = list(list(nb = strat_penalty), 
                          list(penalty = func_penalty)),
                k = c(20, 20)) +
             ti(year, sp_latin_func, 
                bs = c('cr', 'mrf'),
                xt = list(list(penalty = NULL), 
                          list(penalty = func_penalty)),
                k = c(10, 20)) +
             ti(year, ST_12, sp_latin_func, 
                bs = c('cr', 'mrf', 'mrf'),
                xt = list(list(penalty = NULL),
                          list(nb = strat_penalty), 
                          list(penalty = func_penalty)),
                k = c(10, 20, 20)),
           rho = 0.3,
           AR.start = mod_data$year == min(mod_data$year),
           family = poisson(),
           data = mod_data,
           method = 'fREML',
           drop.unused.levels = FALSE,
           discrete = TRUE,
           nthreads = 15)
mod2$coef_posterior <- rmvn(500,
                           mu = coef(mod2),
                           V = mod2$Vp)
saveRDS(mod2, "./models/mod2.rds")
