# Evaluate probabilistic trend forecasts against held-out data
# for competing models
library(mgcv)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(scoringRules)

# Source a few utility functions
source('Functions/utilities.R')

# Source the fitted models and data objects
load("./data/model_objects.rda")
load("./data/training_split.rds")
mod <- readRDS("./models/mod.rds")
mod_bench1 <- readRDS("./models/mod_bench1.rds")

validation_combos <- validation_data %>%
  dplyr::select(sp_latin, strata_name) %>%
  dplyr::distinct()

# Make predictions for the held-out data
newdat <- validation_data
mod_scores <- compute_fcs(model = mod, newdata = newdat)
mod_bench1_scores <- compute_fcs(model = mod_bench1, newdata = newdat)

# Compare multivariate scores from the two models
combo_dat <- data.frame(trait_combo = mod_scores$combo,
                        nontrait_combo = mod_bench1_scores$combo)
combo_mod <- lm(trait_combo ~ 0 + nontrait_combo, 
                 data = combo_dat)
summary(combo_mod)

# Scores from the trait model are between 1% and 5% smaller (better)
# compared to the nontrait model
quantile(mgcv::rmvn(1000, mu = coef(combo_mod),
                    V = vcov(combo_mod))[,1],
         probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

# Proportion of species * region combos in which the trait model
# gives a better prediction
length(which(combo_dat$trait_combo < combo_dat$nontrait_combo)) /
  NROW(combo_dat)

# Also calculate the RMSE of each species * region forecast from 
# competing models. Many ecologists will prefer this point-based evaluation
# to proper scoring rules, mostly because they don't understand proper
# scores
newdat$preds <- predict(mod, 
                        newdata = validation_data, 
                        type = 'response')
validation_data %>%
  dplyr::bind_cols(preds = predict(mod, 
                                   newdata = validation_data, 
                                   type = 'response')) %>%
  dplyr::mutate(sq_resid = (count_sc - preds) ^ 2) %>%
  dplyr::group_by(sp_latin, strata_name) %>%
  dplyr::mutate(rmse = sqrt(mean(sq_resid, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::select(sp_latin, strata_name, rmse) %>%
  dplyr::distinct() -> mod_rmse

validation_data %>%
  dplyr::bind_cols(preds = predict(mod_bench1, 
                                   newdata = validation_data, 
                                   type = 'response')) %>%
  dplyr::mutate(sq_resid = (count_sc - preds) ^ 2) %>%
  dplyr::group_by(sp_latin, strata_name) %>%
  dplyr::mutate(rmse = sqrt(mean(sq_resid, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::select(sp_latin, strata_name, rmse) %>%
  dplyr::distinct() -> mod_bench1_rmse

# Linear regression to ask if the nontrait model's RMSE tends to be
# greater than the trait model's RMSE
rmse_dat <- data.frame(trait_rmse = mod_rmse$rmse,
                       nontrait_rmse = mod_bench1_rmse$rmse)
rmse_mod <- lm(trait_rmse ~ 0 + nontrait_rmse, data = rmse_dat)
summary(rmse_mod)

# Scores from the trait model are between 0.5% and 5% smaller (better)
# compared to the nontrait model
quantile(mgcv::rmvn(1000, mu = coef(rmse_mod),
                    V = vcov(rmse_mod))[,1],
         probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

# Proportion of species * region combos in which the trait model
# gives a better prediction
length(which(mod_rmse$rmse < mod_bench1_rmse$rmse)) /
  NROW(mod_rmse)

# Plot some validation forecasts from the two models
plot_sp_trends(model = mod,
               data = validation_data,
               species = validation_combos$sp_latin[2],
               type = 'response',
               median_records = FALSE) +
  labs(title = 'Trait model')

plot_sp_trends(model = mod_bench1,
               data = validation_data,
               species = validation_combos$sp_latin[2],
               type = 'response',
               median_records = FALSE) +
  labs(title = 'Nontrait model')

plot_sp_trends(model = mod,
               data = validation_data,
               species = validation_combos$sp_latin[5],
               type = 'response',
               median_records = FALSE) +
  labs(title = 'Trait model')

plot_sp_trends(model = mod_bench1,
               data = validation_data,
               species = validation_combos$sp_latin[5],
               type = 'response',
               median_records = FALSE) +
  labs(title = 'Nontrait model')

plot_sp_trends(model = mod,
               data = validation_data,
               species = validation_combos$sp_latin[6],
               type = 'response',
               median_records = FALSE) +
  labs(title = 'Trait model')

plot_sp_trends(model = mod_bench1,
               data = validation_data,
               species = validation_combos$sp_latin[6],
               type = 'response',
               median_records = FALSE) +
  labs(title = 'Nontrait model')
