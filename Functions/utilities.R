#### Utility functions for post-processing and plotting
# from fitted models ####

#'Reduce a large gam object's size and add posterior coefficients to speed up
#'post-processing
#'@param model A fitted \code{gam} object
post_process = function(model){
  
  # Add posterior coefficient draws for faster prediction
  model$coef_posterior <- rmvn(500,
                               mu = coef(model),
                               V = model$Vp)
  
  # Reduce model footprint by removing slots that aren't needed for
  # prediction
  #model$Vc <- NULL
  #model$Ve <- NULL
  #model$Vp <- NULL
  model$db.drho <- NULL
  model$linear.predictors <- NULL
  model$fitted.values <- NULL
  model$control <- NULL
  model$var.summary <- NULL
  model$prior.weights <- NULL
  
  return(model)
}

#'Create a multivariate functional similarity dendrogram for species
#'using lists of trait datasets and use this to make the MRF functional penalty
#'
#'
#'@param trait_dfs A \code{list} containing dataframes of trait variables, where
#'rownames of each dataframe represent the species' names and columns represent trait variables.
#'@param prep_types A \code{list} describing which prepping method to apply to each dataset in
#'\code{trait_dfs}. Each prep method must match one of the following:
#'\code{'prep.fuzzy'} (for proportional variables), \code{'prep.binary'} (for binary variables),
#'or \code{'Q'} (for continuous variables). For more information about prep methods, see
#'\code{\link[ade4]{dist.ktab}}
#'@param ordering A \code{character} string of labels that should be used to reorder
#'the resulting penalty
#'
#'@seealso \code{\link[ade4]{dist.ktab}}
#'
#'@details Multiple functional trait datasets are incorporated with equal weights to
#'create a mixed-variable coefficient of distance, which is a generalized Gower's coefficient
#'of distance. This matrix is then used to create
#'a functional dendrogram with \code{hclust}, which is used to create the MRF penalty
#'
#'
preptrait_penalty = function(trait_dfs, prep_types, ordering){
  
  #### Prep the trait datasets for the Gower's distance calculations ####
  prep_dat <- lapply(seq_along(trait_dfs), function(j){
    
    if(prep_types[j] == 'prep.fuzzy') {
      ade4::prep.fuzzy(trait_dfs[[j]],
                       col.blocks = ncol(trait_dfs[[j]])) }
    
    else if(prep_types[j] == 'prep.binary') {
      ade4::prep.binary(trait_dfs[[j]],
                        col.blocks = ncol(trait_dfs[[j]])) }
    
    else if(prep_types[j] == 'Q') {
      data.frame(trait_dfs[[j]]) }
  })
  ktab1 <- ade4::ktab.list.df(prep_dat)
  
  #### Create a vector of prep types to match ade4 syntax ####
  type <- vector()
  for(i in 1:length(prep_types)){
    type[i] <- if(prep_types[i] == "prep.binary"){'B'}
    else if(prep_types[i] == "prep.fuzzy"){'F'}
    else {'Q'}
  }
  
  # Use the vector of types to create a vector of scaling options
  option <- vector()
  for(i in 1:length(prep_types)){
    option[i] <- if(prep_types[i] == "prep.binary"){'noscale'}
    else if(prep_types[i] == "prep.fuzzy"){'scaledBYrange'}
    else {'scaledBYsd'}
  }
  
  #### Calculate the Gower's distance matrix ####
  distrait <- ade4::dist.ktab(ktab1, type = type, option = option)
  
  # Create the functional dendrogram
  tree_trait <- as.dendrogram(hclust(distrait))
  
  # Create the functional trait MRF penalty and reorder to match
  # the supplied ordering
  penalty <- as.matrix(cophenetic(tree_trait))
  penalty <- penalty - max(penalty)
  diag(penalty) <- -(rowSums(penalty) - diag(penalty))
  if (!is.null(ordering)) {
    if (length(ordering) != nrow(penalty)) {
      stop("'ordering' is not the same length as the number of observations.")
    }
  } else {
    ordering <- rownames(penalty)
  }
  penalty <- penalty[order(match(rownames(penalty), ordering)),
                     order(match(colnames(penalty), ordering))]
  penalty <- MRFtools:::as_mrf_penalty(penalty, 
                                       config = MRFtools:::mrf_config(type = "dendrogram",
                                                                      dendrogram = tree_trait,
                                                                      node_labels = ordering,
                                                                      delta = 0))
  return(penalty)
}

# Compute out of sample forecast scores against new data for fitted
# gam / bam models
compute_fcs <- function(model, newdata){
  
  # Ensure data is arranged by species then by year
  newdata <- newdata %>%
    dplyr::arrange(sp_latin, year)
  
  # Generate probabilistic predictions on the response scale
  preds <- pred_intervals(model, 
                          newdata = newdata, 
                          type = 'response',
                          summarise = FALSE)
  
  es_score = function(truth, idx){
    scoringRules::es_sample(truth, 
                            t(preds[, idx]))
  }
  
  var_score = function(truth, idx){
    scoringRules::vs_sample(truth, 
                            t(preds[, idx]))
  }
  
  # For each species * region combo, compute proper forecast scores
  newdata %>%
    dplyr::mutate(idx = dplyr::row_number()) %>%
    dplyr::group_by(sp_latin, strata_name) %>%
    dplyr::mutate(energy = es_score(count_sc, idx),
                  variogram = var_score(count_sc, idx)) %>%
    dplyr::ungroup() %>%
    # Also compute an evenly-weighted combination score that considers
    # both aspects (calibration and correlation)
    dplyr::mutate(combo = log(energy) * log(variogram)) %>%
    dplyr::select(sp_latin, strata_name, energy, variogram, combo) %>%
    dplyr::distinct() -> out
  
  # For each timepoint, calculate the multivariate energy scores,
  # which penalise multivariate forecasts if they are not well 
  # calibrated (accuracy and sharpness)
  # times <- min(newdata$year) : max(newdata$year)
  # energies <- vector(length = length(times))
  # variograms <- vector(length = length(times))
  # for(t in seq_along(times)){
  #   inds <- which(newdata$year == times[t])
  #   if(length(inds)){
  #     truths <- newdata$count_sc[inds]
  #     fc_sample <- t(preds[,inds])
  #     energies[t] <- scoringRules::es_sample(y = truths,
  #                                            dat = fc_sample)
  #     variograms[t] <- scoringRules::vs_sample(y = truths,
  #                                            dat = fc_sample)
  #   } else {
  #     energies[t] <- NA
  #   }
  # }
  # 
  # Return the scores
  return(out)
}

# Plot an average temporal trend (need to Roxygenise the remaining functions
# once they are a bit further developed)
plot_av_trend = function(model, mod_data, prop_species = 0.5){
  # Subsample the species so we don't have an enormous 
  # datagrid to predict over
  n_species <- length(levels(mod_data$sp_latin))
  n_sp_sample <- ceiling(prop_species * n_species)
  sp_sample <- sample(levels(mod_data$sp_latin), n_sp_sample, replace = FALSE)
  
  # Generate the datagrid
  time_dat <- datagrid(model = mod, 
                       year = unique(mod_data$year),
                       strata_name = unique(mod_data$strata_name),
                       sp_latin = sp_sample) %>%
    dplyr::mutate(sp_latin_phy = sp_latin,
                  sp_latin_func = sp_latin)
  
  # Predict on the link scale
  preds <- pred_intervals(mod, newdata = time_dat, type = 'link')
  
  # Shift predictions by the mean
  time_dat$pred <- preds$mean
  meanpred <- mean(time_dat$pred)
  time_dat$upper <- preds$upper
  time_dat$lower <- preds$lower
  
  # Compute zero-centred average temporal trend
  time_dat %>%
    dplyr::select(year, pred, upper, lower) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise_all(.funs = function(x) 
      mean(x) - meanpred) -> time_dat
  
  # Plot
  ggplot(time_dat, aes(x = year, y = pred)) + 
    geom_line(linewidth = 1, alpha = 0.6) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, col = NA) +
    labs(y = 'Average temporal trend',
         x = 'Year') +
    theme_classic()
}

# Plot estimated trends for a species in a supplied set
# of BCRs
plot_sp_trends = function(model, 
                          data = mod_data,
                          species, 
                          regions,
                          type = 'response',
                          median_records = TRUE){
  
  # Take the first 9 regions in which this species was recorded
  # for default plotting
  if(missing(regions)){
    regions <- as.character((data %>%
      dplyr::filter(sp_latin == !!species) %>%
      dplyr::select(strata_name) %>%
      dplyr::distinct() %>%
      dplyr::pull(strata_name)))
    if(length(regions) > 9){
      regions <- regions[1:9]
    }
  }
  
  type <- match.arg(type, 
                    choices = c('expected', 'response'))
  
  get_offset <- function(model) {
    nm1 <- names(attributes(model$terms)$dataClasses)
    if('(offset)' %in% nm1) {
      deparse(as.list(model$call)$offset)
    } else {
      
      sub("offset\\((.*)\\)$", "\\1", grep('offset', nm1, value = TRUE))
    }
  }
  
  offset_name <- get_offset(model)
  if(!length(offset_name)){
    offset_name <- 'n_records'
  }

  # Use marginaleffects to create the prediction grid; if median_records = TRUE,
  # we want to ask the model what trends it would have expected if the same 
  # number of routes were taken in each region during each year. This gives a 
  # more useful 'relative' trend estimate
  if(median_records){
    newdata <- datagrid(model = model, 
                        year = unique(data$year),
                        sp_latin = species,
                        strata_name = regions) %>%
      dplyr::mutate(sp_latin_phy = species,
                    sp_latin_func = species,
                    sp_strata_name = paste0(sp_latin, '.', strata_name)) %>%
      dplyr::select(-count_sc, -{{offset_name}}) %>%
      dplyr::left_join(data) %>%
      dplyr::group_by(sp_latin, strata_name) %>%
      dplyr::mutate(!!offset_name := median(!!rlang::sym(offset_name), na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(strata_name, year) %>%
      dplyr::filter(!is.na(count))
  } else {
    # If median_records = FALSE, make predictions using the actual
    # recorded number of routes
    newdata <- datagrid(model = model, 
                        year = unique(data$year),
                        sp_latin = species,
                        strata_name = regions) %>%
      dplyr::mutate(sp_latin_phy = species,
                    sp_latin_func = species,
                    sp_strata_name = paste0(sp_latin, '.', strata_name)) %>%
      dplyr::select(-count_sc, -{{offset_name}}) %>%
      dplyr::left_join(data) %>%
      dplyr::group_by(strata_name) %>%
      dplyr::mutate(!!offset_name := ifelse(is.na(!!rlang::sym(offset_name)),
                                       median(!!rlang::sym(offset_name), na.rm = TRUE),
                                       !!rlang::sym(offset_name))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na({{offset_name}})) %>%
      dplyr::arrange(strata_name, year) %>%
      dplyr::filter(!is.na(count))
  }
  
  # Calculate prediction intervals
  preds <- pred_intervals(model, 
                          newdata = newdata, 
                          type = type)
  newdata$pred <- preds$mean
  newdata$upper <- preds$upper
  newdata$lower <- preds$lower
  
  # Make ggplot objects and modify accordingly
  p <- ggplot(newdata, aes(x = year, y = pred)) + 
    geom_line(linewidth = 1, alpha = 0.6) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, col = NA) +
    facet_wrap(~strata_name, scales = 'free') +
    theme_classic()
  
  if(!median_records){
    p <- p + geom_point(aes(y = count_sc)) +
      labs(y = paste0('Predictions for ', 
                      gsub('_', ' ', species)),
           x = 'Year')
  } else {
    p <- p + labs(y = paste0('Expected trend for ', 
                             gsub('_', ' ', species)),
                  x = 'Year')
  }
  p
}

plot_trait_conts = function(trait_derivs,
                            type = 'phylogenetic',
                            ...){
  
  if(missing(trait_derivs)){
    # Compute trait contributions to trend wiggliness
    # for each species x region combination
    trait_derivs <- compute_derivs(type = type, ...)
  }
  

  # Define a colour palatte
  rwb <- colorRampPalette(colors = c("#5c1010", 
                                     "#c30101", 
                                     "#ba3c3c",
                                     "#fbd9d3",
                                     "#e1f1fd",
                                     "darkblue"))
  

  trait_derivs %>%
    # Filter out those with very small wiggliness estimates, as 
    # the contributions don't add any meaningful information here
    # because the estimated trend is flat
    dplyr::filter(wiggliness >= quantile(wiggliness, probs = 0.15,
                                         na.rm = TRUE)) %>%
    
    # Take log of remaining trait contributions and cut into bins
    dplyr::mutate(log_conts = log(trait_conts),
                  grp = cut(log_conts, 
                            seq(min(log_conts), 
                                max(log_conts), by = 0.15), 
                            labels = FALSE, 
                            include.lowest = TRUE)) %>%
    dplyr::mutate(grp = ifelse(is.na(grp), 
                               max(grp, na.rm = TRUE), grp)) -> trait_derivs
  
  trait_derivs %>%
    dplyr::left_join(data.frame(grp = 1:max(trait_derivs$grp),
                                col = rwb(max(trait_derivs$grp)))) -> trait_derivs
  
  trait_derivs %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(brk = min(log_conts)) %>%
    dplyr::pull(brk) -> brks
  
  # Make some nice plots
  col_vec <- trait_derivs$col
  names(col_vec) <- trait_derivs$grp
  ggplot(trait_derivs, aes(log_conts)) +
    geom_histogram(breaks = brks,
                   col = NA, 
                   aes(fill = as.factor(grp)),
                   boundary = min(trait_derivs$trait_conts),
                   show.legend = FALSE) +
    scale_fill_manual(values = col_vec) +
    labs(x = paste0('log(', type, ' contribution to trend wigliness)'),
         y = 'Frequency') +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic()
}

# Use posterior simulation to make predictions from a fitted
# gam / bam model
pred_intervals = function(model, newdata, type = 'link',
                          drop_constants = FALSE,
                          summarise = TRUE){
  
  type <- match.arg(type, 
                    choices = c('expected', 'response', 'link'))
  
  # Generate the linear predictor matrix
  Xp <- predict(model, newdata = newdata,
                type = 'lpmatrix')
  
  # Drop constants if specified; this is useful to ensure the
  # returned predictions are basically zero-centred
  # (giving the types of smooth plots that mgcv tends to return
  # by default)
  if(drop_constants){
    cols_drop <- apply(Xp, 2, sd) == 0
    Xp[,cols_drop] <- 0
  }
  
  # Extract the posterior draws of betas
  betas <- model$coef_posterior
  
  # Calculate linear predictors and flatten to one long vector
  # for faster random number sampling
  all_linpreds <- as.vector(t(apply(betas, 1,
                                    function(row) Xp %*% row +
                                      attr(Xp, 'model.offset'))))
  
  # Calculate the correct type of predictions
  if(type == 'response'){
    if(grepl('Negative Binomial', model$family$family)){
      all_linpreds <- rnbinom(length(all_linpreds),
                              mu = exp(all_linpreds),
                              size = mod$family$getTheta(TRUE))
    } else if(grepl('poisson', model$family$family)){
      all_linpreds <- rpois(length(all_linpreds),
                            lambda = exp(all_linpreds))
    } else {
      all_linpreds <- rnorm(length(all_linpreds),
                            mean = all_linpreds,
                            sd = sqrt(mod$sig2))
    }
  }
  
  if(type == 'expected'){
    if(grepl('Negative Binomial', model$family$family)){
      all_linpreds <- exp(all_linpreds)
    } else if(grepl('poisson', model$family$family)){
      all_linpreds <- exp(all_linpreds)
    } else {
      all_linpreds <- all_linpreds
    }
  }
  
  # Convert predictions back to matrix form, calculate quantiles
  # and return
  all_linpreds <- matrix(all_linpreds, 
                         nrow = nrow(betas),
                         byrow = FALSE)
  if(summarise){
    creds <- apply(all_linpreds, 2, 
                   function(x) 
                     quantile(x, probs = c(0.025, 0.5, 0.975),
                              na.rm = TRUE))
    out <- list(mean = creds[2,],
                upper = creds[3,],
                lower = creds[1,])
  } else {
    out <- all_linpreds
  }

  gc()
  return(out)
}

# Compute variation of estimated trends for each 
# species / region combination and calculate how much the trait information
# contributes to these estimates
compute_derivs = function(data, 
                          model, 
                          type = 'phylogenetic', 
                          n_cores = 4,
                          pred_vars = TRUE){
  n_species <- nlevels(data$sp_latin)
  term_of_interest <- switch(type,
                             'phylogenetic' = 'phy',
                             'functional' = 'func')
  
  # First create the newdata grids so we don't have to send all of the 
  # data to each cluster; we want to predict at the regions where this species
  # has been recorded, and use a fine grid of "years" to allow for more precise estimates
  # of second derivatives with finite differencing
  which_STs <- data %>%
    dplyr::group_by(sp_latin, strata_name) %>%
    dplyr::mutate(recorded = dplyr::case_when(
      any(count > 0) ~ 1,
      TRUE ~ 0
    )) %>%
    dplyr::ungroup() %>%
    dplyr::filter(recorded == 1) %>%
    dplyr::group_by(sp_latin) %>%
    dplyr::summarise(which_ST = list(unique(strata_name)))
  years <- seq(from = min(data$year), 
               to = max(data$year), 
               by = 0.25)
  species <- levels(data$sp_latin)
  
  # Send necessary objects and libraries to each node on the cluster; this function
  # is slow, partly because of the memory overhead to send the huge model object
  # to each node; but I can't think of a better way to do this at present. As a start
  # we can reduce the size of the model object by eliminating things that aren't needed
  # for point-based predictions
  model$residuals <- NULL
  model$coef_posterior <- NULL
  model$linear.predictors <- NULL
  model$fitted.values <- NULL
  model$control <- NULL
  model$var.summary <- NULL
  model$prior.weights <- NULL
  model$R <- NULL
  model$offset <- NULL
  model$Vp <- NULL
  model$Vc <- NULL
  model$Ve <- NULL
  model$edf <- NULL

  cl <- parallel::makePSOCKcluster(n_cores)
  parallel::setDefaultCluster(cl)
  pbapply::pboptions(type = "none")
  parallel::clusterExport(NULL, c('years',
                                  'species',
                                  'model',
                                  'which_STs',
                                  'term_of_interest',
                                  'pred_vars'),
                          envir = environment())
  parallel::clusterEvalQ(cl = cl, library(marginaleffects))
  parallel::clusterEvalQ(cl = cl, library(dplyr))
  parallel::clusterEvalQ(cl = cl, library(mgcv))
  
  # Create newdata grids for each species and send them to each node;
  # using datagrid() from marginaleffects will ensure that the sampling effort
  # (offset variable) will be the same in each year for a given species, allowing
  # us to get a better sense of the overall nonlinear expected trend
  newdat_grids <- pbapply::pblapply(seq_len(n_species), function(x){
    datagrid(sp_latin = species[x],
             strata_name = which_STs$which_ST[x][[1]],
             year = years,
             model = model) %>%
      dplyr::mutate(sp_latin_phy = sp_latin,
                    sp_latin_func = sp_latin)
  }, cl = cl)
  parallel::clusterExport(NULL, c('newdat_grids',
                                  'n_species'),
                          envir = environment())
  
  # Now loop over all species and calculate wiggliness contributions
  deriv_conts <- do.call(rbind, 
                         pbapply::pblapply(
                           seq_len(n_species), 
                           function(x){
                             
                             # Predict from the model and return the specific 
                             # contributions from each term
                             sp_preds <- predict(model, 
                                                 type = 'terms', 
                                                 newdata = newdat_grids[[x]])
                             trait_terms <- grepl(term_of_interest, 
                                                  colnames(sp_preds))
                             
                             # For each region-specific time series, 
                             # compute either the squared second derivative of
                             # trends composed of trait and non-trait term predictions (if pred_vars = FALSE)
                             # or the variances of the trait and non-trait predictions
                             unique_STs <- unique(newdat_grids[[x]]$strata_name)
                             wiggliness <- vector(length = length(unique_STs))
                             trait_conts <- vector(length = length(unique_STs))
                             nontrait_conts <- vector(length = length(unique_STs))
                             for(i in 1:length(unique_STs)){
                               inds <- which(newdat_grids[[x]]$strata_name == 
                                               unique_STs[i])
                               if(!pred_vars){
                                 trait_conts[i] <- sum(diff(diff(rowSums(sp_preds[inds, 
                                                                                  trait_terms,
                                                                                  drop = FALSE]))) ^ 2)
                                 nontrait_conts[i] <- sum(diff(diff(rowSums(sp_preds[inds, 
                                                                                     !trait_terms,
                                                                                     drop = FALSE]))) ^ 2)
                                 wiggliness[i] <- sum(diff(diff(rowSums(sp_preds[inds, ]))) ^ 2)
                               } else {
                                 trait_conts[i] <- var(rowSums(sp_preds[inds, 
                                                                        trait_terms,
                                                                        drop = FALSE]))
                                 nontrait_conts[i] <- var(rowSums(sp_preds[inds, 
                                                                           !trait_terms,
                                                                           drop = FALSE]))
                                 wiggliness[i] <- var(rowSums(sp_preds[inds, ]))
                               }
                               }
                             
                             data.frame(sp_latin = levels(data$sp_latin)[x],
                                        strata_name = unique_STs,
                                        wiggliness = wiggliness,
                                        trait_conts = trait_conts / nontrait_conts)
                           }, cl = cl))
  parallel::stopCluster(cl)
  
  return(deriv_conts)
}
