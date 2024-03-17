# To what extent phylogenetic or functional relationships among species can be used leveraged to inform estimates of nonlinear changes in abundance?

This study aims to use long-term multi-species monitoring data to tackle the above question. 


## Contributors (in no particular order)
[Nicholas Clark](https://github.com/nicholasjclark)
  
[Adam Smith](https://github.com/AdamCSmithCWS)

[Shubhi Sharma](https://github.com/shubhi124081)

[Casey Youngflesh](https://github.com/caseyyoungflesh)
  
[Caleb Robbins](https://github.com/robbinscalebj)

[Hammed Akande](https://github.com/drhammed)

## Proposed methodology
1. Gather multi-species abundance (or relative abundance) measures from long-term monitoring studies
2. Construct phylogenetic and functional trees to represent relationships among species
3. Gather other appropriate information necessary to capture spatial confounding (i.e. coordinates, polygon structures etc...). The script `BBS_trends_data.R` in this repo has some annotated code to walk through such a data gathering / cleaning scheme
4. Build Generalized Additive Models (GAMs) in `mgcv` using tensor product decompositions ([see the help page on tensor products for more information](https://rdrr.io/cran/mgcv/man/te.html)) that can be used to ask how species' relationships inform estimates of nonlinear trend. Make use of the highly flexible `mrf` basis in `mgcv` to incorporate phylogenetic and functional information (see [this post from Cross Validated Gavin Simpson](https://stats.stackexchange.com/questions/638522/gam-model-with-spatial-account-via-mrf) and [this blogpost from myself](https://ecogambler.netlify.app/blog/phylogenetic-smooths-mgcv/) to get a bit more context on how these models work. The script `BBS_trends_models.R` in this repo has some example annotated code to show how these can be fit in `bam()` while also attempting to account for unmodelled temporal autocorrelation
6. Design a model evaluation scheme that allows us to compare fits from phylogenetic, functional and "null" models (that use only the random effect grouping factors of "species", but not their relationships) in a variety of ways (cross-validation by leaving certain species or groups out, with appropriate proper scoring rules; calculating trait contributions to squared second derivatives of trends; comparisons against models that assume trends are linear) 

## Tasks
### Design and justification
- [ ] Review literature to understand approaches that have been used to leverage phylogenetic or functional relationships to inform population estimates
- [ ] Gather information on the types of models / analyses that are commonly used for large multispecies datasets to inform decisions or calculations of indices (for example, how do US and Canadian Governments use NA BBS data? And could the proposed models make any impact on these pipelines?)
- [ ] Also gather information on Gaussian Markov Random Fields and their potential applications in complex nonlinear effect estimates (see for example [this work by Rue and Held](https://www.taylorfrancis.com/books/mono/10.1201/9780203492024/gaussian-markov-random-fields-havard-rue-leonhard-held) and [this post](https://haakonbakkagit.github.io/btopic120.html))

### Methodology
- [ ] Identify appropriate multi-species datasets. There is considerable information (with example code) provided by [this preprint](https://www.biorxiv.org/content/10.1101/2022.11.02.514877v1?rss=1) and the [accompanying Github repo](https://github.com/GitTFJ/correlated_effect_model)
- [ ] For candidate datasets, determine appropriate steps for cleaning and preparing data. We don't want too many shortcuts here (i.e. blindly aggregating with no justification for this), it would be better to think through the data generating process for each dataset
- [ ] Determine appropriate cross-validation schemes, considering blocking over space, time and phylogeny / functional dendrogram, to evaluate candidate models
- [ ] Prepare scoring scripts and justify scoring rules to prioritize; consider CRPS, energy and variogram scores [see this lecture on univariate forecast evaluation](https://nicholasjclark.github.io/physalia-forecasting-course/day3/lecture_4_slidedeck#1) and [this lecture on multivariate forecast evaluation](https://nicholasjclark.github.io/physalia-forecasting-course/day4/lecture_5_slidedeck#1) for context
- [ ] Brainstorm the kinds of outputs that we will need, and make sure we have well-annotated functions that can be applied to any of the models for calculating important metrics (look through the in-development functions in the `Functions/utilities.R` script in this repo for examples; and see the `BBS_trends_analysis.R` script for examples of how these might be used)
