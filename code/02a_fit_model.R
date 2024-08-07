################################################################################
#                                                                              #
#                               MIT License                                    #
#                  Copyright (c) 2024 Christopher Wolfe                        #
#                                                                              #
#     All code herein is subject to the above defined software license.        #
#                                                                              #
################################################################################

# The following script fits a multivariate cumulative probit model assuming a
# continuous predictor variable. This is best used for the H. sapiens dataset.
# This script will take some time to run. Following model fit, the model will 
# saved to the project directory.

## Import Stan model from directory
mod <- cmdstan_model("stan_models/mvp_continuous.stan")

## Prep Stan data
standat <- list(D = 8,
                N = nrow(dat),
                y = as.matrix(dat[,2:9]),
                X = dat$Age,
                cats = 12)

## Fit Model
fit <- mod$sample(data = standat,
                   seed = 1991,
                   refresh = 500,
                   chains = 4,
                   parallel_chains = 4,
                   init = 0.01, 
                   adapt_delta = 0.99, 
                   max_treedepth = 12)

fit$save_object("fitted_models/fit.rds")

######################################END#######################################
