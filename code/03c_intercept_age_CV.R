################################################################################
#                                                                              #
#                               MIT License                                    #
#                  Copyright (c) 2024 Christopher Wolfe                        #
#                                                                              #
#     All code herein is subject to the above defined software license.        #
#                                                                              #
################################################################################

# The following script performs 5-fold cross validation assuming an intercept
# only model (no age). This script outputs a csv to the project directory with 
# log likelihood values for an intercept only model. 

set.seed(855)

## Import Models
full_mod <- cmdstan_model("stan_models/mvp_intercept.stan")
test_mod <- cmdstan_model("stan_models/mvp_intercept_log_lik.stan")

## Prep folds
dat$fold <- kfold_split_random(K = 5, N = nrow(dat))

## Prepare a matrix with the number of post-warmup iterations by observations
log_pd_kfold <- matrix(nrow = 4000, ncol = nrow(dat))

## Perform 5-fold cross-validation storing the log-likelihood values
for(k in 1:5){
  data_train <- list(y = (dat %>% filter(fold != k) %>% select(2:9) %>% 
                            as.matrix),
                     N = nrow(dat[dat$fold != k,]),
                     D = 8,
                     cats = 12)
  data_test <- list(y = (dat %>% filter(fold == k) %>% select(2:9) %>% 
                           as.matrix),
                    N = nrow(dat[dat$fold == k,]),
                    D = 8,
                    cats = 12)
  fit <- full_mod$sample(data = data_train, seed = 855, parallel_chains = 4, 
                         init = 0.01)
  gen_test <- test_mod$generate_quantities(fitted_params = 
                                             fit$draws(c("L","thr", "beta")), 
                                           data= data_test, parallel_chains = 4)
  log_pd_kfold[, dat$fold == k] <- gen_test$draws("lp", format = "draws_matrix")
}

## Save log-likelihood values for model comparison
write.csv(log_pd_kfold, file = "intercept_log_lik.csv", row.names = F)

######################################END#######################################
