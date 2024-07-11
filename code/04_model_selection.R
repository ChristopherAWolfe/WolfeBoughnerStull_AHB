################################################################################
#                                                                              #
#                               MIT License                                    #
#                  Copyright (c) 2024 Christopher Wolfe                        #
#                                                                              #
#     All code herein is subject to the above defined software license.        #
#                                                                              #
################################################################################

# The following script performs model selection between an age model and an 
# intercept-only model. This can be done with ordinal or continuous age. 

## Function to perform ELPD on folds
kfold <- function(log_lik_heldout)  {
  library(matrixStats)
  logColMeansExp <- function(x) {
    # should be more stable than log(colMeans(exp(x)))
    S <- nrow(x)
    colLogSumExps(x) - log(S)
  }
  # See equation (20) of @VehtariEtAl2016
  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)
  colnames(pointwise) <- "elpd"
  # See equation (21) of @VehtariEtAl2016
  elpd_kfold <- sum(pointwise)
  se_elpd_kfold <-  sqrt(ncol(log_lik_heldout) * var(pointwise))
  out <- list(
    pointwise = pointwise,
    elpd_kfold = elpd_kfold,
    se_elpd_kfold = se_elpd_kfold)
  structure(out, class = "loo")
}

## Import Log-Likelihood Values. NA values are set to -Inf -> exp(-Inf) = 0
ll_age <- read.csv("age_log_lik.csv")
ll_age[is.na(ll_age)] <- -Inf
ll_intercept <- read.csv("intercept_log_lik.csv")
ll_intercept[is.na(ll_intercept)] <- -Inf

## Get ELPD
elpd_age <- kfold(ll_age)
elpd_intercept <- kfold(ll_intercept)

## Compare Models
compare(elpd_age, elpd_intercept)

## Print Best
if(elpd_age[c(2,3)] > elpd_intercept[c(2,3)]){
  print("Age Model is Best!")
}else{
  "Intercept-Only Model is Best!"
}

######################################END#######################################
