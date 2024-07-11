################################################################################
#                                                                              #
#                               MIT License                                    #
#                  Copyright (c) 2024 Christopher Wolfe                        #
#                                                                              #
#     All code herein is subject to the above defined software license.        #
#                                                                              #
################################################################################

# The following script performs eigendecomposition of the posterior correlation
# matrix from 05_posterior_corrs.R. This script will show how to create the 
# variable loadings visualizations from the main text. Here we only present
# the 1st 2 dimensions. This can be adjusted based on goals.

## Prep: Necessary function and dimension names
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

dims <- c("I1", "I2", "C", "P3", "P4", "M1", "M2", "M3")

## Posterior Mean Correlation Matrix
fit <- readRDS("fitted_models/fit.rds")

cors <- matrix(data = fit$summary("Omega")$mean, nrow = 8, ncol = 8, 
               dimnames = list(dims, dims))

## Eigendecomposition
eig <- eigen(cors)

## Prepare eigenvectors for plotting
loadings <- as.data.frame(eig$vectors)
loadings$Variables <- dims

## Scaling the arrows for better visualization
loadings$PC1 <- loadings$V1 * sqrt(eig$values[1])
loadings$PC2 <- loadings$V2 * sqrt(eig$values[2])

## Plot
ggplot() + 
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = loadings, aes(x = PC1, y = PC2, label = Variables), 
                  hjust = 1.5, vjust =1.5, color = "red", nudge_y = 0.01, 
                  nudge_x = -0.01) +
  labs(title = "",
       x = paste("PC1 (", round(eig$values[1] / sum(eig$values) * 100, 1), "%)",
                 sep = ""),
       y = paste("PC2 (", round(eig$values[2] / sum(eig$values) * 100, 1), "%)",
                 sep = "")) +
  theme_minimal()

######################################END#######################################
