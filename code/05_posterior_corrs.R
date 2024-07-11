################################################################################
#                                                                              #
#                               MIT License                                    #
#                  Copyright (c) 2024 Christopher Wolfe                        #
#                                                                              #
#     All code herein is subject to the above defined software license.        #
#                                                                              #
################################################################################

# The following script imports the desired model of interest, summarizes the
# posterior information, and demonstrates how to visualize correlation matrices
# from the main text.

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

## Plot
ggplot(melt(get_upper_tri(cors)
            , na.rm = TRUE), aes(Var1, Var2, fill=value)) + 
  geom_tile(height=0.8, width=0.8, color = "black") + 
  scale_fill_gradient2(low="#FFF5F0" , mid="#FCBBA1", high="#99000D", 
                       midpoint = 0.5, limits = c(0,1)) + 
  theme_minimal() + coord_equal() + labs(x="",y="",fill="Corr") + 
  theme(axis.text.x=element_text(size=7, angle=45, vjust=1, hjust=1, 
                                 margin=margin(-3,0,0,0), face = "bold"), 
        axis.text.y=element_text(size=7, margin=margin(0,-3,0,0), 
                                 face = "bold"), 
        panel.grid.major=element_blank(), 
        panel.grid = element_line(color="black")) + 
  geom_text(aes(label=round(value,2)), size=3) + scale_size(guide="none")

######################################END#######################################
