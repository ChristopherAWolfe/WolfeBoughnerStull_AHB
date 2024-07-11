################################################################################
#                                                                              #
#                               MIT License                                    #
#                  Copyright (c) 2024 Christopher Wolfe                        #
#                                                                              #
#     All code herein is subject to the above defined software license.        #
#                                                                              #
################################################################################

# The following script compares the correlation matrices between species and 
# leads to Figure 2 in the main text. This script does assume a user fits and
# saves the output from each species (e.g, fit1, fit2, fit3, fit4).

dims <- c("I1", "I2", "C", "P3", "P4", "M1", "M2", "M3")

## Import Models
### These will have to be adjusted based on a user's names
fit1 <- read.csv("fitted_models/fit1.rds") # Homo
fit2 <- read.csv("fitted_models/fit2.rds") # Pan
fit3 <- read.csv("fitted_models/fit3.rds") # Papio
fit4 <- read.csv("fitted_models/fit4.rds") # Hylobates

## Extract Posterior Samples
fit1_samples <- fit1$draws(variables = "Omega", format = "draws_matrix")
fit2_samples <-fit2$draws(variables = "Omega", format = "draws_matrix")
fit3_samples <-fit3$draws(variables = "Omega", format = "draws_matrix")
fit4_samples <-fit4$draws(variables = "Omega", format = "draws_matrix")

## Prepare data for analysis
fit1_list <- plyr::alply(fit1_samples, 1 ,function(x)
  matrix(data = x, nrow = 8,ncol = 8, dimnames = list(dims, dims)))
fit2_list <- plyr::alply(fit2_samples, 1 ,function(x)
  matrix(data = x, nrow = 8, ncol = 8, dimnames = list(dims, dims)))
fit3_list <- plyr::alply(fit3_samples, 1 ,function(x)
  matrix(data = x, nrow = 8, ncol = 8, dimnames = list(dims, dims)))
fit4_list <- plyr::alply(fit4_samples, 1 ,function(x)
  matrix(data = x, nrow = 8, ncol = 8, dimnames = list(dims, dims)))

### Perform Matrix Comparison in Parallel

## Define the Frobenius norm function
frobenius_norm <- function(mat1, mat2) {
  return(sqrt(sum((mat1 - mat2)^2)))
}

## Combine all lists into one list of lists
all_lists <- list(fit1_list, fit2_list, fit3_list, fit4_list)

## Set up parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

## Export necessary functions to the cluster
clusterExport(cl, varlist = c("frobenius_norm"))

## Function to perform pairwise comparison between two lists
pairwise_comparison <- function(list1, list2) {
  # Generate all combinations of indices
  indices <- expand.grid(i = 1:length(list1), j = 1:length(list2))
  parLapply(cl, 1:nrow(indices), function(k) {
    i <- indices[k, "i"]
    j <- indices[k, "j"]
    frobenius_norm(list1[[i]], list2[[j]])
  })
}

## Perform pairwise comparisons between all pairs of lists
system.time({
  results <- list()
  num_lists <- length(all_lists)
  for (i in 1:(num_lists-1)) {
    for (j in (i+1):num_lists) {
      cat("Comparing list", i, "with list", j, "\n")
      comparison_result <- pairwise_comparison(all_lists[[i]], all_lists[[j]])
      # Reshape the result into a matrix
      results[[paste("List", i, "vs List", j)]] <- matrix(unlist(comparison_result), nrow=length(all_lists[[i]]), ncol=length(all_lists[[j]]))
    }
  }
})

## Stop the cluster
stopCluster(cl)

## Rename to match main text
names(results) <- c("Homo-Pan", "Homo-Papio", "Homo-Hylobates", "Pan-Papio", 
                    "Pan-Hylobates", "Papio-Hylobates")

## Extract unique list names from the comparison names
list_names <- unique(unlist(strsplit(names(results), "-")))

## Create a matrix to store the mean Frobenius norms
mean_results <- matrix(NA, nrow = 4, ncol = 4)
rownames(mean_results) <- colnames(mean_results) <- list_names

## Calculate the mean Frobenius norm for each pairwise comparison
mean_frobenius <- function(result) {
  mean(unlist(result))
}

## Fill the matrix with mean Frobenius norms
for (comparison in names(results)) {
  lists <- strsplit(comparison, "-")[[1]]
  i <- which(list_names == lists[1])
  j <- which(list_names == lists[2])
  mean_results[i, j] <- mean_frobenius(results[[comparison]])
  mean_results[j, i] <- mean_results[i, j] # Symmetric matrix
}

## Convert the matrix to a long-format data frame for ggplot
mean_results_df <- melt(mean_results, na.rm = TRUE)
colnames(mean_results_df) <- c("List1", "List2", "MeanFrobeniusNorm")

## Plot
ggplot(mean_results_df, aes(x = List1, y = List2, fill = MeanFrobeniusNorm)) +
  geom_tile() +
  scale_fill_gradient2( low = "white",
                        high = "steelblue", 
                        limits = c(0, max(mean_results_df$MeanFrobeniusNorm)),
                        midpoint = 1, 
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = "black")) +
  theme_minimal() + labs(x = "", y = "", 
                         fill = "Posterior Mean Frobenius Norm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=round(MeanFrobeniusNorm,3)), size=5) + 
  scale_size(guide="none") + theme(legend.position = "bottom")

######################################END#######################################
