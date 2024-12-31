library('sde')  # For stochastic differential equation simulations

# Function to calculate KS statistic values for uniform distribution
compute_ks_statistic <- function(sample_size, iterations, lower_bound, upper_bound) {
  ks_statistics <- array(0, dim = c(iterations, 2))  # Array to store KS statistics
  
  for (iteration in 1:iterations) {
    # Generate random uniform data and calculate KS statistic for scaled data
    sample_scaled <- runif(sample_size, lower_bound, upper_bound)
    estimated_b <- (sample_size * max(sample_scaled) - min(sample_scaled)) / (sample_size - 1)
    estimated_a <- (sample_size * min(sample_scaled) - max(sample_scaled)) / (sample_size - 1)
    
    ks_statistics[iteration, 1] <- sqrt(sample_size) * 
      (ks.test(sample_scaled, "punif", estimated_a, estimated_b)$statistic)
    
    # Generate another random uniform sample for null distribution
    sample_null <- runif(sample_size)
    estimated_b0 <- (sample_size * max(sample_null) - min(sample_null)) / (sample_size - 1)
    estimated_a0 <- (sample_size * min(sample_null) - max(sample_null)) / (sample_size - 1)
    
    ks_statistics[iteration, 2] <- sqrt(sample_size) * 
      (ks.test(sample_null, "punif", estimated_a0, estimated_b0)$statistic)
  }
  
  return(ks_statistics)
}

###########################################################

# Main program begins:

# Define parameters
lower_bound <- 2                     # Lower bound of the uniform distribution
upper_bound <- 6                     # Upper bound of the uniform distribution
iterations <- 10000                  # Number of iterations per function call
final_iterations <- 50               # Number of final iterations for analysis
quantile_level <- 0.95               # Quantile level to analyze
sample_size <- 200                   # Sample size

# Initialize arrays for results
ks_test_results <- array(0, dim = c(final_iterations, 2))  # KS test statistics and p-values
quantile_differences <- array(0, dim = c(final_iterations, 1))  # Quantile differences
quantile_values <- array(0, dim = c(final_iterations, 2))  # Quantile values for scaled and null data

#################### Quantile Analysis ##########################

for (iteration in 1:final_iterations) {
  # Compute KS statistic values
  ks_statistics <- compute_ks_statistic(sample_size, iterations, lower_bound, upper_bound)
  
  # Calculate quantiles for KS statistics
  quantile_values[iteration, 1] <- quantile(ks_statistics[, 1], quantile_level)  # Scaled data
  quantile_values[iteration, 2] <- quantile(ks_statistics[, 2], quantile_level)  # Null data
  
  # Perform Wilcoxon signed-rank test to compare distributions
  wilcox_test <- wilcox.test(ks_statistics[, 1], ks_statistics[, 2], alternative = "two.sided")
  ks_test_results[iteration, ] <- c(wilcox_test$statistic, wilcox_test$p.value)
  
  # Print progress
  cat("Iteration:", iteration, ";  ")
}

# Calculate differences between quantiles
quantile_differences[, 1] <- quantile_values[, 2] - quantile_values[, 1]

############## Plotting #################

# Plot histogram of p-values from the Wilcoxon test
hist(ks_test_results[, 2], probability = TRUE, 
     main = "Histogram of Wilcoxon Test P-Values", 
     xlab = "P-Value")

# Uncomment the following lines to visualize other results
# hist(ks_test_results[, 1], probability = TRUE, main = "Histogram of KS Statistics", xlab = "KS Statistic")
# hist(quantile_differences[, 1], probability = TRUE, main = "Histogram of Quantile Differences", xlab = "Quantile Difference")
