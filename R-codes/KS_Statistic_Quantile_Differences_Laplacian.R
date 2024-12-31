library('vsgoftest')  # Library for Laplacian distribution functions
library('sde')        # Stochastic Differential Equations library (not used in this code)

# Function to calculate KS statistic values for transformed and null distributions
compute_ks_statistic <- function(sample_size, iterations, mean_param, scale_param) {
  ks_statistics <- array(0, dim = c(iterations, 2))  # Initialize array to store KS statistics
  
  for (iteration in 1:iterations) {
    # Generate unique samples from Laplacian(0, 1)
    base_sample <- unique(rlaplace(sample_size, m = 0, b = 1))
    
    # Transform sample to Laplacian(mean_param, scale_param)
    transformed_sample <- base_sample * scale_param + mean_param
    
    # Estimate parameters from the transformed sample
    estimated_mean <- mean(transformed_sample)
    estimated_scale <- sum(abs(transformed_sample - estimated_mean)) / length(transformed_sample)
    
    # Compute KS statistic for transformed sample
    ks_statistics[iteration, 1] <- sqrt(length(transformed_sample)) * 
      ks.test(transformed_sample, "plaplace", estimated_mean, estimated_scale)$statistic
    
    # Generate another sample from Laplacian(0, 1) for the null distribution
    null_sample <- unique(rlaplace(sample_size, m = 0, b = 1))
    
    # Estimate parameters from the null sample
    null_estimated_mean <- mean(null_sample)
    null_estimated_scale <- sum(abs(null_sample - null_estimated_mean)) / length(null_sample)
    
    # Compute KS statistic for the null sample
    ks_statistics[iteration, 2] <- sqrt(length(null_sample)) * 
      ks.test(null_sample, "plaplace", null_estimated_mean, null_estimated_scale)$statistic
  }
  
  return(ks_statistics)
}

# Main program begins here

# Define parameters for the Laplacian distribution and simulation
mean_param <- 7                          # Mean parameter for Laplacian distribution
scale_param <- 2                         # Scale parameter for Laplacian distribution
iterations <- 10000                      # Number of iterations per function call
final_iterations <- 100                  # Number of iterations for quantile analysis
quantile_level <- 0.95                   # Quantile level to analyze
sample_size <- 500                       # Sample size

# Initialize arrays to store results
ks_test_results <- array(0, dim = c(final_iterations, 2))  # KS test statistics and p-values
quantile_differences <- array(0, dim = c(final_iterations, 1))  # Quantile differences
quantile_values <- array(0, dim = c(final_iterations, 2))  # Quantile values for transformed and null samples

#################### Quantile Analysis ##########################

for (iteration in 1:final_iterations) {
  # Compute KS statistic values for transformed and null distributions
  ks_statistics <- compute_ks_statistic(sample_size, iterations, mean_param, scale_param)
  
  # Calculate quantiles for KS statistics
  quantile_values[iteration, 1] <- quantile(ks_statistics[, 1], quantile_level)  # Transformed sample
  quantile_values[iteration, 2] <- quantile(ks_statistics[, 2], quantile_level)  # Null sample
  
  # Perform Wilcoxon signed-rank test to compare distributions
  wilcox_test <- wilcox.test(ks_statistics[, 1], ks_statistics[, 2], alternative = "two.sided")
  ks_test_results[iteration, ] <- c(wilcox_test$statistic, wilcox_test$p.value)
  
  # Print progress
  cat("Iteration:", iteration, ";  ")
}

# Compute differences between quantiles
quantile_differences[, 1] <- quantile_values[, 2] - quantile_values[, 1]

############## Plotting #################

# Plot histogram of quantile differences
hist(quantile_differences[, 1], probability = TRUE, 
     main = "Histogram of Quantile Differences", 
     xlab = "Quantile Difference")

# Plot histogram of p-values from the Wilcoxon test
hist(ks_test_results[, 2], probability = TRUE, 
     main = "Histogram of Wilcoxon Test P-Values", 
     xlab = "P-Value")
