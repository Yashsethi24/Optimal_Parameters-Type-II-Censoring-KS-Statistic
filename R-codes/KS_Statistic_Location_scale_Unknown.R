library('sde')  # For stochastic differential equation simulations

# Function for calculating KS statistic values
compute_ks_statistic <- function(sample_size, iterations, mean_shift, sd_estimate) {
  ks_statistics <- array(0, dim = c(iterations, 2))  # Array to store KS statistics
  
  for (iteration in 1:iterations) {
    # Generate exponential data and apply a shift
    sample_scaled <- rexp(sample_size) + mean_shift  # Scaled data
    estimated_mean <- mean(sample_scaled)           # Estimate mean
    ks_statistics[iteration, 1] <- sqrt(sample_size) * 
      (ks.test(sample_scaled, "pexp", 1 / estimated_mean)$statistic)
    
    # Generate another exponential sample for null distribution
    sample_null <- rexp(sample_size)
    ks_statistics[iteration, 2] <- sqrt(sample_size) * 
      (ks.test(sample_null, "pexp", 1 / mean(sample_null))$statistic)
  }
  
  return(ks_statistics)
}

###########################################################

# Main program begins:

# Define parameters
mean_shift <- 3                    # Mean shift for the hypothesis
sd_estimate <- 2                   # Standard deviation estimate
iterations <- 10000                # Number of iterations for each set
final_iterations <- 10             # Number of final iterations for analysis
quantile_level <- 0.95             # Quantile level to analyze
sample_size <- 200                 # Sample size
bridge_size <- sample_size * 12    # Bridge size for Brownian bridge

# Initialize arrays for results
ks_test_results <- array(0, dim = c(final_iterations, 2))  # KS test statistics and p-values
quantile_differences <- array(0, dim = c(final_iterations, 1))  # Quantile differences
quantile_values <- array(0, dim = c(final_iterations, 2))  # Quantile values for scaled and null data

#################### Quantile Analysis ##########################

for (iteration in 1:final_iterations) {
  # Compute KS statistic values
  ks_statistics <- compute_ks_statistic(sample_size, iterations, mean_shift, sd_estimate)
  
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

############## Testing for a single sample #################

# Generate a single sample of exponential data
n <- 1000
sample_scaled <- rexp(n) + mean_shift        # Scaled data
estimated_mean <- mean(sample_scaled)       # Estimate mean
ks_test_scaled <- ks.test(sample_scaled, "pexp", 1 / estimated_mean)

# Generate null data
sample_null <- rexp(n)
ks_test_null <- ks.test(sample_null, "pexp", 1 / mean(sample_null))

# Display results of the KS tests
cat("\nKS Test for Scaled Data:\n", ks_test_scaled)
cat("\nKS Test for Null Data:\n", ks_test_null)
