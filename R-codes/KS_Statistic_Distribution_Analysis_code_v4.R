library('sde')  # For stochastic differential equations and Brownian bridges

# Function to calculate KS statistic values
compute_ks_statistic <- function(sample_size, iterations, mean, std_dev) {
  ks_statistics <- array(0, dim = c(iterations, 2))  # Array to store KS statistics
  
  for (iteration in 1:iterations) {
    # Generate scaled data from N(mean, std_dev)
    data_scaled <- rnorm(sample_size) * std_dev + mean
    estimated_variance <- var(data_scaled)         # Estimate variance
    estimated_sd <- sqrt(estimated_variance)       # Estimate standard deviation
    
    # Compute KS statistic for scaled data
    ks_statistics[iteration, 1] <- sqrt(sample_size) * 
      ks.test(data_scaled, "pnorm", mean, estimated_sd)$statistic
    
    # Generate null distribution data from N(0, 1)
    data_null <- rnorm(sample_size)
    ks_statistics[iteration, 2] <- sqrt(sample_size) * 
      ks.test(data_null, "pnorm", 0, sqrt(var(data_null)))$statistic
  }
  
  return(ks_statistics)
}

# Function to calculate Kolmogorov distribution values using Brownian bridge
compute_kolmogorov_distribution <- function(bridge_size, iterations) {
  bb_max_values <- numeric(iterations)  # Initialize vector to store maximum values
  
  for (iteration in 1:iterations) {
    # Simulate a Brownian bridge
    brownian_bridge <- BBridge(x = 0, y = 0, t0 = 0, T = 1, N = bridge_size)
    bb_max_values[iteration] <- max(abs(brownian_bridge))  # Maximum absolute deviation
  }
  
  return(bb_max_values)
}

###########################################################

# Main program begins:

# Define parameters
mean <- 3                    # Mean of the normal distribution
std_dev <- 2                 # Standard deviation
iterations <- 100000         # Number of iterations for each test
final_iterations <- 50       # Number of iterations for quantile analysis
quantile_level <- 0.95       # Quantile to analyze
sample_size <- 200           # Sample size
bridge_size <- sample_size * 12  # Bridge size for Brownian bridge

# Initialize arrays for results
ks_test_results <- array(0, dim = c(final_iterations, 2))  # KS test statistics and p-values
quantile_differences <- array(0, dim = c(final_iterations, 1))  # Quantile differences
quantile_values <- array(0, dim = c(final_iterations, 3))  # Quantile values for scaled and null data

#################### Quantile Analysis ##########################

for (iteration in 1:final_iterations) {
  # Compute KS statistic values
  ks_statistics <- compute_ks_statistic(sample_size, iterations, mean, std_dev)
  
  # Compute quantiles for KS statistics
  quantile_values[iteration, 1] <- quantile(ks_statistics[, 1], quantile_level)  # Scaled data
  quantile_values[iteration, 3] <- quantile(ks_statistics[, 2], quantile_level)  # Null data
  
  # Perform Wilcoxon signed-rank test to compare distributions
  wilcox_test <- wilcox.test(ks_statistics[, 1], ks_statistics[, 2], alternative = "two.sided")
  ks_test_results[iteration, ] <- c(wilcox_test$statistic, wilcox_test$p.value)
  
  # Print progress
  cat("Iteration:", iteration, ";  ")
}

# Calculate differences between quantiles
quantile_differences[, 1] <- quantile_values[, 3] - quantile_values[, 1]

############## Plotting #################

# Plot histogram of p-values
hist(ks_test_results[, 2], probability = TRUE, 
     main = "Histogram of Wilcoxon Test P-Values", 
     xlab = "P-Value")

# Plot histogram of quantile differences
hist(quantile_differences[, 1], probability = TRUE, 
     main = "Histogram of Quantile Differences", 
     xlab = "Quantile Difference")

# Uncomment the following lines to calculate and print variance and median
# print(var(quantile_differences[, 1]))
# print(median(quantile_differences[, 1]))
