library('sde')

# Function to calculate Kolmogorov-Smirnov (KS) statistic values
compute_ks_statistic <- function(sample_size, iterations, mean_value, std_dev) {
  ks_stat_values <- numeric(iterations)  # Initialize vector to store KS statistic values
  
  for (iteration in 1:iterations) {
    
    # Generate random samples from a normal distribution
    random_sample <- rnorm(sample_size, mean = mean_value, sd = std_dev)
    
    # Estimate the mean from the sample
    estimated_mean <- mean(random_sample)
    
    # Transform samples to follow Uniform(0,1) distribution
    uniform_samples <- sort(pnorm(random_sample, mean = estimated_mean, sd = std_dev))
    
    # Calculate empirical cumulative distribution functions (ECDF)
    ecdf_lower <- (0:(sample_size - 1)) / sample_size  # Empirical lower bound
    ecdf_upper <- (1:sample_size) / sample_size        # Empirical upper bound
    
    # Compute KS statistics (D+ and D-)
    d_plus <- abs(ecdf_upper - uniform_samples)       # D+
    d_minus <- abs(ecdf_lower - uniform_samples)      # D-
    
    # Overall KS statistic is the maximum of D+ and D-
    ks_statistic <- max(max(d_plus), max(d_minus))
    
    # Scale KS statistic by sqrt(sample size)
    ks_stat_values[iteration] <- sqrt(sample_size) * ks_statistic
  }
  
  return(ks_stat_values)
}

# Function to simulate Kolmogorov distribution using Brownian bridge
simulate_kolmogorov_distribution <- function(sample_size, bridge_size, iterations) {
  max_deviation <- numeric(iterations)  # Initialize vector to store maximum deviations
  
  for (iteration in 1:iterations) {
    # Generate a Brownian bridge
    brownian_bridge <- BBridge(x = 0, y = 0, t0 = 0, T = 1, N = bridge_size)
    
    # Compute absolute values of the Brownian bridge
    abs_bridge <- abs(brownian_bridge)
    
    # Store the maximum absolute deviation
    max_deviation[iteration] <- max(abs_bridge)
  }
  
  return(max_deviation)
}

# Main program starts here:

# Initialize parameters
mean_hypothesis <- 5                       # Mean under null hypothesis
std_dev <- 0.5                             # Standard deviation (known or estimated)
iterations <- 10000                        # Number of iterations for simulation
final_iterations <- 2500                   # Final number of iterations for quantile analysis
quantile_of_interest <- c(0.95)            # Quantile to analyze
sample_size <- 200                         # Sample size
bridge_size <- sample_size * 9             # Brownian bridge size multiplier

# Initialize arrays to store results
ks_test_results <- array(0, dim = c(final_iterations, 2))  # KS statistic and p-values
quantile_diff <- array(0, dim = c(final_iterations, length(quantile_of_interest)))  # Quantile differences
quantile_values <- array(0, dim = c(final_iterations, 2 * length(quantile_of_interest)))  # KS and Kolmogorov quantiles

# Perform quantile analysis
for (iteration in 1:final_iterations) {
  # Simulate Kolmogorov distribution values
  kolmogorov_values <- simulate_kolmogorov_distribution(sample_size, bridge_size, iterations)
  
  # Compute KS statistic values
  ks_stat_values <- compute_ks_statistic(sample_size, iterations, mean_hypothesis, std_dev)
  
  # Calculate quantiles for KS statistics and Kolmogorov distribution
  quantile_values[iteration, 1] <- quantile(ks_stat_values, quantile_of_interest)
  quantile_values[iteration, 2] <- quantile(kolmogorov_values, quantile_of_interest)
  
  # Perform two-sample KS test between the two distributions
  ks_test <- ks.test(ks_stat_values, kolmogorov_values, alternative = 'two.sided')
  ks_test_results[iteration, ] <- c(ks_test$statistic, ks_test$p.value)
  
  # Print progress
  cat("Iteration:", iteration, "Date:", date(), "\n")
}

# Write quantile values to CSV
#write.csv(quantile_values, "ks_quantile_values.csv", row.names = FALSE)

# Load quantile values for further analysis
quantile_values <- read.csv("ks_quantile_values.csv", header = TRUE)

# Calculate differences between quantiles
quantile_diff[, 1] <- quantile_values[, 2] - quantile_values[, 1]

# Plot histograms for analysis
hist(ks_test_results[, 1], probability = TRUE, main = "Histogram of KS Statistics", xlab = "KS Statistic")
hist(ks_test_results[, 2], probability = TRUE, main = "Histogram of KS Test P-Values", xlab = "P-Value")
hist(quantile_diff[, 1], probability = TRUE, main = "Histogram of Quantile Differences", xlab = "Quantile Difference")

# Uncomment the following lines to print variance and median of quantile differences
# print(var(quantile_diff[, 1]))
# print(median(quantile_diff[, 1]))

# Save KS test results to CSV
#write.csv(ks_test_results, "ks_test_results.csv", row.names = FALSE)
#write.csv(quantile_values, "ks_quantile_differences.csv", row.names = FALSE)
