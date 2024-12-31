library('sde')  # Load library for stochastic differential equations

# Function to compute KS statistic values for complete data
compute_ks_statistic <- function(sample_size, iterations) {
  ks_statistics <- numeric(iterations)  # Initialize vector to store KS statistics
  
  for (iteration in 1:iterations) {
    # Generate sorted uniform random data
    data <- sort(runif(sample_size, 0, 1))
    
    # Empirical CDF bounds
    ecdf_lower <- ((1:sample_size) - 1) / sample_size
    ecdf_upper <- (1:sample_size) / sample_size
    
    # Compute KS statistics
    d_minus <- abs(data - ecdf_lower)  # D-
    d_plus <- abs(data - ecdf_upper)  # D+
    ks_statistics[iteration] <- sqrt(sample_size) * max(max(d_minus), max(d_plus))
  }
  
  return(ks_statistics)
}

# Function to compute Kolmogorov distribution using Brownian bridge
compute_kolmogorov_distribution <- function(sample_size, bridge_size, iterations) {
  kolmogorov_values <- numeric(iterations)  # Initialize vector for Kolmogorov values
  
  for (iteration in 1:iterations) {
    # Simulate a Brownian bridge
    brownian_bridge <- BBridge(x = 0, y = 0, t0 = 0, T = 1, N = bridge_size)
    kolmogorov_values[iteration] <- max(abs(brownian_bridge))  # Max deviation
  }
  
  return(kolmogorov_values)
}

# Main program begins:

# Define parameters
sample_size <- 200          # Sample size
iterations <- 10000         # Number of iterations for each function call
bridge_size <- sample_size * 12  # Grid size for Brownian bridge
final_iterations <- 25      # Number of final iterations
quantile_level <- 0.95      # Quantile to analyze

# Compute Kolmogorov and KS statistic values
kolmogorov_values <- compute_kolmogorov_distribution(sample_size, bridge_size, iterations)
ks_test_values <- compute_ks_statistic(sample_size, iterations)

# Initialize arrays for results
ks_test_results <- array(0, dim = c(final_iterations, 2))  # KS test statistics and p-values
quantile_differences <- array(0, dim = c(final_iterations, 1))  # Quantile differences
quantile_values <- array(0, dim = c(final_iterations, 2))  # Quantile values for KS and Kolmogorov data

for (iteration in 1:final_iterations) {
  # Recompute Kolmogorov and KS statistic values
  kolmogorov_values <- compute_kolmogorov_distribution(sample_size, bridge_size, iterations)
  ks_test_values <- compute_ks_statistic(sample_size, iterations)
  
  # Compute quantiles for both distributions
  quantile_values[iteration, 1] <- quantile(ks_test_values, quantile_level)       # KS statistic quantile
  quantile_values[iteration, 2] <- quantile(kolmogorov_values, quantile_level)   # Kolmogorov quantile
  
  # Perform KS test to compare distributions
  ks_test <- ks.test(ks_test_values, kolmogorov_values, alternative = "two.sided")
  ks_test_results[iteration, ] <- c(ks_test$statistic, ks_test$p.value)
  
  # Print progress
  cat("Iteration:", iteration, "Date:", date(), "\n")
}

# Plotting the results
hist(kolmogorov_values, probability = TRUE, 
     main = "Histogram of Kolmogorov Values", 
     xlab = "Kolmogorov Value")

# Compare one-sided KS statistic between Kolmogorov and KS values
final_ks_test <- ks.test(kolmogorov_values, ks_test_values, alternative = "two.sided")
print(final_ks_test)
