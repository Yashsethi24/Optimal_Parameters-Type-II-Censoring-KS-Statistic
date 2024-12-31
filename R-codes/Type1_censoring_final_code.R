# Function for calculating KS statistic values under Type I censoring
Type1censor <- function(sample_size, iterations) {
  ks_stat_values <- numeric(iterations)  # Initialize vector to store KS statistics
  
  for (iteration in 1:iterations) {
    # Generate data and apply censoring
    data <- runif(sample_size, 0, 1)         # Generate uniform random data
    censor_time <- 0.4                      # Define censoring threshold
    non_censored_data <- sort(data[data < censor_time])  # Retain non-censored data
    
    # Compute empirical CDF and KS statistic components
    n <- sample_size                        # Total sample size
    l <- length(non_censored_data)          # Number of non-censored observations
    ecm <- ((1:l) - 1) / n                  # Empirical CDF lower bound
    ecp <- (1:l) / n                        # Empirical CDF upper bound
    dp <- abs(non_censored_data - ecp)      # D+ values
    dm <- abs(non_censored_data - ecm)      # D- values
    d3 <- abs(censor_time - l / n)          # Difference between censor time and ECDF
    
    # Calculate KS statistic for this iteration
    ks_stat_values[iteration] <- sqrt(n) * max(max(dp), max(dm), d3)
  }
  
  return(ks_stat_values)  # Return KS statistics for all iterations
}

###########################################################

# Main program begins:

# Define parameters
sample_size <- 100           # Sample size
iterations <- 5000           # Number of iterations
bridge_size <- 1000          # Grid size for Brownian bridge simulation

# Compute KS statistic values under Type I censoring
ks_test_values <- Type1censor(sample_size, iterations)

# Compute and display quantiles of the KS statistic values
ks_quantiles <- quantile(ks_test_values)
cat("KS statistic quantiles:\n", ks_quantiles, "\n")

# Calculate the upper 3% quantile using a normal approximation
upper_quantile <- qnorm(0.03, mean(ks_test_values), sqrt(var(ks_test_values)), lower.tail = FALSE)
cat("Upper 3% quantile (normal approximation):", upper_quantile, "\n")

# Simulate Kolmogorov distribution using Brownian bridge
brownian_values <- numeric(iterations)
for (i in 1:iterations) {
  # Generate a Brownian bridge and compute maximum absolute deviation
  brownian_values[i] <- max(abs(BBridge(x = 0, y = 0, t0 = 0, T = 1, N = bridge_size)))
}

# Display Kolmogorov distribution quantiles
brownian_quantiles <- quantile(brownian_values)
cat("Kolmogorov distribution quantiles (Brownian bridge):\n", brownian_quantiles, "\n")
