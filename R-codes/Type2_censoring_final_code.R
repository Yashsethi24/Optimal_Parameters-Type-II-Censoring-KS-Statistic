# Function for calculating KS statistic values under Type II censoring
Type2censor <- function(sample_size, iterations) {
  ks_stat_values <- numeric(iterations)  # Initialize vector to store KS statistics
  
  for (iteration in 1:iterations) {
    # Generate and sort data
    data <- sort(runif(sample_size, 0, 1))  # Generate uniform random data and sort it
    
    # Apply Type II censoring (retain only the first `censor_obj` observations)
    censor_obj <- 40  # Number of non-censored observations
    non_censored_data <- data[1:censor_obj]  # Retain non-censored data
    
    # Compute empirical CDF and KS statistic components
    l <- length(non_censored_data)          # Number of non-censored observations
    ecm <- ((1:l) - 1) / sample_size        # Empirical CDF lower bound
    ecp <- (1:l) / sample_size              # Empirical CDF upper bound
    dp <- abs(non_censored_data - ecp)      # D+ values
    dm <- abs(non_censored_data - ecm)      # D- values
    
    # Calculate KS statistic for this iteration
    ks_stat_values[iteration] <- sqrt(sample_size) * max(max(dp), max(dm))
  }
  
  return(ks_stat_values)  # Return KS statistics for all iterations
}

###########################################################

# Main program begins:

# Define parameters
sample_size <- 100           # Sample size
iterations <- 5000           # Number of iterations
bridge_size <- 1000          # Grid size for Brownian bridge simulation

# Compute KS statistic values under Type II censoring
ks_test_values <- Type2censor(sample_size, iterations)

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
