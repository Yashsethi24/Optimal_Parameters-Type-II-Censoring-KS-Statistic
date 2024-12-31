library('sde')  # Load the sde library for simulating stochastic processes

# Function to calculate KS statistic values
calculate_ks_statistic <- function(sample_size, iterations, true_mean, true_sd) {
  ks_stat_values <- numeric(iterations)  # Initialize vector to store KS statistic values
  
  for (iteration in 1:iterations) {
    # Generate a random sample from a normal distribution
    random_sample <- rnorm(sample_size, mean = true_mean, sd = true_sd)
    
    # Estimate variance and standard deviation from the sample
    sample_variance <- sum((random_sample - true_mean)^2) / sample_size
    estimated_sd <- sqrt(sample_variance)
    
    # Transform samples to follow Uniform(0,1) distribution using the estimated SD
    uniform_samples <- sort(pnorm(random_sample, mean = true_mean, sd = estimated_sd))
    
    # Empirical cumulative distribution function (ECDF)
    ecdf_lower <- (0:(sample_size - 1)) / sample_size  # Empirical lower bound
    ecdf_upper <- (1:sample_size) / sample_size        # Empirical upper bound
    
    # Compute KS statistics (D+ and D-)
    d_plus <- abs(ecdf_upper - uniform_samples)  # D+
    d_minus <- abs(ecdf_lower - uniform_samples)  # D-
    ks_statistic <- max(max(d_plus), max(d_minus))  # Maximum of D+ and D-
    
    # Scale KS statistic by sqrt(sample size)
    ks_stat_values[iteration] <- sqrt(sample_size) * ks_statistic
  }
  
  return(ks_stat_values)
}

# Function to simulate Kolmogorov distribution using Brownian bridge
simulate_kolmogorov_distribution <- function(sample_size, bridge_size, iterations) {
  max_deviations <- numeric(iterations)  # Initialize vector to store maximum deviations
  
  for (iteration in 1:iterations) {
    # Generate a Brownian bridge
    brownian_bridge <- BBridge(x = 0, y = 0, t0 = 0, T = 1, N = bridge_size)
    
    # Compute the maximum absolute deviation
    max_deviations[iteration] <- max(abs(brownian_bridge))
  }
  
  return(max_deviations)
}

###########################################################
# Main program starts here

# Define parameters
mean_hypothesis <- 5                       # Mean under the null hypothesis
std_dev <- 0.5                             # Standard deviation of the sample
iterations <- 10000                        # Number of iterations for simulations
final_iterations <- 2500                   # Number of iterations for quantile analysis
quantile_of_interest <- 0.95               # Quantile to analyze
sample_size <- 200                         # Sample size
bridge_size <- sample_size * 15            # Size of the Brownian bridge

# Initialize arrays to store results
ks_test_results <- array(0, dim = c(final_iterations, 2))  # KS statistic and p-value
quantile_differences <- array(0, dim = c(final_iterations, length(quantile_of_interest)))  # Quantile differences
quantile_values <- array(0, dim = c(final_iterations, 2 * length(quantile_of_interest)))  # Quantile values

#################### Quantile analysis ##########################

for (iteration in 1:final_iterations) {
  # Simulate Kolmogorov distribution values
  kolmogorov_values <- simulate_kolmogorov_distribution(sample_size, bridge_size, iterations)
  
  # Compute KS statistic values
  ks_stat_values <- calculate_ks_statistic(sample_size, iterations, mean_hypothesis, std_dev)
  
  # Calculate quantiles for KS statistics and Kolmogorov distribution
  quantile_values[iteration, 1] <- quantile(ks_stat_values, quantile_of_interest)
  quantile_values[iteration, 2] <- quantile(kolmogorov_values, quantile_of_interest)
  
  # Perform two-sample KS test
  ks_test <- ks.test(ks_stat_values, kolmogorov_values, alternative = "two.sided")
  ks_test_results[iteration, ] <- c(ks_test$statistic, ks_test$p.value)
  
  # Log progress
  cat("Iteration:", iteration, "Date:", date(), "\n")
}

# Save quantile values to a CSV file
write.csv(quantile_values, "ks_quantile_values.csv", row.names = FALSE)

# Load quantile values for further analysis
quantile_values <- read.csv("ks_quantile_values.csv", header = TRUE)

# Calculate differences between quantiles
quantile_differences[, 1] <- quantile_values[, 2] - quantile_values[, 1]

############## Plotting #################

# Histogram of KS statistics
hist(ks_test_results[, 1], probability = TRUE, main = "Histogram of KS Statistics", xlab = "KS Statistic")

# Histogram of p-values
hist(ks_test_results[, 2], probability = TRUE, main = "Histogram of KS Test P-Values", xlab = "P-Value")

# Histogram of quantile differences
hist(quantile_differences[, 1], probability = TRUE, main = "Histogram of Quantile Differences", xlab = "Quantile Difference")

# Print variance and median of quantile differences
print(var(quantile_differences[, 1]))
print(median(quantile_differences[, 1]))

################### KS test for p-values ###################

ks.test(ks_test_results[, 2], "punif")
