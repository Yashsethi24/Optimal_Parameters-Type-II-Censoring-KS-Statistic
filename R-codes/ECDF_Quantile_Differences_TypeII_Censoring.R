# Define the quantile of interest
target_quantile <- c(0.95)

# Set the number of iterations
num_iterations <- 1000

# Initialize arrays to store quantile differences for various methods
quantile_diff_complete <- array(0, dim = c(num_iterations, length(target_quantile))) # Quantile difference for complete data
quantile_diff_typeI <- array(0, dim = c(num_iterations, length(target_quantile)))    # Quantile difference for Type I censoring
quantile_diff_ratio <- array(0, dim = c(num_iterations, length(target_quantile)))    # Quantile difference using ratio (co/N)
quantile_diff_beta <- array(0, dim = c(num_iterations, length(target_quantile)))     # Quantile difference using Beta distribution
quantile_diff_func <- array(0, dim = c(num_iterations, length(target_quantile)))     # Quantile difference using function of Beta

# Initialize an array to store quartile values for all methods
quantile_values <- array(0, dim = c(num_iterations, 8 * length(target_quantile)))

# Read precomputed quantile values from a CSV file
quantile_values <- read.csv("qtval.csv")

# Loop through each iteration to compute quantile differences
for (iteration in 1:num_iterations) {
  # Compute differences for complete data
  quantile_diff_complete[iteration, ] <- quantile_values[iteration, 6] - quantile_values[iteration, 2]
  
  # Compute differences using ratio (co/N)
  quantile_diff_ratio[iteration, ] <- quantile_values[iteration, 7] - quantile_values[iteration, 3]
  
  # Compute differences using Beta distribution
  quantile_diff_beta[iteration, ] <- quantile_values[iteration, 8] - quantile_values[iteration, 4]
  
  # Compute differences using a function of Beta
  quantile_diff_func[iteration, ] <- quantile_values[iteration, 9] - quantile_values[iteration, 5]
}

# Plot the empirical cumulative distribution function (e.c.d.f) for quartile differences
plot(
  ecdf(quantile_diff_ratio[, 1]), 
  col = 2, 
  main = "E.c.d.f of Quartile Differences",
  xlab = "Quartile Difference",
  ylab = "Empirical CDF"
)  # Plot using ratio (co/N) in red

# Add e.c.d.f for Beta distribution in green
lines(ecdf(quantile_diff_beta[, 1]), col = 3)

# Add e.c.d.f for function of Beta in blue
lines(ecdf(quantile_diff_func[, 1]), col = 4)

# Add a vertical line at 0 for reference
abline(v = 0)

# Add a legend to the plot
legend(
  "bottomright", 
  legend = c("Using r/N", "Using Beta", "Using function of Beta"), 
  col = c(2, 3, 4), 
  lty = 1, 
  cex = 0.5
)

