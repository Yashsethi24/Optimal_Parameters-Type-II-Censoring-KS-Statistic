# Set parameters
n <- 1000000  # Number of samples
set.seed(2)   # Set seed for reproducibility

# Generate uniform random samples and sort them
samples <- runif(n, 0, 1)
sorted_samples <- sort(samples)

# Initialize vectors for empirical CDF calculations
ecp <- numeric(n)  # Empirical CDF (upper bound)
ecm <- numeric(n)  # Empirical CDF (lower bound)
cd <- numeric(n)   # Cumulative distribution
dp <- numeric(n)   # D+ values
dm <- numeric(n)   # D- values

# Compute empirical CDF and D+ and D- values
for (i in 1:n) {
  ecp[i] <- i / n                  # ECDF upper bound
  ecm[i] <- (i - 1) / n            # ECDF lower bound
  cd[i] <- sorted_samples[i]       # Sorted cumulative distribution
  dp[i] <- ecp[i] - cd[i]          # D+ calculation
  dm[i] <- ecm[i] - cd[i]          # D- calculation
}

# Compute the Kolmogorov-Smirnov statistic
dm <- abs(dm)                      # Take absolute values of D-
ks_statistic <- max(max(dm), max(dp))  # Maximum of D+ and D-
ks_statistic

# Visualize Brownian bridge and associated statistics
n <- 100000  # Update sample size
y <- runif(n, 0, 1)  # Generate uniform samples
sample_mean <- mean(y)  # Mean of samples
sample_variance <- var(y)  # Variance of samples

# Normalize the samples
x <- (y - sample_mean) / sqrt(sample_variance)
x <- cumsum(x)  # Cumulative sum
x <- x / sqrt(n)

# Generate time points
t <- seq(0, 1, length.out = n)

# Plot the Brownian bridge
plot(t, x, col = 'blue', type = 'l', main = "Brownian Bridge", xlab = "t", ylab = "Bridge")

# Calculate Brownian bridge adjustments
B <- x - t * x[n]
max_B <- max(B)  # Maximum deviation of the Brownian bridge
tmin <- which.min(B)  # Time point of minimum value
new_t <- (t[tmin] - t) %% 1  # Adjust time points
j <- floor(new_t * n)

# Calculate adjusted Brownian bridge
B_min <- min(B)
adjusted_B <- B[j] - B_min

# Plot the adjusted Brownian bridge
plot(1:n, adjusted_B, type = "l", main = "Adjusted Brownian Bridge", xlab = "Index", ylab = "Bridge Adjustment")

# Censoring example
censor_time <- 1.2  # Censoring threshold
data <- rexp(n, rate = 1)  # Generate exponential data
non_censored <- data[data < censor_time]  # Filter non-censored data

# Perform Kolmogorov-Smirnov test
ks_test_result <- ks.test(non_censored, 'pexp', alternative = 'two.sided')
ks_test_result
