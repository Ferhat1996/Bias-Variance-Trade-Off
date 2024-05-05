#********************************Question_0************************************#

# Standard Gaussian kernel function
K <- function(u) 1/sqrt(2*pi) * exp(-u^2/2)

# Estimation routine for kernel density estimator with Gaussian kernel
kernel_density_estimator <- function(X, h, x) {
  n <- length(X)
  # Calculate kernel density estimator at point x with bandwidth h
  density <- sum(K((X - x) / h)) / (n * h)
  return(density)
}

#********************************Question_1************************************#

#************************************_i_***************************************#

# Parameters
set.seed(123) # Setting seed for reproducibility
n <- 1000
S <- 2000
x <- c(0.4, 0.99)
bandwidth_values <- seq(0.05, 3, by = 0.05)
results <- array(NA, dim = c(length(bandwidth_values), length(x), 3))

for (j in 1:length(x)) {
  for (i in 1:length(bandwidth_values)) {
    h <- bandwidth_values[i]
    # Generate S draws of the data and calculate estimator
    estimates <- replicate(S, kernel_density_estimator(rnorm(n), h, x[j]))
    # Calculate bias, standard deviation, and MSE
    bias <- mean(estimates) - dnorm(x[j])
    sd_est <- sd(estimates)
    mse <- mean((estimates - dnorm(x[j]))^2)
    results[i, j, ] <- c(bias, sd_est, mse)
  }
}

# Create data frame from results
result_df <- data.frame(
  Bandwidth = bandwidth_values,
  Bias = results[, 1, 1],
  Standard_Deviation = results[, 1, 2],
  MSE = results[, 1, 3]
)

# Print the resulting data frame
print(result_df)

#************************************_ii_***************************************#

# Plot the results for x = 0.4 and x = 0.99
par(mfrow=c(1,3))
for (j in 1:length(x)) {
  plot(bandwidth_values, results[, j, 1], type = "l", col = "red", ylab = "Bias",
       xlab = "Bandwidth", main = paste("Bias at x =", x[j]))
  plot(bandwidth_values, results[, j, 2], type = "l", col = "blue", ylab = "Standard Deviation",
       xlab = "Bandwidth", main = paste("Standard Deviation at x =", x[j]))
  plot(bandwidth_values, results[, j, 3], type = "l", col = "green", ylab = "MSE",
       xlab = "Bandwidth", main = paste("MSE at x =", x[j]))
}

#***********************************_iii_**************************************#

# Find the MSE-optimal bandwidth for both evaluation points
optimal_bandwidth <- bandwidth_values[apply(results[,,3], 2, which.min)]
cat("MSE-optimal bandwidth for x = 0.4:", optimal_bandwidth[1], "\n")
cat("MSE-optimal bandwidth for x = 0.99:", optimal_bandwidth[2], "\n")

#**********************************_iV_&_v_*************************************#

# Calculate Theoretical AMSE-optimal bandwidth for standard normal distribution at x = 0.4
Theoretical_AMSE_optimal_bandwidth <- function(x, n) {
  kappa_0 <- integrate(function(u) {(K(u))^2}, -Inf, Inf)$value  # Fourth-order kernel constant for Gaussian kernel
  f <- dnorm(x)   # True density function (standard normal in this case)
  f_second_derivative <- -x^2 * f - f  # Second derivative of the true density function
  mu_2 <- integrate(function(u) u^2 * dnorm(u), -Inf, Inf)$value  # Calculate the second raw moment of the true density function
  
  Theoretical_AMSE_bandwidth <- ((kappa_0 * f) / ((mu_2 * f_second_derivative)^2))^(1/5) * n^(-1/5)
  return(Theoretical_AMSE_bandwidth)
}

# Calculate Theoretical AMSE-optimal bandwidth for both evaluation points
Theoretical_AMSE_bandwidth_0.4 <- Theoretical_AMSE_optimal_bandwidth(0.4, n)
Theoretical_AMSE_bandwidth_0.99 <- Theoretical_AMSE_optimal_bandwidth(0.99, n)

cat("Theoretical_AMSE-optimal bandwidth for x = 0.4:", Theoretical_AMSE_bandwidth_0.4, "\n")
cat("Theoretical_AMSE-optimal bandwidth for x = 0.99:", Theoretical_AMSE_bandwidth_0.99, "\n")
