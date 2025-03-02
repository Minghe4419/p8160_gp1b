pacman::p_load(geepack, MASS, Matrix, ggplot2, reshape2, ggcorrplot)

# --------------------------------------------------------------------------#
##### Generating the data #####
# Function to generate Longitudinal Multivariate Poisson Data
generate_longitudinal_poisson <- function(n_subjects, n_time, rho, lambda_base, beta1, sigma_z = 1) {
  # Step 1: AR(1) correlation matrix
  Sigma <- outer(1:n_time, 1:n_time, function(i, j) rho^abs(i - j))
  Sigma <- as.matrix(nearPD(Sigma)$mat)  # Ensure positive definiteness
  
  # Step 2: Generate subject-specific latent normal effects
  subject_random_effects <- mvrnorm(n_subjects, mu = rep(0, n_time), Sigma = sigma_z^2 * Sigma)
  
  # Step 3: Assign a binary predictor (50% probability)
  subject_data <- data.frame(
    Subject = 1:n_subjects,
    X = rbinom(n_subjects, 1, 0.5)  # Binary predictor (0 or 1)
  )
  
  # Step 4: Expand into long format for longitudinal structure
  long_data <- merge(
    subject_data,
    data.frame(Time = rep(1:n_time, times = n_subjects), Subject = rep(1:n_subjects, each = n_time)),
    by = "Subject"
  )
  
  # Step 5: Assign latent normal effects to each subject-time combination
  long_data$Z <- as.vector(subject_random_effects)
  
  # Step 6: Compute Poisson mean (lambda) using a log link with the covariate
  long_data$Lambda <- lambda_base * exp(beta1 * long_data$X + long_data$Z)
  
  # Step 7: Generate Poisson counts
  long_data$Count <- rpois(n_subjects * n_time, lambda = long_data$Lambda)
  
  return(long_data)
}

simulate_ols_analysis <- function(n_subjects, n_time, rho, lambda_base, beta1, sigma_z, n_sim) {
  type1_errors <- numeric(n_sim)
  bias_values <- numeric(n_sim)
  coverage_values <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Generate data with sigma_z controlling non-normality
    sim_data <- generate_longitudinal_poisson(
      n_subjects = n_subjects,
      n_time = n_time,
      rho = rho,
      lambda_base = lambda_base,  # Fixed lambda_base
      beta1 = beta1,
      sigma_z = sigma_z  # Controls non-normality (the larger the value, the more deviation)
    )
    
    # Fit an OLS model
    lm_model <- lm(Count ~ X, data = sim_data)
    
    # Extract p-value for X
    p_value <- summary(lm_model)$coefficients["X", "Pr(>|t|)"]
    type1_errors[i] <- ifelse(p_value < 0.05, 1, 0)
    
    # Compute bias
    beta1_hat <- coef(lm_model)["X"]
    bias_values[i] <- beta1_hat - beta1
    
    # Compute 95% CI for beta1
    beta1_se <- summary(lm_model)$coefficients["X", "Std. Error"]
    beta1_CI <- beta1_hat + c(-1.96, 1.96) * beta1_se
    
    # Check coverage
    coverage_values[i] <- (beta1_CI[1] <= beta1 & beta1_CI[2] >= beta1)
  }
  
  # Return results
  data.frame(
    n_subjects = n_subjects,
    rho = rho,
    sigma_z = sigma_z,
    avg_bias = mean(bias_values),
    type1_error_rate = mean(type1_errors),
    coverage_prob = mean(coverage_values)
  )
}

# Define simulation parameters (removing non_normality_level, replaced with sigma_z)
test_conditions <- expand.grid(
  n_subjects = c(50, 100, 500, 1000),  # Sample size
  rho = c(0, 0.3, 0.6, 0.9),           # Correlation
  sigma_z = c(1, 2, 5, 10)             # Controls non-normality (the larger, the more deviation)
)

# Run the simulation
set.seed(123)
simulation_results_ols <- do.call(rbind, apply(test_conditions, 1, function(row) {
  simulate_ols_analysis(
    n_subjects = row["n_subjects"],
    n_time = 4,
    rho = row["rho"],
    lambda_base = 5,
    beta1 = 0,  # Type I Error scenario
    sigma_z = row["sigma_z"],
    n_sim = 500  # Increase number of simulations
  )
}))

# 1. Bias vs. Correlation
ggplot(simulation_results_ols, aes(x = rho, y = avg_bias, color = as.factor(sigma_z))) +
  geom_line() + 
  geom_point() +
  facet_wrap(~ n_subjects, scales = "fixed") +
  labs(
    title = "Bias in OLS Regression Across Correlation and Non-Normality Levels",
    x = "Correlation Level (rho)", 
    y = "Bias",
    color = "Non-Normality (sigma_z)"
  ) +
  theme_minimal()

# 2. Type I Error Rate vs. Correlation
ggplot(simulation_results_ols, aes(x = rho, y = type1_error_rate, color = as.factor(sigma_z))) +
  geom_line() + 
  geom_point() +
  facet_wrap(~ n_subjects, scales = "fixed") +
  labs(
    title = "Type I Error Rate in OLS Regression",
    x = "Correlation Level (rho)", 
    y = "Type I Error Rate",
    color = "Non-Normality (sigma_z)"
  ) +
  theme_minimal()

# 3. Coverage Probability vs. Correlation
ggplot(simulation_results_ols, aes(x = rho, y = coverage_prob, color = as.factor(sigma_z))) +
  geom_line() + 
  geom_point() +
  facet_wrap(~ n_subjects, scales = "fixed") +
  labs(
    title = "Coverage Probability of Confidence Intervals",
    x = "Correlation Level (rho)", 
    y = "Coverage Probability",
    color = "Non-Normality (sigma_z)"
  ) +
  theme_minimal()
###
#When sigma_z increases, the variance of the random effects increases, 
#which causes the distribution of the Poisson mean (Lambda) 
#to become more right-skewed and produces more extreme high values. 
#OLS assumes constant variance, but in reality, 
#the data exhibit heteroscedasticity (the residual variance is larger for higher Count values). 
#In this situation, the standard errors computed under the homoscedasticity assumption are overestimated.

#Impact:
#Increased Coverage: The confidence intervals become wider (beta1_hat ± 1.96 * the overestimated SE), 
#resulting in a higher probability of covering the true value (beta1 = 0).

#Decreased Type I Error Rate: Due to the overestimated standard errors, 
#the absolute value of the test statistic t = beta1_hat / SE becomes smaller, 
#which leads to larger p-values and a lower frequency of rejecting the null hypothesis.
####