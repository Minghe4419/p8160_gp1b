pacman::p_load(geepack, MASS, Matrix, ggplot2, reshape2, ggcorrplot)

# --------------------------------------------------------------------------#
##### Generating the data #####
generate_longitudinal_poisson <- function(n_subjects, n_time, rho, lambda_base, beta1, sigma_z = 1) {
  # Step 1: AR(1) correlation matrix
  Sigma <- outer(1:n_time, 1:n_time, function(i, j) rho^abs(i - j))
  Sigma <- as.matrix(Matrix::nearPD(Sigma)$mat)
  
  # Step 2: Generate latent effects with adjusted variance
  subject_random_effects <- MASS::mvrnorm(n_subjects, mu = rep(0, n_time), Sigma = Sigma * sigma_z^2)
  
  # Step 3: Assign binary predictor
  subject_data <- data.frame(Subject = 1:n_subjects, X = rbinom(n_subjects, 1, 0.5))
  
  # Step 4: Expand to long format
  long_data <- merge(subject_data, expand.grid(Subject = 1:n_subjects, Time = 1:n_time), by = "Subject")
  
  # Step 5: Assign latent effects
  long_data$Z <- as.vector(subject_random_effects)
  
  # Step 6: Compute Poisson mean (lambda_base固定)
  long_data$Lambda <- lambda_base * exp(beta1 * long_data$X + long_data$Z)
  
  # Step 7: Generate Poisson counts
  long_data$Count <- rpois(nrow(long_data), lambda = long_data$Lambda)
  
  return(long_data)
}



simulate_gee_analysis <- function(n_subjects, n_time, rho, lambda_base, beta1, sigma_z, n_sim) {
  type1_errors <- numeric(n_sim)
  bias_values <- numeric(n_sim)
  coverage_values <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Generate data with sigma_z controlling latent variability
    sim_data <- generate_longitudinal_poisson(
      n_subjects = n_subjects,
      n_time = n_time,
      rho = rho,
      lambda_base = lambda_base,
      beta1 = beta1,
      sigma_z = sigma_z
    )
    
    # Fit GEE model
    gee_model <- geeglm(
      Count ~ X,
      id = Subject,
      data = sim_data,
      family = poisson,
      corstr = "ar1"
    )
    
    # Extract p-value for predictor X
    p_value <- summary(gee_model)$coefficients["X", "Pr(>|W|)"]
    type1_errors[i] <- ifelse(p_value < 0.05, 1, 0)
    
    # Compute bias
    beta1_hat <- coef(gee_model)["X"]
    bias_values[i] <- beta1_hat - beta1
    
    # Compute 95% CI for beta1
    beta1_se <- summary(gee_model)$coefficients["X", "Std.err"]
    beta1_CI <- beta1_hat + c(-1.96, 1.96) * beta1_se
    
    # Check if the CI covers the true beta1
    coverage_values[i] <- (beta1_CI[1] <= beta1 & beta1_CI[2] >= beta1)
  }
  
  # Return results
  return(data.frame(
    n_subjects = n_subjects,
    rho = rho,
    sigma_z = sigma_z,
    avg_bias = mean(bias_values),
    type1_error_rate = mean(type1_errors),
    coverage_prob = mean(coverage_values)
  ))
}

# 定义参数网格
test_conditions <- expand.grid(
  n_subjects = c(50, 100, 500, 1000),
  rho = c(0, 0.3, 0.6, 0.9),
  sigma_z = c(1, 2, 5)  # 控制非正态性
)

# 固定参数
lambda_base <- 5
beta1_null <- 0
n_time <- 4
n_sim <- 500

# 运行GEE模拟
set.seed(123)
simulation_results_gee <- do.call(rbind, apply(test_conditions, 1, function(row) {
  simulate_gee_analysis(
    n_subjects = row["n_subjects"],
    n_time = n_time,
    rho = row["rho"],
    lambda_base = lambda_base,
    beta1 = beta1_null,
    sigma_z = row["sigma_z"],
    n_sim = n_sim
  )
}))

# 运行OLS模拟（使用修改后的函数）
simulation_results_ols <- do.call(rbind, apply(test_conditions, 1, function(row) {
  simulate_ols_analysis(
    n_subjects = row["n_subjects"],
    n_time = n_time,
    rho = row["rho"],
    lambda_base = lambda_base,
    beta1 = beta1_null,
    sigma_z = row["sigma_z"],
    n_sim = n_sim
  )
}))


# 标记模型类型
simulation_results_gee$model <- "GEE"
simulation_results_ols$model <- "OLS"
comparison_results <- rbind(simulation_results_gee, simulation_results_ols)

# 绘制对比图
ggplot(comparison_results, aes(x = rho, y = type1_error_rate, color = model, linetype = model)) +
  geom_line(size = 1) + 
  geom_point(size = 2) +
  facet_grid(n_subjects ~ sigma_z, labeller = label_both) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(
    title = "Type I Error Rate: OLS vs GEE",
    x = "Correlation Level (rho)",
    y = "Type I Error Rate",
    color = "Model", linetype = "Model"
  ) +
  theme_minimal()

ggplot(comparison_results, aes(x = rho, y = coverage_prob, color = model, linetype = model)) +
  geom_line(size = 1) + 
  geom_point(size = 2) +
  facet_grid(n_subjects ~ sigma_z, labeller = label_both) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(
    title = "Coverage Probability: OLS vs GEE",
    x = "Correlation Level (rho)",
    y = "Coverage Probability",
    color = "Model", linetype = "Model"
  ) +
  theme_minimal()