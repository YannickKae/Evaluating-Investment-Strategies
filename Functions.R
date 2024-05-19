### STILL in PROGRESS ###

# Libraries
library(stats)
library(dplyr)

# Computes the Sharpe Ratio from a vector of returns
sharpe_ratio <- function(returns, risk_free_rate = 0) {
  mean_return <- mean(returns) - risk_free_rate
  std_dev <- sd(returns)
  sharpe_ratio <- mean_return / std_dev
  
  return(sharpe_ratio)
}

# Expected maximum Sharpe ratio
expected_max_sharpe_ratio <- function(mean_sharpe, var_sharpe, M) {
  gamma <- 0.5772156649015328606 # Euler-Mascheroni constant
  result <- mean_sharpe + sqrt(var_sharpe) * ((1 - gamma) * qnorm(1 - 1 / M) + gamma * qnorm(1 - 1 / (M * exp(1))))
  
  return(result)
}

# Computes the t-Statistic from a vector of returns
t_statistic <- function(returns, risk_free_rate = 0) {
  N <- length(returns)
  sr <- sharpe_ratio(returns, risk_free_rate)
  t_stat <- sr * sqrt(N)
  
  return(t_stat)
}

# Computes the multiple testing adjusted critical t-values by the Bonferroni method
bonferroni_t_statistic <- function(t_statistics, significance_level = 0.05) {
  num_tests <- length(t_statistics)
  adjusted_alpha <- rep(significance_level / num_tests, num_tests)
  z_critical <- qnorm(1 - adjusted_alpha / 2)
  
  results <- data.frame(
    Test_Number = 1:num_tests,
    t_Statistic = t_statistics,
    Necessary_t_Statistic = z_critical,
    Necessary_p_Value = adjusted_alpha,
    Success = as.integer(t_statistics > z_critical)
  )
  
  return(results)
}

# Computes the multiple testing adjusted critical t-values by the Holm method
holm_t_statistics <- function(t_statistics, significance_level = 0.05) {
  num_tests <- length(t_statistics)
  sorted_indices <- order(t_statistics, decreasing = TRUE)
  sorted_t_statistics <- t_statistics[sorted_indices]
  adjusted_alpha <- significance_level / (num_tests + 1 - 1:num_tests)
  z_critical <- qnorm(1 - adjusted_alpha / 2)
  
  results <- data.frame(
    Test_Number = sorted_indices,
    t_Statistic = sorted_t_statistics,
    Necessary_t_Statistic = z_critical,
    Necessary_p_Value = adjusted_alpha,
    Success = as.integer(sorted_t_statistics > z_critical)
  )
  
  results <- results[order(results$Test_Number), ]
  rownames(results) <- NULL
  
  return(results)
}

# Computes the multiple testing adjusted critical t-values by the BHY method
bhy_t_statistics <- function(t_statistics, significance_level = 0.05) {
  num_tests <- length(t_statistics)
  c_m <- sum(1 / 1:num_tests)
  sorted_indices <- order(t_statistics, decreasing = TRUE)
  sorted_t_statistics <- t_statistics[sorted_indices]
  adjusted_alpha <- (1:num_tests) * significance_level / (num_tests * c_m)
  z_critical <- qnorm(1 - adjusted_alpha / 2)
  
  results <- data.frame(
    Test_Number = sorted_indices,
    t_Statistic = sorted_t_statistics,
    Necessary_t_Statistic = z_critical,
    Necessary_p_Value = adjusted_alpha,
    Success = as.integer(sorted_t_statistics > z_critical)
  )
  
  results <- results[order(results$Test_Number), ]
  rownames(results) <- NULL
  
  return(results)
}

# Bundles all functions to compute the adjusted critical t-values
necessary_t_statistics <- function(t_statistics, significance_level, method = 'bonferroni') {
  if (method == 'bonferroni') {
    return(bonferroni_t_statistic(t_statistics, significance_level))
  } else if (method == 'holm') {
    return(holm_t_statistics(t_statistics, significance_level))
  } else if (method == 'bhy') {
    return(bhy_t_statistics(t_statistics, significance_level))
  } else {
    stop("Method must be 'bonferroni', 'holm', or 'bhy'")
  }
}

# Sharpe Ratio Haircut
haircut_sharpe_ratio <- function(returns, risk_free_rate, num_tests, k = 1, method = 'bonferroni') {
  N <- length(returns)
  t <- t_statistic(returns, risk_free_rate)
  p <- 2 * pnorm(-abs(t))
  min_p_value <- 1e-10
  p <- max(p, min_p_value)
  
  if (method == 'bonferroni') {
    p_adj <- min(p * num_tests, 1)
  } else if (method == 'holm') {
    p_adj <- min(p * (num_tests + 1 - k), 1)
  } else if (method == 'bhy') {
    c_m <- sum(1 / 1:num_tests)
    p_adj <- min(p * num_tests * c_m / k, 1)
  } else {
    stop("Method must be 'bonferroni', 'holm', or 'bhy'")
  }
  
  t_adj <- qnorm(1 - p_adj / 2)
  SR_adj <- t_adj / sqrt(N)
  
  return(SR_adj)
}

# Evaluate all strategies
evaluate_strategies <- function(returns_matrix, risk_free_rate = 0) {
  num_strategies <- ncol(returns_matrix)
  N <- nrow(returns_matrix)
  
  original_sharpe_ratios <- c()
  t_statistics <- c()
  
  for (i in 1:num_strategies) {
    returns <- returns_matrix[, i]
    sr <- sharpe_ratio(returns, risk_free_rate)
    t_stat <- t_statistic(returns, risk_free_rate)
    original_sharpe_ratios <- c(original_sharpe_ratios, sr)
    t_statistics <- c(t_statistics, t_stat)
  }
  
  sorted_indices <- order(t_statistics, decreasing = TRUE)
  
  haircut_sharpe_ratios_bonferroni <- c()
  haircut_sharpe_ratios_holm <- c()
  haircut_sharpe_ratios_bhy <- c()
  
  for (k in 1:length(sorted_indices)) {
    idx <- sorted_indices[k]
    returns <- returns_matrix[, idx]
    haircut_sharpe_ratios_bonferroni <- c(haircut_sharpe_ratios_bonferroni, haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k = k, method = 'bonferroni'))
    haircut_sharpe_ratios_holm <- c(haircut_sharpe_ratios_holm, haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k = k, method = 'holm'))
    haircut_sharpe_ratios_bhy <- c(haircut_sharpe_ratios_bhy, haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k = k, method = 'bhy'))
  }
  
  results <- data.frame(
    Strategy = sorted_indices,
    Original = original_sharpe_ratios[sorted_indices],
    Bonferroni = haircut_sharpe_ratios_bonferroni,
    Holm = haircut_sharpe_ratios_holm,
    BHY = haircut_sharpe_ratios_bhy
  )
  
  return(results)
}
