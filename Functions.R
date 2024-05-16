# Packages
library(dplyr)
library(stats)

# Computes the Sharpe Ratio from a vector of returns
sharpe_ratio <- function(returns, risk_free_rate = 0) {
  mean_return <- mean(returns) - risk_free_rate
  std_dev <- sd(returns)
  sharpe_ratio <- mean_return / std_dev
  
  return(sharpe_ratio)
}

# Computes the t-Statistic from the Sharpe Ratio
t_statistic <- function(sharpe_ratio, N, freq = 'annual') {
  if (freq == 'annual') {
    t_stat <- sharpe_ratio * sqrt(N)
  } else if (freq == 'monthly') {
    t_stat <- sharpe_ratio * sqrt(12 * N)
  } else if (freq == 'daily') {
    t_stat <- sharpe_ratio * sqrt(252 * N)
  } else {
    stop("Frequency must be 'annual', 'monthly', or 'daily'")
  }
  
  return(t_stat)
}

# Computes the multiple testing adjusted critical values
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
  
  results <- results %>% arrange(Test_Number)
  
  return(results)
}

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
  
  results <- results %>% arrange(Test_Number)
  
  return(results)
}

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
haircut_sharpe_ratio <- function(sharpe_ratio, num_tests, N, k = 1, freq = 'annual', method = 'bonferroni') {
  t <- t_statistic(sharpe_ratio, N, freq)
  p <- 2 * pnorm(abs(t), lower.tail = FALSE)
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
