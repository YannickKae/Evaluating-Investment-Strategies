# Packages
library(dplyr)
library(tidyr)
library(purrr)

# Computes the Sharpe Ratio from a vector of returns
sharpe_ratio <- function(returns, risk_free_rate = 0) {
  mean_return <- mean(returns) - risk_free_rate
  std_dev <- sd(returns)
  sharpe_ratio <- mean_return / std_dev
  
  return(sharpe_ratio)
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
  ) %>%
    arrange(Test_Number) %>%
    mutate(Test_Number = row_number())
  
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
  ) %>%
    arrange(Test_Number) %>%
    mutate(Test_Number = row_number())
  
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
  p <- 2 * (1 - pnorm(abs(t)))
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
  
  original_sharpe_ratios <- numeric(num_strategies)
  haircut_sharpe_ratios_bonferroni <- numeric(num_strategies)
  haircut_sharpe_ratios_holm <- numeric(num_strategies)
  haircut_sharpe_ratios_bhy <- numeric(num_strategies)
  
  for (i in 1:num_strategies) {
    returns <- returns_matrix[, i]
    sr <- sharpe_ratio(returns, risk_free_rate)
    original_sharpe_ratios[i] <- sr
    haircut_sharpe_ratios_bonferroni[i] <- haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k = i, method = 'bonferroni')
    haircut_sharpe_ratios_holm[i] <- haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k = i, method = 'holm')
    haircut_sharpe_ratios_bhy[i] <- haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k = i, method = 'bhy')
  }
  
  results <- data.frame(
    Strategy = 1:num_strategies,
    Original_SR = original_sharpe_ratios,
    SR_Bonferroni = haircut_sharpe_ratios_bonferroni,
    SR_Holm = haircut_sharpe_ratios_holm,
    SR_BHY = haircut_sharpe_ratios_bhy
  )
  
  return(results)
}
