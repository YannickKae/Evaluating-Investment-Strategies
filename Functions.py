# Libraries
import numpy as np
import scipy.stats
import pandas as pd

# --- Core Metrics ---

def sharpe_ratio(returns, risk_free_rate=0):
    """
    Computes the Sharpe Ratio from a vector of returns.

    Args:
        returns (np.ndarray or pd.Series): Array of returns.
        risk_free_rate (float, optional): The risk-free rate of return. Defaults to 0.

    Returns:
        float: The calculated Sharpe Ratio.
    """
    if np.std(returns, ddof=1) == 0:
        return np.nan # Avoid division by zero if returns are constant
    mean_return = np.mean(returns) - risk_free_rate
    std_dev = np.std(returns, ddof=1) # Sample standard deviation
    sr = mean_return / std_dev
    return sr

def expected_max_sharpe_ratio(mean_sharpe, var_sharpe, num_tests, num_returns):
    """
    Approximates the expected maximum Sharpe Ratio when testing M strategies.

    Args:
        mean_sharpe (float): The expected Sharpe ratio across tests (E[SR_m]).
        var_sharpe (float): The variance of the Sharpe ratio across tests (Var[SR_m]).
        num_tests (int): The number of independent tests (M).
        num_returns (int): The number of returns observations (N).
                             Note: The original formula provided used N inside Phi^-1,
                             but the description implies M. This implementation uses M
                             as it relates to the number of *tests* being compared.

    Returns:
        float: The approximated expected maximum Sharpe Ratio E[SR_max].
    """
    gamma = 0.5772156649015328606 # Euler-Mascheroni constant
    # Using num_tests (M) based on the context of comparing multiple strategies.
    # If the formula using N (num_returns) was intended, replace num_tests with num_returns below.
    inv_cdf1 = scipy.stats.norm.ppf(1 - 1 / num_tests)
    inv_cdf2 = scipy.stats.norm.ppf(1 - 1 / (num_tests * np.exp(1)))

    expected_max_sr = mean_sharpe + np.sqrt(var_sharpe) * \
                      ((1 - gamma) * inv_cdf1 + gamma * inv_cdf2)
    return expected_max_sr

def t_statistic(returns, risk_free_rate=0):
    """
    Computes the t-Statistic from a vector of returns.

    Args:
        returns (np.ndarray or pd.Series): Array of returns.
        risk_free_rate (float, optional): The risk-free rate of return. Defaults to 0.

    Returns:
        float: The calculated t-Statistic. Returns np.nan if Sharpe Ratio is nan.
    """
    N = len(returns)
    if N == 0:
        return np.nan
    sr = sharpe_ratio(returns, risk_free_rate)
    if np.isnan(sr):
        return np.nan
    t_stat = sr * np.sqrt(N)
    return t_stat

# --- Multiple Testing Adjustments ---

def bonferroni_critical_t(t_statistics, significance_level=0.05):
    """
    Computes the multiple testing adjusted critical t-value using the Bonferroni method.
    Compares each t-statistic against this single critical value.

    Args:
        t_statistics (list or np.ndarray): List or array of observed t-statistics.
        significance_level (float, optional): The desired family-wise error rate (alpha). Defaults to 0.05.

    Returns:
        pd.DataFrame: DataFrame with test number, original t-statistic,
                      the necessary critical t-statistic (same for all tests),
                      the adjusted p-value threshold, and a success flag.
    """
    num_tests = len(t_statistics)
    if num_tests == 0:
        return pd.DataFrame(columns=['Test Number', 't-Statistic', 'Necessary t-Statistic', 'Necessary p-Value', 'Success'])

    # Adjusted significance level per test (two-sided)
    adjusted_alpha_per_test = significance_level / num_tests
    # Critical t-value (or z-value for large N approximation)
    # We look for t > critical_value, so we use 1 - alpha_adj/2 for the upper tail
    z_critical = scipy.stats.norm.ppf(1 - adjusted_alpha_per_test / 2)

    results = pd.DataFrame({
        'Test Number': np.arange(1, num_tests + 1),
        't-Statistic': t_statistics,
        'Necessary t-Statistic': z_critical, # Same critical value for all tests
        'Necessary p-Value': adjusted_alpha_per_test,
        'Success': (np.abs(np.array(t_statistics)) > z_critical).astype(int) # Check if abs(t) > critical
    })
    return results

def holm_critical_t(t_statistics, significance_level=0.05):
    """
    Computes the multiple testing adjusted critical t-values using the Holm method.
    Compares ordered t-statistics against sequentially adjusted critical values.

    Args:
        t_statistics (list or np.ndarray): List or array of observed t-statistics.
        significance_level (float, optional): The desired family-wise error rate (alpha). Defaults to 0.05.

    Returns:
        pd.DataFrame: DataFrame with original test number, t-statistic,
                      the necessary critical t-statistic for its rank,
                      the adjusted p-value threshold for its rank, and a success flag.
                      The DataFrame is sorted by the original test number.
    """
    num_tests = len(t_statistics)
    if num_tests == 0:
        return pd.DataFrame(columns=['Test Number', 't-Statistic', 'Necessary t-Statistic', 'Necessary p-Value', 'Success'])

    # Sort t-statistics descending by absolute value to find most significant first
    # Keep track of original indices
    abs_t_statistics = np.abs(np.array(t_statistics))
    sorted_indices = np.argsort(abs_t_statistics)[::-1]
    sorted_t_statistics = np.array(t_statistics)[sorted_indices]
    sorted_abs_t_statistics = abs_t_statistics[sorted_indices]

    # Calculate Holm-adjusted alpha levels and corresponding critical t-values
    adjusted_alpha = np.array([significance_level / (num_tests + 1 - k) for k in range(1, num_tests + 1)])
    z_critical = scipy.stats.norm.ppf(1 - adjusted_alpha / 2) # Upper tail critical value

    # Determine success based on Holm's sequential rejection principle
    success = np.zeros(num_tests, dtype=int)
    reject_hypothesis = True
    for k in range(num_tests):
        if reject_hypothesis and sorted_abs_t_statistics[k] > z_critical[k]:
            success[k] = 1
        else:
            reject_hypothesis = False # Stop rejecting once a test fails
            success[k] = 0 # Mark this and subsequent tests as failed


    results = pd.DataFrame({
        'Original Index': sorted_indices,
        'Test Number': sorted_indices + 1,
        't-Statistic': sorted_t_statistics,
        'Necessary t-Statistic': z_critical,
        'Necessary p-Value': adjusted_alpha,
        'Success': success
    })

    # Sort back to original order based on 'Test Number'
    results = results.sort_values(by='Original Index').reset_index(drop=True)
    results = results.drop(columns=['Original Index']) # Clean up helper column


    return results


def bhy_critical_t(t_statistics, significance_level=0.05):
    """
    Computes the multiple testing adjusted critical t-values using the BHY method (controls FDR).
    Compares ordered t-statistics against sequentially adjusted critical values.

    Args:
        t_statistics (list or np.ndarray): List or array of observed t-statistics.
        significance_level (float, optional): The desired False Discovery Rate (FDR) level (alpha). Defaults to 0.05.

    Returns:
        pd.DataFrame: DataFrame with original test number, t-statistic,
                      the necessary critical t-statistic for its rank,
                      the adjusted p-value threshold for its rank, and a success flag.
                      The DataFrame is sorted by the original test number.
    """
    num_tests = len(t_statistics)
    if num_tests == 0:
        return pd.DataFrame(columns=['Test Number', 't-Statistic', 'Necessary t-Statistic', 'Necessary p-Value', 'Success'])

    # Constant c(m) for BHY
    c_m = np.sum([1.0 / i for i in range(1, num_tests + 1)])

    # Sort t-statistics descending by absolute value
    # Keep track of original indices
    abs_t_statistics = np.abs(np.array(t_statistics))
    sorted_indices = np.argsort(abs_t_statistics)[::-1]
    sorted_t_statistics = np.array(t_statistics)[sorted_indices]
    sorted_abs_t_statistics = abs_t_statistics[sorted_indices]


    # Calculate BHY-adjusted alpha levels and corresponding critical t-values
    adjusted_alpha = np.array([k * significance_level / (num_tests * c_m) for k in range(1, num_tests + 1)])
    z_critical = scipy.stats.norm.ppf(1 - adjusted_alpha / 2) # Upper tail critical value

    # Determine success based on BHY procedure: find the largest k such that |t|_(k) > critical_t(k)
    # All tests with rank j <= k are considered successful (reject null hypothesis)
    success = np.zeros(num_tests, dtype=int)
    max_k = 0
    for k in range(num_tests):
        if sorted_abs_t_statistics[k] > z_critical[k]:
             max_k = k + 1 # rank is k+1

    if max_k > 0:
        success[:max_k] = 1 # Mark all up to rank max_k as successful

    results = pd.DataFrame({
        'Original Index': sorted_indices,
        'Test Number': sorted_indices + 1,
        't-Statistic': sorted_t_statistics,
        'Necessary t-Statistic': z_critical,
        'Necessary p-Value': adjusted_alpha,
        'Success': success
    })

    # Sort back to original order
    results = results.sort_values(by='Original Index').reset_index(drop=True)
    results = results.drop(columns=['Original Index'])

    return results


def necessary_critical_t(t_statistics, significance_level=0.05, method='bonferroni'):
    """
    Wrapper function to compute multiple testing adjusted critical t-values
    using the specified method.

    Args:
        t_statistics (list or np.ndarray): List or array of observed t-statistics.
        significance_level (float, optional): The significance or FDR level (alpha). Defaults to 0.05.
        method (str, optional): The adjustment method ('bonferroni', 'holm', 'bhy'). Defaults to 'bonferroni'.

    Returns:
        pd.DataFrame: DataFrame with results from the chosen method.

    Raises:
        ValueError: If the method is not one of the allowed options.
    """
    if method == 'bonferroni':
        return bonferroni_critical_t(t_statistics, significance_level)
    elif method == 'holm':
        return holm_critical_t(t_statistics, significance_level)
    elif method == 'bhy':
        return bhy_critical_t(t_statistics, significance_level)
    else:
        raise ValueError("Method must be 'bonferroni', 'holm', or 'bhy'")


# --- Strategy Evaluation ---

def haircut_sharpe_ratio(returns, risk_free_rate, num_tests, rank=1, method='bonferroni'):
    """
    Adjusts the Sharpe Ratio based on multiple testing correction using p-value adjustment.

    Args:
        returns (np.ndarray or pd.Series): Array of returns for a single strategy.
        risk_free_rate (float): The risk-free rate of return.
        num_tests (int): The total number of tests (M) being conducted.
        rank (int, optional): The rank of this test when sorted by significance (p-value ascending, or t-stat descending).
                               Required for 'holm' and 'bhy'. Defaults to 1.
        method (str, optional): The adjustment method ('bonferroni', 'holm', 'bhy'). Defaults to 'bonferroni'.

    Returns:
        float: The adjusted Sharpe Ratio (SR_cor). Returns np.nan if original t-stat is nan.

    Raises:
        ValueError: If the method is not one of the allowed options.
    """
    N = len(returns)
    if N == 0:
        return np.nan

    t = t_statistic(returns, risk_free_rate)
    if np.isnan(t):
        return np.nan

    # Calculate original two-sided p-value
    # Use survival function (1 - cdf) for potentially better precision in the tail
    p = 2 * scipy.stats.norm.sf(abs(t))
    min_p_value = 1e-10 # Floor p-value to avoid issues with p=0
    p = max(p, min_p_value)

    # Adjust p-value based on the method
    if method == 'bonferroni':
        p_adj = min(p * num_tests, 1.0)
    elif method == 'holm':
        # Holm adjustment factor depends on rank k (here called 'rank')
        p_adj = min(p * (num_tests + 1 - rank), 1.0)
    elif method == 'bhy':
        # BHY adjustment factor depends on rank k and c(m)
        c_m = np.sum([1.0 / i for i in range(1, num_tests + 1)])
        p_adj = min(p * num_tests * c_m / rank, 1.0)
    else:
        raise ValueError("Method must be 'bonferroni', 'holm', or 'bhy'")

    # Calculate the t-statistic corresponding to the adjusted p-value
    # This is the 'corrected' t-statistic (t_cor)
    # ppf(1 - p_adj / 2) gives the upper critical value for a two-sided test
    t_adj = scipy.stats.norm.ppf(1 - p_adj / 2) if p_adj < 1.0 else 0.0 # If p_adj=1, t_adj=0

    # If original t was negative, the adjusted t should also be negative (or zero)
    # The haircut should shrink the SR towards zero, not flip its sign due to p-value adjustment.
    t_adj = np.sign(t) * t_adj if t != 0 else 0.0


    # Calculate the adjusted Sharpe Ratio
    SR_adj = t_adj / np.sqrt(N) if N > 0 else np.nan

    return SR_adj

# Evaluate all strategies
def evaluate_strategies(returns_matrix, risk_free_rate=0):
    """
    Evaluates multiple strategies, calculating original and haircut Sharpe Ratios.

    Args:
        returns_matrix (np.ndarray): 2D array where each column represents the returns of one strategy.
        risk_free_rate (float, optional): The risk-free rate. Defaults to 0.

    Returns:
        pd.DataFrame: DataFrame containing the original and adjusted Sharpe Ratios
                      for each strategy, sorted by the original t-statistic descending.
                      Columns: 'Strategy', 't-Statistic', 'Original SR', 'Bonferroni SR', 'Holm SR', 'BHY SR'.
    """
    num_strategies = returns_matrix.shape[1]
    num_returns = returns_matrix.shape[0]

    if num_strategies == 0 or num_returns == 0:
        return pd.DataFrame(columns=['Strategy', 't-Statistic', 'Original SR', 'Bonferroni SR', 'Holm SR', 'BHY SR'])

    original_sharpe_ratios = []
    t_statistics_list = []
    for i in range(num_strategies):
        returns = returns_matrix[:, i]
        sr = sharpe_ratio(returns, risk_free_rate)
        t_stat = t_statistic(returns, risk_free_rate)
        original_sharpe_ratios.append(sr)
        t_statistics_list.append(t_stat)

    # Sort indices based on t_statistics in descending order (most significant first)
    # Handle NaNs in t-statistics: place them at the end
    t_stats_array = np.array(t_statistics_list)
    valid_indices = np.where(~np.isnan(t_stats_array))[0]
    nan_indices = np.where(np.isnan(t_stats_array))[0]
    sorted_valid_indices = valid_indices[np.argsort(t_stats_array[valid_indices])[::-1]]
    # Combine sorted valid indices with NaN indices
    sorted_indices = np.concatenate((sorted_valid_indices, nan_indices))


    # Compute adjusted Sharpe Ratios using the haircut function
    haircut_sr_bonferroni = []
    haircut_sr_holm = []
    haircut_sr_bhy = []

    num_valid_tests = len(valid_indices) # Use number of valid tests for adjustments

    rank_map = {index: rank + 1 for rank, index in enumerate(sorted_valid_indices)}

    for i in range(num_strategies):
        idx = i # Use the original index 'i' to fetch returns
        returns = returns_matrix[:, idx]
        original_t = t_statistics_list[idx]

        if np.isnan(original_t):
             # Assign NaN for adjusted SR if original t is NaN
             haircut_sr_bonferroni.append(np.nan)
             haircut_sr_holm.append(np.nan)
             haircut_sr_bhy.append(np.nan)
        else:
            # Determine rank for Holm and BHY based on the sorted valid indices
            current_rank = rank_map[idx]
            # Calculate haircuts using the number of valid tests
            haircut_sr_bonferroni.append(haircut_sharpe_ratio(returns, risk_free_rate, num_valid_tests, rank=current_rank, method='bonferroni'))
            haircut_sr_holm.append(haircut_sharpe_ratio(returns, risk_free_rate, num_valid_tests, rank=current_rank, method='holm'))
            haircut_sr_bhy.append(haircut_sharpe_ratio(returns, risk_free_rate, num_valid_tests, rank=current_rank, method='bhy'))


    # Reorder results according to the sorted t-statistics
    results = pd.DataFrame({
        'Strategy': sorted_indices + 1, # Strategy numbers based on sorted significance
        't-Statistic': np.array(t_statistics_list)[sorted_indices],
        'Original SR': np.array(original_sharpe_ratios)[sorted_indices],
        'Bonferroni SR': np.array(haircut_sr_bonferroni)[sorted_indices],
        'Holm SR': np.array(haircut_sr_holm)[sorted_indices],
        'BHY SR': np.array(haircut_sr_bhy)[sorted_indices]
    })

    return results
