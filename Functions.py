# Libraries
import numpy as np
import scipy.stats
import pandas as pd

# Computes the Sharpe Ratio rom a vector of returns
def sharpe_ratio(returns, risk_free_rate=0):
    mean_return = np.mean(returns) - risk_free_rate
    std_dev = np.std(returns, ddof=1)
    sharpe_ratio = mean_return / std_dev

    return sharpe_ratio

# Computes the t-Statistic rom a vector of returns
def t_statistic(returns, risk_free_rate=0):
    N = len(returns)
    sr = sharpe_ratio(returns, risk_free_rate)
    t_stat = sr * np.sqrt(N)

    return t_stat

# Computes the multiple testing adjusted critical t-values by the Bonferroni method
def bonferroni_t_statistic(t_statistics, significance_level = 0.05):
    num_tests = len(t_statistics)
    adjusted_alpha = np.array([significance_level / num_tests for _ in range(num_tests)])
    z_critical = scipy.stats.norm.ppf(1 - adjusted_alpha / 2)

    results = pd.DataFrame({
        'Test Number': np.arange(1, num_tests + 1),
        't-Statistic': t_statistics,
        'Necessary t-Statistic': z_critical,
        'Necessary p-Value': adjusted_alpha,
        'Success': (t_statistics > z_critical).astype(int)
    })

    return results

# Computes the multiple testing adjusted critical t-values by the Holm method
def holm_t_statistics(t_statistics, significance_level = 0.05):
    num_tests = len(t_statistics)
    sorted_indices = np.argsort(t_statistics)[::-1]  # Sort in descending order, because we use t-Statistics
    sorted_t_statistics = np.array(t_statistics)[sorted_indices]
    adjusted_alpha = np.array([significance_level / (num_tests + 1 - k) for k in range(1, num_tests + 1)])
    z_critical = scipy.stats.norm.ppf(1 - adjusted_alpha / 2)

    results = pd.DataFrame({
        'Test Number': sorted_indices + 1,
        't-Statistic': sorted_t_statistics,
        'Necessary t-Statistic': z_critical,
        'Necessary p-Value': adjusted_alpha,
        'Success': (sorted_t_statistics > z_critical).astype(int)
    })

    results = results.sort_values(by='Test Number').reset_index(drop=True)

    return results

# Computes the multiple testing adjusted critical t-values by the BHY method
def bhy_t_statistics(t_statistics, significance_level = 0.05):
    num_tests = len(t_statistics)
    c_m = np.sum([1.0 / i for i in range(1, num_tests + 1)])
    sorted_indices = np.argsort(t_statistics)[::-1] # Sort in descending order, because we use t-Statistics
    sorted_t_statistics = np.array(t_statistics)[sorted_indices]
    adjusted_alpha = np.array([k * significance_level / (num_tests * c_m) for k in range(1, num_tests + 1)])
    z_critical = scipy.stats.norm.ppf(1 - adjusted_alpha / 2)

    results = pd.DataFrame({
        'Test Number': sorted_indices + 1,
        't-Statistic': sorted_t_statistics,
        'Necessary t-Statistic': z_critical,
        'Necessary p-Value': adjusted_alpha,
        'Success': (sorted_t_statistics > z_critical).astype(int)
    })

    results = results.sort_values(by='Test Number').reset_index(drop=True)

    return results

# Bundles all functions to compute the adjusted critical t-values
def necessary_t_statistics(t_statistics, significance_level, method='bonferroni'):
    if method == 'bonferroni':
        return bonferroni_t_statistic(t_statistics, significance_level)
    elif method == 'holm':
        return holm_t_statistics(t_statistics, significance_level)
    elif method == 'bhy':
        return bhy_t_statistics(t_statistics, significance_level)
    else:
        raise ValueError("Method must be 'bonferroni', 'holm', or 'bhy'")


# Sharpe Ratio Haircut
def haircut_sharpe_ratio(returns, risk_free_rate, num_tests, k=1, method='bonferroni'):
    N = len(returns)
    t = t_statistic(returns, risk_free_rate)
    p = 2 * scipy.stats.norm.sf(abs(t))
    min_p_value = 1e-10
    p = max(p, min_p_value)
    if method == 'bonferroni':
        p_adj = min(p * num_tests, 1)
    elif method == 'holm':
        p_adj = min(p * (num_tests + 1 - k), 1)
    elif method == 'bhy':
        c_m = np.sum([1.0 / i for i in range(1, num_tests + 1)])
        p_adj = min(p * num_tests * c_m / k, 1)
    else:
        raise ValueError("Method must be 'bonferroni', 'holm', or 'bhy'")
    t_adj = scipy.stats.norm.ppf(1 - p_adj / 2)
    SR_adj = t_adj / np.sqrt(N)

    return SR_adj

# Evaluate all strategies
def evaluate_strategies(returns_matrix, risk_free_rate=0):
    num_strategies = returns_matrix.shape[1]
    N = returns_matrix.shape[0]

    original_sharpe_ratios = []
    t_statistics = []
    for i in range(num_strategies):
        returns = returns_matrix[:, i]
        sr = sharpe_ratio(returns, risk_free_rate)
        t_stat = t_statistic(returns, risk_free_rate)
        original_sharpe_ratios.append(sr)
        t_statistics.append(t_stat)

    # Sort indices based on t_statistics in descending order
    sorted_indices = np.argsort(t_statistics)[::-1]

    # Compute adjusted Sharpe Ratios
    haircut_sharpe_ratios_bonferroni = []
    haircut_sharpe_ratios_holm = []
    haircut_sharpe_ratios_bhy = []

    for k, idx in enumerate(sorted_indices):
        returns = returns_matrix[:, idx]
        haircut_sharpe_ratios_bonferroni.append(haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k=k+1, method='bonferroni'))
        haircut_sharpe_ratios_holm.append(haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k=k+1, method='holm'))
        haircut_sharpe_ratios_bhy.append(haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k=k+1, method='bhy'))

    # Arrange the results in the sorted order
    results = pd.DataFrame({
        'Strategy': sorted_indices + 1,
        'Original': np.array(original_sharpe_ratios)[sorted_indices],
        'Bonferroni': haircut_sharpe_ratios_bonferroni,
        'Holm': haircut_sharpe_ratios_holm,
        'BHY': haircut_sharpe_ratios_bhy
    })

    return results
