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


# Computes the t-Statistic from the Sharpe Ratio
def t_statistic(sharpe_ratio, N, freq='annual'):
    if freq == 'annual':
        t_stat = sharpe_ratio * np.sqrt(N)
    elif freq == 'monthly':
        t_stat = sharpe_ratio * np.sqrt(12 * N)
    elif freq == 'daily':
        t_stat = sharpe_ratio * np.sqrt(252 * N)
    else:
        raise ValueError("Frequency must be 'annual', 'monthly', or 'daily'")

    return t_stat


from scipy.stats import norm
import numpy as np
import pandas as pd


def bonferroni_t_statistic(t_statistics, significance_level = 0.05):
    num_tests = len(t_statistics)
    adjusted_alpha = np.array([significance_level / num_tests for _ in range(num_tests)])
    z_critical = norm.ppf(1 - adjusted_alpha / 2)

    results = pd.DataFrame({
        'Test Number': np.arange(1, num_tests + 1),
        't-Statistic': t_statistics,
        'Necessary t-Statistic': z_critical,
        'Necessary p-Value': adjusted_alpha,
        'Success': (t_statistics > z_critical).astype(int)
    })

    return results


def holm_t_statistics(t_statistics, significance_level = 0.05):
    num_tests = len(t_statistics)
    sorted_indices = np.argsort(t_statistics)[::-1]  # Sort in descending order
    sorted_t_statistics = np.array(t_statistics)[sorted_indices]
    adjusted_alpha = np.array([significance_level / (num_tests + 1 - k) for k in range(1, num_tests + 1)])
    z_critical = norm.ppf(1 - adjusted_alpha / 2)

    results = pd.DataFrame({
        'Test Number': sorted_indices + 1,
        't-Statistic': sorted_t_statistics,
        'Necessary t-Statistic': z_critical,
        'Necessary p-Value': adjusted_alpha,
        'Success': (sorted_t_statistics > z_critical).astype(int)
    })

    results = results.sort_values(by='Test Number').reset_index(drop=True)

    return results


def bhy_t_statistics(t_statistics, significance_level = 0.05):
    num_tests = len(t_statistics)
    c_m = np.sum([1.0 / i for i in range(1, num_tests + 1)])
    sorted_indices = np.argsort(t_statistics)  # Sort in ascending order
    sorted_t_statistics = np.array(t_statistics)[sorted_indices]
    adjusted_alpha = np.array([k * significance_level / (num_tests * c_m) for k in range(1, num_tests + 1)])
    z_critical = norm.ppf(1 - adjusted_alpha / 2)

    results = pd.DataFrame({
        'Test Number': sorted_indices + 1,
        't-Statistic': sorted_t_statistics,
        'Necessary t-Statistic': z_critical,
        'Necessary p-Value': adjusted_alpha,
        'Success': (sorted_t_statistics > z_critical).astype(int)
    })

    results = results.sort_values(by='Test Number').reset_index(drop=True)

    return results


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
def haircut_sharpe_ratio(sharpe_ratio, num_tests, N, k=1, freq='annual', method='bonferroni'):
    t = t_statistic(sharpe_ratio, N, freq='annual')
    p = 2 * norm.sf(abs(t))
    if method == 'bonferroni':
        p_adj = min(p * num_tests, 1)
    elif method == 'holm':
        p_adj = min(p * (num_tests + 1 - k), 1)
    elif method == 'bhy':
        c_m = np.sum([1.0 / i for i in range(1, num_tests + 1)])
        p_adj = min(p * num_tests * c_m / k, 1)
    else:
        raise ValueError("Method must be 'bonferroni', 'holm', or 'bhy'")
    t_adj = norm.ppf(1 - p_adj / 2)
    SR_adj = t_adj / np.sqrt(N)

    return SR_adj
