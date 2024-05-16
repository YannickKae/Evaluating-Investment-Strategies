# Trading Strategy Evaluation Functions

This repository contains functions for evaluating trading strategies as described in the paper "Evaluating Trading Strategies". The functions implement various statistical techniques to assess the significance of trading strategies considering multiple testing.

## 1. Calculating the Sharpe Ratio and T-Statistic

### Sharpe Ratio

The Sharpe Ratio measures the average return that exceeds the risk-free rate, relative to the volatility of the return. It is a commonly used metric to understand the risk-adjusted return of an investment.

$$
SR = \frac{\mu - r_f}{\sigma}
$$

- $\mu$: Mean return
- $r_f$: Risk-free rate
- $\sigma$: Standard deviation of the return

### $t$-Statistic

The $t$-Statistic is used to assess the significance of the Sharpe Ratio. It scales the Sharpe Ratio by the square root of the number of observations, which helps determine if the observed Sharpe Ratio is statistically significant.

$$
t = SR \times \sqrt{N}
$$

- $SR$: Sharpe Ratio
- $N$: Number of returns

## 2. Multiple Testing Adjustments

### Bonferroni Method

The Bonferroni Method is a conservative approach for multiple testing correction. It reduces the chance of type I errors (false positives) by dividing the significance level by the number of tests.

$$
p = \frac{\alpha}{m}
$$

- $p$: unadjusted $p$-value
- $m$: Number of tests

### Holm Method

The Holm Method is a stepwise correction that is less conservative than the Bonferroni Method. It adjusts the $p$-values sequentially, starting from the smallest p-value, and ensures that the significance level is maintained across multiple tests. The Holm method outputs adjusted p-values that control the family-wise error rate (FWER), which is the probability of making one or more false discoveries among all the tests.

$$
p_{k} = \frac{\alpha}{m + 1 - k}
$$

- $\alpha$: Significance level
- $m$: Number of tests
- $k$: Index of the test sorted by ascending $p$-value

### Benjamini-Hochberg-Yekutieli Method

The BHY Method controls the False Discovery Rate (FDR) and is less conservative than Family-wise Error Rate (FWER) methods like Bonferroni and Holm. FDR is the expected proportion of false discoveries among the rejected hypotheses.

$$
p_{k} = \frac{k \times \alpha}{m \times \left(\frac{1}{1} + \frac{1}{2} + \ldots + \frac{1}{m}\right)}
$$

- $\alpha$: Significance level
- $m$: Number of tests
- $k$: Rank of the test by descending p-value

## 4. Adjusting the Sharpe Ratio by Haircut

The Sharpe Ratio is adjusted (haircut) considering multiple testing adjustments by mapping the adjusted test statistic back to the Sharpe Ratio.
