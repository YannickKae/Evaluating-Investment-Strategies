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
p_{\text{adjusted}} = \min(p \times m, 1.0)
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
p_{k} = \frac{k \times \alpha}{M \times \sum_{i=1}^{M} \frac{1}{i}}
$$

- $\alpha$: Significance level
- $M$: Number of tests
- $k$: Rank of the test by ascending p-value

## 3. In-Sample (IS) and Out-of-Sample (OOS) Tests

### In-Sample and Out-of-Sample Split

The data is split into In-Sample (IS) and Out-of-Sample (OOS) periods to test the robustness of the trading strategy. This helps determine if the strategy performs well on unseen data and avoids overfitting.

### Probability of Backtest Overfitting (PBO)

The Probability of Backtest Overfitting (PBO) measures how likely it is that an IS strategy will underperform in the OOS period. It is calculated as the probability that the out-of-sample Sharpe Ratio is less than the median in-sample Sharpe Ratio.

$$
PBO = \mathbb{P}(SR_{\text{OOS}} < \text{median}(SR_{\text{IS}}))
$$

- $SR_{\text{OOS}}$: Sharpe Ratio of the out-of-sample period
- $SR_{\text{IS}}$: Sharpe Ratio of the in-sample period

## 4. Adjusting the Sharpe Ratio by Haircut

### Haircutting Sharpe Ratio

The Sharpe Ratio is adjusted (haircut) considering multiple testing adjustments. This adjustment accounts for the increased likelihood of finding a spurious result when multiple tests are conducted.

## 5. Function to Evaluate Trading Strategies

The main function integrates all the above methods to comprehensively evaluate trading strategies.
