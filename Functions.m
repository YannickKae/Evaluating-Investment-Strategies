% Computes the Sharpe Ratio from a vector of returns
function sr = sharpe_ratio(returns, risk_free_rate)
    if nargin < 2
        risk_free_rate = 0;
    end
    mean_return = mean(returns) - risk_free_rate;
    std_dev = std(returns, 1);
    sr = mean_return / std_dev;
end

% Expected maximum Sharpe ratio
function result = expected_max_sharpe_ratio(mean_sharpe, var_sharpe, M)
    gamma = 0.5772156649015328606; % Euler-Mascheroni constant
    result = mean_sharpe + sqrt(var_sharpe) * ((1 - gamma) * norminv(1 - 1 / M) + gamma * norminv(1 - 1 / (M * exp(1))));
end

% Computes the t-Statistic from a vector of returns
function t_stat = t_statistic(returns, risk_free_rate)
    if nargin < 2
        risk_free_rate = 0;
    end
    N = length(returns);
    sr = sharpe_ratio(returns, risk_free_rate);
    t_stat = sr * sqrt(N);
end

% Computes the multiple testing adjusted critical t-values by the Bonferroni method
function results = bonferroni_t_statistic(t_statistics, significance_level)
    if nargin < 2
        significance_level = 0.05;
    end
    num_tests = length(t_statistics);
    adjusted_alpha = significance_level / num_tests * ones(num_tests, 1);
    z_critical = norminv(1 - adjusted_alpha / 2);
    
    results = table((1:num_tests)', t_statistics(:), z_critical, adjusted_alpha, t_statistics(:) > z_critical, ...
                    'VariableNames', {'TestNumber', 'tStatistic', 'NecessarytStatistic', 'NecessarypValue', 'Success'});
end

% Computes the multiple testing adjusted critical t-values by the Holm method
function results = holm_t_statistics(t_statistics, significance_level)
    if nargin < 2
        significance_level = 0.05;
    end
    num_tests = length(t_statistics);
    [sorted_t_statistics, sorted_indices] = sort(t_statistics, 'descend');
    adjusted_alpha = significance_level ./ (num_tests + 1 - (1:num_tests)');
    z_critical = norminv(1 - adjusted_alpha / 2);
    
    results = table(sorted_indices, sorted_t_statistics, z_critical, adjusted_alpha, sorted_t_statistics > z_critical, ...
                    'VariableNames', {'TestNumber', 'tStatistic', 'NecessarytStatistic', 'NecessarypValue', 'Success'});
    
    results = sortrows(results, 'TestNumber');
end

% Computes the multiple testing adjusted critical t-values by the BHY method
function results = bhy_t_statistics(t_statistics, significance_level)
    if nargin < 2
        significance_level = 0.05;
    end
    num_tests = length(t_statistics);
    c_m = sum(1 ./ (1:num_tests));
    [sorted_t_statistics, sorted_indices] = sort(t_statistics, 'descend');
    adjusted_alpha = (1:num_tests)' * significance_level / (num_tests * c_m);
    z_critical = norminv(1 - adjusted_alpha / 2);
    
    results = table(sorted_indices, sorted_t_statistics, z_critical, adjusted_alpha, sorted_t_statistics > z_critical, ...
                    'VariableNames', {'TestNumber', 'tStatistic', 'NecessarytStatistic', 'NecessarypValue', 'Success'});
    
    results = sortrows(results, 'TestNumber');
end

% Bundles all functions to compute the adjusted critical t-values
function results = necessary_t_statistics(t_statistics, significance_level, method)
    if nargin < 3
        method = 'bonferroni';
    end
    switch method
        case 'bonferroni'
            results = bonferroni_t_statistic(t_statistics, significance_level);
        case 'holm'
            results = holm_t_statistics(t_statistics, significance_level);
        case 'bhy'
            results = bhy_t_statistics(t_statistics, significance_level);
        otherwise
            error("Method must be 'bonferroni', 'holm', or 'bhy'");
    end
end

% Sharpe Ratio Haircut
function SR_adj = haircut_sharpe_ratio(returns, risk_free_rate, num_tests, k, method)
    if nargin < 4
        k = 1;
    end
    if nargin < 5
        method = 'bonferroni';
    end
    N = length(returns);
    t = t_statistic(returns, risk_free_rate);
    p = 2 * normcdf(-abs(t));
    min_p_value = 1e-10;
    p = max(p, min_p_value);
    switch method
        case 'bonferroni'
            p_adj = min(p * num_tests, 1);
        case 'holm'
            p_adj = min(p * (num_tests + 1 - k), 1);
        case 'bhy'
            c_m = sum(1 ./ (1:num_tests));
            p_adj = min(p * num_tests * c_m / k, 1);
        otherwise
            error("Method must be 'bonferroni', 'holm', or 'bhy'");
    end
    t_adj = norminv(1 - p_adj / 2);
    SR_adj = t_adj / sqrt(N);
end

% Evaluate all strategies
function results = evaluate_strategies(returns_matrix, risk_free_rate)
    if nargin < 2
        risk_free_rate = 0;
    end
    [N, num_strategies] = size(returns_matrix);
    original_sharpe_ratios = zeros(num_strategies, 1);
    t_statistics = zeros(num_strategies, 1);
    for i = 1:num_strategies
        returns = returns_matrix(:, i);
        original_sharpe_ratios(i) = sharpe_ratio(returns, risk_free_rate);
        t_statistics(i) = t_statistic(returns, risk_free_rate);
    end
    
    [sorted_t_statistics, sorted_indices] = sort(t_statistics, 'descend');
    
    haircut_sharpe_ratios_bonferroni = zeros(num_strategies, 1);
    haircut_sharpe_ratios_holm = zeros(num_strategies, 1);
    haircut_sharpe_ratios_bhy = zeros(num_strategies, 1);
    
    for k = 1:num_strategies
        idx = sorted_indices(k);
        returns = returns_matrix(:, idx);
        haircut_sharpe_ratios_bonferroni(k) = haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k, 'bonferroni');
        haircut_sharpe_ratios_holm(k) = haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k, 'holm');
        haircut_sharpe_ratios_bhy(k) = haircut_sharpe_ratio(returns, risk_free_rate, num_strategies, k, 'bhy');
    end
    
    results = table(sorted_indices + 1, original_sharpe_ratios(sorted_indices), haircut_sharpe_ratios_bonferroni, ...
                    haircut_sharpe_ratios_holm, haircut_sharpe_ratios_bhy, ...
                    'VariableNames', {'Strategy', 'Original', 'Bonferroni', 'Holm', 'BHY'});
end
