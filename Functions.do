* Computes the Sharpe Ratio from a vector of returns
program define sharpe_ratio
    syntax varlist(min=1 max=1) [, RiskFreeRate(real 0)]
    
    qui sum `varlist'
    local mean_return = r(mean) - `RiskFreeRate'
    qui sum `varlist', detail
    local std_dev = r(sd)
    local sharpe_ratio = `mean_return' / `std_dev'
    
    return scalar sharpe_ratio = `sharpe_ratio'
end

* Expected maximum Sharpe ratio
program define expected_max_sharpe_ratio
    syntax, MeanSharpe(real) VarSharpe(real) M(real)
    
    local gamma = 0.5772156649015328606
    local term1 = invnormal(1 - 1/`M')
    local term2 = invnormal(1 - 1/(`M' * exp(1)))
    local result = `MeanSharpe' + sqrt(`VarSharpe') * ((1 - `gamma') * `term1' + `gamma' * `term2')
    
    return scalar result = `result'
end

* Computes the t-Statistic from a vector of returns
program define t_statistic
    syntax varlist(min=1 max=1) [, RiskFreeRate(real 0)]
    
    qui count
    local N = r(N)
    sharpe_ratio `varlist', RiskFreeRate(`RiskFreeRate')
    local sr = r(sharpe_ratio)
    local t_stat = `sr' * sqrt(`N')
    
    return scalar t_stat = `t_stat'
end

* Computes the multiple testing adjusted critical t-values by the Bonferroni method
program define bonferroni_t_statistic
    syntax varlist(min=1) [, SignificanceLevel(real 0.05)]
    
    local num_tests = wordcount("`varlist'")
    tempname results
    tempvar t_stat necessary_t_stat necessary_p_value success
    
    gen `results' = .
    gen `t_stat' = .
    gen `necessary_t_stat' = .
    gen `necessary_p_value' = .
    gen `success' = .
    
    forvalues i = 1/`num_tests' {
        local test_var : word `i' of `varlist'
        t_statistic `test_var'
        local t = r(t_stat)
        local adjusted_alpha = `SignificanceLevel' / `num_tests'
        local z_critical = invnormal(1 - `adjusted_alpha' / 2)
        
        replace `results' = `i' in `i'
        replace `t_stat' = `t' in `i'
        replace `necessary_t_stat' = `z_critical' in `i'
        replace `necessary_p_value' = `adjusted_alpha' in `i'
        replace `success' = (`t' > `z_critical') in `i'
    }
    
    list `results' `t_stat' `necessary_t_stat' `necessary_p_value' `success'
end

* Computes the multiple testing adjusted critical t-values by the Holm method
program define holm_t_statistics
    syntax varlist(min=1) [, SignificanceLevel(real 0.05)]
    
    local num_tests = wordcount("`varlist'")
    tempname results
    tempvar t_stat necessary_t_stat necessary_p_value success
    
    gen `results' = .
    gen `t_stat' = .
    gen `necessary_t_stat' = .
    gen `necessary_p_value' = .
    gen `success' = .
    
    local t_values = ""
    forvalues i = 1/`num_tests' {
        local test_var : word `i' of `varlist'
        t_statistic `test_var'
        local t = r(t_stat)
        local t_values = "`t_values' `t'"
    }
    
    * Sort t-values in descending order
    local sorted_t_values = `: list sort descending local t_values'
    
    forvalues k = 1/`num_tests' {
        local t = word(`k') of `sorted_t_values'
        local adjusted_alpha = `SignificanceLevel' / (`num_tests' + 1 - `k')
        local z_critical = invnormal(1 - `adjusted_alpha' / 2)
        
        replace `results' = `k' in `k'
        replace `t_stat' = `t' in `k'
        replace `necessary_t_stat' = `z_critical' in `k'
        replace `necessary_p_value' = `adjusted_alpha' in `k'
        replace `success' = (`t' > `z_critical') in `k'
    }
    
    list `results' `t_stat' `necessary_t_stat' `necessary_p_value' `success'
end

* Computes the multiple testing adjusted critical t-values by the BHY method
program define bhy_t_statistics
    syntax varlist(min=1) [, SignificanceLevel(real 0.05)]
    
    local num_tests = wordcount("`varlist'")
    local c_m = 0
    forvalues i = 1/`num_tests' {
        local c_m = `c_m' + 1/`i'
    }
    
    tempname results
    tempvar t_stat necessary_t_stat necessary_p_value success
    
    gen `results' = .
    gen `t_stat' = .
    gen `necessary_t_stat' = .
    gen `necessary_p_value' = .
    gen `success' = .
    
    local t_values = ""
    forvalues i = 1/`num_tests' {
        local test_var : word `i' of `varlist'
        t_statistic `test_var'
        local t = r(t_stat)
        local t_values = "`t_values' `t'"
    }
    
    * Sort t-values in descending order
    local sorted_t_values = `: list sort descending local t_values'
    
    forvalues k = 1/`num_tests' {
        local t = word(`k') of `sorted_t_values'
        local adjusted_alpha = `k' * `SignificanceLevel' / (`num_tests' * `c_m')
        local z_critical = invnormal(1 - `adjusted_alpha' / 2)
        
        replace `results' = `k' in `k'
        replace `t_stat' = `t' in `k'
        replace `necessary_t_stat' = `z_critical' in `k'
        replace `necessary_p_value' = `adjusted_alpha' in `k'
        replace `success' = (`t' > `z_critical') in `k'
    }
    
    list `results' `t_stat' `necessary_t_stat' `necessary_p_value' `success'
end

* Bundles all functions to compute the adjusted critical t-values
program define necessary_t_statistics
    syntax varlist(min=1) [, SignificanceLevel(real 0.05) Method(string "bonferroni")]
    
    if "`Method'" == "bonferroni" {
        bonferroni_t_statistic `varlist', SignificanceLevel(`SignificanceLevel')
    }
    else if "`Method'" == "holm" {
        holm_t_statistics `varlist', SignificanceLevel(`SignificanceLevel')
    }
    else if "`Method'" == "bhy" {
        bhy_t_statistics `varlist', SignificanceLevel(`SignificanceLevel')
    }
    else {
        di as error "Method must be 'bonferroni', 'holm', or 'bhy'"
    }
end

* Sharpe Ratio Haircut
program define haircut_sharpe_ratio
    syntax varlist(min=1 max=1) [, RiskFreeRate(real 0) NumTests(real) k(real 1) Method(string "bonferroni")]
    
    qui count
    local N = r(N)
    t_statistic `varlist', RiskFreeRate(`RiskFreeRate')
    local t = r(t_stat)
    local p = 2 * (1 - normal(abs(`t')))
    local min_p_value = 1e-10
    local p = max(`p', `min_p_value')
    
    if "`Method'" == "bonferroni" {
        local p_adj = min(`p' * `NumTests', 1)
    }
    else if "`Method'" == "holm" {
        local p_adj = min(`p' * (`NumTests' + 1 - `k'), 1)
    }
    else if "`Method'" == "bhy" {
        local c_m = 0
        forvalues i = 1/`NumTests' {
            local c_m = `c_m' + 1/`i'
        }
        local p_adj = min(`p' * `NumTests' * `c_m' / `k', 1)
    }
    else {
        di as error "Method must be 'bonferroni', 'holm', or 'bhy'"
    }
    
    local t_adj = invnormal(1 - `p_adj' / 2)
    local SR_adj = `t_adj' / sqrt(`N')
    
    return scalar SR_adj = `SR_adj'
end

* Evaluate all strategies
program define evaluate_strategies
    syntax varlist(min=1) [, RiskFreeRate(real 0)]
    
    local num_strategies = wordcount("`varlist'")
    qui count
    local N = r(N)
    
    tempname original_sharpe_ratios t_statistics sorted_indices
    
    gen `original_sharpe_ratios' = .
    gen `t_statistics' = .
    
    forvalues i = 1/`num_strategies' {
        local strategy_var : word `i' of `varlist'
        sharpe_ratio `strategy_var', RiskFreeRate(`RiskFreeRate')
        local sr = r(sharpe_ratio)
        t_statistic `strategy_var', RiskFreeRate(`RiskFreeRate')
        local t = r(t_stat)
        
        replace `original_sharpe_ratios' = `sr' in `i'
        replace `t_statistics' = `t' in `i'
    }
    
    * Sort indices based on t_statistics in descending order
    sort `t_statistics' (descending)
    
    tempname haircut_sharpe_ratios_bonferroni haircut_sharpe_ratios_holm haircut_sharpe_ratios_bhy
    
    gen `haircut_sharpe_ratios_bonferroni' = .
    gen `haircut_sharpe_ratios_holm' = .
    gen `haircut_sharpe_ratios_bhy' = .
    
    forvalues k = 1/`num_strategies' {
        local strategy_var : word `k' of `varlist'
        haircut_sharpe_ratio `strategy_var', RiskFreeRate(`RiskFreeRate') NumTests(`num_strategies') k(`k') Method("bonferroni")
        local sr_bonf = r(SR_adj)
        haircut_sharpe_ratio `strategy_var', RiskFreeRate(`RiskFreeRate') NumTests(`num_strategies') k(`k') Method("holm")
        local sr_holm = r(SR_adj)
        haircut_sharpe_ratio `strategy_var', RiskFreeRate(`RiskFreeRate') NumTests(`num_strategies') k(`k') Method("bhy")
        local sr_bhy = r(SR_adj)
        
        replace `haircut_sharpe_ratios_bonferroni' = `sr_bonf' in `k'
        replace `haircut_sharpe_ratios_holm' = `sr_holm' in `k'
        replace `haircut_sharpe_ratios_bhy' = `sr_bhy' in `k'
    }
    
    * Arrange the results in the sorted order
    list `original_sharpe_ratios' `haircut_sharpe_ratios_bonferroni' `haircut_sharpe_ratios_holm' `haircut_sharpe_ratios_bhy'
end
