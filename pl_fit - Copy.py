# ============================================================================
# Copyright (c) 2025 Thomas Broggini and Liu Xiao. All rights reserved.
#
# This code was jointly written by:
# Thomas Broggini (Frankfurt University, Germany)
# Liu Xiao (Xiangyang First Hospital, China)
#
# The way how the mathematical models and methodologies implemented in this code have been
# individually customized and are not intended for generic use. 
#
# For permission requests or inquiries, please contact:
# Thomas Broggini  : broggini@med.uni-frankfurt.de
# Liu Xiao         : careyneurosurgery@gmail.com  /  carey-lau@foxmail.com
#
# Unauthorized use will be considered a violation of intellectual property rights.
# ============================================================================



# maximum-likelihood fitting methods with goodness-of-fit tests based on the Kolmogorov-Smirnov statistic and likelihood ratios as described by Clauset et al.
# Clauset, A., Shalizi, C. R. & Newman, M. E. J. Power-Law Distributions in Empirical Data. Siam Rev 51, 661-703 (2009). https://doi.org:10.1137/070710111
# Xiao Liu
# 2025/07/26
# copyright reserved

import numpy as np
import powerlaw
from scipy.stats import poisson, chi2
from scipy.special import gammaln


def powerlaw_analysis(k_list, freq_list):
    result = {
        'xmin': np.nan,
        'alpha': np.nan,
        'll_powerlaw': np.nan,
        'll_poisson': np.nan,
        'likelihood_ratio': np.nan,
        'p_value_lr': np.nan,
        'p_value_ks': np.nan,
        'powerlaw_plausible': False,
        'powerlaw_better': False
    }

    try:
        k = np.array(k_list, dtype=int)
        freq = np.array(freq_list, dtype=int)
        data = np.repeat(k, freq)
	
        fit = powerlaw.Fit(data, discrete=True, verbose=False)
        xmin = float(fit.xmin)
        alpha = float(fit.alpha)

        data_cut = data[data >= xmin]
        xmin_adjusted = False

        if len(data_cut) < 50:
            sorted_data = np.sort(np.unique(data))
            for x in sorted_data:
                if np.sum(data >= x) >= 50:
                    xmin = x
                    xmin_adjusted = True
                    break
            fit = powerlaw.Fit(data, discrete=True, xmin=xmin, verbose=False)
            alpha = float(fit.alpha)
            data_cut = data[data >= xmin]

        # 打印最终使用的 xmin
        if xmin_adjusted:
            print(f"Adjusted xmin: {xmin}, with {len(data_cut)} data points ≤ xmin.")
        else:
            print(f"Original xmin is sufficient: {xmin}, with {len(data_cut)} data points ≥ xmin.")

        n = len(data_cut)

        result['xmin'] = xmin
        result['alpha'] = alpha

        # log-likelihood for power-law
        ll_powerlaw = n * np.log((alpha - 1) / xmin) - alpha * np.sum(np.log(data_cut / xmin))
        result['ll_powerlaw'] = ll_powerlaw

        # poisson log-likelihood
        lambda_poisson = np.mean(data_cut)
        ll_poisson = np.sum(data_cut * np.log(lambda_poisson) - lambda_poisson - gammaln(data_cut + 1))
        result['ll_poisson'] = ll_poisson

        # likelihood ratio
        LR = 2 * (ll_powerlaw - ll_poisson)
        result['likelihood_ratio'] = LR

        # p value for LR
        p_value_lr = chi2.sf(LR, df=1)
        result['p_value_lr'] = p_value_lr

        # KS test 
        R, p_ks = fit.distribution_compare('power_law', 'poisson')
        result['p_value_ks'] = p_ks

        # decision
        result['powerlaw_plausible'] = bool(p_ks > 0.1)
        result['powerlaw_better'] = bool(p_value_lr < 0.025)

    except Exception as e:
        result['error'] = str(e)
        print(f"Error in powerlaw_analysis: {str(e)}")

    return result