
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

import numpy as np
from scipy.stats import poisson, kstest
from scipy.special import gammaln

def poisson_log_likelihood(data, lam):
    # 泊松对数似然
    return np.sum(data * np.log(lam) - lam - gammaln(data + 1))

def powerlaw_log_likelihood(data, alpha, xmin):
    # 幂律分布的对数似然，离散幂律近似公式
    # p(x) = x^{-alpha} / Z, Z = sum_{x=xmin}^\infty x^{-alpha}
    # 这里用连续近似Zeta函数： Z = zeta(alpha, xmin)
    # 为简单起见，用近似式：
    # ll = n * log((alpha -1)/xmin) - alpha * sum(log(x/xmin))
    n = len(data)
    return n * np.log((alpha - 1) / xmin) - alpha * np.sum(np.log(data / xmin))

def discrete_ks_statistic(data, cdf_func):
    # 计算离散数据的KS统计量
    data_sorted = np.sort(data)
    n = len(data)
    empirical_cdf = np.arange(1, n+1) / n
    model_cdf = np.array([cdf_func(x) for x in data_sorted])
    ks = np.max(np.abs(empirical_cdf - model_cdf))
    return ks

def powerlaw_cdf(x, alpha, xmin):
    # 幂律离散分布的累计分布函数，近似计算
    # CDF(x) = (sum_{k=xmin}^x k^{-alpha}) / (sum_{k=xmin}^\infty k^{-alpha})
    # 这里使用Hurwitz zeta函数做归一化
    from scipy.special import zeta
    if x < xmin:
        return 0.0
    else:
        numerator = np.sum([k**(-alpha) for k in range(xmin, int(x)+1)])
        denominator = zeta(alpha, xmin)
        return numerator / denominator

def poisson_cdf(x, lam):
    return poisson.cdf(x, lam)

def powerlaw_analysis(k_list, freq_list):
    result = {
        'xmin': np.nan,
        'alpha': np.nan,
        'll_powerlaw': np.nan,
        'll_poisson': np.nan,
        'likelihood_ratio': np.nan,
        'p_value_lr': np.nan,
        'ks_powerlaw': np.nan,
        'ks_poisson': np.nan,
        'powerlaw_plausible': False,
        'powerlaw_better': False
    }

    try:
        k = np.array(k_list, dtype=int)
        freq = np.array(freq_list, dtype=int)
        data = np.repeat(k, freq)

        # 初步确定 xmin 为最小值，或者可以优化
        xmin = np.min(data)
        data_cut = data[data >= xmin]

        # 使用 Clauset 方法估计 alpha
        # alpha = 1 + n / sum(log(x/xmin))
        n = len(data_cut)
        if n == 0:
            raise ValueError("No data >= xmin")
        alpha = 1 + n / np.sum(np.log(data_cut / xmin))

        # 计算幂律对数似然
        ll_powerlaw = powerlaw_log_likelihood(data_cut, alpha, xmin)

        # 泊松参数估计 lambda 为样本均值
        lambda_poisson = np.mean(data_cut)
        ll_poisson = poisson_log_likelihood(data_cut, lambda_poisson)

        # 似然比检验
        LR = 2 * (ll_powerlaw - ll_poisson)
        from scipy.stats import chi2
        p_value_lr = chi2.sf(LR, df=1)

        # 计算KS统计量，检验幂律拟合
        ks_powerlaw = discrete_ks_statistic(data_cut, lambda x: powerlaw_cdf(x, alpha, xmin))

        # KS统计量，泊松拟合
        ks_poisson = discrete_ks_statistic(data_cut, lambda x: poisson_cdf(x, lambda_poisson))

        # 结果赋值
        result['xmin'] = xmin
        result['alpha'] = alpha
        result['ll_powerlaw'] = ll_powerlaw
        result['ll_poisson'] = ll_poisson
        result['likelihood_ratio'] = LR
        result['p_value_lr'] = p_value_lr
        result['ks_powerlaw'] = ks_powerlaw
        result['ks_poisson'] = ks_poisson
        # 判断幂律是否合理（KS距离小）
        result['powerlaw_plausible'] = (ks_powerlaw < 0.1)
        # 判断幂律是否比泊松好（似然比检验p值小，且KS更优）
        result['powerlaw_better'] = (p_value_lr < 0.05) and (ks_powerlaw < ks_poisson)

    except Exception as e:
        result['error'] = str(e)
        print(f"Error in powerlaw_analysis: {str(e)}")

    return result
