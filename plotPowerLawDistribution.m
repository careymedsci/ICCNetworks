
 
% ============================================================================
% Copyright (c) 2025 Thomas Broggini and Liu Xiao. All rights reserved.
%
% This code was jointly written by:
% Thomas Broggini (Frankfurt University, Germany)
% Liu Xiao (Xiangyang First Hospital, China)
%
% The way how the mathematical models and methodologies implemented in this code have been
% individually customized and are not intended for generic use. 
%
% For permission requests or inquiries, please contact:
% Thomas Broggini  : broggini@med.uni-frankfurt.de
% Liu Xiao         : careyneurosurgery@gmail.com  /  carey-lau@foxmail.com
%
% Unauthorized use will be considered a violation of intellectual property rights.
% ============================================================================


% 作者：Thomas Broggini 和 刘晓；法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Thomas Broggini (Supervisor) 
% Xiao Liu (Student)
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany

% Xiangyang No.1 people's Hospital, China
% Also an python version coded 
% xiao.liu@stud.uni-frankfurt.de
% 2025.07.06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ========================================================================
% Function: plotPowerLawDistribution
% Purpose:  Analyze degree distribution using power law and alternative models
%           (Poisson, lognormal, truncated power law) with statistical validation
%           功能：使用幂律及替代模型（泊松、对数正态、截断幂律）分析度分布并进行统计验证
%
% Method:   Implements maximum-likelihood fitting with goodness-of-fit tests
%           based on Kolmogorov-Smirnov statistic and likelihood ratios as 
%           described by Clauset et al.
%           方法：基于Clauset等人提出的KS检验和似然比检验的最大似然拟合法
%
% Reference: 
% Clauset, A., Shalizi, C. R. & Newman, M. E. J. Power-Law Distributions in 
% Empirical Data. Siam Rev 51, 661-703 (2009). 
% https://doi.org:10.1137/070710111
%
% Author:   Xiao Liu and Thomas Broggini
% Date:     2025/07/26
% Copyright reserved
% ========================================================================

function[params, analysis_summary] = plotPowerLawDistribution(params, analysis_summary)
    % params.probabilityDistribution assumed format:
    % columns: [k, N_k, ..., P(k)]
    % where k = degree (number of co-active cells), N_k = number of nodes with degree k,
    % P(k) = probability distribution
    % 参数 params.probabilityDistribution 假定格式：
    % 列：[k, N_k, ..., P(k)]
    % 其中 k 为度（共活跃细胞数），N_k 为度为 k 的节点数，P(k) 为概率分布
    PD = params.probabilityDistribution;
    % 去除 k=0 节点且概率为0的行
    % Remove rows where k=0 or probability P(k)=0
    valid_idx = (PD(:,1) > 0) & (PD(:,4) > 0) & (PD(:,2) > 0);
    k_all = PD(valid_idx, 1);
    N_k_all = PD(valid_idx, 2);

    % 重新计算概率分布，确保规范（概率总和为1）
    % Recalculate probability distribution to normalize sum to 1
    N_total = sum(N_k_all);
    P_k_all = N_k_all / N_total;
    k = k_all;
    P_k = P_k_all;

    % ---------- 幂律拟合过滤 xmin ----------
    % Filter xmin threshold
    answer = inputdlg('input xmin value (Nature accepts xmin ≤ 2):', ...
                  'set xmin', [1 50], {'1'});
    xmin = str2double(answer{1});
    valid_fit = k >= xmin;
    k_fit = k(valid_fit);
    P_k_fit = P_k(valid_fit);

    % ---------- 幂律拟合 (log-log线性拟合) ----------
    % Power law fit (log-log linear fitting)
    logk = log10(k_fit);
    logP = log10(P_k_fit);
    p = polyfit(logk, logP, 1);
    gamma = -p(1);
    a_power = 10^p(2);

    P_power = a_power * k_fit.^(-gamma);

    % 计算拟合优度 R^2
    % Calculate goodness of fit R^2
    r = corr(logk, logP);
    R_squared_power = r^2;

    % ---------- 截断幂律拟合 ----------
    % Truncated power law fit
    trunc_model = fittype('C * x.^(-gamma) .* exp(-lambda * x)', ...
        'independent', 'x', 'coefficients', {'C', 'gamma', 'lambda'});
    opts = fitoptions(trunc_model);
    opts.StartPoint = [max(P_k_fit), gamma, 0.01];
    opts.Lower = [0, 0, 0];
    opts.Upper = [Inf, Inf, Inf];
    try
        [fit_trunc, gof_trunc] = fit(k_fit, P_k_fit, trunc_model, opts);
    catch
        warning('Truncated power law fit failed. Using defaults.');
        fit_trunc.C = NaN; fit_trunc.gamma = NaN; fit_trunc.lambda = NaN;
        gof_trunc.rsquare = NaN;
    end
    if ~isnan(fit_trunc.C)
        P_trunc = fit_trunc.C * k_fit.^(-fit_trunc.gamma) .* exp(-fit_trunc.lambda * k_fit);
    else
        P_trunc = nan(size(k_fit));
    end

   % ---------- 泊松拟合 ----------
    % Poisson Distribution Fit
    % Fit the data to a Poisson distribution
    lambda = sum(k .* N_k_all) / sum(N_k_all); 
    k_max = max([max(k_all), ceil(lambda + 8*sqrt(lambda)), max(k_all)*2]);
    k_range = 0:k_max;
    P_poisson_full = poisspdf(k_range, lambda); 
    P_poisson_at_k = poisspdf(k, lambda);
    min_nonzero_Pk = min(P_k(P_k > 0));
    if isempty(min_nonzero_Pk)
        min_nonzero_Pk = max(P_k); 
    end
    threshold = max(1e-12, min_nonzero_Pk/100); 
    P_poisson_full_plot = P_poisson_full;
    P_poisson_full_plot(P_poisson_full_plot < threshold) = NaN;
    
    P_poisson_plot = P_poisson_at_k;
    P_poisson_plot(P_poisson_plot < threshold) = NaN;
    
    valid_idx_poiss = (P_poisson_at_k > 0 & P_k > 0 & P_poisson_at_k >= threshold);
    if sum(valid_idx_poiss) >= 2
        r_poisson = corr(log10(k(valid_idx_poiss)), log10(P_poisson_at_k(valid_idx_poiss)));
        R_squared_poisson = r_poisson^2;
    else
        R_squared_poisson = NaN;
    end
    % ---------- 对数正态拟合 ----------
    % Lognormal fit
    mu = mean(log(k));
    sigma = std(log(k));
    P_logn = lognpdf(k, mu, sigma);
    P_logn = P_logn / max(P_logn) * max(P_k);
    r_logn = corr(log10(k), log10(P_logn), 'Type', 'Pearson');
    R_squared_logn = r_logn^2;

    % ---------- 绘制图形 ----------
    % Plot figures
    figure; set(gcf,'Color','w');
    loglog(k, P_k, 'ko-', 'MarkerFaceColor','k', 'DisplayName', 'Empirical data');
    hold on;
    loglog(k_fit, P_power, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Power law fit (γ=%.3f, R²=%.3f)', gamma, R_squared_power));
    if ~isnan(fit_trunc.C)
        loglog(k_fit, P_trunc, 'm:', 'LineWidth', 2, 'DisplayName', sprintf('Trunc. Power law (R²=%.3f)', gof_trunc.rsquare));
    end
loglog(k_range, P_poisson_full_plot, 'g--', 'LineWidth', 2, ...
    'DisplayName', sprintf('Poisson fit (\\lambda=%.2f, R^2=%.3g)', lambda, R_squared_poisson));
loglog(k, P_logn, 'b-.', 'LineWidth', 2, 'DisplayName', sprintf('Lognormal fit (R²=%.3f)', R_squared_logn));

    % 标出网络 hubs 区域 
    % Highlight network hubs area
    hub_thresh = prctile(k, 60);
    y_bottom = 1e-4;
    y_top = max(P_k(k >= hub_thresh));
    fill([hub_thresh max(k) max(k) hub_thresh], [y_bottom y_bottom y_top y_top], [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    text(mean([hub_thresh max(k)]), y_top*1.5, 'Network hubs', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    xlabel('Degree k (number of co-active cells per cell)');
    ylabel('Probability P(k)');
    title('Degree distribution and model fits');
    legend('Location','SouthWest');
    grid on;
    hold off;

    % ---------- 统计量计算 ----------
    % Statistical calculations
    k_mean = sum(k .* N_k_all) / N_total;
    k_var = sum(((k - k_mean).^2) .* N_k_all) / N_total;
    k_std = sqrt(k_var);

    % ER 网络标准差理论曲线
    % Erdős-Rényi theoretical std dev curve
    figure; set(gcf,'Color','w');
    k_range = linspace(0, max(k_mean*1.5,4), 100);
    sigma_er = sqrt(k_range);
    fill([k_range fliplr(k_range)], [zeros(size(k_range)) fliplr(sigma_er)], 'c', 'FaceAlpha', 0.2, 'EdgeColor','none');
    hold on;
    plot(k_range, sigma_er, 'r--', 'LineWidth', 2);
    scatter(k_mean, k_std, 100, 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'w');
    text(k_mean + 0.1, k_std, 'Observed network', 'FontSize', 12);
    xlabel('<k> (mean degree)');
    ylabel('\sigma (std dev of degree)');
    title('\sigma vs <k> with Erdős-Rényi comparison');
    grid on; box off;
    legend({'Erdős-Rényi expectation \sigma = \surd<k>', 'Observed network'}, 'Location', 'NorthWest');

    % Output print
    fprintf('Estimated gamma (power law): %.3f\n', gamma);
    if gamma > 1 && gamma < 3
        fprintf('Network is scale-free according to gamma criterion.\n');
    else
        fprintf('Network is NOT scale-free according to gamma criterion.\n');
    end

    fprintf('<k> = %.3f\n', k_mean);
    fprintf('σ = %.3f\n', k_std);
    fprintf('R2 poisson = %.3f\n', R_squared_poisson);
    fprintf('R2 power = %.3f\n', R_squared_power);
    fprintf('R2 truc power = %.3f\n', gof_trunc.rsquare);

    analysis_summary(end+1,5:6) = {"gamma", gamma};
    analysis_summary(end+1,5:6) = {"k_mean", k_mean};
    analysis_summary(end+1,5:6) = {"k_std", k_std};
    save('_all_result.mat', 'params','analysis_summary' );


    
    % ----------------------- R setting ----------------------------------

    % Here we need Rs to process the powerlaw test
    % setting the environment
    input_csv = [tempname, '.csv'];
    output_json = [tempname, '.json'];
    
    %  CSV data for R
    T = table(k_all, N_k_all, 'VariableNames', {'k', 'freq'});
    writetable(T, input_csv);

    %  Rscript.exe path
    rscript = '"C:\Program Files\R\R-4.5.0\bin\x64\Rscript.exe"';
    
    %  R code path
    rfile = '"D:\BaiduSyncdisk\Side Project\Code for getting cell corr\current_pipline_for_blinking_analysis\pipline\function_module\function_module\powerlaw_analysis.R"';
    
    % cmd
    cmd = sprintf('%s %s "%s" "%s"', rscript, rfile, input_csv, output_json);
    
    % -----R analysisi  using maximum-likelihood fitting----
    status = system(cmd);
    if status ~= 0
        error('Power-law hypothesis rejected (p2 <= 0.1).');
    end
    
    % get the analysis results from  R JSON 
    json_data = jsondecode(fileread(output_json));
    
    % output
    % disp(['xmin: ', num2str(json_data.xmin_initial)]);
    % disp(['alpha: ', num2str(json_data.alpha_initial)]);
    disp(['p-value (KS): ', num2str(json_data.p_value_ks)]);
    disp(['p-value (LR): ', num2str(json_data.p_value_lr)]);
    disp(['alpha: ', num2str(json_data.alpha_final)]);
    % final xmin
    xmin_final = json_data.xmin_final;
    fprintf('Final xmin: %d (Data points: %d)\n', xmin_final, sum(k_all >= xmin_final));

    % decision on powerlaw
        powerlaw_plausible = (json_data.p_value_ks > 0.1);
    if powerlaw_plausible
        disp('Power-law is a plausible hypothesis for the data (p2 > 0.1).');
    else
        disp('Power-law hypothesis rejected (p2 <= 0.1).');
    end
    
    % decision on powerlaw better or not
        powerlaw_better = (json_data.p_value_lr < 0.05) && (json_data.LR > 0);
    if powerlaw_better
        disp('Power-law is a significantly different fit with Poisson (p1 < 0.025).');
    else
        disp('Power-law is not significantly different with Poisson (p1 >= 0.025).');
    end


    
    % traunc powerlaw results
    % disp(['Truncated alpha: ', num2str(json_data.truncated_alpha)]);
    disp(['Truncated lambda: ', num2str(json_data.truncated_lambda)]);
    disp(['p-value (KS Trunc): ', num2str(json_data.p_value_ks_trunc)]);
    disp(['p-value (LR Trunc): ', num2str(json_data.p_value_lr_trunc)]);
    disp(['p-value (poneside): ', num2str(json_data.p_one_sided)]);
    disp(['p-value (ptwoside): ', num2str(json_data.p_two_sided)]);
    disp(['LR: ', num2str(json_data.LR)]);
    powerlaw_traunc_plausible = (json_data.p_value_ks_trunc > 0.1);
    if powerlaw_traunc_plausible 
        disp('Traunc-power-law is a plausible hypothesis for the data (p2 > 0.1).');
    else
        disp('Power-law hypothesis rejected (p2 <= 0.1).');
    end
     powerlaw_traunc_better = (json_data.p_value_lr_trunc < 0.025);
    if powerlaw_traunc_better
        disp('Power-law is a significantly better fit than Poisson (p1 < 0.025).');
    else
        disp('Power-law is not significantly better than Poisson (p1 >= 0.025).');
    end    % Save results

    analysis_summary(end+1,5:6) = {"R^2_powerlaw", R_squared_power};
    analysis_summary(end+1,5:6) = {"R^2_truncated_powerlaw", gof_trunc.rsquare};
    analysis_summary(end+1,5:6) = {"R^2_poisson", R_squared_poisson};
    analysis_summary(end+1,5:6) = {"R^2_lognormal", R_squared_logn};
    analysis_summary(end+1,5:6) = {"truncated_gamma", fit_trunc.gamma};
    analysis_summary(end+1,5:6) = {"truncated_lambda", fit_trunc.lambda};
    analysis_summary(end+1,5:6) = {"p_LR(p1)", json_data.p_value_lr};
    analysis_summary(end+1,5:6) = {"p_KS(p2)", json_data.p_value_ks};
    % analysis_summary(end+1,5:6) = {"xmin", json_data.xmin_initial};
    analysis_summary(end+1,5:6) = {"xmin_final", xmin_final};
    analysis_summary(end+1,5:6) = {"p_LR_traunc(p1)", json_data.p_value_lr_trunc};
    analysis_summary(end+1,5:6) = {"p_KS_traunc(p2)", json_data.p_value_ks_trunc};
    analysis_summary(end+1,5:6) = {"decision_on_power", powerlaw_plausible};
    analysis_summary(end+1,5:6) = {"power_better_poision", powerlaw_better};
    analysis_summary(end+1,5:6) = {"decision_on_traunc_power", powerlaw_traunc_plausible};
    analysis_summary(end+1,5:6) = {"traunc_power_better_power", powerlaw_traunc_better};

    params.powerlaw_fit = P_power;
    params.truncated_powerlaw_fit = P_trunc;
    params.truncated_gamma = fit_trunc.gamma;
    params.truncated_lambda = fit_trunc.lambda;
    params.poisson_fit = P_poisson_plot;
    params.lognormal_fit = P_logn;
end