
 
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


function results = analyze_co_culture()

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

% === 文件选择 / File selection ===
[filename, pathname] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel File');
if isequal(filename, 0)
    disp('文件选择已取消 / File selection canceled');
    return;
end
file = fullfile(pathname, filename);
sheetname = 1;

% === 指标 / Metrics ===
metrics = {'Peaks per 100 cells in 60s', ...
           'Mean peak width of all cells (s)'};

% === 数据位置 / Data check ===
data_positions = struct( ...
    'Peaks', struct( ...
        'Original', [4, 'A', 'E'], ...
        'F10_N2', [4, 'F', 'J'], ...
        'F10_N3', [12, 'A', 'E'], ...
        'Stable_N3', [12, 'F', 'J']), ...
    'Width', struct( ...
        'Original', [7, 'A', 'E'], ...
        'F10_N2', [7, 'F', 'J'], ...
        'F10_N3', [15, 'A', 'E'], ...
        'Stable_N3', [15, 'F', 'J']) ...
);

group_names = {'Original F10','F10.N2','F10.N3','Stable F10.N3'};

% === 颜色 / Colors ===
colors = [
    0.2 0.2 0.2;
    0.2 0.5 0.9;
    0.9 0.3 0.3;
    0.2 0.7 0.4
];

results = struct();

for m = 1:length(metrics)
    metric_name = metrics{m};
    fprintf('\n正在处理指标 / Processing metric: %s\n', metric_name);

    % ==== 读取位置 / Select data positions ====
    if contains(metric_name, 'Peaks')
        pos = data_positions.Peaks;
    else
        pos = data_positions.Width;
    end

    % ==== 读取数据 / Read data ====
    Original = readmatrix(file, 'Sheet', sheetname, ...
        'Range', sprintf('%s%d:%s%d', pos.Original(2), pos.Original(1), pos.Original(3), pos.Original(1)));
    F10_N2 = readmatrix(file, 'Sheet', sheetname, ...
        'Range', sprintf('%s%d:%s%d', pos.F10_N2(2), pos.F10_N2(1), pos.F10_N2(3), pos.F10_N2(1)));
    F10_N3 = readmatrix(file, 'Sheet', sheetname, ...
        'Range', sprintf('%s%d:%s%d', pos.F10_N3(2), pos.F10_N3(1), pos.F10_N3(3), pos.F10_N3(1)));
    Stable_N3 = readmatrix(file, 'Sheet', sheetname, ...
        'Range', sprintf('%s%d:%s%d', pos.Stable_N3(2), pos.Stable_N3(1), pos.Stable_N3(3), pos.Stable_N3(1)));

    groups = {Original, F10_N2, F10_N3, Stable_N3};
    all_data = [Original(:); F10_N2(:); F10_N3(:); Stable_N3(:)];
    group_idx = [repmat(1, numel(Original), 1);
                 repmat(2, numel(F10_N2), 1);
                 repmat(3, numel(F10_N3), 1);
                 repmat(4, numel(Stable_N3), 1)];

    % ==== 正态性检验 / Normality test ====
    normal_flags = false(1,4);
    for g = 1:4
        data_g = groups{g};
        data_g = data_g(~isnan(data_g));
        if numel(data_g) < 3
            normal_flags(g) = true;
        else
            try
                [h, ~] = swtest(data_g, 0.05);
            catch
                [h, ~] = kstest((data_g - mean(data_g)) / std(data_g));
            end
            normal_flags(g) = (h == 0);
        end
    end

    %% ==== 一、ANOVA + Tukey ====
    [p_anova, ~, stats_anova] = anova1(all_data, group_idx, 'off');
    comp_tukey = multcompare(stats_anova, 'CType', 'tukey-kramer', 'Display', 'off');
    tukey_results = {};
    for i = 1:size(comp_tukey,1)
        g1 = group_names{comp_tukey(i,1)};
        g2 = group_names{comp_tukey(i,2)};
        tukey_results{i,1} = sprintf('%s vs %s', g1, g2);
        tukey_results{i,2} = comp_tukey(i,6); % p值
        tukey_results{i,3} = comp_tukey(i,4); % 均值差
        tukey_results{i,4} = [comp_tukey(i,3), comp_tukey(i,5)]; % 置信区间
    end
    tukey_table = cell2table(tukey_results, ...
        'VariableNames', {'Comparison','pValue','MeanDiff','95CI'});

    %% ==== 二、Kruskal-Wallis + Dunn ====
    [p_kw, ~, stats_kw] = kruskalwallis(all_data, group_idx, 'off');
    comp_dunn = multcompare(stats_kw, 'CType', 'dunn-sidak', 'Display', 'off');
    dunn_results = {};
    for i = 1:size(comp_dunn,1)
        g1 = group_names{comp_dunn(i,1)};
        g2 = group_names{comp_dunn(i,2)};
        dunn_results{i,1} = sprintf('%s vs %s', g1, g2);
        dunn_results{i,2} = comp_dunn(i,6); % p
        dunn_results{i,3} = NaN; % 
        dunn_results{i,4} = NaN; % 
    end
    dunn_table = cell2table(dunn_results, ...
        'VariableNames', {'Comparison','pValue','MeanDiff','95CI'});

    %% ==== 整理输出 / Save to struct ====
    results(m).Metric = metric_name;
    results(m).GroupNames = group_names;
    results(m).Normality = normal_flags;

    results(m).ANOVA_Tukey.p_main = p_anova;
    results(m).ANOVA_Tukey.Pairwise = tukey_table;

    results(m).Kruskal_Dunn.p_main = p_kw;
    results(m).Kruskal_Dunn.Pairwise = dunn_table;

    %% ==== 绘图 / Plot ====
    figure('Color','w','Position',[200,200,700,550]); hold on;
    for g = 1:4
        data_g = groups{g};
        data_g = data_g(~isnan(data_g));
        x_jitter = g + (rand(size(data_g)) - 0.5)*0.15;
        scatter(x_jitter, data_g, 70, 'filled', 'MarkerFaceColor', colors(g,:));

        % 小样本 (<=5) 用 SD，大样本 (>5) 用 SEM
        if ~isempty(data_g)
            n = numel(data_g);
            if n <= 5
                err = std(data_g);  % SD
            else
                err = std(data_g)/sqrt(n);  % SEM
            end
            errorbar(g, mean(data_g), err, 'k', 'CapSize', 15, 'LineWidth', 1.5);
        end
    end
    xlim([0.5 4.5]);
    xticks(1:4);
    xticklabels(group_names);
    ylabel(metric_name, 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('%s\nANOVA p = %.4g, KW p = %.4g', metric_name, p_anova, p_kw), 'FontSize', 14);
    box on;
end

% === 保存最终结果 / Save results ===
save('stat_results_full.mat', 'results');
fprintf('\n所有分析完成并保存至 stat_results_full_wnt.mat / All analyses completed and saved to stat_results_full.mat\n');

end
