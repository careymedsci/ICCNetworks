
 
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


function results_freq_stats = plot_freq_Min_Max_group_comparison()
   
% 作者：刘晓； 法兰克福大学医院；湖北医药学院襄阳市第一人民医院
    % By Xiao Liu
    % Department of Neurosurgery
    % Neuroscience Centre
    % University Hospital Frankfurt
    % Goethe University Frankfurt, Germany
    % Frankfurt Cancer Institute, Germany
    % Xiangyang No.1 people's Hospital, China
    % Also an python version coded 
    % Xiao.Liu@stud.uni-frankfurt.de
    % 2024.10.09
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot violin plots for multiple group comparisons and automatically select 
% Welch t test or Kruskal-Wallis tests, plus mixed ANOVA with post-hoc if parametric

% ==== 选择 Sheet name / Select Sheet ====
clear;
sheet_options = {'cbx', 'gap26', 'suramin', 'apyrase', 'EGTA', 'Vera', 'skf', 'dan','AP5',...
    'NieCl2', 'fgf', 'tram34'};
[selection, ok] = listdlg('PromptString', 'choose（Sheet Name）:', ...
                          'SelectionMode', 'single', ...
                          'ListString', sheet_options);
if ~ok
    disp('Sheet cancelled');
    results_freq_stats = [];
    return;
end
sheetname = sheet_options{selection};
fprintf('Sheet name: %s\n', sheetname);

[filename, pathname] = uigetfile({'*.xls;*.xlsx'}, 'choose the Excel file');
if isequal(filename, 0)
    disp('cancel');
    results_freq_stats = [];
    return;
end
file = fullfile(pathname, filename);

% loading the right Sheet data
data_all = readmatrix(file, 'Sheet', sheetname);
rows = [28; 24; 25; 29]; % 对应4组 choose the cells string the data 

% ==== of there is 3 groups ====
if any(~isnan(data_all(rows,13:end)),'all')
    % if there is 3 groups, do the analysis by another function
    results_freq_stats = plot_freq_Min_Max_group_comparison_3groups(file, sheetname);
    return;
end

% or we will do the analysis by this function
data = data_all(rows, 1:12); % control + drug treating

% 颜色定义 | Color definitions
color_control = [0, 113, 188]/255;  
color_cbx = [163, 30, 50]/255;      
violin_color_control = [174, 215, 255]/255; 
violin_color_cbx = [255, 188, 188]/255;  

figure; 
hold on;

% =============================
% 准备数据 | Prepare data
% =============================
group_data = cell(1, 8); 
group_means = zeros(1, 8); 
group_stderr = zeros(1, 8); 
normality_p = zeros(1, 8); 

% Control (1,3,5,7)
for i = 1:2:7
    row_idx = ceil(i/2);
    raw_data = data(row_idx, 1:6);
    valid_data = raw_data(~isnan(raw_data));
    group_data{i} = valid_data;
    
    group_means(i) = mean(valid_data);
    group_stderr(i) = std(valid_data) / sqrt(numel(valid_data));
    
    if numel(valid_data) >= 3
        [~, normality_p(i)] = lillietest(valid_data);
    else
        normality_p(i) = NaN;
    end
end

% Drug (2,4,6,8)
for i = 2:2:8
    row_idx = i/2;
    raw_data = data(row_idx, 7:12);
    valid_data = raw_data(~isnan(raw_data));
    group_data{i} = valid_data;
    
    group_means(i) = mean(valid_data);
    group_stderr(i) = std(valid_data) / sqrt(numel(valid_data));
    
    if numel(valid_data) >= 3
        [~, normality_p(i)] = lillietest(valid_data);
    else
        normality_p(i) = NaN;
    end
end

fprintf('\n==== Normality test ====\n');
for i = 1:8
    if isnan(normality_p(i))
        fprintf('Group %d: Insufficient data for normality test (n=%d)\n', i, numel(group_data{i}));
    else
        fprintf('Group %d: p = %.4f %s\n', i, normality_p(i), ...
            ternary(normality_p(i) > 0.05, '(Normal)', '(Non-normal)'));
    end
end

non_normal_count = sum(normality_p < 0.05 & ~isnan(normality_p));
use_nonparametric = non_normal_count >= 4;

fprintf('\nNon-normal groups: %d/%d, using %s test\n', non_normal_count, 8, ...
    ternary(use_nonparametric, 'Mann-Whitney U test', 'Welch’s t-test'));
fprintf ('FDR（False Discovery Rate）by Benjamini-Hochberg');

% =============================
% 绘制图形 | Plot visualization
% =============================
violin_width = 0.6;

for i = 1:8
    current_data = group_data{i};
    if isempty(current_data)
        continue;
    end
    
    if mod(i, 2) == 1
        violin_color = violin_color_control;
        scatter_color = color_control;
    else
        violin_color = violin_color_cbx;
        scatter_color = color_cbx;
    end
    
    [f, xi] = ksdensity(current_data);
    if max(f) > 0
        f = f/max(f)*violin_width/2;
    else
        f = zeros(size(f));
    end
    fill([i+f, i-fliplr(f)], [xi, fliplr(xi)], violin_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    
    jitter = (rand(size(current_data)) - 0.5) * 0.3;
    scatter(i + jitter, current_data, 40, scatter_color, 'filled', 'MarkerEdgeColor', 'none');
    
    plot(i, group_means(i), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
    line([i, i], [group_means(i) - group_stderr(i), group_means(i) + group_stderr(i)], 'Color', 'k', 'LineWidth', 1.5);
    line([i-0.1, i+0.1], [group_means(i) - group_stderr(i), group_means(i) - group_stderr(i)], 'Color', 'k', 'LineWidth', 1.5);
    line([i-0.1, i+0.1], [group_means(i) + group_stderr(i), group_means(i) + group_stderr(i)], 'Color', 'k', 'LineWidth', 1.5);
end

% =============================
% 原有统计分析 | Original pairwise tests
% =============================
group_labels = {
    'C1 (control global cells‘ min frequency)';
    'D1 (treatment global cells‘ min frequency)';
    'C2 (control periodical cells‘ min frequency)';
    'D2 (treatment periodical cells‘ min frequency)';
    'C3 (control periodical cells‘ max frequency)';
    'D3 (treatment periodical cells‘ max frequency)';
    'C4 (control global cells‘ max frequency)';
    'D4 (treatment global cells‘ max frequency)';
};

comparisons = {
    [1, 2], 'C1 vs D1'; 
    [3, 4], 'C2 vs D2'; 
    [5, 6], 'C3 vs D3'; 
    [7, 8], 'C4 vs D4';
    [1, 3], 'C1 vs C2';
    [2, 4], 'D1 vs D2';
    [5, 7], 'C3 vs C4';
    [6, 8], 'D3 vs D4'
};

num_comparisons = size(comparisons,1);
raw_p = zeros(1,num_comparisons);
effect_size = zeros(1,num_comparisons); 
effect_ci = zeros(num_comparisons,2);  

for i = 1:num_comparisons
    idx1 = comparisons{i,1}(1);
    idx2 = comparisons{i,1}(2);
    data1 = group_data{idx1};
    data2 = group_data{idx2};
    
    if mod(idx1,2) ~= mod(idx2,2)
        comp_type = 'independent';
    else
        comp_type = 'paired';
    end
    
    if use_nonparametric
        [p, ~] = ranksum(data1, data2);
        raw_p(i) = p;
    else
        [~, p] = ttest2(data1, data2, 'Vartype', 'unequal');
        raw_p(i) = p;
    end
    
    if strcmp(comp_type,'independent')
        [d, ci] = compute_cohen_d(data1, data2);
        effect_size(i) = d;
        effect_ci(i,:) = ci;
    else
        [dz, ci] = compute_cohen_dz(data1, data2);
        effect_size(i) = dz;
        effect_ci(i,:) = ci;
    end
end

% FDR校正
[sorted_p, sort_idx] = sort(raw_p);
m = length(raw_p);
adjusted_p = zeros(size(raw_p));
for i = 1:m
    adjusted_p(sort_idx(i)) = min(sorted_p(i) * m / i, 1);
end

% =============================
% 混合设计 ANOVA
% =============================
anova_results = struct();
if ~use_nonparametric
    % set 1：Min frequency
    n_c = min(numel(group_data{1}), numel(group_data{3}));
    n_d = min(numel(group_data{2}), numel(group_data{4}));
    data_control = [group_data{1}(1:n_c)', group_data{3}(1:n_c)'];
    data_drug =    [group_data{2}(1:n_d)', group_data{4}(1:n_d)'];
    tbl1 = array2table([data_control; data_drug], 'VariableNames', {'Global', 'Periodical'});
    tbl1.Treatment = [repmat({'Control'}, n_c, 1); repmat({'Drug'}, n_d, 1)];
    WithinDesign1 = table(categorical({'Global';'Periodical'}), 'VariableNames', {'CellType'});
    rm1 = fitrm(tbl1, 'Global-Periodical ~ Treatment', 'WithinDesign', WithinDesign1);
    ranovatbl1 = ranova(rm1, 'WithinModel', 'CellType');
    mc1_celltype = multcompare(rm1, 'CellType', 'By', 'Treatment', 'ComparisonType', 'bonferroni');
    mc1_treatment = multcompare(rm1, 'Treatment', 'By', 'CellType', 'ComparisonType', 'bonferroni');
    anova_results.freq_min.ranova = ranovatbl1;
    anova_results.freq_min.posthoc_celltype = mc1_celltype;
    anova_results.freq_min.posthoc_treatment = mc1_treatment;
    
    % set 2：Max frequency
    n_c = min(numel(group_data{5}), numel(group_data{7}));
    n_d = min(numel(group_data{6}), numel(group_data{8}));
    data_control = [group_data{5}(1:n_c)', group_data{7}(1:n_c)'];
    data_drug =    [group_data{6}(1:n_d)', group_data{8}(1:n_d)'];
    tbl2 = array2table([data_control; data_drug], 'VariableNames', {'Periodical', 'Global'});
    tbl2.Treatment = [repmat({'Control'}, n_c, 1); repmat({'Drug'}, n_d, 1)];
    WithinDesign2 = table(categorical({'Periodical';'Global'}), 'VariableNames', {'CellType'});
    rm2 = fitrm(tbl2, 'Periodical-Global ~ Treatment', 'WithinDesign', WithinDesign2);
    ranovatbl2 = ranova(rm2, 'WithinModel', 'CellType');
    mc2_celltype = multcompare(rm2, 'CellType', 'By', 'Treatment', 'ComparisonType', 'bonferroni');
    mc2_treatment = multcompare(rm2, 'Treatment', 'By', 'CellType', 'ComparisonType', 'bonferroni');
    anova_results.freq_max.ranova = ranovatbl2;
    anova_results.freq_max.posthoc_celltype = mc2_celltype;
    anova_results.freq_max.posthoc_treatment = mc2_treatment;
end

hold off;

% =============================
% 存储统计结果
% =============================
results_freq_stats = struct();
results_freq_stats.group_data = group_data;
results_freq_stats.group_means = group_means;
results_freq_stats.group_stderr = group_stderr;
results_freq_stats.normality_p = normality_p;

results_freq_stats.raw_p = cell(num_comparisons,3);
results_freq_stats.adjusted_p = cell(num_comparisons,3);
results_freq_stats.effect_size = cell(num_comparisons,3);
results_freq_stats.effect_ci = cell(num_comparisons,4);

for i = 1:num_comparisons
    idx1 = comparisons{i,1}(1);
    idx2 = comparisons{i,1}(2);
    
    name1 = group_labels{idx1};
    name2 = group_labels{idx2};
    
    results_freq_stats.raw_p{i,1} = raw_p(i);
    results_freq_stats.raw_p{i,2} = name1;
    results_freq_stats.raw_p{i,3} = name2;
    
    results_freq_stats.adjusted_p{i,1} = adjusted_p(i);
    results_freq_stats.adjusted_p{i,2} = name1;
    results_freq_stats.adjusted_p{i,3} = name2;
    
    results_freq_stats.effect_size{i,1} = effect_size(i);
    results_freq_stats.effect_size{i,2} = name1;
    results_freq_stats.effect_size{i,3} = name2;
    
    results_freq_stats.effect_ci{i,1} = effect_ci(i,1);
    results_freq_stats.effect_ci{i,2} = effect_ci(i,2);
    results_freq_stats.effect_ci{i,3} = name1;
    results_freq_stats.effect_ci{i,4} = name2;
end

results_freq_stats.anova_results = anova_results;
save('freq_results_stats.mat', 'results_freq_stats');

end

%% ======= Helper functions =======
function out = ternary(cond, trueVal, falseVal)
    if cond
        out = trueVal;
    else
        out = falseVal;
    end
end

function [d, ci] = compute_cohen_d(x1, x2)
    n1 = length(x1);
    n2 = length(x2);
    s1 = var(x1); 
    s2 = var(x2);
    pooled_std = sqrt(((n1-1)*s1 + (n2-1)*s2) / (n1+n2-2));
    d = (mean(x1) - mean(x2)) / pooled_std;
    se = sqrt((n1+n2)/(n1*n2) + d^2/(2*(n1+n2)));
    ci = [d - 1.96*se, d + 1.96*se];
end

function [dz, ci] = compute_cohen_dz(x1, x2)
    D = x1 - x2;
    dz = mean(D) / std(D); 
    n = length(D);
    t_val = tinv(0.975, n-1); 
    se = sqrt((1/n) + (dz^2)/(2*n)); 
    ci = [dz - t_val*se, dz + t_val*se];
end
