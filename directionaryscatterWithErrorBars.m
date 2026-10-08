
 
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
% 2025.03.09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function directionaryscatterWithErrorBars()
% directionaryScatterWithErrorBars_KW_SD: Scatter + mean±SD plot for three cell groups
%
% Reads Excel file with sheet 'direction' (columns B-D: all cells, periodic, non-periodic)
% Removes rows containing outliers in any column
% Plots scatter points with mean ± SD as "Wang-style" error bars
% Performs Kruskal-Wallis test with Dunn's post-hoc
% Saves statistical results to 'KW_results.mat', including all pairwise p-values
%
% Usage:
%   directionaryscatterWithErrorBars_KW_SD()
%
% Output:
%   - Figure with scatter points + mean ± SD
%   - Kruskal-Wallis test and Dunn's post-hoc results
%   - Results saved as 'KW_results.mat'
%
% Author: Thomas Broggini & Xiao Liu
% Date: 2025.03.09

%% 1. Select Excel file
[filename, pathname] = uigetfile('*.xls;*.xlsx', 'Select Excel file');
if isequal(filename,0)
    disp('User canceled file selection.');
    return;
end
fullpath = fullfile(pathname, filename);

%% 2. Read 'direction' sheet and remove outliers
sheetname = 'direction';

% Read numeric data from columns B-D (2-4), starting from row 2
data_raw = xlsread(fullpath, sheetname, 'B2:D12'); 

% Remove NaN rows first
data_noNaN = data_raw(~any(isnan(data_raw),2), :);

% Detect outliers for each column
outlier_mask = any(isoutlier(data_noNaN), 2); % any row contains outlier

% Remove rows containing outliers
data = data_noNaN(~outlier_mask, :);

labels = {'All cells', 'Periodic cells', 'Non-periodic cells'};
nGroups = size(data,2);

%% 3. Plotting scatter + mean ± SD ("Wang-style")
figure; hold on;
groupColors = lines(nGroups);
scatter_width = 0.15; % x-axis jitter
bar_width = 0.15;     % half-width for error lines

for i = 1:nGroups
    % Scatter points
    x_scatter = (rand(size(data,1),1)*scatter_width) + i - 0.3;
    scatter(x_scatter, data(:,i), 30, ...
        'MarkerEdgeColor', groupColors(i,:), ...
        'MarkerFaceColor', groupColors(i,:));
    
    % Mean ± SD
    mean_val = mean(data(:,i));
    sd_val = std(data(:,i));
    
    x_center = i + 0.1;
    y_top = mean_val + sd_val;
    y_bottom = mean_val - sd_val;
    x_left = x_center - bar_width;
    x_right = x_center + bar_width;
    
    % Vertical line
    plot([x_center, x_center], [y_bottom, y_top], 'Color', groupColors(i,:), 'LineWidth', 1.5);
    % Top horizontal
    plot([x_left, x_right], [y_top, y_top], 'Color', groupColors(i,:), 'LineWidth', 1.5);
    % Bottom horizontal
    plot([x_left, x_right], [y_bottom, y_bottom], 'Color', groupColors(i,:), 'LineWidth', 1.5);
    % Mean line
    plot([x_left, x_right], [mean_val, mean_val], 'Color', groupColors(i,:), 'LineWidth', 2);
end

xlim([0.5, nGroups+0.5]);
set(gca, 'XTick', 1:nGroups, 'XTickLabel', labels);
ylabel('Value');
title('Scatter + Mean ± SD ("Wang-style")');
box off;

%% 4. Statistical analysis: Kruskal-Wallis + Dunn's post-hoc
[p_kw, tbl, stats] = kruskalwallis(data, labels, 'off');

% Dunn's test (pairwise comparison)
results_dunn = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');

% results_dunn columns: [group1, group2, lowerCI, diff, upperCI, p-value]
% Save all pairwise comparisons
pairwise_info = table(results_dunn(:,1), results_dunn(:,2), results_dunn(:,4), results_dunn(:,6), ...
    'VariableNames', {'Group1','Group2','Difference','pValue'});

% Display Kruskal-Wallis p-value on plot
y_max = max(data, [], 'all')*1.1;
text(2, y_max, sprintf('Kruskal-Wallis p=%.4f', p_kw), ...
    'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');

hold off;

%% 5. Save all results
save('KW_results.mat', 'p_kw', 'tbl', 'stats', 'results_dunn', 'pairwise_info', 'data', 'labels');

end







% % 选择xls文件
% [filename, pathname] = uigetfile('*.xls;*.xlsx', 'Select Excel file');
% if isequal(filename,0)
%     disp('User canceled file selection.');
%     return;
% end
% fullpath = fullfile(pathname, filename);
% 
% % 读取'sheet'名字为'direction'的数据
% sheetname = 'direction';
% 
% % 先读第一行，获取前三列标签
% [~, txt] = xlsread(fullpath, sheetname, 'A1:C1');
% labels = txt(1,:);
% 
% % 读第二行开始的前三列数据
% data = xlsread(fullpath, sheetname, 'A2:C1000'); % 假设最多读取1000行
% 
% % 去除NaN的行
% data = data(~any(isnan(data(:,1:3)),2), :);
% 
% figure; hold on;
% groupColors = lines(3);
% nGroups = 3;
% scatter_width = 0.15; % 散点半边宽度
% bar_width = 0.15;     % 误差棒半边宽度
% 
% for i = 1:nGroups
%     % 散点：x轴左半边
% x_scatter = (rand(size(data,1),1)*scatter_width) + i - 0.3;
% scatter(x_scatter, data(:,i), 30, ...
%     'MarkerEdgeColor', groupColors(i,:), ...
%     'MarkerFaceColor', groupColors(i,:));  % 实心点
%     % 计算均值和标准差
%     mean_val = mean(data(:,i));
%     std_val = std(data(:,i));
% 
%     % “王”字形误差棒（无柱子，只线）
%     x_center = i + 0.1; % 误差棒偏右半边，离散点近
%     y_top = mean_val + std_val;
%     y_bottom = mean_val - std_val;
%     x_left = x_center - bar_width;
%     x_right = x_center + bar_width;
% 
%     % 竖线
%     plot([x_center, x_center], [y_bottom, y_top], 'Color', groupColors(i,:), 'LineWidth', 1.5);
%     % 上横线
%     plot([x_left, x_right], [y_top, y_top], 'Color', groupColors(i,:), 'LineWidth', 1.5);
%     % 下横线
%     plot([x_left, x_right], [y_bottom, y_bottom], 'Color', groupColors(i,:), 'LineWidth', 1.5);
%     % 中间均值横线
%     plot([x_left, x_right], [mean_val, mean_val], 'Color', groupColors(i,:), 'LineWidth', 2);
% end
% 
% xlim([0.5, nGroups+0.5]);
% set(gca, 'XTick', 1:nGroups, 'XTickLabel', labels);
% ylabel('Value');
% title('Scatter + Mean with "Wang" style Error Bars (No bars)');
% 
% box off;
% 
% % 统计分析
% % 1. ANOVA
% [p_anova, ~, ~] = anova1(data, labels, 'off');
% 
% % 2. 后两组t检验
% [~, p_ttest] = ttest2(data(:,2), data(:,3));
% 
% 
% 
% % 显示ANOVA星号在所有柱子正上方中间
% y_max = max(data, [], 'all') * 1.1;
% text(2, y_max, stars(p_anova), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
% 
% % 后两组比较的星号及连线
% y_line_height = y_max * 0.95;
% plot([2.1, 3.1], [y_line_height, y_line_height], 'k-', 'LineWidth', 1.5);
% plot([2.1, 2.1], [y_line_height, y_line_height-0.02*y_max], 'k-', 'LineWidth', 1.5);
% plot([3.1, 3.1], [y_line_height, y_line_height-0.02*y_max], 'k-', 'LineWidth', 1.5);
% text(2.6, y_line_height*1.01, stars(p_ttest), 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
% 
% hold off;
% 
% 
% function s = stars(p)
%     if p < 0.0001
%         s = '****';
%     elseif p < 0.001
%         s = '***';
%     elseif p < 0.01
%         s = '**';
%     elseif p < 0.05
%         s = '*';
%     else
%         s = 'ns';
%     end
% end
