
 
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

function analyze_periodic_power_distribution(toplot, params, x_coords, y_coords)
% 析Periodic细胞在功率分布中的富集情况
% Analyze Periodic cell enrichment in power distribution

% 作者：Thomas Broggini 和 刘晓；法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Thomas Broggini (Supervisor) and Xiao Liu (Student)
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany

% Xiangyang No.1 people's Hospital, China
% Also an python version coded 
% xiao.liu@stud.uni-frankfurt.de
% 2025.07.06

% Function description:
% This function analyzes Periodic cell enrichment in high power regions, including:
% 1. Visualization of Periodic cell colocalization with high power cells
% 2. Calculation of Periodic cell proportion at different power percentiles
% 3. Comparison of power distribution between Periodic and non-Periodic cells
% 4. Validation of enrichment significance through randomization tests
%
% Input parameters:
%   toplot    - Multitaper analysis result structure with fields: mask, mask_ind, extpwr
%   params    - Parameter structure with field: periodicity
%   x_coords  - X coordinates of all cells (from mFFT analysis)
%   y_coords  - Y coordinates of all cells (from mFFT analysis)
%
% Output:
%   Multiple figure windows displaying analysis results

%
% 功能描述：
% 此函数分析周期性振荡细胞(Periodic cells)在高功率区域的富集情况，包括：
% 1. 可视化Periodic细胞与高功率细胞的共定位
% 2. 计算Periodic细胞在不同功率分位中的比例
% 3. Periodic细胞与非Periodic细胞的功率分布比较
% 4. 通过随机化检验验证富集显著性
%
% 输入参数：
%   toplot    - 多锥谱分析结果结构体，包含mask, mask_ind, extpwr等字段
%   params    - 参数结构体，包含periodicity字段
%   x_coords  - 所有细胞的x坐标（来自mFFT分析）
%   y_coords  - 所有细胞的y坐标（来自mFFT分析）
%
% 输出：
%   生成多个图形窗口展示分析结果

%% 1. 参数验证与初始化 / Parameter validation and initialization
% 检查必要的输入参数 / Check necessary input parameters
if nargin < 4
    error('Insufficient input arguments. Required: toplot, params, x_coords, y_coords');
end

% 验证toplot结构体字段 / Validate toplot structure fields
required_toplot_fields = {'mask', 'mask_ind', 'extpwr'};
for i = 1:length(required_toplot_fields)
    if ~isfield(toplot, required_toplot_fields{i})
        error('toplot structure missing required field: %s', required_toplot_fields{i});
    end
end

% 验证params结构体字段 / Validate params structure fields
if ~isfield(params, 'periodicity')
    error('params structure missing required field: periodicity');
end

fprintf('\n========================================\n');
fprintf('开始分析Periodic细胞功率分布...\n');
fprintf('Starting Periodic cell power distribution analysis...\n');
fprintf('========================================\n\n');

%% 2. 提取Periodic细胞坐标 / Extract Periodic cell coordinates
% 设置缩放因子（预留，当前未使用）/ Scaling factor (reserved, currently not used)
scaling = 1;

% 从params获取Periodic细胞索引 / Get Periodic cell indices from params
periodic_indices = find(params.periodicity == 1);

% 提取Periodic细胞坐标 / Extract Periodic cell coordinates
all_x = x_coords;
all_y = y_coords;
periodic_coords = [all_x(periodic_indices); all_y(periodic_indices)]';

fprintf('检测到 %d 个Periodic细胞（周期性振荡细胞）\n', length(periodic_indices));
fprintf('Detected %d Periodic cells\n', length(periodic_indices));

%% 3. 识别高功率细胞 / Identify high power cells
% 获取功率值 / Get power values
power_values = toplot.extpwr;

% 按功率值降序排序 / Sort power values in descending order
[sorted_power, sorted_indices_number_lable] = sort(power_values, 'descend');

% 选择前50%的细胞作为高功率细胞 / Select top 50% cells as high power cells
num_high_power = ceil(numel(sorted_power) * 0.5);

% 获取高功率细胞的索引 / Get indices of high power cells
sorted_mask_inds = toplot.mask_ind(sorted_indices_number_lable);
top_inds = sorted_mask_inds(1:num_high_power);

% 将线性索引转换为行列坐标 / Convert linear indices to row-column coordinates
[row, col] = ind2sub(size(toplot.mask), top_inds);
high_power_coords = [col(:), row(:)];

fprintf('识别出 %d 个高功率细胞（前50%%）\n', num_high_power);
fprintf('Identified %d high power cells (top 50%%)\n\n', num_high_power);

%% 4. 计算Periodic细胞与高功率细胞的共定位 / Calculate colocalization of Periodic and high power cells
% 设置精度以避免浮点误差 / Set precision to avoid floating point errors
precision = 4;
periodic_coords_round = unique(round(periodic_coords, precision), 'rows');
high_power_coords_round = unique(round(high_power_coords, precision), 'rows');

% 查找共同的坐标 / Find common coordinates
[common_coords, idx_periodic, ~] = intersect(...
    periodic_coords_round, high_power_coords_round, 'rows');

% 计算共定位比例 / Calculate colocalization proportion
num_common = size(common_coords, 1);
num_periodic_total = size(periodic_coords, 1);
percentage = (num_common / num_periodic_total) * 100;

fprintf('共定位统计：\n');
fprintf('Colocalization statistics:\n');
fprintf('  Periodic细胞总数: %d\n', num_periodic_total);
fprintf('  Total Periodic cells: %d\n', num_periodic_total);
fprintf('  共定位细胞数: %d\n', num_common);
fprintf('  Colocalized cells: %d\n', num_common);
fprintf('  共定位比例: %.1f%%\n', percentage);
fprintf('  Colocalization proportion: %.1f%%\n\n', percentage);

%% 5. 可视化Periodic细胞与高功率细胞共定位 / Visualize Periodic and high power cell colocalization
figure('Name', 'Periodic细胞与高功率细胞共定位', 'Position', [300 300 800 600]);
% Periodic Cells and High Power Cells Colocalization

% 绘制所有细胞（浅灰色）/ Plot all cells (light gray)
scatter(all_x, all_y, 20, [0.7 0.7 0.7], 'filled');
set(gca, 'YDir', 'reverse');  % 反转Y轴方向匹配图像坐标 / Reverse Y-axis for image coordinates
hold on;

% 绘制Periodic细胞（深红色）/ Plot Periodic cells (dark red)
scatter(periodic_coords(:,1), periodic_coords(:,2), ...
        20, [0.6431, 0.1176, 0.1961], 'filled');
set(gca, 'YDir', 'reverse');

% 绘制高功率细胞（蓝色）/ Plot high power cells (blue)
scatter(high_power_coords(:,1), high_power_coords(:,2), ...
        20, [0, 0.4431, 0.7373], 'filled');
set(gca, 'YDir', 'reverse');

% 绘制共定位细胞（黄色）/ Plot colocalized cells (yellow)
if ~isempty(common_coords)
    scatter(common_coords(:,1), common_coords(:,2), ...
            35, [0.9290 0.6940 0.1250], 'filled', ...
            'MarkerEdgeColor', [0.9290 0.6940 0.1250]);
    set(gca, 'YDir', 'reverse');
end

% 美化图形 / Beautify figure
axis image;
title(sprintf('Periodic细胞在高功率区域的富集 (%.1f%%)', percentage));
% Periodic Cell Enrichment in High Power Regions
xlabel('X (μm)');
ylabel('Y (μm)');
legend({'All Cells', 'Periodic Cells', 'High Power Cells', 'Colocalized Cells'}, ...
       'Location', 'bestoutside');
grid on;

%% 6. 分析不同功率分位中Periodic细胞的比例 / Analyze Periodic cell proportion at different power percentiles
percentiles = 0.1:0.1:1.0;  % 从10%到100% / From 10% to 100%
prop_periodic = zeros(1, length(percentiles));  % 预分配数组 / Preallocate array

fprintf('计算不同功率分位数的富集曲线...\n');
fprintf('Calculating enrichment curve at different power percentiles...\n');

for i = 1:length(percentiles)
    % 计算当前分位数对应的细胞数 / Calculate number of cells for current percentile
    num_high_power = ceil(numel(sorted_power) * percentiles(i));
    
    % 获取高功率细胞索引和坐标 / Get high power cell indices and coordinates
    top_indices = sorted_indices_number_lable(1:num_high_power);
    top_inds = sorted_mask_inds(1:num_high_power);
    [row, col] = ind2sub(size(toplot.mask), top_inds);
    high_power_coords = [col(:), row(:)];
    
    % 坐标去重并四舍五入 / Remove duplicates and round coordinates
    precision = 4;
    periodic_coords_round = unique(round(periodic_coords, precision), 'rows');
    high_power_coords_round = unique(round(high_power_coords, precision), 'rows');
    
    % 查找共定位细胞 / Find colocalized cells
    [common_coords, ~, ~] = intersect(periodic_coords_round, high_power_coords_round, 'rows');
    
    % 计算比例 / Calculate proportion
    num_common = size(common_coords, 1);
    num_periodic_total = size(periodic_coords_round, 1);
    prop_periodic(i) = num_common / num_periodic_total;
    
    % 显示进度 / Show progress
    fprintf('  分位数 %.0f%%: 共定位比例 = %.3f\n', percentiles(i)*100, prop_periodic(i));
end

% 绘制富集曲线 / Plot enrichment curve
figure('Name', 'Periodic细胞富集曲线', 'Position', [100 100 800 600]);
% Periodic Cell Enrichment Curve
plot(percentiles * 100, prop_periodic, 'o-', 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]);
xlabel('Top X% High-Power Cells');
ylabel('Proportion of Periodic Cells');
title('Enrichment of Periodic Cells in High-Power Regions');
grid on;

fprintf('富集曲线计算完成\n\n');
fprintf('Enrichment curve calculation completed\n\n');

%% 7. Periodic细胞与非Periodic细胞的功率分布比较 / Compare power distribution between Periodic and non-Periodic cells
fprintf('比较Periodic细胞与非Periodic细胞的功率分布...\n');
fprintf('Comparing power distribution between Periodic and non-Periodic cells...\n');

% 获取所有细胞的坐标（去重）/ Get coordinates of all cells (unique)
all_coords = [all_x(:), all_y(:)];
all_coords_round = unique(round(all_coords, precision), 'rows');

% 识别Periodic细胞在去重坐标中的索引 / Identify Periodic cell indices in unique coordinates
[~, idx_periodic_all] = intersect(all_coords_round, periodic_coords_round, 'rows');
periodic_mask = false(size(all_coords_round, 1), 1);
periodic_mask(idx_periodic_all) = true;

% 提取功率值 / Extract power values
power_all = power_values;
periodic_power = power_all(periodic_mask);
non_periodic_power = power_all(~periodic_mask);

% Kolmogorov-Smirnov检验 / Kolmogorov-Smirnov test
[~, p_ks, ks_stat] = kstest2(periodic_power, non_periodic_power);

fprintf('\n功率分布比较结果：\n');
fprintf('Power distribution comparison results:\n');
fprintf('  Periodic细胞数: %d\n', length(periodic_power));
fprintf('  Periodic cells: %d\n', length(periodic_power));
fprintf('  非Periodic细胞数: %d\n', length(non_periodic_power));
fprintf('  Non-Periodic cells: %d\n', length(non_periodic_power));
fprintf('  KS检验统计量: %.3f\n', ks_stat);
fprintf('  KS test statistic: %.3f\n', ks_stat);
fprintf('  p值: %.2e\n', p_ks);
fprintf('  p-value: %.2e\n', p_ks);

% 绘制累积分布函数图 / Plot cumulative distribution function
figure('Name', 'Periodic细胞功率分布CDF比较', 'Position', [100 100 800 600]);
% Power Distribution CDF Comparison for Periodic Cells

h1 = cdfplot(periodic_power);
hold on;
h2 = cdfplot(non_periodic_power);

% 设置线型属性 / Set line properties
set(h1, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]);
set(h2, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);

% 添加统计信息注释 / Add statistical information annotation
annotation_text = {
    sprintf('KS检验: D = %.3f, p = %.2e', ks_stat, p_ks)
    sprintf('KS test: D = %.3f, p = %.2e', ks_stat, p_ks)
};

text(0.05, 0.6, annotation_text, 'Units', 'normalized', ...
     'FontSize', 10, 'BackgroundColor', [1 1 1 0.7]);

% 添加图例和标题 / Add legend and title
legend({'Periodic Cells', 'Non-Periodic Cells'}, 'Location', 'best');
title(sprintf('功率分布比较 (KS检验 p = %.2e)', p_ks));
% Power Distribution Comparison (KS test p = ...)
xlabel('Power');
ylabel('Cumulative Probability');
grid on;
set(gca, 'YColor', 'k');

%% 8. 随机化检验计算95%置信区间 / Randomization test to calculate 95% confidence interval
fprintf('\n进行随机化检验...\n');
fprintf('Performing randomization test...\n');

% 重新获取所有细胞的坐标和功率值 / Re-get coordinates and power values of all cells
sorted_mask_inds = toplot.mask_ind;
[row_all, col_all] = ind2sub(size(toplot.mask), sorted_mask_inds);
coords_all = round([col_all(:), row_all(:)], 4);

% 获取Periodic细胞索引和数量 / Get Periodic cell indices and count
[~, periodic_idx, ~] = intersect(coords_all, periodic_coords, 'rows');
n_periodic = length(periodic_idx);

% 设置随机化参数 / Set randomization parameters
n_iter = 500;  % 随机化迭代次数 / Number of randomization iterations

fprintf('  随机化迭代次数: %d\n', n_iter);
fprintf('  Number of randomization iterations: %d\n', n_iter);

% 重新计算富集曲线 / Recalculate enrichment curve
percentiles = 0.1:0.1:1.0;
prop_periodic_actual = zeros(length(percentiles), 1);
periodic_coords_round = unique(periodic_coords, 'rows');

[sorted_power, sorted_indices_number_lable] = sort(power_values, 'descend');
sorted_mask_inds = toplot.mask_ind(sorted_indices_number_lable);

for i = 1:length(percentiles)
    num_top = ceil(numel(sorted_power) * percentiles(i));
    top_inds = sorted_mask_inds(1:num_top);
    [row, col] = ind2sub(size(toplot.mask), top_inds);
    high_power_coords = unique(round([col(:), row(:)], 4), 'rows');
    [common_coords, ~, ~] = intersect(periodic_coords_round, high_power_coords, 'rows');
    prop_periodic_actual(i) = size(common_coords, 1) / size(periodic_coords_round, 1);
end

% 随机化检验 / Randomization test
prop_random_all = zeros(length(percentiles), n_iter);

for iter = 1:n_iter
    % 随机选择与Periodic细胞数量相同的细胞 / Randomly select same number of cells as Periodic cells
    rand_idx = randsample(size(coords_all, 1), n_periodic);
    rand_coords = coords_all(rand_idx, :);
    
    % 计算每个分位数的比例 / Calculate proportion for each percentile
    for j = 1:length(percentiles)
        num_top = ceil(numel(sorted_power) * percentiles(j));
        top_inds = sorted_mask_inds(1:num_top);
        [row, col] = ind2sub(size(toplot.mask), top_inds);
        high_power_coords = unique(round([col(:), row(:)], 4), 'rows');
        [common_coords, ~, ~] = intersect(round(rand_coords, 4), high_power_coords, 'rows');
        prop_random_all(j, iter) = size(common_coords, 1) / n_periodic;
    end
    
    % 显示进度 / Show progress
    if mod(iter, 100) == 0
        fprintf('    已完成 %d/%d 次迭代\n', iter, n_iter);
        fprintf('    Completed %d/%d iterations\n', iter, n_iter);
    end
end

% 计算随机化的统计量 / Calculate randomization statistics
mean_prop_rand = mean(prop_random_all, 2);
std_prop_rand = std(prop_random_all, 0, 2);
ci_low = mean_prop_rand - 1.96 * std_prop_rand;
ci_high = mean_prop_rand + 1.96 * std_prop_rand;

% 绘制富集曲线与置信区间 / Plot enrichment curve with confidence interval
figure('Name', 'Periodic细胞随机化检验富集曲线', 'Position', [100 100 900 700]);
% Randomization Test Enrichment Curve for Periodic Cells
hold on;

% 绘制95%置信区间 / Plot 95% confidence interval
fill([percentiles * 100, fliplr(percentiles * 100)], ...
     [ci_low', fliplr(ci_high')], [0.9 0.9 0.9], 'EdgeColor', 'none', ...
     'DisplayName', '95% CI (Random)');

% 绘制随机化平均值 / Plot randomization mean
plot(percentiles * 100, mean_prop_rand, 'k--', 'LineWidth', 1.5, ...
     'DisplayName', 'Mean Random');

% 绘制实际Periodic细胞富集曲线 / Plot actual Periodic cell enrichment curve
plot(percentiles * 100, prop_periodic_actual, 'r-', 'LineWidth', 2, ...
     'DisplayName', 'Real Periodic');

% 美化图形 / Beautify figure
xlabel('Top X% High-Power Cells');
ylabel('Proportion of Periodic Cells');
title('富集曲线：Periodic细胞 vs 随机基线');
% Enrichment Curve: Periodic Cells vs Random Baseline
legend('Location', 'best');
grid on;
box on;

% 添加显著性标注 / Add significance annotation
for i = 1:length(percentiles)
    if prop_periodic_actual(i) > ci_high(i)
        text(percentiles(i) * 100, prop_periodic_actual(i) + 0.02, '*', ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    end
end

fprintf('\n========================================\n');
fprintf('Periodic细胞功率分布分析完成！\n');
fprintf('Periodic cell power distribution analysis completed!\n');
fprintf('生成了 %d 个图形窗口\n', 4);
fprintf('Generated %d figure windows\n', 4);
fprintf('========================================\n');

end  % 函数结束 / End of function