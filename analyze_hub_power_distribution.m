
 
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

function analyze_hub_power_distribution(toplot, params, x_coords, y_coords)

% Analyze hub cell in power distribution
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function description:
% This function analyzes Hub/perodoci cell enrichment in high power regions, including:
% 1. Visualization of Hub cell colocalization with high power cells
% 2. Calculation of Hub cell proportion at different power percentiles
% 3. Comparison of power distribution between Hub and non-Hub cells
% 4. Validation of enrichment significance through randomization tests

% Input parameters:
%   toplot - Multitaper analysis result structure with fields: mask, mask_ind, extpwr
%   params - Parameter structure with fields: hubX, hubY, hubs
%   x_coords - X coordinates of all cells
%   y_coords - Y coordinates of all cells

% Output:
%   Multiple figure windows displaying analysis results

% 功能描述：
% 此函数分析Hub细胞在高功率区域的富集情况，包括：
% 1. 可视化Hub细胞与高功率细胞的共定位
% 2. 计算Hub细胞在不同功率分位中的比例
% 3. Hub细胞与非Hub细胞的功率分布比较
% 4. 通过随机化检验验证富集显著性
%
% 输入参数：
%   toplot - 多锥谱分析结果结构体，包含mask, mask_ind, extpwr等字段
%   params - 参数结构体，包含hubX, hubY, hubs等字段
%   x_coords - 所有细胞的x坐标
%   y_coords - 所有细胞的y坐标
%
% 输出：
%   生成多个图形窗口展示分析结果
%


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
required_params_fields = {'hubX', 'hubY', 'hubs'};
for i = 1:length(required_params_fields)
    if ~isfield(params, required_params_fields{i})
        error('params structure missing required field: %s', required_params_fields{i});
    end
end

fprintf('开始分析Hub细胞功率分布...\n');
fprintf('Starting Hub cell power distribution analysis...\n');

%% 2. 提取Hub细胞坐标 / Extract Hub cell coordinates
% 从params获取Hub细胞坐标 / Get Hub cell coordinates from params
all_x = params.hubX;
all_y = params.hubY;

fprintf('检测到 %d 个Hub细胞\n', length(all_x));
fprintf('Detected %d Hub cells\n', length(all_x));

%% 3. 可视化Hub细胞分布 / Visualize Hub cell distribution
% 创建图形窗口 / Create figure window
figure('Name', 'Hub细胞分布可视化', 'Position', [100 100 800 600]);
% Hub Cells Distribution Visualization

% 绘制所有细胞（浅灰色） / Plot all cells (light gray)
scatter(x_coords, y_coords, 20, [0.7 0.7 0.7], 'filled');
set(gca, 'YDir', 'reverse');  % 反转Y轴方向匹配图像坐标 / Reverse Y-axis for image coordinates
hold on;

% 叠加绘制Hub细胞（深红色） / Overlay Hub cells (dark red)
scatter(all_x, all_y, 20, [0.6431, 0.1176, 0.1961], 'filled');
axis image;  % 保持坐标轴比例一致 / Maintain equal axis scaling

% 添加标题和标签 / Add title and labels
title('Hub Cells Distribution');
xlabel('X (μm)');
ylabel('Y (μm)');
legend({'All Cells', 'Hub Cells'}, 'Location', 'bestoutside');
grid on;

% 保存Hub细胞坐标 / Save Hub cell coordinates
hub_coords = [all_x(:), all_y(:)];

%% 4. 识别高功率细胞 / Identify high power cells
% 获取功率值 / Get power values
power_values = toplot.extpwr;

% 按功率值降序排序 / Sort power values in descending order
[sorted_power, sorted_indices_number_lable] = sort(power_values, 'descend');

% 选择前50%的细胞作为高功率细胞 / Select top 50% cells as high power cells
num_high_power = ceil(numel(sorted_power) * 0.50);

% 获取高功率细胞的索引 / Get indices of high power cells
sorted_mask_inds = toplot.mask_ind(sorted_indices_number_lable);
top_inds = sorted_mask_inds(1:num_high_power);

% 将线性索引转换为行列坐标 / Convert linear indices to row-column coordinates
[row, col] = ind2sub(size(toplot.mask), top_inds);
high_power_coords = [col(:), row(:)];

fprintf('识别出 %d 个高功率细胞（前50%%）\n', num_high_power);
fprintf('Identified %d high power cells (top 50%%)\n', num_high_power);

%% 5. 计算Hub细胞与高功率细胞的共定位 / Calculate colocalization of Hub and high power cells
% 设置精度以避免浮点误差 / Set precision to avoid floating point errors
precision = 4;
hub_coords_round = round(hub_coords, precision);
high_power_coords_round = round(high_power_coords, precision);

% 查找共同的坐标 / Find common coordinates
[common_coords, idx_hub, ~] = intersect(hub_coords_round, high_power_coords_round, 'rows');

% 计算共定位比例 / Calculate colocalization proportion
num_common = size(common_coords, 1);
num_hub_total = size(hub_coords, 1);
percentage = (num_common / num_hub_total) * 100;

fprintf('共定位统计：\n');
fprintf('Colocalization statistics:\n');
fprintf('  Hub细胞总数: %d\n', num_hub_total);
fprintf('  Total Hub cells: %d\n', num_hub_total);
fprintf('  共定位细胞数: %d\n', num_common);
fprintf('  Colocalized cells: %d\n', num_common);
fprintf('  共定位比例: %.1f%%\n', percentage);
fprintf('  Colocalization proportion: %.1f%%\n', percentage);

%% 6. 可视化共定位结果 / Visualize colocalization results
figure('Name', 'Hub细胞与高功率细胞共定位', 'Position', [100 100 900 700]);
% Hub Cells and High Power Cells Colocalization

% 绘制所有细胞 / Plot all cells
scatter(x_coords, y_coords, 20, [0.7 0.7 0.7], 'filled');
set(gca, 'YDir', 'reverse');
hold on;

% 绘制Hub细胞 / Plot Hub cells
scatter(hub_coords(:, 1), hub_coords(:, 2), 20, [0.6431, 0.1176, 0.1961], 'filled');
set(gca, 'YDir', 'reverse');

% 绘制高功率细胞 / Plot high power cells
scatter(high_power_coords(:, 1), high_power_coords(:, 2), ...
        20, [0, 0.4431, 0.7373], 'filled');
set(gca, 'YDir', 'reverse');

% 如果存在共定位细胞，用特殊颜色标记 / If common cells exist, mark with special color
if ~isempty(common_coords)
    scatter(common_coords(:, 1), common_coords(:, 2), ...
            35, [0.9290 0.6940 0.1250], 'filled', ...
            'MarkerEdgeColor', [0.9290 0.6940 0.1250]);
    set(gca, 'YDir', 'reverse');
end

% 美化图形 / Beautify figure
axis image;
title(sprintf('Hub细胞在高功率区域的富集 (%.1f%%)', percentage));
% Hub Cell Enrichment in High Power Regions
xlabel('X (μm)');
ylabel('Y (μm)');
legend({'All Cells', 'Hub Cells', 'High Power Cells', 'Colocalized Cells'}, ...
       'Location', 'bestoutside');
grid on;

%% 7. 分析不同功率分位中Hub细胞的比例 / Analyze Hub cell proportion at different power percentiles
percentiles = 0.1:0.1:1.0;  % 从10%到100% / From 10% to 100%
prop_hub = zeros(size(percentiles));  % 预分配数组 / Preallocate array

% 获取唯一坐标（四舍五入）以避免重复 / Get unique coordinates (rounded) to avoid duplicates
hub_coords_round = unique(round(hub_coords, 4), 'rows');

% 重新获取排序后的索引 / Get sorted indices again
[sorted_power, sorted_indices_number_lable] = sort(power_values, 'descend');
sorted_mask_inds = toplot.mask_ind(sorted_indices_number_lable);

% 遍历每个分位数 / Iterate through each percentile
for i = 1:length(percentiles)
    % 计算当前分位数对应的细胞数 / Calculate number of cells for current percentile
    num_high_power = ceil(numel(sorted_power) * percentiles(i));
    
    % 获取高功率细胞索引 / Get high power cell indices
    top_inds = sorted_mask_inds(1:num_high_power);
    [row, col] = ind2sub(size(toplot.mask), top_inds);
    
    % 获取坐标并去重 / Get coordinates and remove duplicates
    high_power_coords = [col(:), row(:)];
    high_power_coords_round = unique(round(high_power_coords, 4), 'rows');
    
    % 查找共定位细胞 / Find colocalized cells
    [common_coords, ~, ~] = intersect(hub_coords_round, high_power_coords_round, 'rows');
    
    % 计算比例 / Calculate proportion
    num_common = size(common_coords, 1);
    num_hub_total = size(hub_coords_round, 1);
    prop_hub(i) = num_common / num_hub_total;
end

% 绘制富集曲线 / Plot enrichment curve
figure('Name', 'Hub细胞富集曲线', 'Position', [100 100 800 600]);
% Hub Cell Enrichment Curve
plot(percentiles * 100, prop_hub, 'o-', 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]);
xlabel('Top X% High-Power Cells');
ylabel('Proportion of Hub Cells');
title('Enrichment of Hub Cells in High-Power Regions');
grid on;

%% 8. Hub细胞与非Hub细胞的功率分布比较 / Compare power distribution between Hub and non-Hub cells
% 获取所有细胞的坐标 / Get coordinates of all cells
row_all = row;
col_all = col;
coords_all = round([col_all(:), row_all(:)], 4);

% 识别Hub细胞索引 / Identify Hub cell indices
[~, hub_idx, ~] = intersect(coords_all, hub_coords, 'rows');
all_idx = 1:size(coords_all, 1);
non_hub_idx = setdiff(all_idx, hub_idx);

% 提取功率值 / Extract power values
hub_power = power_values(hub_idx);
non_hub_power = power_values(non_hub_idx);

% Kolmogorov-Smirnov检验 / Kolmogorov-Smirnov test
[~, p_ks, ks_stat] = kstest2(hub_power, non_hub_power);

fprintf('\n功率分布比较：\n');
fprintf('Power distribution comparison:\n');
fprintf('  Hub细胞数: %d\n', length(hub_power));
fprintf('  Hub cells: %d\n', length(hub_power));
fprintf('  非Hub细胞数: %d\n', length(non_hub_power));
fprintf('  Non-Hub cells: %d\n', length(non_hub_power));
fprintf('  KS检验统计量: %.3f\n', ks_stat);
fprintf('  KS test statistic: %.3f\n', ks_stat);
fprintf('  p值: %.2e\n', p_ks);
fprintf('  p-value: %.2e\n', p_ks);

% 绘制累积分布函数图 / Plot cumulative distribution function
figure('Name', '功率分布CDF比较', 'Position', [100 100 800 600]);
% Power Distribution CDF Comparison
h1 = cdfplot(hub_power);
hold on;
h2 = cdfplot(non_hub_power);

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
legend({'Hub Cells', 'Non-Hub Cells'}, 'Location', 'best');
title(sprintf('功率分布比较 (KS检验 p = %.2e)', p_ks));
% Power Distribution Comparison (KS test p = ...)
xlabel('Power');
ylabel('Cumulative Probability');
grid on;
set(gca, 'YColor', 'k');

%% 9. 随机化检验计算95%置信区间 / Randomization test to calculate 95% confidence interval
% 重新获取所有细胞的坐标和功率值 / Re-get coordinates and power values of all cells
sorted_mask_inds = toplot.mask_ind;
[row_all, col_all] = ind2sub(size(toplot.mask), sorted_mask_inds);
coords_all = round([col_all(:), row_all(:)], 4);

% 获取Hub细胞索引和功率值 / Get Hub cell indices and power values
[~, hub_idx, ~] = intersect(coords_all, hub_coords, 'rows');
all_idx = 1:size(coords_all, 1);
non_hub_idx = setdiff(all_idx, hub_idx);
hub_power = power_values(hub_idx);
n_hub = length(hub_idx);

% 设置随机化参数 / Set randomization parameters
n_iter = 1000;  % 随机化迭代次数 / Number of randomization iterations

% 重新计算富集曲线 / Recalculate enrichment curve
percentiles = 0.1:0.1:1.0;
prop_hub = zeros(length(percentiles), 1);
hub_coords_round = unique(hub_coords, 'rows');

[sorted_power, sorted_indices_number_lable] = sort(power_values, 'descend');
sorted_mask_inds = toplot.mask_ind(sorted_indices_number_lable);

for i = 1:length(percentiles)
    num_top = ceil(numel(sorted_power) * percentiles(i));
    top_inds = sorted_mask_inds(1:num_top);
    [row, col] = ind2sub(size(toplot.mask), top_inds);
    high_power_coords = unique(round([col(:), row(:)], 4), 'rows');
    [common_coords, ~, ~] = intersect(hub_coords_round, high_power_coords, 'rows');
    prop_hub(i) = size(common_coords, 1) / size(hub_coords_round, 1);
end

fprintf('\n进行随机化检验（%d次迭代）...\n', n_iter);
fprintf('Performing randomization test (%d iterations)...\n', n_iter);

% 随机化检验 / Randomization test
prop_random_all = zeros(length(percentiles), n_iter);

for iter = 1:n_iter
    % 随机选择与Hub细胞数量相同的细胞 / Randomly select same number of cells as Hub cells
    rand_idx = randsample(size(coords_all, 1), n_hub);
    rand_coords = coords_all(rand_idx, :);
    
    % 计算每个分位数的比例 / Calculate proportion for each percentile
    for j = 1:length(percentiles)
        num_top = ceil(numel(sorted_power) * percentiles(j));
        top_inds = sorted_mask_inds(1:num_top);
        [row, col] = ind2sub(size(toplot.mask), top_inds);
        high_power_coords = unique(round([col(:), row(:)], 4), 'rows');
        [common_coords, ~, ~] = intersect(round(rand_coords, 4), high_power_coords, 'rows');
        prop_random_all(j, iter) = size(common_coords, 1) / n_hub;
    end
    
    % 显示进度 / Show progress
    if mod(iter, 100) == 0
        fprintf('  已完成 %d/%d 次迭代\n', iter, n_iter);
        fprintf('  Completed %d/%d iterations\n', iter, n_iter);
    end
end

% 计算随机化的统计量 / Calculate randomization statistics
mean_prop_rand = mean(prop_random_all, 2);
std_prop_rand = std(prop_random_all, 0, 2);
ci_low = mean_prop_rand - 1.96 * std_prop_rand;
ci_high = mean_prop_rand + 1.96 * std_prop_rand;

% 绘制富集曲线与置信区间 / Plot enrichment curve with confidence interval
figure('Name', '随机化检验富集曲线', 'Position', [100 100 900 700]);
% Randomization Test Enrichment Curve
hold on;

% 绘制95%置信区间 / Plot 95% confidence interval
fill([percentiles * 100, fliplr(percentiles * 100)], ...
     [ci_low', fliplr(ci_high')], [0.9 0.9 0.9], 'EdgeColor', 'none', ...
     'DisplayName', '95% CI (Random)');

% 绘制随机化平均值 / Plot randomization mean
plot(percentiles * 100, mean_prop_rand, 'k--', 'LineWidth', 1.5, ...
     'DisplayName', 'Mean Random');

% 绘制实际Hub细胞富集曲线 / Plot actual Hub cell enrichment curve
plot(percentiles * 100, prop_hub, 'r-', 'LineWidth', 2, ...
     'DisplayName', 'Real Hub');

% 美化图形 / Beautify figure
xlabel('Top X% High-Power Cells');
ylabel('Proportion of Hub Cells');
title('富集曲线：Hub细胞 vs 随机基线');
% Enrichment Curve: Hub Cells vs Random Baseline
legend('Location', 'best');
grid on;
box on;

% 添加显著性标注 / Add significance annotation
for i = 1:length(percentiles)
    if prop_hub(i) > ci_high(i)
        text(percentiles(i) * 100, prop_hub(i) + 0.02, '*', ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    end
end

fprintf('\n分析完成！\n');
fprintf('Analysis completed!\n');
fprintf('生成了 %d 个图形窗口\n', 5);
fprintf('Generated %d figure windows\n', 5);

end  % 函数结束 / End of function