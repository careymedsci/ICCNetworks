
 
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


function [params, analysis_summary] = plotFrequencyHistogram(params, analysis_summary)

% 输入参数结构体 Input struct:
%   params.secondsPerFrame  - 每帧持续时间（秒） Time per frame (s)
%   params.numFrames        - 总帧数 Total number of frames
%   params.peakLocations    - [numFrames x numCells]，细胞峰值位置的二值矩阵 Binary matrix of peak positions
%   params.pk_num_define    - 峰值数量的最小阈值 Minimum required number of peaks
%   params.periodicity      - 每个细胞是否具有周期性，1=周期性，0=非周期性 1=periodic, 0=non-periodic
%   params.tracesleft       - 有效细胞数量 Number of valid cells

% === Step 1: 初始化与预处理 Initialization ===
Fs = 1 / params.secondsPerFrame;  % 采样率 Sampling frequency (Hz)
numFrames = params.numFrames;
tracesleft = params.tracesleft;

peakCounts = squeeze(sum(params.peakLocations, 1));  % 每个细胞的峰数量 Peak count per cell

% 选取峰值数量大于等于阈值的细胞 Valid cells with enough peaks
minPeaks = params.pk_num_define;
validCells = find(peakCounts >= minPeaks);

freqAll = [];       % 所有合格细胞的频率 Frequencies of all valid cells
freqPeriodic = [];  % 周期细胞的频率 Frequencies of periodic cells

% === Step 2: 计算每个细胞的平均频率 Compute frequency for each valid cell ===
for i = validCells
    peakTimes = find(params.peakLocations(:, i)) / Fs;  % 峰的时间点 Peak time points (s)
    if numel(peakTimes) >= minPeaks
        peakIntervals = diff(peakTimes);       % 峰间间隔 Interval between peaks
        meanFreq = 1 / mean(peakIntervals);    % 平均频率 Mean frequency (Hz)
        freqAll = [freqAll, meanFreq];         % 添加至总频率 Add to all frequencies
        if params.periodicity(i) == 1
            freqPeriodic = [freqPeriodic, meanFreq];  % 添加至周期频率 Add to periodic frequencies
        end
    end
end

% === Step 3: 获取最大/最小频率及其细胞索引 Find extreme frequency values ===
[minFreq, minIdx] = min(freqPeriodic);
[maxFreq, maxIdx] = max(freqPeriodic);
[minFreq1, minIdx1] = min(freqAll);
[maxFreq1, maxIdx1] = max(freqAll);

% === Step 4: 汇总结果 Summary results ===
analysis_summary(1,5) = "periodic_min_frq";     % 周期细胞最小频率
analysis_summary(1,6) = minFreq;
analysis_summary(2,5) = "cell label of 'periodic_min_frq'";
analysis_summary(2,6) = minIdx;

analysis_summary(3,5) = "periodic_max_frq";     % 周期细胞最大频率
analysis_summary(3,6) = maxFreq;
analysis_summary(4,5) = "cell label of 'periodic_max_frq'";
analysis_summary(4,6) = maxIdx;

analysis_summary(5,5) = "general_min_frq";      % 所有细胞最小频率
analysis_summary(5,6) = minFreq1;
analysis_summary(6,5) = "cell label of 'general_min_frq'";
analysis_summary(6,6) = minIdx1;

analysis_summary(7,5) = "general_max_frq";      % 所有细胞最大频率
analysis_summary(7,6) = maxFreq1;
analysis_summary(8,5) = "cell label of 'general_max_frq'";
analysis_summary(8,6) = maxIdx1;

% === Step 5: 构建频率分布直方图 Frequency histogram ===
[countsAll, edges] = histcounts(freqAll);                  % 所有细胞频率直方图 Histogram for all
countsPeriodic = histcounts(freqPeriodic, edges);          % 周期细胞频率 Histogram for periodic

percentAll = 100 * countsAll / tracesleft;                 % 所有细胞百分比百分比转换
percentPeriodic = 100 * countsPeriodic / tracesleft;       % 周期细胞百分比

% === Step 6: 绘制直方图 Histogram plotting ===
figure; hold on;
centers = edges(1:end-1) + diff(edges)/2;  % 计算每个bin的中心 Bin centers

% 绘制条形图 Bar plots
bar(centers, percentAll, 'FaceAlpha', 0.5, 'BarWidth', 1);      % 所有细胞 All
bar(centers, percentPeriodic, 'FaceAlpha', 0.5, 'BarWidth', 1); % 周期细胞 Periodic

xlabel('Frequency (Hz)');               
ylabel('Percentage of Total Cells (%)');
title('Frequency Distribution of Cells with ≥4 Peaks'); 
legend('All (≥4 peaks)', 'Periodic Cells');
box on;
grid on;

% 设置y轴范围 Set Y-axis limits
ymax = max([percentAll, percentPeriodic]) + 2;
ylim([0 ymax]);

% 添加数值标签 Add value labels above bars
for i = 1:length(centers)
    if percentAll(i) > 0
        text(centers(i), percentAll(i) + 0.5, ...
             sprintf('%.1f%%', percentAll(i)), ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 8);
    end
    if percentPeriodic(i) > 0
        text(centers(i), percentPeriodic(i) + 0.5, ...
             sprintf('%.1f%%', percentPeriodic(i)), ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 8, 'Color', 'red');
    end
end

% 设置x轴刻度与标签 X-axis ticks and labels
xticks(edges);
xticklabels(arrayfun(@(x) sprintf('%.3f', x), edges, 'UniformOutput', false));

% === Step 7: 保存图形与数据 Save figure and results ===
savefig('frequency_vs_all_percentage.fig');
save('_all_result', 'params', 'analysis_summary');
end
