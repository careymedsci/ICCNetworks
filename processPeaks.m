

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
% 2024.10.09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [params, analysis_summary] = processPeaks(params)
% PROCESSPEAKS Detects and quantifies calcium transient peaks from ROI traces
% 峰值检测与统计分析函数：用于提取ROI钙信号瞬变事件并计算统计指标
%
% Description:
% This function preprocesses calcium fluorescence traces, detects transient
% peaks using prominence and width thresholds, and computes event statistics
% including peak rate, amplitude, and duration.
%
% 本函数对钙成像时间序列进行预处理（去均值和平滑），使用峰值检测算法
% 提取Ca2+瞬变事件，并计算峰频率、峰值幅度和宽度等统计特征。
%
% Input:
%   params - structure containing analysis parameters and data
%            输入结构体，包含数据和分析参数
%
% Required fields in params:
%   allTraces         - fluorescence traces (time × ROI)
%                       钙信号矩阵（时间 × ROI）
%   numFrames         - number of frames
%                       时间帧数量
%   numRois           - number of ROIs
%                       ROI数量
%   secondsPerFrame   - acquisition interval (s)
%                       采样时间间隔（秒）
%   minamplitude      - minimum peak prominence threshold
%                       峰值最小显著性阈值
%   maxwidth          - maximum allowed peak width (seconds)
%                       最大允许峰宽（秒）
%   output            - output filename prefix
%                       输出文件名前缀
%   plot_cells_traces - index for file naming
%                       输出文件编号标识
%   mode              - processing mode ID
%                       模式编号
%
% Output:
%   params            - updated parameter structure
%                       更新后的参数结构体
%   analysis_summary  - summary table of peak statistics
%                       峰值统计摘要表

%% Initialization 初始化变量

% Initialize full-size peak location matrix
% 初始化完整尺寸峰值矩阵
peakLocationsfullsize = zeros(params.numFrames, params.numRois);

% Convert sampling interval to frames per minute
% 将采样时间转换为“每分钟帧数”
framesPerMinute = 60 / params.secondsPerFrame;

%% Preprocessing 数据预处理

% Remove mean baseline for each ROI trace
% 对每个ROI时间序列进行去均值处理（基线校正）
params.allTraces = params.allTraces - ...
    repmat(mean(params.allTraces,2), 1, params.numRois);

%% Temporal smoothing 时间平滑

% Apply Gaussian smoothing to reduce noise
% 使用高斯滤波平滑时间信号以抑制高频噪声
% sigma = 7 seconds

for g = 1:params.numRois
    params.allTraces(:,g) = imgaussfilt( ...
        params.allTraces(:,g), ...
        7/params.secondsPerFrame);
end

%% Peak storage initialization 峰值存储初始化

% Structure for storing individual peak info
% 峰值结构体（位置和幅度）
peaksInfo = struct('peakLocations', [], 'peakHeight', []);

% Matrices for amplitude and width storage
% 峰值幅度和宽度矩阵
amplitudes = zeros(params.numFrames, params.numRois);
widths = zeros(params.numFrames, params.numRois);

%% Peak detection 峰值检测

% Loop through each ROI trace
% 遍历所有ROI时间序列
for trace = 1:params.numRois
    
    % Detect peaks using prominence and width constraints
    % 使用显著性和宽度约束进行峰值检测
    [height, locs, width, amp] = findpeaks( ...
        params.allTraces(:,trace), ...
        'MinPeakProminence', params.minamplitude, ...
        'MaxPeakWidth', params.maxwidth/params.secondsPerFrame);
    
    % Store peak locations
    % 存储峰位置
    peaksInfo(trace).peakLocations = locs;
    
    % Store raw peak intensities
    % 存储真实信号强度峰值
    peaksInfo(trace).peakHeight = params.allTraces(locs,trace);
    
    % Store peak prominence (relative amplitude)
    % 存储相对振幅（峰-谷差）
    amplitudes(locs,trace) = amp;
    
    % Store peak width (in frames)
    % 存储峰宽（单位：帧）
    widths(locs,trace) = width;
end

%% Peak location matrix generation 峰位置矩阵生成

% Binary peak occurrence matrix
% 峰值出现的二值矩阵
peakLocations = amplitudes > 0;

% Resize to full frame length
% 对齐完整时间长度
s = size(peakLocations,1);
if s < params.numFrames
    peakLocationsfullsize(1:s,:) = peakLocations;
end

%% Save peak location matrices 保存峰位置文件

outputValue = sprintf('%s_peakLocations%s.txt', ...
    params.output, num2str(params.plot_cells_traces));
writematrix(peakLocations, outputValue);

outputValue = sprintf('%s_peakLocations_full%s.txt', ...
    params.output, num2str(params.plot_cells_traces));
writematrix(peakLocationsfullsize, outputValue);

%% Peak statistics calculation 峰值统计计算

% Total number of detected peaks
% 峰值总数量
totalPeaks = sum(sum(amplitudes > 0));

% Convert to peaks per 10 minutes
% 计算10分钟内峰事件数量
peaksPerTenMinutes = totalPeaks / params.numFrames ...
                     * framesPerMinute * 10;

% Normalize by ROI count (per 100 cells)
% 归一化为每100个细胞的峰事件数
normalizedPeaks = peaksPerTenMinutes / params.numRois * 100;

%% Peak amplitude statistics 峰幅度统计

% Extract all peak amplitudes
% 提取所有峰值振幅
peakAmplitudes = reshape(amplitudes(peakLocations), 1, []).';

meanAmplitude = mean(peakAmplitudes);

outputValue = sprintf('%s_peakamplitudes%s.txt', ...
    params.output, num2str(params.plot_cells_traces));
writematrix(peakAmplitudes, outputValue);

%% Peak width statistics 峰宽统计

% Extract all peak widths
% 提取峰宽
peakWidths = reshape(widths(peakLocations), 1, []);

meanWidth = mean(peakWidths, 'all');

% Convert width from frames to seconds
% 峰宽由帧转换为秒
peakWidths = peakWidths.' * params.secondsPerFrame;

outputValue = sprintf('%s_peakwidths%s.txt', ...
    params.output, num2str(params.plot_cells_traces));
writematrix(peakWidths, outputValue);

%% Summary table generation 统计摘要表

analysis_summary(1,1) = "numFrames";
analysis_summary(1,2) = params.numFrames;

analysis_summary(2,1) = "numRois";
analysis_summary(2,2) = params.numRois;

analysis_summary(3,1) = "peaks / 100 cells in 600s";
analysis_summary(3,2) = normalizedPeaks;

analysis_summary(4,1) = "Mean peak amplitude";
analysis_summary(4,2) = meanAmplitude;

analysis_summary(5,1) = "Mean peak width (seconds)";
analysis_summary(5,2) = meanWidth;

%% Update params structure 更新参数结构体

params.peakLocations = peakLocations;
params.peakLocationsfullsize = peakLocationsfullsize;
params.peaksInfo = peaksInfo;
params.amplitudes = amplitudes;
params.widths = widths;

%% Save all variables 保存分析结果

filename = sprintf('all_ver_%d.mat', params.mode);
save(filename);

end










































% backup 

% function [params, analysis_summary] = processPeaks(params)
%     % 处理和分析峰值数据
%     % 输入: params - 包含所有参数和数据的结构体
%     % 输出: 
%     %   params - 更新后的参数结构体
%     %   analysis_summary - 分析结果摘要矩阵
% 
%     % 初始化矩阵
%     peakLocationsfullsize = zeros(params.numFrames, params.numRois);
%     framesPerMinute = 60 / params.secondsPerFrame;
% 
%     % 数据预处理
%     params.allTraces = params.allTraces - repmat(mean(params.allTraces,2), 1, params.numRois);
% 
%     % 平滑处理
%     for g = 1:params.numRois
%         params.allTraces(:,g) = imgaussfilt(params.allTraces(:,g), 7/params.secondsPerFrame); 
%     end
% 
%     % 初始化峰值存储
%     peaksInfo = struct('peakLocations', [], 'peakHeight', []);
%     amplitudes = zeros(params.numFrames, params.numRois);
%     widths = zeros(params.numFrames, params.numRois);
% 
%     % 寻找峰值
%     for trace = 1:params.numRois
%         [height, locs, width, amp] = findpeaks(params.allTraces(:,trace), ...
%             'MinPeakProminence', params.minamplitude, ...
%             'MaxPeakWidth', params.maxwidth/params.secondsPerFrame);
% 
%         % 存储峰值信息
%         peaksInfo(trace).peakLocations = locs;
%         peaksInfo(trace).peakHeight = params.allTraces(locs,trace);
%         % 这个是真实的intensities值的高度集合
%         amplitudes(locs,trace) = amp;
%         % 这个是相对振幅（最高点与两边最低点的差值）集合
%         widths(locs,trace) = width;
%         % 这个是真实的波峰宽度集合
%     end
% 
%     % 生成峰值位置矩阵
%     peakLocations = amplitudes > 0;
%     s = size(peakLocations,1);
%     if s < params.numFrames
%         peakLocationsfullsize(1:s,:) = peakLocations;
%     end
% 
% 
%     outputValue = sprintf('%s_peakLocations%s.txt', params.output, num2str(params.plot_cells_traces));
%     writematrix(peakLocations, outputValue);
%     outputValue = sprintf('%s_peakLocations_full%s.txt', params.output, num2str(params.plot_cells_traces));
%     writematrix(peakLocationsfullsize, outputValue);
% 
% 
%     totalPeaks = sum(sum(amplitudes > 0));
%     peaksPerTenMinutes = totalPeaks / params.numFrames * framesPerMinute * 10;
%     normalizedPeaks = peaksPerTenMinutes / params.numRois * 100;
% 
% 
%     peakAmplitudes = reshape(amplitudes(peakLocations), 1, []).';
%     meanAmplitude = mean(peakAmplitudes);
%     outputValue = sprintf('%s_peakamplitudes%s.txt', params.output, num2str(params.plot_cells_traces));
%     writematrix(peakAmplitudes, outputValue);
% 
% 
%     peakWidths = reshape(widths(peakLocations), 1, []);
%     meanWidth = mean(peakWidths, 'all');
%     peakWidths = peakWidths.' * params.secondsPerFrame;
%     outputValue = sprintf('%s_peakwidths%s.txt', params.output, num2str(params.plot_cells_traces));
%     writematrix(peakWidths, outputValue);
% 
%     analysis_summary(1,1) = "numFrames";
%     analysis_summary(1,2) = params.numFrames;
%     analysis_summary(2,1) = "numRois";
%     analysis_summary(2,2) = params.numRois;
%     analysis_summary(3,1) = "peaks / 100 cells in 600s";
%     analysis_summary(3,2) = normalizedPeaks;
%     analysis_summary(4,1) = "Mean peak amplitude";
%     analysis_summary(4,2) = meanAmplitude;
%     analysis_summary(5,1) = "Mean peak width (seconds)";
%     analysis_summary(5,2) = meanWidth;
% 
%     % 更新params结构体中的必要字段
%     params.peakLocations = peakLocations;
%     params.peakLocationsfullsize = peakLocationsfullsize;
%     params.peaksInfo = peaksInfo;
%     params.amplitudes = amplitudes;
%     params.widths = widths;
% 
%     % 保存所有变量
%     filename = sprintf('all_ver_%d.mat', params.mode);
%     save(filename);
% end