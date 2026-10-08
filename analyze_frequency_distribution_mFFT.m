

 
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


function analyze_frequency_distribution_mFFT()
% ANALYZE_FREQUENCY_DISTRIBUTION 分析钙信号频率分布及细胞类型特征
% Analyze calcium signal frequency distribution and cell type characteristics
%
% Xiao Liu 
% Spuervisor: Thomas Broggini
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany
% Xiangyang No.1 people's Hospital, China
% Modified from  UCSD
% Xiao.Liu@stud.uni-frankfurt.de
% 2025-04-19

% Function description:
% This function loads multitaper spectral analysis results, calculates the
% maximum power frequency distribution, and specifically analyzes frequency
% characteristics of Hub and Periodic cells, generating segmented frequency
% distribution plots.

%% 用户输入参数设置 / User parameter input
% 设置频率显示范围 / Set frequency display range
user_decision = true;

while user_decision
    % 输入最高频率阈值 / Input maximum frequency threshold
    prompt_max = {'For the Max Freq distribution, what phase parameter you want to put in<the highest frequency you choose,like 0.2?>'};
    answer_max = inputdlg(prompt_max);
    value_max_freq = str2double(answer_max{1});
    disp(['Max frequency value you input: ', num2str(value_max_freq)]);
    
    % 输入最低频率阈值 / Input minimum frequency threshold
    prompt_min = {'For the Max Freq distribution, what phase parameter you want to put in<the lowest frequency you choose,like 0.04?>'};
    answer_min = inputdlg(prompt_min);
    value_min_freq = str2double(answer_min{1});
    disp(['Min frequency value you input: ', num2str(value_min_freq)]);

    %% 加载分析结果 / Load analysis results
    % 加载之前保存的多锥谱分析结果 / Load previously saved multitaper analysis results
    load('dfof4_coordintes.txt.mat');  % 确保此文件包含toplot结构 / Ensure this file contains toplot structure
    
    %% 初始化参数和变量 / Initialize parameters and variables
    % 加载ROI掩膜 / Load ROI mask
    im_mask = toplot.mask;  
    
    % 设置模式参数（通常为1） / Set mode parameter (typically 1)
    tmp_mode = 1;  
    
    % 初始化映射数组 / Initialize mapping arrays
    map = zeros(toplot.mask_size);  
    im_size = size(toplot.mask);  
    im_phase = zeros(size(im_mask));  
    im_mag = zeros(size(im_mask));  
    
    % 设置默认显示参数 / Set default display parameters
    tmp_max_frm = 10;  
    tmp_min_frm = -10;  
    tmp_max_phase = pi/2;  
    tmp_min_phase = -pi/2;  
    tmp_max_mag = 0.003;  
    tmp_min_mag = 0;  
    tmp_max_pwr = log10(130);  
    tmp_min_pwr = log10(0.1);  
    tmp_max_f0 = 0.13;  
    tmp_min_f0 = 0.0;  
    
    % 存储参数到结构体 / Store parameters to structure
    toplot.prams.max_frm = tmp_max_frm;
    toplot.prams.min_frm = -tmp_min_frm;
    toplot.prams.max_phase = tmp_max_phase;
    toplot.prams.min_phase = tmp_min_phase;
    toplot.prams.max_mag = tmp_max_mag;
    toplot.prams.min_mag = tmp_min_mag;
    toplot.prams.max_pwr = tmp_max_pwr;
    toplot.prams.min_pwr = tmp_min_pwr;
    toplot.prams.binsize = linspace(-pi, pi, 360);  
    ext = 1;  
    
    % 频率插值初始化 / Frequency interpolation initialization
    % 注意：此处假设findx, tapers, nfft, Fs已存在于工作空间
    % Note: Assumes findx, tapers, nfft, Fs already exist in workspace
    toplot.findx = findx;
    tapers_FT = fft(tapers, nfft) / Fs;  
    tapers_FT = tapers_FT(1, :);  
    t_norm = sum(tapers_FT.^2);  
    
    % 清理临时变量 / Clear temporary variables
    clear Un Sn i J
    
    %% 频率峰值处理 / Frequency peak processing
    % 构建扩展的频率数组 / Build extended frequency array
    if size(toplot.f_peak, 1) > 3
        % 处理4个以上峰值的情况 / Handle case with more than 3 peaks
        f0 = horzcat(toplot.f_peak(1), ...
            linspace(toplot.f_peak(2), toplot.f_peak(2) * fix(1 / toplot.f_peak(2)), fix(1 / toplot.f_peak(2))), ...
            linspace(toplot.f_peak(3), toplot.f_peak(3) * fix(1 / toplot.f_peak(3)), fix(1 / toplot.f_peak(3))), ...
            linspace(toplot.f_peak(4), toplot.f_peak(4) * fix(1 / toplot.f_peak(4)), fix(1 / toplot.f_peak(4))));
    elseif size(toplot.f_peak, 1) > 2
        % 处理3个峰值的情况 / Handle case with 3 peaks
        f0 = horzcat(toplot.f_peak(1), ...
            linspace(toplot.f_peak(2), toplot.f_peak(2) * fix(1 / toplot.f_peak(2)), fix(1 / toplot.f_peak(2))), ...
            linspace(toplot.f_peak(3), toplot.f_peak(3) * fix(1 / toplot.f_peak(3)), fix(1 / toplot.f_peak(3))));
    else
        % 处理少于3个峰值的情况 / Handle case with less than 3 peaks
        f0 = horzcat(toplot.f_peak(1), toplot.f_peak(2), toplot.f_peak(3));
    end
    
    %% FFT插值和功率计算 / FFT interpolation and power calculation
    % 为每个模式插值FFT / Interpolate FFT for each mode
    for mode = 1:size(scores, 2)
        for k = 1:ntapers(2)
            interp_FFT(:, k, mode) = interp1(f, taperedFFT(:, k, mode), f0);
        end
    end
    
    % 初始化数组存储相干性和功率结果 / Initialize array for coherence and power results
    A = zeros(length(f0), size(interp_FFT, 3), 'single');
    A = complex(A);
    
    % 累加所有taper的贡献 / Sum contributions from all tapers
    for k = 1:ntapers(2)
        A = A + tapers_FT(k) .* squeeze(interp_FFT(:, k, :));
    end
    A = A ./ t_norm;  % 归一化结果 / Normalize result
    
    % 计算每个模式的平均功率 / Calculate mean power for each mode
    mu = scores * A';  
    clear A k
    
    %% 频率窗口分析 / Frequency window analysis
    % 设置频率窗口 / Set frequency window
    toplot.wind(2) = round(length(f) / 2);
    A = zeros(length(f(toplot.wind(1):toplot.wind(2))), size(taperedFFT, 3), 'single');
    A = complex(A);
    
    for k = 1:ntapers(2)
        A = A + tapers_FT(k) .* squeeze(taperedFFT(toplot.wind(1):toplot.wind(2), k, :));
    end
    A = A ./ t_norm;
    toplot.ampwind = scores * A';  % 存储窗口频率的振幅 / Store amplitude for windowed frequencies
    clear A k
    
    %% 频率分解计算 / Frequency decomposition calculation
    tic
    Fde = zeros(size(mu), 'single');
    for k = 1:ntapers(2)
        % 减去平均值 / Subtract mean
        z = scores * squeeze(interp_FFT(:, k, :))' - mu * tapers_FT(k);  
        Fde = Fde + conj(z) .* z;  % 计算分解 / Calculate decomposition
        tmp_tic = tic;
        fprintf('Finish calculating %d/%d taper. Elapsed time is %f seconds\n', k, ntapers(2), toc(tmp_tic));
    end
    toc
    
    %% 最终频率计算 / Final frequency calculation
    F = (size(taperedFFT, 2) - 1) * (conj(mu) .* mu) * t_norm ./ Fde;
    clear Fde nsvd
    
    % 存储分析结果 / Store analysis results
    toplot.amps = mu;  % 振幅谱 / Amplitude spectrum
    toplot.fval = F;  % 频率值 / Frequency values
    
    % 移除chronux工具箱路径 / Remove chronux toolbox path
    rmpath(genpath('C:\chronux_2_12'));  
    
    toplot.f0 = f0;  % 频率峰值 / Frequency peaks
    dfact = toplot.Delta_f * 2;  % 缩放因子 / Scaling factor
    toplot.maxf0 = zeros(1, length(toplot.amps));  % 初始化最大频率 / Initialize max frequencies
    
    % 找到每个模式的最大频率 / Find maximum frequency for each mode
    for n = 1:length(toplot.amps)
        [~, idx] = max(abs(toplot.amps(n, :)).^2);
        toplot.maxf0(n) = toplot.f0(idx);  
    end
    
    %% 功率结果初始化 / Power results initialization
    toplot.maxf = zeros(1, length(toplot.ampwind));
    toplot.extpwr = zeros(1, length(toplot.ampwind));
    
    % 找到最大功率及对应频率 / Find maximum power and corresponding frequencies
    for n = 1:length(toplot.ampwind)
        [val, idx] = max(abs(toplot.ampwind(n, :)).^2);
        toplot.maxf(n) = toplot.f(toplot.wind(1) + idx);  
        toplot.extpwr(n) = val / dfact;  
    end
    
    %% 特定频率峰值功率存储 / Store power at specific frequency peaks
    % 注意：峰值的命名是项目特定的，可以忽略 / Note: Peak names are project-specific, can be ignored
    toplot.resonancepwr = abs(toplot.amps(:,1)).^2 ./ dfact;  % 第一峰值功率 / First peak power
    toplot.puffpwr = abs(toplot.amps(:,2)).^2 ./ dfact;  % 第二峰值功率 / Second peak power
    toplot.puffharmpwr = abs(toplot.amps(:,3)).^2 ./ dfact;  % 第三峰值功率 / Third peak power
    
    % 基于骨骼标签调整结果 / Adjust results based on skeletal label
    toplot.maxf = toplot.maxf(:, toplot.skel_label);
    toplot.maxf0 = toplot.maxf0(:, toplot.skel_label);
    toplot.extpwr = toplot.extpwr(:, toplot.skel_label);
    
    toplot.pwrext = zeros(size(toplot.skel_label));
    toplot.pwrext(1,:) = toplot.resonancepwr(toplot.skel_label);
    toplot.pwrext(2,:) = toplot.puffpwr(toplot.skel_label);
    toplot.pwrext(3,:) = toplot.puffharmpwr(toplot.skel_label);
    
    %% 处理额外峰值频率 / Handle additional peak frequencies
    % 5个峰值的情况 / Case with 5 peaks
    if size(toplot.f_peak, 2) == 5
        audiopeak = 1 + length(horzcat(toplot.f_peak(1), ...
            linspace(toplot.f_peak(2), toplot.f_peak(2) * fix(1 / toplot.f_peak(2)), fix(1 / toplot.f_peak(2))), ...
            linspace(toplot.f_peak(3), toplot.f_peak(3) * fix(1 / toplot.f_peak(3)), fix(1 / toplot.f_peak(3)))));
        
        toplot.audpwr = abs(toplot.amps(:, audiopeak)).^2 ./ dfact;
        
        vispeak = 2 + length(linspace(toplot.f_peak(2), toplot.f_peak(2) * fix(1 / toplot.f_peak(2)), fix(1 / toplot.f_peak(2))));
        toplot.vispwr = abs(toplot.amps(:, vispeak)).^2 ./ dfact;
        
        toplot.pwrext(3,:) = toplot.vispwr(toplot.skel_label);
        toplot.pwrext(4,:) = toplot.audpwr(toplot.skel_label);
        toplot.pwrext(5,:) = toplot.puffharmpwr(toplot.skel_label);
    end
    
    % 4个峰值的情况 / Case with 4 peaks
    if size(toplot.f_peak, 2) == 4
        vispeak = 2 + length(linspace(toplot.f_peak(2), toplot.f_peak(2) * fix(1 / toplot.f_peak(2)), fix(1 / toplot.f_peak(2))));
        toplot.vispwr = abs(toplot.amps(:, vispeak)).^2 ./ dfact;
        
        toplot.pwrext(3,:) = toplot.vispwr(toplot.skel_label);
        toplot.pwrext(4,:) = toplot.puffharmpwr(toplot.skel_label);
    end
    
    % 最终频率峰值结果 / Final frequency peak results
    toplot.f_peak = toplot.f_peak';
    toplot.maxf = toplot.maxf(:, toplot.skel_label);
    
    %% 设置频率分析参数 / Set frequency analysis parameters
    tmp_max_f0 = value_max_freq;  
    tmp_min_f0 = value_min_freq; 
    frequency_values = toplot.maxf; 
    all_x = x_coords;  % 假设x_coords已定义 / Assume x_coords is defined
    all_y = y_coords;  % 假设y_coords已定义 / Assume y_coords is defined
    totalnumber = size(y_coords, 2);
    
    %% 创建频率分箱 / Create frequency bins
    num_bins = 20; 
    bin_edges = linspace(tmp_min_f0, tmp_max_f0, num_bins); 
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    % 获取掩膜坐标 / Get mask coordinates
    [mask_y, mask_x] = ind2sub(size(toplot.mask), toplot.mask_ind);
    
    %% 过滤数值伪影 / Filter numerical artifacts
    % 移除小于0.005Hz的频率值（预处理后可能残留的数值误差）
    % Remove frequencies < 0.005 Hz (possible numerical artifacts after preprocessing)
    frequency_values = frequency_values(frequency_values >= 0.005);
    
    %% ----------------------- Hub细胞分析 -----------------------
    %% ----------------------- Hub Cell Analysis -----------------------
    hub_indices = find(params.hubs == 1);
    hub_freq_values = zeros(length(hub_indices), 1);
    
    % 找到Hub细胞对应的频率值 / Find frequency values for Hub cells
    for i = 1:length(hub_indices)
        x = all_x(hub_indices(i));
        y = all_y(hub_indices(i));
        dist = sqrt((mask_x - x).^2 + (mask_y - y).^2);
        [~, min_idx] = min(dist);
        hub_freq_values(i) = toplot.maxf(min_idx);
    end
    hub_freq_values = hub_freq_values(hub_freq_values >= 0.005);
    
    % 计算频率分布 / Calculate frequency distribution
    frequency_counts = histcounts(frequency_values, bin_edges) / totalnumber * 100;
    hub_counts = histcounts(hub_freq_values, bin_edges) / totalnumber * 100;
    
    %% 创建Hub细胞频率分布图 / Create Hub cell frequency distribution plot
    figure;
    t = tiledlayout(2, 1, 'TileSpacing', 'none', 'Padding', 'compact');
    
    % 上子图：显示>15%的部分 / Upper subplot: Show >15% portion
    ax1 = nexttile(t);
    bar(ax1, bin_centers, frequency_counts, 'FaceColor', [0.6 0.6 0.6]); 
    hold(ax1, 'on');
    bar(ax1, bin_centers, hub_counts, 'FaceColor', 'r');
    ylim(ax1, [15 max(frequency_counts) + 5]);
    xlim(ax1, [0.005 max(bin_centers)]);
    xticks(ax1, []);
    ylabel(ax1, 'Cell Percentage [%]');
    title(ax1, 'Frequency Distribution with Hub Cells');
    legend(ax1, {'All Cells', 'Hub Cells'}, 'Location', 'best');
    % 添加参考线 / Add reference line
    line(ax1, [0.005 max(bin_centers)], [15 15], 'Color', 'k', 'LineStyle', '--');
    
    % 下子图：显示≤5%的部分 / Lower subplot: Show ≤5% portion
    ax2 = nexttile(t);
    bar(ax2, bin_centers, frequency_counts, 'FaceColor', [0.6 0.6 0.6]); 
    hold(ax2, 'on');
    bar(ax2, bin_centers, hub_counts, 'FaceColor', 'r');
    ylim(ax2, [0 5]);
    xlim(ax2, [0.005 max(bin_centers)]);
    xticks(ax2, 0.005:0.004:max(bin_centers));
    xlabel(ax2, 'Frequency [Hz]');
    ylabel(ax2, 'Cell Percentage [%]');
    line(ax2, [0.005 max(bin_centers)], [5 5], 'Color', 'k', 'LineStyle', '--');
    
    % 添加文本标签 / Add text labels
    for i = 1:length(frequency_counts)
        if frequency_counts(i) <= 5
            text(ax2, bin_centers(i), frequency_counts(i) + 0.3, ...
                sprintf('%.1f', frequency_counts(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
        end
        if hub_counts(i) <= 5 && hub_counts(i) > 0
            text(ax2, bin_centers(i), hub_counts(i) + 0.3, ...
                sprintf('%.1f', hub_counts(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'r');
        end
    end
    
    % 保存图像 / Save figure
    saveas(gcf, 'Freq_distri_with_hub_broken_percent.tif');
    
    %% ----------------------- Periodic细胞分析 -----------------------
    %% ----------------------- Periodic Cell Analysis -----------------------
    periodic_indices = find(params.periodicity == 1);
    periodic_freq_values = zeros(length(periodic_indices), 1);
    
    % 找到Periodic细胞对应的频率值 / Find frequency values for Periodic cells
    for i = 1:length(periodic_indices)
        x = all_x(periodic_indices(i));
        y = all_y(periodic_indices(i));
        dist = sqrt((mask_x - x).^2 + (mask_y - y).^2);
        [~, min_idx] = min(dist);
        periodic_freq_values(i) = toplot.maxf(min_idx);
    end
    periodic_freq_values = periodic_freq_values(periodic_freq_values >= 0.005);
    
    % 计算频率分布 / Calculate frequency distribution
    frequency_counts = histcounts(frequency_values, bin_edges) / totalnumber * 100;
    periodic_counts = histcounts(periodic_freq_values, bin_edges) / totalnumber * 100;
    
    %% 创建Periodic细胞频率分布图 / Create Periodic cell frequency distribution plot
    figure;
    t = tiledlayout(2, 1, 'TileSpacing', 'none', 'Padding', 'compact');
    
    % 上子图：显示>15%的部分 / Upper subplot: Show >15% portion
    ax1 = nexttile(t);
    bar(ax1, bin_centers, frequency_counts, 'FaceColor', [0.6 0.6 0.6]); 
    hold(ax1, 'on');
    bar(ax1, bin_centers, periodic_counts, 'FaceColor', 'r');
    ylim(ax1, [15 max(frequency_counts) + 5]);
    xlim(ax1, [0.005 max(bin_centers)]);
    xticks(ax1, []);
    ylabel(ax1, 'Cell Percentage [%]');
    title(ax1, 'Frequency Distribution with Periodic Cells');
    legend(ax1, {'All Cells', 'Periodic Cells'}, 'Location', 'best');
    line(ax1, [0.005 max(bin_centers)], [15 15], 'Color', 'k', 'LineStyle', '--');
    
    % 下子图：显示≤5%的部分 / Lower subplot: Show ≤5% portion
    ax2 = nexttile(t);
    bar(ax2, bin_centers, frequency_counts, 'FaceColor', [0.6 0.6 0.6]); 
    hold(ax2, 'on');
    bar(ax2, bin_centers, periodic_counts, 'FaceColor', 'r');
    ylim(ax2, [0 5]);
    xlim(ax2, [0.005 max(bin_centers)]);
    xticks(ax2, 0.005:0.004:max(bin_centers));
    xlabel(ax2, 'Frequency [Hz]');
    ylabel(ax2, 'Cell Percentage [%]');
    line(ax2, [0.005 max(bin_centers)], [5 5], 'Color', 'k', 'LineStyle', '--');
    
    % 添加文本标签 / Add text labels
    for i = 1:length(frequency_counts)
        if frequency_counts(i) <= 5
            text(ax2, bin_centers(i), frequency_counts(i) + 0.3, ...
                sprintf('%.1f', frequency_counts(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
        end
        if periodic_counts(i) <= 5 && periodic_counts(i) > 0
            text(ax2, bin_centers(i), periodic_counts(i) + 0.3, ...
                sprintf('%.1f', periodic_counts(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'r');
        end
    end
    
    % 保存图像 / Save figure
    saveas(gcf, 'Freq_distri_with_periodic_broken_percent.tif');
    
    %% ----------------------- 频率调试信息 -----------------------
    %% ----------------------- Frequency Debug Information -----------------------
    [minFreq, minIdx] = min(hub_freq_values);
    [maxFreq, maxIdx] = max(hub_freq_values);
    [minFreq1, minIdx1] = min(periodic_freq_values);
    [maxFreq1, maxIdx1] = max(periodic_freq_values);
    [minFreq2, minIdx2] = min(frequency_values);
    [maxFreq2, maxIdx2] = max(frequency_values);
    
    % 显示统计信息 / Display statistical information
    fprintf('Periodic_min_frq %.4f Hz cellnumber %d\n', minFreq1, minIdx1);
    fprintf('Periodic_max_frq %.4f Hz cellnumber %d\n', maxFreq1, maxIdx1);
    fprintf('Hub_min_frq %.4f Hz cellnumber %d\n', minFreq, minIdx);
    fprintf('Hub_max_frq %.4f Hz cellnumber %d\n', maxFreq, maxIdx);
    fprintf('All_min_frq %.4f Hz cellnumber %d\n', minFreq2, minIdx2);
    fprintf('All_max_frq %.4f Hz cellnumber %d\n', maxFreq2, maxIdx2);
    
    %% 询问用户是否继续 / Ask user if continue
    choice = inputdlg('End and go to the next plotting? Input y for Yes; any other key for No');
    if strcmp(choice, 'y')
        user_decision = false;
    else
        user_decision = true;
    end
end

end  % 函数结束 / End of function