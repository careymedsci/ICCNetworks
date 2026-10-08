

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
% 2026.06.09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [params, analysis_summary] = analyzeNetworkConnectivity(params, analysis_summary)
% analyzeNetworkConnectivity - Analyzes network connectivity and properties
%
% Inputs:
%   params - Structure containing analysis parameters and data
%   analysis_summary - Structure containing analysis results
%
% Outputs:
%   params - Updated parameters structure with network metrics
%   analysis_summary - Updated analysis summary with connectivity results
%   Also creates and saves connectivity plots and data

%% SECTION 1: Initialize variables
% Extract variables from params struct
cacul_corr = params.cacul_corr;
plot_connectivity = params.plot_connectivity;
C = params.C;
corrleft = params.corrleft;
output = params.output;
mode = params.mode;
threshcorr = params.threshcorr;
plot_C_vs_D = params.plot_C_vs_D;
D = params.D;
dcutoff = params.dcutoff;


%% SECTION 2: Connectivity Analysis
if cacul_corr > 0
    if plot_connectivity > 0
        % Calculate connectivity for different cutoffs
        conn = zeros(11,2);
        for i = 0:10
            cutoff = i/10;
            conn(i+1,1) = cutoff;
            v = size(C(C>=cutoff))/2;
            if v(1)>0
                conn(i+1,2) = v(1) / corrleft;
            else
                conn(i+1,2) = 0;
            end
        end

        % Plot connectivity vs cutoff
        figure('Name','connectivity vs. cut-off')
        semilogy(conn(:,1),conn(:,2));
        hold on
        a = semilogy(conn(:,1),conn(:,2),'k.');
        set(a,'MarkerSize',10);
        xlim([0 1]);
        ylim([0.00001 1]);
        xlabel("cut-off",'FontSize', 14);
        ylabel("connectivity",'FontSize', 14);
        strValues = strtrim(cellstr(num2str(conn(:,2))));
        text(conn(:,1),conn(:,2),strValues,'VerticalAlignment','bottom');
        hold off
        drawnow

        % Save connectivity data
        outputValue = sprintf('%s_connectivity%s.fig',output,num2str(mode));
        savefig(outputValue);
        outputValue = sprintf('%s_connectivity%s.txt',output,num2str(mode));
        writematrix(conn, outputValue);

        % Update analysis summary
        allcorr = size(C(C>threshcorr))/2;
        analysis_summary(21,1) = "all correlations above threshcorr";
        analysis_summary(21,2) = allcorr(1);
        analysis_summary(22,1) = "connectivity (all correlations above threshcorr / corrleft)";
        analysis_summary(22,2) = allcorr(1) / corrleft;
    end
end

%% SECTION 3: Correlation vs Distance Analysis
if cacul_corr > 0
    if plot_C_vs_D > 0
        figure('Name','Correlation vs. Distance');
        hold on;

        % ---- 1. 散点图（左侧 y 轴） ----
        yyaxis left;
        a = plot(D, C, 'k.');
        set(a, 'MarkerSize', 10);
        ylabel("correlation", 'FontSize', 14);
        xlabel("distance (um)", 'FontSize', 14);
        if mode ~= 4
            xlim([0 dcutoff]);
        end
        ylim([0 1]);

        % ---- 2. 百分位截止线（绿色虚线） ----
        if params.mode ~= 0
            nCells = size(C,1);
            mask = triu(true(nCells), 1);
            corr_vals = C(mask);
            corr_vals = corr_vals(~isnan(corr_vals));
            if ~isempty(corr_vals)
                sortedCorr = sort(corr_vals, 'ascend');
                Npairs = numel(sortedCorr);
                ce = params.im_percentiles/100;
                idxcutoff = ceil(ce * Npairs);
                if idxcutoff < 1
                    idxcutoff = 1;
                elseif idxcutoff > Npairs
                    idxcutoff = Npairs;
                end
                corrcutoff = sortedCorr(idxcutoff);
                yline(corrcutoff, 'g--', 'LineWidth', 1.5);
                txt_x = params.dcutoff * 0.02;  
                txt_y = corrcutoff + 0.02;   
                text(txt_x, txt_y, sprintf('%d%% cutoff: %.2f', params.im_percentiles, corrcutoff), ...
                    'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
            end
        end

        % ---- 3. 固定阈值线（红色实线） ----
        line([0, dcutoff], [threshcorr threshcorr], 'Color','red', 'LineWidth', 1.5);

        % ---- 4. 计算高于阈值的点的距离分布 ----
        idx_above = (C >= threshcorr);
        D_above = D(idx_above);
        binEdges = 0:1:dcutoff;
        counts = histcounts(D_above, binEdges);
        binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
        maxCount = max(counts);
        if maxCount == 0
            maxCount = 1;   % 防止除零
        end

        % 归一化计数到 [threshcorr, 1]
        normCounts = threshcorr + (counts / maxCount) * (1 - threshcorr);

        % ---- 5. 绘制从阈值线开始的矩形（条形） ----
        for i = 1:length(binCenters)
            x = binCenters(i) - 0.5;   % 宽度为1
            width = 1;
            y_bottom = threshcorr;
            y_top = normCounts(i);
            rectangle('Position', [x, y_bottom, width, y_top - y_bottom], ...
                      'FaceColor', [0.6392, 0.1176, 0.1961], ...
                      'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        % ---- 6. 添加右侧 y 轴（显示真实计数） ----
        yyaxis right;
        ylim([0 1]);   % 与左侧轴范围一致
        % 设置几个刻度位置（例如 5 个等分点），并对应计数
        tickPos = linspace(0, 1, 6);   % 0, 0.2, 0.4, 0.6, 0.8, 1.0
        tickLabels = round(linspace(0, maxCount, 6));
        yticks(tickPos);
        yticklabels(tickLabels);
        ylabel('Number of pairs (corr ≥ threshcorr)', 'FontSize', 14);
        % 右侧轴颜色与条形匹配（可选）
        ax = gca;
        ax.YAxis(2).Color = [0.6392, 0.1176, 0.1961];

        % ---- 7. 添加说明文字 ----
        text(0.02, 0.95, 'Bar height normalized to [thresh, 1]', ...
             'Units', 'normalized', 'FontSize', 10, 'Color', 'blue');

        hold off;
        drawnow;
        outputValue = sprintf('%s_corr_vs_distance%s.fig', output, num2str(mode));
        savefig(outputValue);
    end
end



%% SECTION 5: Update params and export variables
% Store network metrics in params
params.conn = conn;
if params.mode ~= 0
mode_label = sprintf('Mode %d %d%% cutoff evaluation', params.mode, params.im_percentiles);
analysis_summary(params.mode,5) = mode_label;
analysis_summary(params.mode,6) =  corrcutoff;
end

end