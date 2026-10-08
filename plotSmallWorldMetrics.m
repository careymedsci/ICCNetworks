
 
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

function S = plotSmallWorldMetrics()
% plotSmallWorldMetrics
% Reads small-world network metrics from Excel, plots λ and σ with scatter + SEM, 
% performs paired t-test, and calculates small-worldness index S.
%
% Output:
%   S - Small-worldness metric

    % ===== Select Excel File =====
    [filename, pathname] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel File');
    if isequal(filename,0)
        error('No file selected.');
    end
    filepath = fullfile(pathname, filename);

    % ===== Fixed sheet name =====
    sheetname = 'small-world-test';
    fprintf('Reading sheet: %s\n', sheetname);

    % ===== Read only 20 rows starting from row 3 =====
    row_start = 3;
    nrows = 20;

    % λ data (shortest path)
    dataL1 = readmatrix(filepath, 'Sheet', sheetname, ...
        'Range', sprintf('B%d:B%d', row_start, row_start+nrows-1)); % empirical
    dataL2 = readmatrix(filepath, 'Sheet', sheetname, ...
        'Range', sprintf('C%d:C%d', row_start, row_start+nrows-1)); % random

    % σ data (clustering coefficient)
    dataR1 = readmatrix(filepath, 'Sheet', sheetname, ...
        'Range', sprintf('D%d:D%d', row_start, row_start+nrows-1)); % empirical
    dataR2 = readmatrix(filepath, 'Sheet', sheetname, ...
        'Range', sprintf('E%d:E%d', row_start, row_start+nrows-1)); % random

    % Remove NaNs
    dataL1 = dataL1(~isnan(dataL1));
    dataL2 = dataL2(~isnan(dataL2));
    dataR1 = dataR1(~isnan(dataR1));
    dataR2 = dataR2(~isnan(dataR2));

    % ===== Colors & jitter =====
    color_emp = [0 0.4 1]; % blue
    color_rnd = [0 0 0];   % black
    scatter_jitter = 0.05;

    figure('Color','w'); hold on;

    % ===== Left Y-axis: λ =====
    yyaxis left
    xL_emp = 1 + scatter_jitter*randn(size(dataL1));
    xL_rnd = 1.4 + scatter_jitter*randn(size(dataL2));
    scatter(xL_emp, dataL1, 30, color_emp, 'filled');
    scatter(xL_rnd, dataL2, 30, color_rnd, 'filled');

    % Add SEM
    meanL = [mean(dataL1), mean(dataL2)];
    semL  = [std(dataL1)/sqrt(numel(dataL1)), std(dataL2)/sqrt(numel(dataL2))];
    errorbar([1, 1.4], meanL, semL, 'k', 'LineStyle','none','LineWidth',1.5);

    ylabel('Mean Shortest Path Length λ');

    % ===== Right Y-axis: σ =====
    yyaxis right
    xR_emp = 2.5 + scatter_jitter*randn(size(dataR1));
    xR_rnd = 2.9 + scatter_jitter*randn(size(dataR2));
    scatter(xR_emp, dataR1, 30, color_emp, 'filled');
    scatter(xR_rnd, dataR2, 30, color_rnd, 'filled');

    % Add SEM
    meanR = [mean(dataR1), mean(dataR2)];
    semR  = [std(dataR1)/sqrt(numel(dataR1)), std(dataR2)/sqrt(numel(dataR2))];
    errorbar([2.5, 2.9], meanR, semR, 'k', 'LineStyle','none','LineWidth',1.5);

    ylabel('Clustering Coefficient σ');

    % ===== X-axis =====
    xticks([1.2 2.7]);
    xticklabels({'λ','σ'});
    title('Small-world Network Metrics');

    % ===== Paired t-test & significance annotation =====
    yyaxis left
    valid_idx = ~isnan(dataL1) & ~isnan(dataL2);
    [~, pL] = ttest(dataL1(valid_idx), dataL2(valid_idx));
    ytxt = max([dataL1; dataL2])*1.05;
    text(1.2, ytxt, significance_text(pL), 'HorizontalAlignment','center','FontSize',10);

    yyaxis right
    valid_idx = ~isnan(dataR1) & ~isnan(dataR2);
    [~, pR] = ttest(dataR1(valid_idx), dataR2(valid_idx));
    ytxtR = max([dataR1; dataR2])*1.05;
    text(2.7, ytxtR, significance_text(pR), 'HorizontalAlignment','center','FontSize',10);

    grid on; box off; set(gca,'FontSize',12);

    % ===== Small-worldness index =====
    L_emp = mean(dataL1); L_rnd = mean(dataL2);
    C_emp = mean(dataR1); C_rnd = mean(dataR2);
    S = (C_emp/C_rnd)/(L_emp/L_rnd);
    fprintf('Small-worldness S = %.4f\n', S);
    str_S = sprintf('S = %.2f', S);
    annotation('textbox',[0.45,0.4,0.1,0.1],'String',str_S,...
        'FitBoxToText','on','EdgeColor','none','FontSize',12,'FontWeight','bold');
end

%% ===== Supporting function for significance =====
function txt = significance_text(p)
    if p < 0.0001
        txt = 'P < 0.0001';
    elseif p < 0.001
        txt = sprintf('P = %.4f', p);
    elseif p < 0.05
        txt = sprintf('P = %.3f', p);
    else
        txt = 'ns';
    end
end






% function S = plotSmallWorldMetrics()
% % plotSmallWorldMetrics
% % This function reads an Excel file containing network metrics data
% % (mean shortest path length and clustering coefficient), plots them on 
% % dual Y-axes with scatter plots and background highlighting, and calculates 
% % the small-worldness index S.
% %
% % 本函数读取包含网络指标的 Excel 文件（如平均最短路径长度和聚类系数），
% % 并在双 Y 轴上绘图，最后计算小世界性指数 S。
% %
% % Output:
% %   S - Small-worldness metric 小世界性指标
% 
%     % ------------------  Excel  ------------------
%     [filename, pathname] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel File');
%     if isequal(filename,0)
%         disp('No file selected.');
%         return;
%     end
%     filepath = fullfile(pathname, filename);
% 
%     % 
%     opts = detectImportOptions(filepath, 'VariableNamingRule','preserve');
%     opts.DataRange = '3:10000';
%     opts.SelectedVariableNames = opts.VariableNames(2:end);  
%     T = readtable(filepath, opts);
% 
% 
%     [~, ~, raw2] = xlsread(filepath);
%     group_labels = string(raw2(2,2:5));
% 
% 
%     dataL1 = T{:,1}; dataL2 = T{:,2};  
%     dataR1 = T{:,3}; dataR2 = T{:,4};  
% 
% 
%     dataL1 = dataL1(~isnan(dataL1));
%     dataL2 = dataL2(~isnan(dataL2));
%     dataR1 = dataR1(~isnan(dataR1));
%     dataR2 = dataR2(~isnan(dataR2));
% 
% 
%     colors = {[0 0.4 1], [1 0.7 0]};   
%     scatter_jitter = 0.05;  
%     violin_width = 0.3;
% 
%     figure('Color','w'); hold on;
% 
% 
%     yyaxis left;
%     xL1_center = 1;  xL2_center = 1 + violin_width;
%     xL1 = xL1_center + scatter_jitter*randn(size(dataL1));
%     xL2 = xL2_center + scatter_jitter*randn(size(dataL2));
%     scatter(xL1, dataL1, 20, colors{1}, 'filled');
%     scatter(xL2, dataL2, 20, colors{2}, 'filled');
% 
%     yl_left = ylim; h_left = diff(yl_left); y_base_left = yl_left(1);
% 
% 
%     fill([xL1_center-violin_width/2, xL1_center+violin_width/2, ...
%           xL1_center+violin_width/2, xL1_center-violin_width/2], ...
%          [y_base_left, y_base_left, y_base_left+h_left, y_base_left+h_left], ...
%          colors{1}*0.3 + [1 1 1]*0.7, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
%     fill([xL2_center-violin_width/2, xL2_center+violin_width/2, ...
%           xL2_center+violin_width/2, xL2_center-violin_width/2], ...
%          [y_base_left, y_base_left, y_base_left+h_left, y_base_left+h_left], ...
%          colors{2}*0.3 + [1 1 1]*0.7, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
%     ylabel('Mean Shortest Path Length 平均最短路径长度');
% 
% 
%     yyaxis right;
%     xR1_center = 3; xR2_center = 3 + violin_width;
%     xR1 = xR1_center + scatter_jitter*randn(size(dataR1));
%     xR2 = xR2_center + scatter_jitter*randn(size(dataR2));
%     scatter(xR1, dataR1, 20, colors{1}, 'filled');
%     scatter(xR2, dataR2, 20, colors{2}, 'filled');
% 
%     yl_right = ylim; h_right = diff(yl_right); y_base_right = yl_right(1);
% 
% 
%     fill([xR1_center-violin_width/2, xR1_center+violin_width/2, ...
%           xR1_center+violin_width/2, xR1_center-violin_width/2], ...
%          [y_base_right, y_base_right, y_base_right+h_right, y_base_right+h_right], ...
%          colors{1}*0.3 + [1 1 1]*0.7, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
%     fill([xR2_center-violin_width/2, xR2_center+violin_width/2, ...
%           xR2_center+violin_width/2, xR2_center-violin_width/2], ...
%          [y_base_right, y_base_right, y_base_right+h_right, y_base_right+h_right], ...
%          colors{2}*0.3 + [1 1 1]*0.7, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
%     ylabel('Clustering Coefficient 聚类系数');
% 
% 
%     xticks([1.25 3.25]);
%     xticklabels({'Mean Shortest Path', 'Clustering Coefficient'});
%     title('Small-world Network Metrics 小世界网络指标');
%     legend({'Empirical Data 实测', 'Random Data 随机'}, 'Location','northeast');
%     grid on; box off;
% 
% 
%     % Ensure NaN removal
%     L_empirical = mean(dataL1, 'omitnan');
%     L_random    = mean(dataL2, 'omitnan');
%     C_empirical = mean(dataR1, 'omitnan');
%     C_random    = mean(dataR2, 'omitnan');
% 
%     % Calculate S
%     S = (C_empirical / C_random) / (L_empirical / L_random);
%     fprintf('Small-worldness S = %.4f\n', S);
% end
