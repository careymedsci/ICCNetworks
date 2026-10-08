
 
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


function layer_analysis()
% ==== excel to be loaded ====
[filename, pathname] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel File');
if isequal(filename,0)
    disp('File selection canceled');
    return;
end
file = fullfile(pathname, filename);

% ==== detect sheet names ====
[~, sheet_names] = xlsfinfo(file);
fprintf('Detected sheets: %s\n', strjoin(sheet_names, ', '));

% ==== color map ====
colors = {
    [0.3 0.3 0.3];  % Ca2+
    [0.0 0.6 0.0];  % Active
    [0.0 0.4 0.7];  % Co-active
    [1.0 0.5 0.0];  % Periodic
    [0.8 0.1 0.2];  % Hub
    [0.1 0.3 0.8]   % Coordinating rhythmic
    };

% ==== mode 1: only one sheet (original mode) ====
if numel(sheet_names) == 1
    sheetname = sheet_names{1};
    fprintf('Only one sheet detected → using %s\n', sheetname);

    data_rows = [62,10,13,16,19,58]; % include line 58
    x_labels = {'Ca^{2+} cells (%)','Active cells (%)','Co-active cells (%)',...
        'Periodic cells (%)','Hub cells (%)','Coordinating rhythmic cells ratio'};

    figure('Color', 'w', 'Position', [100, 100, 550, 500]); hold on;

    for i = 1:length(data_rows)
        row_idx = data_rows(i);
        row_data = readmatrix(file, 'Sheet', sheetname, 'Range', sprintf("A%d:F%d", row_idx, row_idx));
        if isempty(row_data) || all(isnan(row_data))
            continue;
        end
        y = row_data(~isnan(row_data));
        x = i * ones(size(y)) + (rand(size(y)) - 0.5) * 0.2;
        scatter(x, y, 60, 'filled', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
        errorbar(i, mean(y), std(y)/sqrt(length(y)), 'k', 'CapSize', 12, 'LineWidth', 1.2);
    end

% ==== mode 2: multiple sheets (read from "layer" and "periodic rate") ====
else
    if any(strcmpi(sheet_names, 'layer'))
        sheet_layer = 'layer';
    else
        sheet_layer = sheet_names{1};
        fprintf('Warning: No "layer" sheet found, using first sheet: %s\n', sheet_layer);
    end

    if any(strcmpi(sheet_names, 'periodic rate'))
        sheet_rate = 'periodic rate';
    else
        sheet_rate = [];
        fprintf('Warning: No "periodic rate" sheet found, skipping rhythmic ratio.\n');
    end

    % data definition (columns)
    col_idx = [2,3,4,5,6];  % for layer sheet
    metric_names = {'Ca^{2+} cells (%)','Active cells (%)','Co-active cells (%)',...
        'Periodic cells (%)','Hub cells (%)'};
    if ~isempty(sheet_rate)
        metric_names{6} = 'Coordinating rhythmic cells ratio';
    end

    figure('Color', 'w', 'Position', [100, 100, 550, 500]); hold on;

    % read columns from "layer"
    for i = 1:length(col_idx)
        data = readmatrix(file, 'Sheet', sheet_layer, 'Range', sprintf("%s3:%s100", ...
            excelColumnName(col_idx(i)), excelColumnName(col_idx(i))));
        data = data(~isnan(data));
        x = i * ones(size(data)) + (rand(size(data)) - 0.5) * 0.2;
        scatter(x, data, 60, 'filled', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
        errorbar(i, mean(data), std(data)/sqrt(length(data)), 'k', 'CapSize', 12, 'LineWidth', 1.2);
    end

    % read rhythmic ratio if available
    if ~isempty(sheet_rate)
        data_rhythm = readmatrix(file, 'Sheet', sheet_rate, 'Range', 'C2:C100');
        data_rhythm = data_rhythm(~isnan(data_rhythm));
        i = length(col_idx) + 1;
        x = i * ones(size(data_rhythm)) + (rand(size(data_rhythm)) - 0.5) * 0.2;
        scatter(x, data_rhythm, 60, 'filled', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
        errorbar(i, mean(data_rhythm), std(data_rhythm)/sqrt(length(data_rhythm)), 'k', 'CapSize', 12, 'LineWidth', 1.2);
    end
end

% ==== formatting (common for both modes) ====
num_metrics = 6;
xlim([0.5, num_metrics + 0.5]);
ylabel('Percentage (%)', 'FontSize', 12, 'FontWeight', 'bold');

xticks(1:num_metrics);
xticklabels({'Ca^{2+} cells','Active cells','Co-active cells','Periodic cells','Hub cells','Coord. rhythmic'});
xtickangle(30);

% background
patch([0.5 3.5 3.5 0.5], [0 0 100 100], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
patch([3.5 6.5 6.5 3.5], [0 0 100 100], [1.0 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% titles
text(2, -8, 'of All Cells', 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
text(5, -8, 'of Active Cells', 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.5 0 0]);

ylim([0, 100]);
box off;
set(gca, 'FontSize', 11, 'LineWidth', 1.2);
title('Combined Ca^{2+} cell classification', 'FontWeight', 'bold');

results_treatment = [];
disp('Plot done. No stats performed.');
end



% ==== helper: convert numeric col index to Excel column name ====
function name = excelColumnName(colNumber)
letters = '';
while colNumber > 0
    remainder = mod(colNumber - 1, 26);
    letters = [char(65 + remainder), letters];
    colNumber = floor((colNumber - 1) / 26);
end
name = letters;
end
