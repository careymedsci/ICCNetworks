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


function Fig_2A_summary_cells_postive_ca_events_matalab()

% 作者：Thomas Broggini 和 刘晓；法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Thomas Broggini and Xiao Liu
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

% Fig_2A_summary_cells_postive_ca_events_matalab
% This function reads an Excel file containing calcium event counts,
% processes the data, and visualizes cell population distribution using
% transformed axes with density, scatter, and bar overlay.
%
% No input/output needed. Prompts to select Excel file.

    clear; clc; close all;

    % Load Excel data
    [file, path] = uigetfile({'*.xlsx;*.xls'}, 'Load the Excel data');
    if isequal(file, 0)
        disp('File selection canceled.');
        return;
    end
    filename = fullfile(path, file);

    data = readtable(filename);
    group_names = data.Properties.VariableNames;
    data_array = table2array(data);

    group = [];
    value = [];
    for i = 1:size(data_array,2)
        valid_idx = ~isnan(data_array(:,i));
        group = [group; repmat(group_names(i), sum(valid_idx), 1)];
        value = [value; data_array(valid_idx,i)];
    end

    groups = unique(group);
    n_groups = length(groups);

    pos_rate = zeros(n_groups,1);
    mean_val = zeros(n_groups,1);
    sem_val = zeros(n_groups,1);
    for i = 1:n_groups
        idx = strcmp(group, groups{i});
        pos_rate(i) = mean(value(idx) ~= 0) * 100;
        mean_val(i) = mean(value(idx));
        sem_val(i) = std(value(idx)) / sqrt(sum(idx));
    end

    % Axis transformation functions
    stretch_1_15 = 10;
    compress_above500 = 0.5;
    transform_value = @(y) (y<=15).*y*stretch_1_15 + ...
                          (y>15 & y<=500).*(y-15+15*stretch_1_15) + ...
                          (y>500).*((y-500)*compress_above500 + (500-15+15*stretch_1_15));
    inverse_value = @(y) (y<=15*stretch_1_15).*(y/stretch_1_15) + ...
                         (y>15*stretch_1_15 & y<= (500-15+15*stretch_1_15)).*(y - 15*stretch_1_15 +15) + ...
                         (y>(500-15+15*stretch_1_15)).*((y - (500-15+15*stretch_1_15))/compress_above500 + 500);

    % Define custom group colors
    custom_colors = [
        64, 128, 128;
        199, 121, 134;
        25, 25, 112;
        150, 150, 150;
        128, 128, 0;
        102, 205, 170
    ] / 255;

    % Begin plotting
    figure('Color','w','Position',[100 100 1200 600]);
    hold on;

    % Plot data for each group
    for i = 1:n_groups
        idx = strcmp(group, groups{i});
        y_raw = value(idx);
        y = transform_value(y_raw);

        [f, xi] = ksdensity(y, 'Bandwidth', [], 'Function','pdf');
        f = f / max(f) * 0.2;
        fill(i + [0 f], [xi xi(end)], custom_colors(i,:), ...
            'FaceAlpha',0.5, 'EdgeColor','none');

        jitterX = i - 0.12*randn(size(y));
        scatter(jitterX, y, 35, custom_colors(i,:), 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',0.2, 'MarkerFaceAlpha',0.8);

        med = median(y);
        plot([i-0.1, i+0.1], [med med], 'k-', 'LineWidth',2);
    end

    % Plot mean ± SEM
    for i = 1:n_groups
        mean_y = transform_value(mean_val(i));
        sem_y = sem_val(i);
        sem_upper = transform_value(mean_val(i) + sem_y);
        sem_lower = transform_value(max(mean_val(i) - sem_y, 0));

        plot([i-0.07, i+0.07], [mean_y mean_y], 'k-', 'LineWidth',2);
        plot([i i], [sem_lower sem_upper], 'k-', 'LineWidth',1.5);

        cap_width = 0.07;
        plot([i-cap_width i+cap_width], [sem_upper sem_upper], 'k-', 'LineWidth',1.5);
        plot([i-cap_width i+cap_width], [sem_lower sem_lower], 'k-', 'LineWidth',1.5);
    end

    % Customize x-axis
    xlim([0.5 n_groups+0.5]);
    xticks(1:n_groups);
    xticklabels(groups);
    set(gca, 'XTickLabelRotation', 0, 'FontSize',16, 'FontWeight','bold');

    % Customize y-axis
    yticks_raw = [0 1 5 10 15 500 600 700 800 900 1000];
    yticks_trans = transform_value(yticks_raw);
    yticklabels_text = string(yticks_raw);
    set(gca, 'YTick', yticks_trans, 'YTickLabel', yticklabels_text, 'FontSize',16);

    ylim([0 max(transform_value([value; 1000]))*1.1]);

    % Add y-axis break lines
    yline(transform_value(15), 'k--', 'LineWidth', 1.5);
    yline(transform_value(500), 'k--', 'LineWidth', 1.5);

    for i = 1:n_groups
        plot([i-0.2, i+0.2], [transform_value(15) transform_value(15)], 'k', 'LineWidth',1.5)
        plot([i-0.2, i+0.2], [transform_value(500) transform_value(500)], 'k', 'LineWidth',1.5)
    end

    % Annotate bars
    for i = 1:n_groups
        str = sprintf('%.1f%%\n%.1f±%.1f', pos_rate(i), mean_val(i), sem_val(i));
        text(i+0.35, max(transform_value(value))*1.05, str, ...
            'Rotation',90, 'FontSize',11, 'FontWeight','bold', 'HorizontalAlignment','left', 'Color','k');
    end

    ylabel('Cell Number', 'FontSize',14);
    title('Calcium Event Counts Across Cell Lines', ...
        'FontWeight','bold', 'FontSize',18);

    box on;
    hold off;
end
