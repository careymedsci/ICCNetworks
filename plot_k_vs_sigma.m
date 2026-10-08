
 
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

function plot_k_vs_sigma_drug()
% Xiao Liu - 2025-7-29 (modified 2026-02-16)
% Plots σ vs <k> for cells, compared against theoretical curve σ = sqrt(k)
% Data are read from the 'Scale-free-test' sheet of an Excel file.
% All points are plotted in black.

    % ======= Select Excel file =======
    [file, path] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel file');
    if isequal(file, 0)
        error('No file selected.');
    end
    filepath = fullfile(path, file);

    % ======= Sheet selection =======
    sheetname = 'Scale-free-test';
    fprintf('Selected sheet: %s\n', sheetname);

    % ======= Read data (rows 2-21, columns 2 & 3) =======
    data = readmatrix(filepath, 'Sheet', sheetname, 'Range', 'B2:C21');
    k_vals = data(:,1);      % 第2列为k
    sigma_vals = data(:,2);  % 第3列为σ

    % Remove NaN
    valid_idx = ~isnan(k_vals) & ~isnan(sigma_vals);
    k_vals = k_vals(valid_idx);
    sigma_vals = sigma_vals(valid_idx);

    % ======= Create Figure =======
    figure;
    hold on;
    set(gcf, 'Color', 'w');

    % ======= Plot theoretical shaded area σ = sqrt(k) =======
    k_range = linspace(0, max(k_vals)*1.2, 200);
    sigma_theory = sqrt(k_range);
    fill([k_range, fliplr(k_range)], [sigma_theory, zeros(size(sigma_theory))], ...
        [0.85 0.85 0.85], 'EdgeColor', 'none'); % light grey

    % ======= Plot theoretical curve =======
    plot(k_range, sigma_theory, 'k--', 'LineWidth', 2, ...
        'DisplayName', '\sigma = \surd<k>');

    % ======= Plot data points (all black) =======
    scatter(k_vals, sigma_vals, 50, 'k', 'filled');

    % ======= Axis formatting =======
    xlabel('<k> (mean number of co-active cells per cell)', 'FontSize', 12);
    ylabel('\sigma (std. dev. of connections per cell)', 'FontSize', 12);
    title(['\sigma vs <k> — Sheet: ' sheetname], 'FontSize', 14);
    legend('Location', 'NorthWest', 'Box', 'off');
    grid on;
    box off;

    % ======= Axis limits =======
    xlim([0, max(k_vals) * 1.2]);
    ylim([0, max(sigma_vals) * 1.2]);
end
