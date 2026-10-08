
 
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

function peaks_and_coordinating()
    % ==============================================================
    % Compare two sets of paired variables:
    % (1) Hub vs Non-Hub
    % (2) Periodic vs Non-Periodic
    % ==============================================================

    % =======================
    % 1. Select Excel file
    % =======================
    [filename, pathname] = uigetfile({'*.xls;*.xlsx', 'Excel Files (*.xls, *.xlsx)'}, ...
                                     'Select Excel File');
    if isequal(filename, 0)
        disp('User canceled file selection.');
        return;
    end
    fullpath = fullfile(pathname, filename);

    % =======================
    % 2. Select sheet name
    % =======================
% ======= Sheet selection / 选择工作表名称 =======
sheet_options = {'B16-F10.N3','B16-BrM.3'};
[selection, ok] = listdlg('PromptString', 'Select sheet name:', ...
                          'SelectionMode', 'single', ...
                          'ListString', sheet_options);
if ~ok
    disp('Sheet selection cancelled.');
    return;
end
sheetname = sheet_options{selection};
fprintf('Selected sheet: %s\n', sheetname);

    % =======================
    % 3. Load data from selected sheet
    % =======================
    data = readmatrix(fullpath, 'Sheet', sheetname);

    % =======================
    % 4. Define colors
    % =======================
    color_ctrl = [0 113 188] / 255;  % Blue (control)
    color_exp  = [163 30 50] / 255;  % Red (experimental)

    % =======================
    % 5. Extract datasets
    % =======================
    hub     = data(57, :);
    nonhub  = data(58, :);
    periodic     = data(53, :);
    nonperiodic  = data(54, :);

    % Split control and experiment groups
    hub_ctrl        = hub(1:6);
    nonhub_ctrl     = nonhub(1:6);
    hub_exp         = hub(7:end);
    nonhub_exp      = nonhub(7:end);
    periodic_ctrl   = periodic(1:6);
    nonperiodic_ctrl= nonperiodic(1:6);
    periodic_exp    = periodic(7:end);
    nonperiodic_exp = nonperiodic(7:end);

    % =======================
    % 6. Plot: Hub vs Non-Hub
    % =======================
    figure('Color', 'w', 'Name', 'Hub vs Non-hub');
    hold on;

    % Scatter plots
    scatter(nonhub_ctrl, hub_ctrl, 80, 'filled', 'MarkerFaceColor', color_ctrl);
    scatter(nonhub_exp, hub_exp, 80, 'filled', 'MarkerFaceColor', color_exp);

    % Ratio lines
    max_val = max([hub, nonhub]) * 1.1;
    x_line = linspace(0, max_val, 200);
    plot(x_line, x_line, 'k--', 'LineWidth', 1.5);
    plot(x_line, 2*x_line, 'k:', 'LineWidth', 1.2);
    plot(x_line, 0.5*x_line, 'k-.', 'LineWidth', 1.2);

    xlabel('Non-hub (%)');
    ylabel('Hub (%)');
    title('Hub vs Non-hub periodicity ratio');
    axis equal;
    xlim([0 max_val]);
    ylim([0 max_val]);
    grid on; box off;

    % Statistics
    [~, pval1, ~, stats1] = ttest(hub, nonhub);
    sig_label1 = getSigLabel(pval1);

    text(0.65*max_val, 0.15*max_val, ...
        sprintf('Paired t-test\np = %.4f\n%s', pval1, sig_label1), ...
        'FontSize', 11, 'Color', 'k', ...
        'BackgroundColor', [1 1 1 0.6], 'Margin', 4);

    hold off;

    % =======================
    % 7. Plot: Periodic vs Non-Periodic
    % =======================
    figure('Color', 'w', 'Name', 'Mean number of peaks');
    hold on;

    scatter(nonperiodic_ctrl, periodic_ctrl, 80, 'filled', 'MarkerFaceColor', color_ctrl);
    scatter(nonperiodic_exp, periodic_exp, 80, 'filled', 'MarkerFaceColor', color_exp);

    max_val2 = max([periodic, nonperiodic]) * 1.1;
    x_line2 = linspace(0, max_val2, 200);
    plot(x_line2, x_line2, 'k--', 'LineWidth', 1.5);
    plot(x_line2, 2*x_line2, 'k:', 'LineWidth', 1.2);
    plot(x_line2, 0.5*x_line2, 'k-.', 'LineWidth', 1.2);

    xlabel('Non-periodic (%)');
    ylabel('Periodic (%)');
    title('Mean number of peaks');
    axis equal;
    xlim([0 max_val2]);
    ylim([0 max_val2]);
    grid on; box off;

    [~, pval2, ~, stats2] = ttest(periodic, nonperiodic);
    sig_label2 = getSigLabel(pval2);

    text(0.65*max_val2, 0.15*max_val2, ...
        sprintf('Paired t-test\np = %.4f\n%s', pval2, sig_label2), ...
        'FontSize', 11, 'Color', 'k', ...
        'BackgroundColor', [1 1 1 0.6], 'Margin', 4);

    hold off;

    % =======================
    % 8. Print + save results
    % =======================
    fprintf('\nHub vs Non-hub: t(%d)=%.3f, p=%.4f (%s)\n', stats1.df, stats1.tstat, pval1, sig_label1);
    fprintf('Periodic vs Non-periodic: t(%d)=%.3f, p=%.4f (%s)\n', stats2.df, stats2.tstat, pval2, sig_label2);

    analysis_summary(end+1, 1:3) = {"Hub vs Non-hub", stats1.tstat, pval1};
    analysis_summary(end+1, 1:3) = {"Periodic vs Non-periodic", stats2.tstat, pval2};

end


% =======================
%  Helper function
% =======================
function sig_label = getSigLabel(p)
    if p < 0.001
        sig_label = '***';
    elseif p < 0.01
        sig_label = '**';
    elseif p < 0.05
        sig_label = '*';
    else
        sig_label = 'n.s.';
    end
end
