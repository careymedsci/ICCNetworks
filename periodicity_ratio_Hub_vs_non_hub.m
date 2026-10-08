

 
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


function periodicity_ratio_Hub_vs_non_hub()
    
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
    % 2025.03.09
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ==========================
    % 1. read Excel data file name: summary_on_all_math_fit.xlsx
    % ==========================
    [filename, pathname] = uigetfile({'*.xls;*.xlsx', 'Excel Files (*.xls, *.xlsx)'}, ...
                                     'Choose Excel file');
    if isequal(filename,0)
        disp('User canceled file selection.');
        return;
    end

    fullpath = fullfile(pathname, filename);
    opts = detectImportOptions(fullpath, 'Sheet', 'periodic rate');
    data = readtable(fullpath, opts);

    nonhub = data{:,2}';
    hub = data{:,3}';

    fprintf('Loaded %d paired samples.\n', numel(nonhub));

    % ==========================
    % 2. pair t test
    % ==========================
    [~, pval, ~, stats] = ttest(hub, nonhub);

    % significance label
    if pval < 0.001
        sig_label = '***';
    elseif pval < 0.01
        sig_label = '**';
    elseif pval < 0.05
        sig_label = '*';
    else
        sig_label = 'n.s.';
    end

    % ==========================
    % 3. plotting
    % ==========================
    figure('Color','w','Name','Hub vs Non-hub Periodic Cells');
    hold on;

    % 散点
    scatter(nonhub, hub, 60, 'filled', ...
            'MarkerFaceColor', [0.2 0.45 0.8], 'MarkerEdgeColor', 'none');

    max_val = max([nonhub hub]) * 1.1;
    x_line = linspace(0, max_val, 200);

    %  y = x, y = 2x, y = 0.5x lines
    h1 = plot(x_line, x_line, 'r--', 'LineWidth', 1.5);       % y = x
    h2 = plot(x_line, 2*x_line, 'k:', 'LineWidth', 1.5);      % y = 2x
    h3 = plot(x_line, 0.5*x_line, 'k-.', 'LineWidth', 1.5);   % y = 0.5x

    xlabel('Non-hub (%)', 'FontSize', 12);
    ylabel('Hub (%)', 'FontSize', 12);
    title('Hub vs Non-hub percentage of periodic cells', 'FontSize', 13);
    axis equal;
    xlim([0 max_val]);
    ylim([0 max_val]);
    grid on;
    box off;

    % legend
    legend([h1 h2 h3], {'y = x', 'y = 2x', 'y = 0.5x'}, 'Location', 'northwest', 'FontSize', 10);

    % p
    text(0.65*max_val, 0.15*max_val, ...
        sprintf('Paired t-test\np = %.4f\n%s', pval, sig_label), ...
        'FontSize', 11, 'Color', 'k', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'BackgroundColor', [1 1 1 0.6], ...
        'Margin', 4);

    % ==========================
    % 4. mean ± SEM
    % ==========================

    mean_nonhub = mean(nonhub);
    sem_nonhub = std(nonhub) / sqrt(length(nonhub));
    mean_hub = mean(hub);
    sem_hub = std(hub) / sqrt(length(hub));


    offset = 0.02 * max_val;   % 可调整


    errorbar_x = errorbar(mean_nonhub, offset, sem_nonhub, 'horizontal', ...
        'Color', [0.2 0.45 0.8], 'LineWidth', 1.5, 'CapSize', 6);

    plot(mean_nonhub, offset, 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', [0.2 0.45 0.8], 'MarkerEdgeColor', 'none');

    errorbar_y = errorbar(offset, mean_hub, sem_hub, 'vertical', ...
        'Color', [0.2 0.45 0.8], 'LineWidth', 1.5, 'CapSize', 6);
    plot(offset, mean_hub, 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', [0.2 0.45 0.8], 'MarkerEdgeColor', 'none');

    hold off;

    % ==========================
    % 5. results
    % ==========================
    mean_ratio = mean(hub ./ nonhub);
    mean_diff = mean(hub - nonhub);
    fprintf('Mean ratio (Hub / Non-hub) = %.3f\n', mean_ratio);
    fprintf('Mean difference = %.3f\n', mean_diff);
    fprintf('Paired t-test: t(%d) = %.3f, p = %.4f (%s)\n', ...
            stats.df, stats.tstat, pval, sig_label);

end

    % % ==========================
    % % 1. read Excel data file name: summary_on_all_math_fit.xlsx
    % % ==========================
    % [filename, pathname] = uigetfile({'*.xls;*.xlsx', 'Excel Files (*.xls, *.xlsx)'}, ...
    %                                  'Choose Excel file');
    % if isequal(filename,0)
    %     disp('User canceled file selection.');
    %     return;
    % end
    % 
    % fullpath = fullfile(pathname, filename);
    % opts = detectImportOptions(fullpath, 'Sheet', 'periodic rate');
    % data = readtable(fullpath, opts);
    % 
    % nonhub = data{:,2}';
    % hub = data{:,3}';
    % 
    % fprintf('Loaded %d paired samples.\n', numel(nonhub));
    % 
    % % ==========================
    % % 2. pair t test
    % % ==========================
    % [~, pval, ~, stats] = ttest(hub, nonhub);
    % 
    % % significance label
    % if pval < 0.001
    %     sig_label = '***';
    % elseif pval < 0.01
    %     sig_label = '**';
    % elseif pval < 0.05
    %     sig_label = '*';
    % else
    %     sig_label = 'n.s.';
    % end
    % 
    % % ==========================
    % % 3. plotting
    % % ==========================
    % figure('Color','w','Name','Hub vs Non-hub Periodic Cells');
    % hold on;
    % 
    % 
    % scatter(nonhub, hub, 60, 'filled', ...
    %         'MarkerFaceColor', [0.2 0.45 0.8], 'MarkerEdgeColor', 'none');
    % 
    % 
    % max_val = max([nonhub hub]) * 1.1;
    % x_line = linspace(0, max_val, 200);
    % 
    % %  y = x, y = 2x, y = 0.5x 
    % h1 = plot(x_line, x_line, 'r--', 'LineWidth', 1.5);       % y = x
    % h2 = plot(x_line, 2*x_line, 'k:', 'LineWidth', 1.5);      % y = 2x
    % h3 = plot(x_line, 0.5*x_line, 'k-.', 'LineWidth', 1.5);   % y = 0.5x
    % 
    % % 
    % xlabel('Non-hub (%)', 'FontSize', 12);
    % ylabel('Hub (%)', 'FontSize', 12);
    % title('Hub vs Non-hub percentage of periodic cells', 'FontSize', 13);
    % axis equal;
    % xlim([0 max_val]);
    % ylim([0 max_val]);
    % grid on;
    % box off;
    % 
    % % 
    % legend([h1 h2 h3], {'y = x', 'y = 2x', 'y = 0.5x'}, 'Location', 'northwest', 'FontSize', 10);
    % 
    % % P value
    % text(0.65*max_val, 0.15*max_val, ...
    %     sprintf('Paired t-test\np = %.4f\n%s', pval, sig_label), ...
    %     'FontSize', 11, 'Color', 'k', ...
    %     'HorizontalAlignment', 'left', ...
    %     'VerticalAlignment', 'middle', ...
    %     'BackgroundColor', [1 1 1 0.6], ...
    %     'Margin', 4);
    % 
    % hold off;
    % 
    % % ==========================
    % % 5. results
    % % ==========================
    % mean_ratio = mean(hub ./ nonhub);
    % mean_diff = mean(hub - nonhub);
    % fprintf('Mean ratio (Hub / Non-hub) = %.3f\n', mean_ratio);
    % fprintf('Mean difference = %.3f\n', mean_diff);
    % fprintf('Paired t-test: t(%d) = %.3f, p = %.4f (%s)\n', ...
    %         stats.df, stats.tstat, pval, sig_label);


