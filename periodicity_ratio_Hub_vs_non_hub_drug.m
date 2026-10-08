

 
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



function periodicity_ratio_Hub_vs_non_hub_drug()

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

    % ==========================
    % 1. Excel + Sheet selection
    % ==========================
    [filename, pathname] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel File');
    if isequal(filename,0)
        disp('file selection canceled');
        return;
    end
    file = fullfile(pathname, filename);

    sheet_options = {'cbx','gap26','suramin','apyrase','EGTA','Vera', ...
                     'skf','dan','AP5','NieCl2','fgf','tram34'};
    [selection, ok] = listdlg('PromptString','choose Sheet:', ...
                              'SelectionMode','single', ...
                              'ListString', sheet_options);
    if ~ok
        disp('Sheet cancelled');
        return;
    end
    sheetname = sheet_options{selection};
    fprintf('Sheet name: %s\n', sheetname);

    data = readmatrix(file, 'Sheet', sheetname);

    % ==========================
    % color definition
    % ==========================
    color_ctrl = [0 113 188] / 255;
    color_drug = [163 30 50] / 255;

    %% =========================================================
    % 2. Hub vs Non-hub (periodic ratio)
    %% =========================================================
    hub     = data(52, :);
    nonhub  = data(53, :);

    hub_ctrl      = hub(1:6);
    hub_drug      = hub(7:12);
    nonhub_ctrl   = nonhub(1:6);
    nonhub_drug   = nonhub(7:12);

    % paired t-test
    [~, pval, ~, stats] = ttest(hub, nonhub);

    % ===== plotting =====
    figure('Color','w','Name','Hub vs Non-hub periodic ratio');
    hold on;

    % scatter
    scatter(nonhub_ctrl, hub_ctrl, 70, 'filled', ...
        'MarkerFaceColor', color_ctrl, 'MarkerEdgeColor','none');
    scatter(nonhub_drug, hub_drug, 70, 'filled', ...
        'MarkerFaceColor', color_drug, 'MarkerEdgeColor','none');

    % axis range
    max_val = max([hub nonhub], [], 'omitnan') * 1.1;
    x = linspace(0, max_val, 400);

    % reference lines
    plot(x, x,     'k--', 'LineWidth', 1.6);      % y = x
    plot(x, 2*x,   'k:',  'LineWidth', 1.2);      % y = 2x
    plot(x, 0.5*x, 'k:',  'LineWidth', 1.2);      % y = 1/2 x

    axis equal
    xlim([0 max_val]); ylim([0 max_val]);

    xlabel('Non-hub cells (%)');
    ylabel('Hub cells (%)');
    title('Percentage of periodic cells');

    grid on;
    box off;

    % annotation
    text(0.05*max_val, 0.85*max_val, ...
        {'y = x','y = 2x','y = 0.5x'}, ...
        'FontSize',10,'Color',[0.2 0.2 0.2]);

    text(0.6*max_val, 0.15*max_val, ...
        sprintf('Paired t-test\np = %.4f', pval), ...
        'FontSize',11,'BackgroundColor',[1 1 1 0.75], ...
        'EdgeColor',[0.8 0.8 0.8]);

    hold off;

    %% =========================================================
    % 3. Mean number of peaks: periodic vs non-periodic
    %% =========================================================
    periodic     = data(56, :);
    nonperiodic  = data(57, :);

    per_ctrl     = periodic(1:6);
    per_drug     = periodic(7:12);
    nonper_ctrl  = nonperiodic(1:6);
    nonper_drug  = nonperiodic(7:12);

    % paired t-test
    [~, pval2, ~, stats2] = ttest(periodic, nonperiodic);

    % ===== plotting =====
    figure('Color','w','Name','Mean number of peaks');
    hold on;

    scatter(nonper_ctrl, per_ctrl, 70, 'filled', ...
        'MarkerFaceColor', color_ctrl, 'MarkerEdgeColor','none');
    scatter(nonper_drug, per_drug, 70, 'filled', ...
        'MarkerFaceColor', color_drug, 'MarkerEdgeColor','none');

    max_val2 = max([periodic nonperiodic], [], 'omitnan') * 1.1;
    x2 = linspace(0, max_val2, 400);

    plot(x2, x2,     'k--', 'LineWidth', 1.6);    % y = x
    plot(x2, 2*x2,   'k:',  'LineWidth', 1.2);    % y = 2x
    plot(x2, 0.5*x2, 'k:',  'LineWidth', 1.2);    % y = 1/2 x

    axis equal
    xlim([0 max_val2]); ylim([0 max_val2]);

    xlabel('Non-periodic cells');
    ylabel('Periodic cells');
    title('Mean number of peaks');

    grid on;
    box off;

    text(0.05*max_val2, 0.85*max_val2, ...
        {'y = x','y = 2x','y = 0.5x'}, ...
        'FontSize',10,'Color',[0.2 0.2 0.2]);

    text(0.6*max_val2, 0.15*max_val2, ...
        sprintf('Paired t-test\np = %.4f', pval2), ...
        'FontSize',11,'BackgroundColor',[1 1 1 0.75], ...
        'EdgeColor',[0.8 0.8 0.8]);

    hold off;

    %% ==========================
    % 4. console output
    %% ==========================
    fprintf('\n[Hub vs Non-hub]\n');
    fprintf('Mean ratio = %.3f\n', mean(hub ./ nonhub, 'omitnan'));
    fprintf('t(%d)=%.3f, p=%.4f\n', stats.df, stats.tstat, pval);

    fprintf('\n[Periodic vs Non-periodic]\n');
    fprintf('Mean ratio = %.3f\n', mean(periodic ./ nonperiodic, 'omitnan'));
    fprintf('t(%d)=%.3f, p=%.4f\n', stats2.df, stats2.tstat, pval2);

end
