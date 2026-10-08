
 
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
    % 2025-7-29 (modified 2025-8-23)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function plots σ vs <k> for control and drug groups, 
    % compared against the theoretical curve σ = sqrt(k)


    % ======= Select Excel file / 选择Excel文件 =======
    [file, path] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel file');
    if isequal(file, 0)
        error('No file selected.');
    end
    filepath = fullfile(path, file);

    % ======= Sheet selection / 选择工作表名称 =======
    sheet_options = {'cbx', 'gap26', 'Suramin', 'apyrase', 'EGTA', 'Vera', 'skf', 'dan','AP5','B16-F10.N3','B16-BrM.3',...
       'NieCl2', 'fgf', 'tram34' };
    [selection, ok] = listdlg('PromptString', 'Select sheet name:', ...
                              'SelectionMode', 'single', ...
                              'ListString', sheet_options);
    if ~ok
        disp('Sheet selection cancelled.');
        return;
    end
    sheetname = sheet_options{selection};
    fprintf('Selected sheet: %s\n', sheetname);

    % ======= Read data (up to 18 columns) =======
    k_vals     = readmatrix(filepath, 'Sheet', sheetname, 'Range', 'A41:R41');
    sigma_vals = readmatrix(filepath, 'Sheet', sheetname, 'Range', 'A42:R42');

    ncols = max(find(~isnan(k_vals) | ~isnan(sigma_vals))); 

    % ======= Colors =======
    color_blue   = [0, 113, 188] / 255;   % blue色
    color_red    = [163, 30, 50] / 255;   % red色
    color_black  = [0, 0, 0];             % black色
    color_fill   = [215, 215, 215] / 255; % light grey浅灰

    % ======= Create Figure =======
    figure;
    hold on;
    set(gcf, 'Color', 'w');

    % ======= Plot theoretical shaded area σ = sqrt(k) =======
    valid_k = k_vals(~isnan(k_vals));
    k_range = linspace(0, max(valid_k) * 1.2, 200);
    sigma_theory = sqrt(k_range);
    fill([k_range, fliplr(k_range)], [sigma_theory, zeros(size(sigma_theory))], ...
        color_fill, 'EdgeColor', 'none');

    % ======= Plot theoretical curve =======
    plot(k_range, sigma_theory, 'k--', 'LineWidth', 2, ...
        'DisplayName', '\sigma = \surd<k>');

    % ======= Plot groups =======
    if ncols <= 12
        % ---- 原规则 ----
        for i = 1:min(6, ncols)
            if isnan(k_vals(i)) || isnan(sigma_vals(i)), continue; end
            scatter(k_vals(i), sigma_vals(i), 80, ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', color_blue);
            text(k_vals(i)+0.1, sigma_vals(i), ['C' num2str(i)], ...
                'FontSize', 10, 'Color', color_blue);
        end
        for i = 7:min(12, ncols)
            if isnan(k_vals(i)) || isnan(sigma_vals(i)), continue; end
            scatter(k_vals(i), sigma_vals(i), 80, ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', color_red);
            text(k_vals(i)+0.1, sigma_vals(i), ['E' num2str(i-6)], ...
                'FontSize', 10, 'Color', color_red);
        end
    else
        % ---- ＞12 col ----
        for i = 1:min(6, ncols)
            if isnan(k_vals(i)) || isnan(sigma_vals(i)), continue; end
            scatter(k_vals(i), sigma_vals(i), 80, ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', color_blue);
            text(k_vals(i)+0.1, sigma_vals(i), ['C' num2str(i)], ...
                'FontSize', 10, 'Color', color_blue);
        end
        for i = 7:min(12, ncols)
            if isnan(k_vals(i)) || isnan(sigma_vals(i)), continue; end
            scatter(k_vals(i), sigma_vals(i), 80, ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', color_red);
            text(k_vals(i)+0.1, sigma_vals(i), ['E' num2str(i)], ...
                'FontSize', 10, 'Color', color_red);
        end
        for i = 13:min(18, ncols)
            if isnan(k_vals(i)) || isnan(sigma_vals(i)), continue; end
            scatter(k_vals(i), sigma_vals(i), 80, ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', color_black);
            text(k_vals(i)+0.1, sigma_vals(i), ['E' num2str(i)], ...
                'FontSize', 10, 'Color', color_black);
        end
    end

    % ======= Axis formatting =======
    xlabel('<k> (mean number of co-active cells per cell)', 'FontSize', 12);
    ylabel('\sigma (std. dev. of connections per cell)', 'FontSize', 12);
    title(['\sigma vs <k> — Sheet: ' sheetname], 'FontSize', 14);
    legend('Location', 'NorthWest', 'Box', 'off');
    grid on;
    box off;

    % ======= Axis limits =======
    all_k = k_vals(~isnan(k_vals));
    all_sigma = sigma_vals(~isnan(sigma_vals));
    xlim([0, max(all_k) * 1.2]);
    ylim([0, max(all_sigma) * 1.2]);
end
