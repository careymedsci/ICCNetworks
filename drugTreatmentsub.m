
 
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


function results_treatment = drugTreatmentsub()

% ==== excel to be laoded====
[filename, pathname] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel File');
if isequal(filename,0)
    disp('file selection canceled');
    return;
end
file = fullfile(pathname, filename);

% ==== Sheet  ====
sheet_options = {'cbx','gap26','suramin','apyrase','EGTA','Vera', 'skf', 'dan', 'AP5',...
    'NieCl2', 'fgf', 'tram34'};
[selection, ok] = listdlg('PromptString','choose Sheet:',...
                          'SelectionMode','single',...
                          'ListString',sheet_options);
if ~ok
    disp('Sheet cancelled');
    return;
end
sheetname = sheet_options{selection};
fprintf('Sheet name: %s\n', sheetname);

% ==== data row & lables ====
data_rows = [4,10,7,13,16,19,48,45,22,51];
x_labels = {'Ca^{2+} peaks (counts)','Active cells (%)','Mean peak width (all cell)',...
    'Co-active cells (%)','Hub cells (%)','Periodic cells (% of all cells)',...
    'Mean peak width (active cell, s)','Mean peak width (periodic cells, s)',...
    'Mean correlation value (of all cells)','connectivity'};


% detect the data to figure out 2 groups or 3 groups
row_data = readmatrix(file,'Sheet',sheetname,'Range',sprintf("A%d:R%d",data_rows(1),data_rows(1)));
num_nonNaN_cols = sum(~isnan(row_data));

if num_nonNaN_cols <= 12
    fprintf('2 groups data detected，analysis by t test or mw...\n');
    
    % ==================== analysis on 2 groups data ====================
    ctrl_color = [0, 113, 188]/255;
    treatment_color  = [164, 30, 50]/255;
    
    p_values = zeros(1, length(data_rows));
    results_treatment = struct('Metric', {}, 'p_value', {}, 'CI_low', {}, 'CI_high', {}, ...
                     'cliffs_d', {}, 'TestUsed', {}, 'p_mannwhitney', {}, 'cohen_d', {}, 'power', {}, ...
                     'normality_p_ctrl', {}, 'normality_p_treat', {}, 'test_used_final', {}, ...
                     'p_value_final', {}, 'cohen_d_final', {}, 'CI_low_final', {}, 'CI_high_final', {});
    
    for i = 1:length(data_rows)
        row_idx = data_rows(i);
        row_data = readmatrix(file, 'Sheet', sheetname, 'Range', sprintf("A%d:L%d", row_idx, row_idx));
        if isempty(row_data) || all(isnan(row_data))
            continue;
        end
        ctrl = row_data(1:6); ctrl = ctrl(~isnan(ctrl));
        treatment = row_data(7:12); treatment = treatment(~isnan(treatment));
        if isempty(ctrl) || isempty(treatment)
            continue;
        end
    
        ctrl_mean = mean(ctrl); ctrl_sem = std(ctrl)/sqrt(length(ctrl));
        treatment_mean = mean(treatment); treatment_sem = std(treatment)/sqrt(length(treatment));
    
        % === Mann-Whitney U test ===
        [p_wilcoxon, ~] = ranksum(ctrl, treatment);
        cliffs_d = cliffsDelta(ctrl, treatment);
        combined_data = [ctrl'; treatment'];
        group = [ones(length(ctrl),1); 2*ones(length(treatment),1)];
        diff_fun = @(data) median(data(group == 1)) - median(data(group == 2));
        ci_mw = bootci(2000, {diff_fun, combined_data}, 'type', 'percentile');
    
        % G*Power 
        cohen_d = cliffs_d * pi / sqrt(3);
        alpha = 0.05;
        try
            power_val = sampsizepwr('t2', [0 1], cohen_d, [], min(length(ctrl), length(treatment)), 'Alpha', alpha);
        catch
            power_val = NaN;
        end
    
        % normality
        [h_ctrl, p_sw_ctrl] = lillietest(ctrl, 'Alpha', 0.05);
        [h_treat, p_sw_treat] = lillietest(treatment, 'Alpha', 0.05);
        normal_ctrl = (h_ctrl == 0);
        normal_treat = (h_treat == 0);
    
        if normal_ctrl && normal_treat
            test_used_final = 't-test';
            [~, p_ttest, ci_ttest, ~] = ttest2(ctrl, treatment, 'Vartype', 'unequal');
            main_p_final = p_ttest;
            pooled_std = sqrt(((length(ctrl)-1)*var(ctrl) + (length(treatment)-1)*var(treatment)) / (length(ctrl)+length(treatment)-2));
            cohen_d_final = (mean(ctrl) - mean(treatment)) / pooled_std;
            CI_low_final = ci_ttest(1);
            CI_high_final = ci_ttest(2);
        else
            test_used_final = 'Mann-Whitney U';
            main_p_final = p_wilcoxon;
            cohen_d_final = cohen_d;
            CI_low_final = ci_mw(1);
            CI_high_final = ci_mw(2);
        end
    
        % plotting
        figure('Color', 'w', 'Position', [100, 100, 400, 500]); hold on;
        scatter(ones(size(ctrl)) + (rand(1, length(ctrl)) - 0.5)*0.1, ctrl, 70, 'filled', 'MarkerFaceColor', ctrl_color);
        scatter(2*ones(size(treatment)) + (rand(1, length(treatment)) - 0.5)*0.1, treatment, 70, 'filled', 'MarkerFaceColor', treatment_color);
        errorbar(1, ctrl_mean, ctrl_sem, 'k', 'CapSize', 15, 'LineWidth', 1.5);
        errorbar(2, treatment_mean, treatment_sem, 'k', 'CapSize', 15, 'LineWidth', 1.5);
        y_range = max([ctrl, treatment]) - min([ctrl, treatment]);
        maxY = max([ctrl_mean + ctrl_sem, treatment_mean + treatment_sem]) + 0.2*y_range;
        plot([1 2], [maxY maxY], 'k-', 'LineWidth', 1.5);
        if main_p_final < 0.001, sig = '***';
        elseif main_p_final < 0.01, sig = '**';
        elseif main_p_final < 0.05, sig = '*';
        else, sig = 'ns'; end
        text(1.5, maxY * 1.05, sig, 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        xticks([1 2]); xticklabels({'Control', 'treatment'});
        ylabel(x_labels{i}, 'FontSize', 12, 'FontWeight', 'bold');
    
        % save
        results_treatment(i).Metric = x_labels{i};
        results_treatment(i).p_value = p_wilcoxon;
        results_treatment(i).CI_low = ci_mw(1);
        results_treatment(i).CI_high = ci_mw(2);
        results_treatment(i).cliffs_d = cliffs_d;
        results_treatment(i).TestUsed = 'Mann-Whitney U';
        results_treatment(i).p_mannwhitney = p_wilcoxon;
        results_treatment(i).cohen_d = cohen_d;
        results_treatment(i).power = power_val;
        results_treatment(i).normality_p_ctrl = p_sw_ctrl;
        results_treatment(i).normality_p_treat = p_sw_treat;
        results_treatment(i).test_used_final = test_used_final;
        results_treatment(i).p_value_final = main_p_final;
        results_treatment(i).cohen_d_final = cohen_d_final;
        results_treatment(i).CI_low_final = CI_low_final;
        results_treatment(i).CI_high_final = CI_high_final;
    end
    
    % to Excel
    T = struct2table(results_treatment);
    writetable(T, 'treatment_comparison_stats.xlsx');
    fprintf('\nSaved to treatment_comparison_stats.xlsx\n');

else
    fprintf('3 groups data detected，anaysis by kw and anova...\n');
    
    % ==================== 3 goups analysis ====================
    ctrl_color = [0, 113, 188]/255;
    treat1_color1 = [163, 30, 50]/255;  % Treat1 7-9
    treat1_color2 = [0, 0, 0]/255;  % Treat1 10-12
    treat2_color1 = [163, 30, 50]/255;  % Treat2 13-15
    treat2_color2 = [0, 0, 0]/255;  % Treat2 16-18

    results_treatment = struct();
    
    for i = 1:length(data_rows)
        row_idx = data_rows(i);
        row_data = readmatrix(file,'Sheet',sheetname,'Range',sprintf("A%d:R%d",row_idx,row_idx));
        if isempty(row_data) || all(isnan(row_data))
            continue;
        end
        
        ctrl = row_data(1:6); ctrl = ctrl(~isnan(ctrl));
        treat1 = row_data(7:12); treat1 = treat1(~isnan(treat1));
        treat2 = row_data(13:18); treat2 = treat2(~isnan(treat2));

        all_data = [ctrl(:); treat1(:); treat2(:)];
        group = [repmat({'Control'}, numel(ctrl),1); repmat({'Conc1'}, numel(treat1),1); repmat({'Conc2'}, numel(treat2),1)];

        % ANOVA + Dunnett
        [p_anova,~,stats] = anova1(all_data,group,'off');
        c_dunnett = multcompare(stats,'CType','dunnett','Display','off');

        % Kruskal-Wallis + Dunn
        p_kw = kruskalwallis(all_data,group,'off');
        c_dunn = multcompare(stats,'CType','dunn-sidak','Display','off');

        % plotting
        figure('Color','w','Position',[100,100,500,500]); hold on;

        % Control 1-6
        for k = 1:length(ctrl)
            scatter(1+(rand-0.5)*0.1, ctrl(k), 70, 'filled', 'MarkerFaceColor', ctrl_color);
        end

        % Treat1
        for k = 1:length(treat1)
            col_idx = k + 6; 
            if col_idx >= 7 && col_idx <= 9
                color = treat1_color1;
            else
                color = treat1_color2;
            end
            scatter(2+(rand-0.5)*0.1, treat1(k), 70, 'filled', 'MarkerFaceColor', color);
        end

        % Treat2
        for k = 1:length(treat2)
            col_idx = k + 12; 
            if col_idx >= 13 && col_idx <= 15
                color = treat2_color1;
            else
                color = treat2_color2;
            end
            scatter(3+(rand-0.5)*0.1, treat2(k), 70, 'filled', 'MarkerFaceColor', color);
        end

        % Mean+SEM
        means = [mean(ctrl), mean(treat1), mean(treat2)];
        sems  = [std(ctrl)/sqrt(length(ctrl)), std(treat1)/sqrt(length(treat1)), std(treat2)/sqrt(length(treat2))];
        errorbar(1:3, means, sems, 'k', 'CapSize', 15, 'LineWidth', 1.5, 'LineStyle', 'none');

        % Dunnett significance
        y_range = max(all_data)-min(all_data);
        y_base = max(means+sems)+0.2*y_range;
        offset = 0.15*y_range;
        for cc = 1:size(c_dunnett,1)
            g1 = c_dunnett(cc,1); g2 = c_dunnett(cc,2); pval = c_dunnett(cc,6);
            if pval<0.001, sig='***';
            elseif pval<0.01, sig='**';
            elseif pval<0.05, sig='*';
            else, sig='ns'; end
            plot([g1 g2], [y_base y_base], 'k-', 'LineWidth', 1.5);
            text(mean([g1 g2]), y_base+0.05*y_range, sig, 'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
            y_base = y_base + offset;
        end

        xticks([1 2 3]); xticklabels({'Control','Conc1','Conc2'});
        ylabel(x_labels{i}, 'FontSize', 12, 'FontWeight', 'bold');

        % save
        results_treatment(i).Metric = x_labels{i};
        results_treatment(i).anova_p = p_anova;
        results_treatment(i).dunnett = c_dunnett;
        results_treatment(i).kw_p = p_kw;
        results_treatment(i).dunn = c_dunn;
    end

    %  to Excel
    results_excel = results_treatment;
    for i = 1:length(results_treatment)
        if isfield(results_treatment(i),'dunnett')
            results_excel(i).dunnett = mat2str(results_treatment(i).dunnett);
        end
        if isfield(results_treatment(i),'dunn')
            results_excel(i).dunn = mat2str(results_treatment(i).dunn);
        end
    end

    T = struct2table(results_excel);
    writetable(T, 'treatment_comparison_stats.xlsx');
    fprintf('\nSaved to treatment_comparison_stats.xlsx\n');
end

save('treatment_comparison_stats.mat')
end
