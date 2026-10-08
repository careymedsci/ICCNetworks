
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

function results = analyze_animal_part()
% === 文件选择 / File selection ===
[filename, pathname] = uigetfile({'*.xls;*.xlsx'}, 'Select Excel File');
if isequal(filename, 0)
    disp('文件选择已取消 / File selection canceled');
    return;
end
file = fullfile(pathname, filename);
sheetname = 1;

% === 指标 / Metrics ===
metrics = {'Peaks per 100 cells in 60s', ...
           'Mean peak width of all cells (s)'};

% === 数据位置 / Data positions ===
data_positions = struct( ...
    'Peaks', struct( ...
        'Original', [4, 'A', 'E'], ...
        'BrM1', [4, 'F', 'J'], ...
        'BrM2', [12, 'A', 'E'], ...
        'BrM3', [12, 'F', 'J'], ...
        'IC_BrM', [21, 'A', 'E'], ...
        'TV_BrM', [21, 'F', 'J']), ...
    'Width', struct( ...
        'Original', [7, 'A', 'E'], ...
        'BrM1', [7, 'F', 'J'], ...
        'BrM2', [15, 'A', 'E'], ...
        'BrM3', [15, 'F', 'J'], ...
        'IC_BrM', [24, 'A', 'E'], ...
        'TV_BrM', [24, 'F', 'J']) ...
);

% === 组别名称更新 / Group labels ===
group_names = {'Original F10','B16-BrM.1','B16-BrM.2','B16-BrM.3','B16-IC-BrM','B16-TV-BrM'};

% === Colors for plotting ===
colors = [
    0.2 0.2 0.2;
    0.2 0.5 0.9;
    0.9 0.3 0.3;
    0.2 0.7 0.4;
    0.7 0.4 0.9;
    0.3 0.7 0.7
];

results = struct();

for m = 1:length(metrics)
    metric_name = metrics{m};
    fprintf('\nProcessing metric: %s\n', metric_name);

    
    if contains(metric_name, 'Peaks')
        pos = data_positions.Peaks;
    else
        pos = data_positions.Width;
    end

    % ==== read data ====
    Original = readmatrix(file,'Sheet',sheetname,'Range',sprintf('%s%d:%s%d',pos.Original(2),pos.Original(1),pos.Original(3),pos.Original(1)));
    BrM1 = readmatrix(file,'Sheet',sheetname,'Range',sprintf('%s%d:%s%d',pos.BrM1(2),pos.BrM1(1),pos.BrM1(3),pos.BrM1(1)));
    BrM2 = readmatrix(file,'Sheet',sheetname,'Range',sprintf('%s%d:%s%d',pos.BrM2(2),pos.BrM2(1),pos.BrM2(3),pos.BrM2(1)));
    BrM3 = readmatrix(file,'Sheet',sheetname,'Range',sprintf('%s%d:%s%d',pos.BrM3(2),pos.BrM3(1),pos.BrM3(3),pos.BrM3(1)));
    IC_BrM = readmatrix(file,'Sheet',sheetname,'Range',sprintf('%s%d:%s%d',pos.IC_BrM(2),pos.IC_BrM(1),pos.IC_BrM(3),pos.IC_BrM(1)));
    TV_BrM = readmatrix(file,'Sheet',sheetname,'Range',sprintf('%s%d:%s%d',pos.TV_BrM(2),pos.TV_BrM(1),pos.TV_BrM(3),pos.TV_BrM(1)));

    groups = {Original, BrM1, BrM2, BrM3, IC_BrM, TV_BrM};
    all_data = [Original(:); BrM1(:); BrM2(:); BrM3(:); IC_BrM(:); TV_BrM(:)];
    group_idx = [repmat(1,numel(Original),1);
                 repmat(2,numel(BrM1),1);
                 repmat(3,numel(BrM2),1);
                 repmat(4,numel(BrM3),1);
                 repmat(5,numel(IC_BrM),1);
                 repmat(6,numel(TV_BrM),1)];

    % ==== normality test ====
    normal_flags = false(1,6);
    for g = 1:6
        data_g = groups{g};
        data_g = data_g(~isnan(data_g));
        if numel(data_g) < 3
            normal_flags(g) = true;
        else
            try
                [h, ~] = swtest(data_g, 0.05);
            catch
                [h, ~] = kstest((data_g-mean(data_g))/std(data_g));
            end
            normal_flags(g) = (h == 0);
        end
    end

    %% ==== ANOVA + Tukey ====
    [p_anova, ~, stats_anova] = anova1(all_data, group_idx, 'off');
    comp_tukey = multcompare(stats_anova,'CType','tukey-kramer','Display','off');
    tukey_results = cell(size(comp_tukey,1),4);
    for i = 1:size(comp_tukey,1)
        tukey_results{i,1} = sprintf('%s vs %s',group_names{comp_tukey(i,1)},group_names{comp_tukey(i,2)});
        tukey_results{i,2} = comp_tukey(i,6);
        tukey_results{i,3} = comp_tukey(i,4);
        tukey_results{i,4} = [comp_tukey(i,3),comp_tukey(i,5)];
    end
    tukey_table = cell2table(tukey_results,'VariableNames',{'Comparison','pValue','MeanDiff','95CI'});

    %% ==== Kruskal-Wallis + Dunn ====
    [p_kw, ~, stats_kw] = kruskalwallis(all_data, group_idx, 'off');
    comp_dunn = multcompare(stats_kw,'CType','dunn-sidak','Display','off');
    dunn_results = cell(size(comp_dunn,1),4);
    for i = 1:size(comp_dunn,1)
        dunn_results{i,1} = sprintf('%s vs %s',group_names{comp_dunn(i,1)},group_names{comp_dunn(i,2)});
        dunn_results{i,2} = comp_dunn(i,6);
    end
    dunn_table = cell2table(dunn_results,'VariableNames',{'Comparison','pValue','MeanDiff','95CI'});

    %% ==== save ====
    results(m).Metric = metric_name;
    results(m).GroupNames = group_names;
    results(m).Normality = normal_flags;
    results(m).ANOVA_Tukey.p_main = p_anova;
    results(m).ANOVA_Tukey.Pairwise = tukey_table;
    results(m).Kruskal_Dunn.p_main = p_kw;
    results(m).Kruskal_Dunn.Pairwise = dunn_table;

    %% ==== plot ====
    figure('Color','w','Position',[200,200,760,550]); hold on;
    for g = 1:6
        data_g = groups{g};
        data_g = data_g(~isnan(data_g));
        x_jitter = g + (rand(size(data_g)) - 0.5)*0.15;
        scatter(x_jitter, data_g, 70, 'filled', 'MarkerFaceColor', colors(g,:));
        errorbar(g, mean(data_g), std(data_g)/sqrt(length(data_g)), ...
            'k', 'CapSize', 15, 'LineWidth', 1.5);
    end
    xlim([0.5 6.5]);
    xticks(1:6);
    xticklabels(group_names);
    ylabel(metric_name,'FontSize',12,'FontWeight','bold');
    title(sprintf('%s\nANOVA p=%.4g | KW p=%.4g',metric_name,p_anova,p_kw),'FontSize',14);
    box on;
end

% === save results ===
save('stat_results_full.mat','results');
fprintf('\n✅ All analyses completed and saved to stat_results_full.mat\n');

end
