function results_stats = plot_freq_Min_Max_group_comparison_3groups(file, sheetname)
% Xiao Liu
% 2025-8-23
% 3-group comparison with manual violin plots for freq_min (left y-axis) and freq_max (right y-axis)

% ==== Load data ====
data_all = readmatrix(file, 'Sheet', sheetname);

% =============================
% 分组 | Group assignment
% =============================
% periodic
periodic.control_min  = data_all(23,1:6);
periodic.control_max  = data_all(24,1:6);
periodic.conc1_min    = data_all(23,7:12);
periodic.conc1_max    = data_all(24,7:12);
periodic.conc2_min    = data_all(23,13:18);
periodic.conc2_max    = data_all(24,13:18);

% all cells
globalcells.control_min  = data_all(27,1:6);
globalcells.control_max  = data_all(28,1:6);
globalcells.conc1_min    = data_all(27,7:12);
globalcells.conc1_max    = data_all(28,7:12);
globalcells.conc2_min    = data_all(27,13:18);
globalcells.conc2_max    = data_all(28,13:18);

% =============================
% Normality & homogeneity tests
% =============================
all_groups = {periodic.control_min, periodic.conc1_min, periodic.conc2_min, ...
              periodic.control_max, periodic.conc1_max, periodic.conc2_max, ...
              globalcells.control_min, globalcells.conc1_min, globalcells.conc2_min, ...
              globalcells.control_max, globalcells.conc1_max, globalcells.conc2_max};

normality_p = nan(size(all_groups));
for i = 1:numel(all_groups)
    x = all_groups{i}; x = x(~isnan(x));
    if numel(x) >= 3
        [~,p] = lillietest(x);
        normality_p(i) = p;
    end
end

all_data = []; all_group_labels = [];
for i = 1:numel(all_groups)
    g = all_groups{i};
    all_data = [all_data; g(:)];
    all_group_labels = [all_group_labels; repmat(i, numel(g), 1)];
end
p_levene = vartestn(all_data, all_group_labels, 'TestType','LeveneAbsolute','Display','off');

% =============================
% Parametric: Mixed ANOVA
% =============================
anova_results = struct();

% ===== freq_min =====
n_c = min(numel(globalcells.control_min), numel(periodic.control_min));
n_1 = min(numel(globalcells.conc1_min), numel(periodic.conc1_min));
n_2 = min(numel(globalcells.conc2_min), numel(periodic.conc2_min));
data_control = [globalcells.control_min(1:n_c)', periodic.control_min(1:n_c)'];
data_conc1   = [globalcells.conc1_min(1:n_1)',   periodic.conc1_min(1:n_1)'];
data_conc2   = [globalcells.conc2_min(1:n_2)',   periodic.conc2_min(1:n_2)'];

tbl1 = array2table([data_control; data_conc1; data_conc2], ...
                   'VariableNames', {'Global','Periodical'});
tbl1.Treatment = [repmat({'Control'}, n_c, 1);
                  repmat({'Concent1'}, n_1, 1);
                  repmat({'Concent2'}, n_2, 1)];
WithinDesign1 = table(categorical({'Global';'Periodical'}), 'VariableNames', {'CellType'});
rm1 = fitrm(tbl1, 'Global-Periodical ~ Treatment', 'WithinDesign', WithinDesign1);
ranovatbl1 = ranova(rm1, 'WithinModel', 'CellType');
mc1_celltype   = multcompare(rm1, 'CellType', 'By', 'Treatment', 'ComparisonType','bonferroni');
mc1_treatment  = multcompare(rm1, 'Treatment', 'By', 'CellType', 'ComparisonType','bonferroni');
anova_results.freq_min.ranova = ranovatbl1;
anova_results.freq_min.posthoc_celltype = mc1_celltype;
anova_results.freq_min.posthoc_treatment = mc1_treatment;

% ===== freq_max =====
n_c = min(numel(globalcells.control_max), numel(periodic.control_max));
n_1 = min(numel(globalcells.conc1_max), numel(periodic.conc1_max));
n_2 = min(numel(globalcells.conc2_max), numel(periodic.conc2_max));
data_control = [globalcells.control_max(1:n_c)', periodic.control_max(1:n_c)'];
data_conc1   = [globalcells.conc1_max(1:n_1)',   periodic.conc1_max(1:n_1)'];
data_conc2   = [globalcells.conc2_max(1:n_2)',   periodic.conc2_max(1:n_2)'];

tbl2 = array2table([data_control; data_conc1; data_conc2], ...
                   'VariableNames', {'Global','Periodical'});
tbl2.Treatment = [repmat({'Control'}, n_c, 1);
                  repmat({'Concent1'}, n_1, 1);
                  repmat({'Concent2'}, n_2, 1)];
WithinDesign2 = table(categorical({'Global';'Periodical'}), 'VariableNames', {'CellType'});
rm2 = fitrm(tbl2, 'Global-Periodical ~ Treatment', 'WithinDesign', WithinDesign2);
ranovatbl2 = ranova(rm2, 'WithinModel', 'CellType');
mc2_celltype   = multcompare(rm2, 'CellType', 'By', 'Treatment', 'ComparisonType','bonferroni');
mc2_treatment  = multcompare(rm2, 'Treatment', 'By', 'CellType', 'ComparisonType','bonferroni');
anova_results.freq_max.ranova = ranovatbl2;
anova_results.freq_max.posthoc_celltype = mc2_celltype;
anova_results.freq_max.posthoc_treatment = mc2_treatment;

% =============================
% Nonparametric: Kruskal-Wallis + Dunn’s test with effect sizes + CI
% =============================
nonparam_results = struct();
exp_names = {'periodic_min','periodic_max','global_min','global_max'};
data_sets = { ...
    {periodic.control_min, periodic.conc1_min, periodic.conc2_min}, ...
    {periodic.control_max, periodic.conc1_max, periodic.conc2_max}, ...
    {globalcells.control_min, globalcells.conc1_min, globalcells.conc2_min}, ...
    {globalcells.control_max, globalcells.conc1_max, globalcells.conc2_max} };

group_labels = {'Control','Concent1','Concent2'};
n_boot = 5000;  % bootstrap 次数

for e = 1:4
    data_e = data_sets{e};
    data_all = []; group_id = [];
    n_total = 0;
    for g = 1:3
        d = data_e{g}(:);
        data_all = [data_all; d];
        group_id = [group_id; repmat(g, numel(d), 1)];
        n_total = n_total + numel(d);
    end
    
    [p_kw, tbl, stats] = kruskalwallis(data_all, group_id, 'off');
    H = tbl{2,3};
    k = numel(data_e);
    epsilon2 = (H - k + 1)/(n_total - k);
    
    c = multcompare(stats, 'CType','dunn-sidak','Display','off');
    
    posthoc_named = struct();
    for kidx = 1:size(c,1)
        g1 = group_labels{c(kidx,1)};
        g2 = group_labels{c(kidx,2)};
        x1 = data_e{c(kidx,1)};
        x2 = data_e{c(kidx,2)};
        n1 = numel(x1); n2 = numel(x2);

        r = tiedrank([x1(:);x2(:)]);
        R1 = sum(r(1:n1));
        U = R1 - n1*(n1+1)/2;
        r_rb = 1 - 2*U/(n1*n2);

        rng('default');  
        r_boot = zeros(n_boot,1);
        for b = 1:n_boot
            samp1 = x1(randi(n1,n1,1));
            samp2 = x2(randi(n2,n2,1));
            r_tmp = tiedrank([samp1(:);samp2(:)]);
            R1b = sum(r_tmp(1:n1));
            U_b = R1b - n1*(n1+1)/2;
            r_boot(b) = 1 - 2*U_b/(n1*n2);
        end
        CI = prctile(r_boot,[2.5 97.5]);
        
        posthoc_named(kidx).comparison = [g1 ' vs ' g2];
        posthoc_named(kidx).r_rb_effect_size = r_rb;
        posthoc_named(kidx).row_p = c(kidx,6);
        posthoc_named(kidx).adjust_p = c(kidx,6);
        posthoc_named(kidx).CI = CI;
    end
    
    nonparam_results.(exp_names{e}).kw_p = p_kw;
    nonparam_results.(exp_names{e}).epsilon2 = epsilon2;
    nonparam_results.(exp_names{e}).posthoc = posthoc_named;
end

% =============================
% Plotting: manual violin + scatter + mean + SEM
% =============================
figure; hold on;
violin_width = 0.6;
blue  = [0,113,188]/255;
red   = [163,30,50]/255;
black = [0,0,0];

% --- freq_min 左轴 ---
freq_min_data = {globalcells.control_min, globalcells.conc1_min, globalcells.conc2_min, ...
                 periodic.control_min, periodic.conc1_min, periodic.conc2_min};
freq_min_colors = {blue, red, black, blue, red, black};
freq_min_labels = {'G-Ctrl','G-C1','G-C2','P-Ctrl','P-C1','P-C2'};
x_pos_min = 1:numel(freq_min_data);

ax1 = axes; hold(ax1,'on'); set(ax1,'Box','off');
ylabel(ax1,'freq_{min}');
for i = 1:numel(freq_min_data)
    data_i = freq_min_data{i}(:);
    if isempty(data_i), continue; end
    [f, xi] = ksdensity(data_i);
    if max(f)>0, f=f/max(f)*violin_width/2; else f=zeros(size(f)); end
    fill([x_pos_min(i)+f, x_pos_min(i)-fliplr(f)], [xi, fliplr(xi)], freq_min_colors{i}, ...
         'EdgeColor','none','FaceAlpha',0.6);
    jitter = (rand(size(data_i))-0.5)*0.3;
    scatter(x_pos_min(i)+jitter, data_i, 40, freq_min_colors{i}, 'filled','MarkerEdgeColor','none');
    mean_i = mean(data_i);
    sem_i = std(data_i,'omitnan')/sqrt(numel(data_i));
    plot(x_pos_min(i), mean_i,'ks','MarkerSize',10,'MarkerFaceColor','w','LineWidth',1.5);
    line([x_pos_min(i),x_pos_min(i)],[mean_i-sem_i, mean_i+sem_i],'Color','k','LineWidth',1.5);
    line([x_pos_min(i)-0.1,x_pos_min(i)+0.1],[mean_i-sem_i, mean_i-sem_i],'Color','k','LineWidth',1.5);
    line([x_pos_min(i)-0.1,x_pos_min(i)+0.1],[mean_i+sem_i, mean_i+sem_i],'Color','k','LineWidth',1.5);
end
xticks(x_pos_min); xticklabels(freq_min_labels); xtickangle(45);

% --- freq_max 右轴 ---
freq_max_data = {periodic.control_max, periodic.conc1_max, periodic.conc2_max, ...
                 globalcells.control_max, globalcells.conc1_max, globalcells.conc2_max};
freq_max_colors = {blue, red, black, blue, red, black};
freq_max_labels = {'P-Ctrl','P-C1','P-C2','G-Ctrl','G-C1','G-C2'};
x_pos_max = numel(freq_min_data) + (1:numel(freq_max_data));

ax2 = axes; hold(ax2,'on'); set(ax2,'YAxisLocation','right','Color','none','XTick',[]);
ylabel(ax2,'freq_{max}');
for i = 1:numel(freq_max_data)
    data_i = freq_max_data{i}(:);
    if isempty(data_i), continue; end
    [f, xi] = ksdensity(data_i);
    if max(f)>0, f=f/max(f)*violin_width/2; else f=zeros(size(f)); end
    fill([x_pos_max(i)+f, x_pos_max(i)-fliplr(f)], [xi, fliplr(xi)], freq_max_colors{i}, ...
         'EdgeColor','none','FaceAlpha',0.6);
    jitter = (rand(size(data_i))-0.5)*0.3;
    scatter(x_pos_max(i)+jitter, data_i, 40, freq_max_colors{i}, 'filled','MarkerEdgeColor','none');
    mean_i = mean(data_i);
    sem_i = std(data_i,'omitnan')/sqrt(numel(data_i));
    plot(x_pos_max(i), mean_i,'ks','MarkerSize',10,'MarkerFaceColor','w','LineWidth',1.5);
    line([x_pos_max(i),x_pos_max(i)],[mean_i-sem_i, mean_i+sem_i],'Color','k','LineWidth',1.5);
    line([x_pos_max(i)-0.1,x_pos_max(i)+0.1],[mean_i-sem_i, mean_i-sem_i],'Color','k','LineWidth',1.5);
    line([x_pos_max(i)-0.1,x_pos_max(i)+0.1],[mean_i+sem_i, mean_i+sem_i],'Color','k','LineWidth',1.5);
end

linkaxes([ax1,ax2],'x');

% =============================
% Save results
% =============================
results_stats = struct();
results_stats.normality_p = normality_p;
results_stats.p_levene = p_levene;
results_stats.anova_results = anova_results;
results_stats.nonparam_results = nonparam_results;
save('freq_results_stats.mat', 'results_stats');

end
