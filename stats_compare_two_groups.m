

function [report] = stats_compare_two_groups(ctrl, cbx)
% stats_compare_two_groups - Performs multiple statistical tests on two groups
% 输入: ctrl - control group (vector), cbx - treatment group (vector)
% 输出: report - struct 包含所有检验与效应量结果

% 中英文报告结构 / Initialize report
report = struct();

% 检查数据有效性 / Data check
if length(ctrl) < 3 || length(cbx) < 3
    error('每组数据至少需要3个样本 / At least 3 samples per group are required');
end

% 合并数据
all_data = [ctrl(:); cbx(:)];
labels = [repmat({'Control'}, length(ctrl), 1); repmat({'CBX'}, length(cbx), 1)];

% 正态性检验 / Normality Test
[~, p_ctrl_norm] = lillietest(ctrl);
[~, p_cbx_norm]  = lillietest(cbx);
report.normality.p_ctrl = p_ctrl_norm;
report.normality.p_cbx  = p_cbx_norm;

% 方差齐性检验 / Variance Equality Test
[p_levene, ~] = vartestn([ctrl(:); cbx(:)], [ones(length(ctrl),1); 2*ones(length(cbx),1)], 'TestType','LeveneAbsolute','Display','off');
report.variance_homogeneity.p = p_levene;

% 选择检验方法 / Choose Test Type
useParametric = (p_ctrl_norm > 0.05 && p_cbx_norm > 0.05);

% 非参数检验（Wilcoxon） / Non-parametric: ranksum + Cliff's delta
[p_ranksum, ~] = ranksum(ctrl, cbx);
d = cliffs_delta(ctrl, cbx);
report.ranksum.p = p_ranksum;
report.cliffs_delta.value = d;
report.cliffs_delta.magnitude = interpret_cliffs_delta(d);

%% t检验 + Cohen's d
if useParametric
    [~, p_t] = ttest2(ctrl, cbx, 'Vartype','unequal');
    report.ttest.p = p_t;
    report.cohens_d.value = compute_cohens_d(ctrl, cbx);
end

% Welch ANOVA + Brown-Forsythe
group = [ones(length(ctrl),1); 2*ones(length(cbx),1)];
[p_anova, tbl, stats] = anova1([ctrl(:); cbx(:)], group, 'off');
report.welch_anova.p = p_anova;
report.welch_anova.df = tbl{2,3};
report.welch_anova.F = tbl{2,5};

end

function d = compute_cohens_d(a, b)
% 计算Cohen's d
na = length(a);
nb = length(b);
sd_pooled = sqrt(((na-1)*var(a) + (nb-1)*var(b)) / (na+nb-2));
d = (mean(a) - mean(b)) / sd_pooled;
end

function d = cliffs_delta(a, b)
% 计算Cliff's delta
n = 0;
s = 0;
for i = 1:length(a)
    for j = 1:length(b)
        n = n + 1;
        if a(i) > b(j)
            s = s + 1;
        elseif a(i) < b(j)
            s = s - 1;
        end
    end
end
d = s / n;
end

function label = interpret_cliffs_delta(d)
% 解释Cliff's delta效应量 / Interpret magnitude
abs_d = abs(d);
if abs_d < 0.147
    label = 'negligible';
elif abs_d < 0.33
    label = 'small';
elif abs_d < 0.474
    label = 'medium';
else
    label = 'large';
end
end