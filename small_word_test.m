
 
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


function small_word_test()

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
% 2025.07.06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ==== 文件选择 / Select Excel File ====
[filename, pathname] = uigetfile({'*.xls;*.xlsx'}, '选择 Excel 文件 / Select Excel file');
if isequal(filename, 0)
    error('取消选择 / selection cancelled');
end
filepath = fullfile(pathname, filename);

% ======= Sheet selection / 选择工作表名称 =======
sheet_options = {'cbx', 'gap26', 'suramin', 'apyrase', 'EGTA', 'Vera','skf','dan','AP5','B16-F10.N3','B16-BrM.3',...
    'NieCl2', 'fgf', 'tram34'};
[selection, ok] = listdlg('PromptString', 'Select sheet name:', ...
                          'SelectionMode', 'single', ...
                          'ListString', sheet_options);
if ~ok
    disp('Sheet selection cancelled.');
    return;
end
sheetname = sheet_options{selection};
fprintf('Selected sheet: %s\n', sheetname);

% ==== 读取数据 / Read Data ====
emp_lambda = readmatrix(filepath, 'Sheet', sheetname, 'Range', 'A33:R33'); % 一次读到18列
rnd_lambda = readmatrix(filepath, 'Sheet', sheetname, 'Range', 'A34:R34');
emp_sigma  = readmatrix(filepath, 'Sheet', sheetname, 'Range', 'A37:R37');
rnd_sigma  = readmatrix(filepath, 'Sheet', sheetname, 'Range', 'A38:R38');

% ==== 原始分析 (1–12列) ====
do_analysis(emp_lambda(1:12), rnd_lambda(1:12), emp_sigma(1:12), rnd_sigma(1:12), sheetname, 'Cols 1–12');

% ==== 如果有13–18列，再做一次分析 (1–6 + 13–18) ====
if any(~isnan(emp_lambda(13:18))) || any(~isnan(rnd_lambda(13:18))) || ...
   any(~isnan(emp_sigma(13:18)))  || any(~isnan(rnd_sigma(13:18)))
   
   emp_lambda2 = [emp_lambda(1:6), emp_lambda(13:18)];
   rnd_lambda2 = [rnd_lambda(1:6), rnd_lambda(13:18)];
   emp_sigma2  = [emp_sigma(1:6),  emp_sigma(13:18)];
   rnd_sigma2  = [rnd_sigma(1:6),  rnd_sigma(13:18)];
   
   do_analysis(emp_lambda2, rnd_lambda2, emp_sigma2, rnd_sigma2, sheetname, 'Cols 1–6 + 13–18');
end

end % function small_word_test

%% ===== 内部函数，跑一次完整分析 =====
function do_analysis(emp_lambda, rnd_lambda, emp_sigma, rnd_sigma, sheetname, tag)

% ==== color setting ====
color_ctrl = [0, 113, 188] / 255;  % 蓝色
color_cbx  = [163, 30, 50] / 255;  % 红色

% x axis
x_lambda = 1; 
x_sigma  = 2;
jitter = 0.01;

figure('Color', 'w'); hold on;
sgtitle(sprintf('%s - %s', sheetname, tag))

% ===== λ - left y axis (配对散点+连线) =====
yyaxis left
for i = 1:numel(emp_lambda)
    if isnan(emp_lambda(i)) || isnan(rnd_lambda(i)), continue; end
    if i <= 6
        clr_emp = color_ctrl; 
        clr_rnd = color_ctrl;
    else
        clr_emp = color_cbx; 
        clr_rnd = color_cbx;
    end
    plot([x_lambda-0.1, x_lambda+0.1], [emp_lambda(i), rnd_lambda(i)], '-', ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
    scatter(x_lambda-0.1, emp_lambda(i), 40, clr_emp, 'filled');
    scatter(x_lambda+0.1, rnd_lambda(i), 40, clr_rnd, 'filled');
end
mean_lambda = [nanmean(emp_lambda), nanmean(rnd_lambda)];
sem_lambda  = [nanstd(emp_lambda)/sqrt(sum(~isnan(emp_lambda))), ...
               nanstd(rnd_lambda)/sqrt(sum(~isnan(rnd_lambda)))];
errorbar([x_lambda-0.1, x_lambda+0.1], mean_lambda, sem_lambda, ...
    'k', 'LineStyle', 'none', 'LineWidth', 2);
ylabel('\lambda (Mean Shortest Path Length)', 'FontSize', 12)
ylim([0, max([emp_lambda, rnd_lambda], [], 'omitnan') * 1.2]);

% ===== σ - right y axis (原散点抖动) =====
yyaxis right
scatter_jittered(x_sigma, emp_sigma, color_ctrl, color_cbx, jitter);
scatter_jittered(x_sigma+0.4, rnd_sigma, color_ctrl, color_cbx, jitter);
mean_sigma = [nanmean(emp_sigma), nanmean(rnd_sigma)];
sem_sigma  = [nanstd(emp_sigma)/sqrt(sum(~isnan(emp_sigma))), ...
              nanstd(rnd_sigma)/sqrt(sum(~isnan(rnd_sigma)))];
errorbar([x_sigma, x_sigma+0.4], mean_sigma, sem_sigma, ...
    'k', 'LineStyle', 'none', 'LineWidth', 2);
ylabel('\sigma (Clustering Coefficient)', 'FontSize', 12)
ylim([0, max([emp_sigma, rnd_sigma], [], 'omitnan') * 1.2]);

% ==== significance test (配对t检验) ====
yyaxis left
valid_idx = ~isnan(emp_lambda) & ~isnan(rnd_lambda);
[~, p_lambda] = ttest(emp_lambda(valid_idx), rnd_lambda(valid_idx));
y1 = max([emp_lambda, rnd_lambda], [], 'omitnan') * 1.15;
text(x_lambda, y1, significance_text(p_lambda), 'HorizontalAlignment', 'center', 'FontSize', 10);
plot([x_lambda-0.1, x_lambda+0.1], [y1, y1], 'k', 'LineWidth', 1.2)

yyaxis right
valid_idx = ~isnan(emp_sigma) & ~isnan(rnd_sigma);
[~, p_sigma] = ttest(emp_sigma(valid_idx), rnd_sigma(valid_idx));
y2 = max([emp_sigma, rnd_sigma], [], 'omitnan') * 1.15;
text(x_sigma, y2, significance_text(p_sigma), 'HorizontalAlignment', 'center', 'FontSize', 10);
plot([x_sigma, x_sigma+0.4], [y2, y2], 'k', 'LineWidth', 1.2)

% ==== coordinates setting ====
xlim([0.5 2.8])
xticks([1.1 2.2])
xticklabels({'\lambda', '\sigma'})
set(gca, 'FontSize', 12)
box off

% ==== small-worldness ====
L_empirical = nanmean(emp_lambda);
L_random    = nanmean(rnd_lambda);
C_empirical = nanmean(emp_sigma);
C_random    = nanmean(rnd_sigma);
if any(isnan([L_empirical, L_random, C_empirical, C_random]))
    warning('存在 NaN 值，无法计算小世界性 S / NaNs detected. Cannot compute S.');
    S = NaN;
else
    S = (C_empirical / C_random) / (L_empirical / L_random);
    fprintf('[%s] Small-worldness S = %.4f\n', tag, S);
end
str_S = sprintf('S = %.2f', S);
annotation('textbox', [0.46, 0.4, 0.1, 0.2], ...
    'String', str_S, 'FitBoxToText', 'on', 'EdgeColor', 'none', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'Rotation', 90, 'HorizontalAlignment', 'center');

end % do_analysis

%% supporting functions
function scatter_jittered(x_center, y_vals, clr_ctrl, clr_cbx, jitter)
    for i = 1:numel(y_vals)
        if isnan(y_vals(i)), continue; end
        if i <= 6
            clr = clr_ctrl;
        else
            clr = clr_cbx;
        end
        x_jittered = x_center + (rand - 0.5)*2*jitter;
        scatter(x_jittered, y_vals(i), 40, clr, 'filled');
    end
end

function txt = significance_text(p)
    if p < 0.0001
        txt = 'P < 0.0001';
    elseif p < 0.001
        txt = sprintf('P = %.4f', p);
    elseif p < 0.05
        txt = sprintf('P = %.3f', p);
    else
        txt = 'ns';
    end
end
