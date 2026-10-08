

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
% 2024.10.09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function KS_CDFs_test(params)

    % PLOTKS_ECDF_COMPARISON - Plot empirical CDFs of correlation values across real and surrogate modes.
    %
    % This function compares the distribution of pairwise correlation coefficients
    % between real data and three control models (mode 1~3), and visualizes the
    % empirical CDFs along with KS statistics.
    %
    % Input:
    %   params - structure containing outputDiary and base path info
    %
    % Author: ChatGPT + User Collaboration
    % Date: 2025

    figure('Color', 'w', 'Position', [100, 100, 900, 600]);

    % Load correlation matrices
    % baseName = extractBefore(params.outputDiary, strlength(params.outputDiary) - 23);
   baseName = params.outputname;
    C_real  = load(sprintf('%s_C0.txt', baseName));
    C_mode1 = load(sprintf('%s_circular_shift_C1.txt', baseName));
    C_mode2 = load(sprintf('%s_scrambled_C2.txt', baseName));
    C_mode3 = load(sprintf('%s_linear_shift_C3.txt', baseName));

    % Extract upper triangle (excluding diagonal)
    get_upper = @(C) C(triu(true(size(C)), 1));
    real_corr = double(get_upper(C_real));
    mode1_corr = double(get_upper(C_mode1));
    mode2_corr = get_upper(C_mode2);
    mode3_corr = get_upper(C_mode3);

    % Compute ECDFs
    [f_real, x_real] = ecdf(real_corr);
    [f_m1, x_m1] = ecdf(mode1_corr);
    [f_m2, x_m2] = ecdf(mode2_corr);
    [f_m3, x_m3] = ecdf(mode3_corr);

    % Interpolate ECDFs to x_real for KS comparison
    f_m1_interp = interp_ecdf_safe(x_real, x_m1, f_m1);
    f_m2_interp = interp_ecdf_safe(x_real, x_m2, f_m2);
    f_m3_interp = interp_ecdf_safe(x_real, x_m3, f_m3);

    % KS statistics
    [ks1, idx1] = max(abs(f_real - f_m1_interp));
    [ks2, idx2] = max(abs(f_real - f_m2_interp));
    [ks3, idx3] = max(abs(f_real - f_m3_interp));

    % Corresponding points
    x_d1 = x_real(idx1); y1 = f_real(idx1); y2 = f_m1_interp(idx1);
    x_d2 = x_real(idx2); y3 = f_real(idx2); y4 = f_m2_interp(idx2);
    x_d3 = x_real(idx3); y5 = f_real(idx3); y6 = f_m3_interp(idx3);

    % Plot ECDFs
    hold on;
    plot(x_real, f_real, 'k-', 'LineWidth', 2.5);
    plot(x_m1, f_m1, 'r--', 'LineWidth', 2);
    plot(x_m2, f_m2, 'b--', 'LineWidth', 2);
    plot(x_m3, f_m3, 'g--', 'LineWidth', 2);

    % KS distance lines
    plot([x_d1 x_d1], [y1 y2], 'r:', 'LineWidth', 1.5);
    plot([x_d2 x_d2], [y3 y4], 'b:', 'LineWidth', 1.5);
    plot([x_d3 x_d3], [y5 y6], 'g:', 'LineWidth', 1.5);

    % Fill area between curves (optional visualization)
    fill([x_real(idx1:end); flipud(x_real(idx1:end))], ...
         [f_real(idx1:end); flipud(f_m1_interp(idx1:end))], ...
         'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    % KS tests
    [~, p1] = kstest2(real_corr, mode1_corr);
    [~, p2] = kstest2(real_corr, mode2_corr);
    [~, p3] = kstest2(real_corr, mode3_corr);

    % Legend and annotation
    legend({
        'Real', ...
        'Mode 1 (circular shift)', ...
        'Mode 2 (scrambled)', ...
        'Mode 3 (linear shift)', ...
        sprintf('KS D1 = %.3f, p = %.1e', ks1, p1), ...
        sprintf('KS D2 = %.3f, p = %.1e', ks2, p2), ...
        sprintf('KS D3 = %.3f, p = %.1e', ks3, p3)
    }, 'Location', 'southeast', 'FontSize', 10, 'Box', 'off');

    xlabel('Correlation Coefficient', 'FontSize', 12, 'FontName', 'Arial');
    ylabel('Cumulative Probability', 'FontSize', 12, 'FontName', 'Arial');
    title('Empirical CDFs of Pairwise Correlations', 'FontSize', 14, 'FontName', 'Arial');
    set(gca, 'FontSize', 12, 'Box', 'off', 'LineWidth', 1.5);
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridColor', [0.7 0.7 0.7]);

    % Annotate KS points
    text(x_d1, (y1 + y2)/2, sprintf('D1 = %.2f', ks1), 'Color', 'r', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(x_d2, (y3 + y4)/2, sprintf('D2 = %.2f', ks2), 'Color', 'b', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(x_d3, (y5 + y6)/2, sprintf('D3 = %.2f', ks3), 'Color', 'g', 'HorizontalAlignment', 'left', 'FontSize', 10);

end

% Helper function
function f_interp = interp_ecdf_safe(xq, x, f)
    [x_unique, idx] = unique(x);
    f_unique = f(idx);
    f_interp = interp1(x_unique, f_unique, xq, 'linear', 'extrap');
end
