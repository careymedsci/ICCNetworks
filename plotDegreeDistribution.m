
 
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

function plotDegreeDistribution(params)
    % This function takes the input 'params' structure containing the
    % co-activation data and degree distribution. It performs data cleaning, 
    % fitting of Power Law, Poisson, and Lognormal distributions and plots the 
    % results.

    % Extract co-activation data
    data = params.CvsD;
    C = data(:,1);
    D = data(:,2);

    % Remove cell pairs without co-activation (C > 0 and D > 0)
    valid_indices = (C > 0) & (D > 0);
    C = C(valid_indices);
    D = D(valid_indices);

    % Number of nodes
    N = size(params.D, 1);

    % Generate ER random graph
    p = 0.01;  % Probability of edge creation
    ER_adjacency = rand(N) < p;
    ER_adjacency = triu(ER_adjacency, 1) + triu(ER_adjacency, 1)';  % Symmetric matrix

    % Compute ER network degree distribution and remove nodes with degree 0
    ER_degrees = sum(ER_adjacency, 2);
    ER_degrees = ER_degrees(ER_degrees > 0);

    % Construct the adjacency matrix of the real network from the correlation matrix
    threshold = params.threshcorr;
    adjacency = params.C >= threshold;

    % Compute real network node degrees and remove nodes with degree 0
    degrees = sum(adjacency, 2);

    % Plot degree distributions (linear scale)
    figure;
    set(gcf, 'Color', 'w');
    subplot(2,2,1);
    hist_real = histogram(degrees, 'Normalization', 'probability');
    hold on;
    xlabel('Degree');
    ylabel('Probability');
    title('Degree Distribution: Real Network');
    grid on;

    subplot(2,2,2);
    hist_er = histogram(ER_degrees, 'Normalization', 'probability');
    hold on;
    xlabel('Degree');
    ylabel('Probability');
    title('Degree Distribution: ER Network');
    grid on;

    % Power-law fitting for real network: log(P(k)) = a * log(k) + b
    figure;
    subplot(2,2,3);
    hist_real_log = histogram(degrees, 'Normalization', 'probability');
    hold on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Degree (log scale)');
    ylabel('Probability (log scale)');
    title('Real Network: Degree Distribution (log-log)');
    grid on;

    deg_vals = hist_real_log.BinEdges(1:end-1);
    prob_vals = hist_real_log.Values;
    valid = (deg_vals > 0) & (prob_vals > 0);
    coeffs_real = polyfit(log(deg_vals(valid)), log(prob_vals(valid)), 1);
    fitted_line_real = exp(polyval(coeffs_real, log(deg_vals(valid))));
    plot(deg_vals(valid), fitted_line_real, 'r--', 'LineWidth', 2);

    % ER network degree distribution (log-log scale)
    subplot(2,2,4);
    hist_er_log = histogram(ER_degrees, 'Normalization', 'probability');
    hold on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Degree (log scale)');
    ylabel('Probability (log scale)');
    title('ER Network: Degree Distribution (log-log)');
    grid on;

    % Add normal distribution fit to ER network (green dashed line)
    mu_er = mean(ER_degrees);   
    sigma_er = std(ER_degrees);
    x_vals_er = linspace(min(ER_degrees), max(ER_degrees), 100);
    normal_fit_er = normpdf(x_vals_er, mu_er, sigma_er);  
    plot(x_vals_er, normal_fit_er, 'g--', 'LineWidth', 2);

    % Update legend with both fits
    legend('Real Network (Empirical)', ...
           sprintf('Power-law fit (Real, slope = %.2f)', coeffs_real(1)), ...
           'ER Network (Empirical)', ...
           'Normal fit (ER)');

    % Disable scientific notation
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;

end
