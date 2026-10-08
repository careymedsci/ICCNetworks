function dcutoff_final = compute_local_dcutoff_with_dps(rparams, params)
% ============================================
% compute_local_dcutoff_with_dps
% ============================================
% Purpose:
% Determine the local connectivity distance threshold (dcutoff) based on
% real data, shuffled data, and distance-preserving shuffled data, in order
% to exclude spurious long-range correlations.
%
% 使用真实数据、完全打乱数据和距离保持打乱数据确定局部连接阈值 dcutoff
% 目的是排除长距离随机耦合造成的伪相关
%
% Inputs:
%   rparams.CvsD - [correlation, distance] for real data
%   params.CvsD  - [correlation, distance] for fully shuffled/random data
%
% Outputs:
%   dcutoff_final - final selected distance threshold (μm)
%
% ============================================

% Parameters
threshcorr = params.threshcorr;  % correlation threshold
d_vals = 20:1:200;               % candidate distance thresholds
k = 2;                            % minimum ratio threshold

% Preallocate arrays
conn_real = zeros(size(d_vals));
conn_shuf = zeros(size(d_vals));
conn_dps = zeros(size(d_vals));

% ========== 1. Distance-preserving shuffle (DPS) ==========
% Keep distances, shuffle correlations
corr_vals = rparams.CvsD(:,1);
d_vals_all = rparams.CvsD(:,2);
rng('shuffle');  % random seed
corr_dps = corr_vals(randperm(length(corr_vals)));
rparams_dps.CvsD = [corr_dps, d_vals_all];

% ========== 2. Compute connection densities ==========
for i = 1:length(d_vals)
    idx_real = rparams.CvsD(:,2) <= d_vals(i);
    idx_shuf = params.CvsD(:,2) <= d_vals(i);
    idx_dps  = rparams_dps.CvsD(:,2) <= d_vals(i);
    
    % Real data
    conn_real(i) = sum(rparams.CvsD(idx_real,1) >= threshcorr) / max(sum(idx_real),1);
    
    % Fully shuffled
    conn_shuf(i) = sum(params.CvsD(idx_shuf,1) >= threshcorr) / max(sum(idx_shuf),1);
    
    % Distance-preserving shuffled
    conn_dps(i) = sum(rparams_dps.CvsD(idx_dps,1) >= threshcorr) / max(sum(idx_dps),1);
end

% 3. Compute ratios
ratio_shuf = conn_real ./ conn_shuf;
ratio_dps  = conn_real ./ conn_dps;

valid_shuf = isfinite(ratio_shuf) & ratio_shuf >= k;
valid_dps  = isfinite(ratio_dps) & ratio_dps  >= k;

if any(valid_shuf) && any(valid_dps)
    % Conservative: take min of two maximum valid distances
    dcutoff_shuf = max(d_vals(valid_shuf));
    dcutoff_dps  = max(d_vals(valid_dps));
    dcutoff = min(dcutoff_shuf, dcutoff_dps);
elseif any(valid_shuf)
    dcutoff = max(d_vals(valid_shuf));
elseif any(valid_dps)
    dcutoff = max(d_vals(valid_dps));
else
    dcutoff = d_vals(1);
    warning('No distance satisfies ratio threshold, using minimum distance %d μm', dcutoff);
end

% 4. Exponential and power-law fitting (real data)
valid_fit = isfinite(conn_real) & conn_real > 0;
d_fit = d_vals(valid_fit);
y_fit = conn_real(valid_fit);

if length(d_fit) >= 3
    ft = fittype('a*exp(-x/lambda)+c','independent','x');
    opts = fitoptions('Method','NonlinearLeastSquares');
    opts.StartPoint = [max(y_fit),50,min(y_fit)];
    [fit_exp, ~] = fit(d_fit(:),y_fit(:),ft,opts);
    lambda_exp = fit_exp.lambda;
else
    lambda_exp = NaN;
end

if length(d_fit) >= 3
    log_x = log(d_fit); log_y = log(y_fit);
    valid_power = y_fit>0 & isfinite(log_y);
    if sum(valid_power) >= 3
        p_poly = polyfit(log_x(valid_power), log_y(valid_power),1);
        b_power = -p_poly(1);
    else
        b_power = NaN;
    end
else
    b_power = NaN;
end

% 5. Plateau point and final dcutoff
diff_conn = diff(conn_real);
[~, idx_plateau] = min(abs(diff_conn));
plateau_dist = d_vals(idx_plateau);
dcutoff_final = min(dcutoff, plateau_dist);

% 6. Plotting
figure;
plot(d_vals, conn_real,'r','LineWidth',1.5); hold on;
plot(d_vals, conn_dps,'yo','LineWidth',1.5);
plot(d_vals, conn_shuf,'b','LineWidth',1.5); 

if exist('fit_exp','var') && ~isnan(lambda_exp)
    x_fit_curve = linspace(min(d_fit),max(d_fit),100);
    y_exp_fit = fit_exp.a * exp(-x_fit_curve/lambda_exp) + fit_exp.c;
    plot(x_fit_curve, y_exp_fit,'r','LineWidth',2);
end

if exist('b_power','var') && ~isnan(b_power)
    y_power_fit = exp(p_poly(2)) * x_fit_curve.^(-b_power);
    plot(x_fit_curve, y_power_fit,'g','LineWidth',2);
end

xline(dcutoff_final,'k',sprintf('dcutoff = %d μm',dcutoff_final));
xlabel('Distance threshold d_{cutoff} (μm)');
ylabel('Connection density');
legend({'Real','Shuffled','Distance-preserving shuffled','Exp fit','Power-law fit','Selected dcutoff'},'Location','best');
grid on;

% 7. Print summary
fprintf('dcutoff (final) = %d μm\n', dcutoff_final);
fprintf('Exponential decay length λ = %.2f μm\n', lambda_exp);
fprintf('Power-law exponent b = %.3f\n', b_power);
end