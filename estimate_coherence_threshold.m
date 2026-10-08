
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

function coh_threshold = estimate_coherence_threshold(Vn, scores, tapers, Fs, nfft, f, findx, toplot_f_vector, n_iter, alpha, block_size)

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

% Estimate coherence threshold using bootstrap null distribution
%
% Inputs:
%   Vn                - data matrix (time x mode/channel)
%   scores            - spatial mode weights (e.g., SVD left/right singular vectors)
%   tapers            - tapers for multitaper method
%   Fs                - sampling frequency
%   nfft              - number of FFT points
%   f                 - frequency vector from multitaper FFT
%   findx             - frequency indices of interest
%   toplot_f_vector   - frequency points for final coherence curve (sparse for speed)
%   n_iter            - number of bootstrap iterations (e.g., 500)
%   alpha             - significance level (e.g., 0.05)
%   block_size        - block size for memory-efficient parallel bootstrap
%
% Output:
%   coh_threshold     - threshold of coherence under null hypothesis

    if nargin < 11
        block_size = 10;
    end

    null_coherence = zeros(n_iter, 1);
    f_vector_sparse = toplot_f_vector(1:4:end);  % downsample for speed

    total_start_time = tic;
    for block_start = 1:block_size:n_iter
        current_block_size = min(block_size, n_iter - block_start + 1);
        block_results = zeros(current_block_size, 1);

        parfor kk = 1:current_block_size
            ii = block_start + kk - 1;

            Vn_shuffled = Vn(randperm(size(Vn,1)), :);
            max_coh = 0;

            for ff = 1:length(f_vector_sparse)
                interp_FFT = zeros(size(scores,2), size(Vn_shuffled,2));
                for ch = 1:size(Vn_shuffled,2)
                    J = mtfftc(Vn_shuffled(:,ch), tapers, nfft, Fs);
                    J = J(findx,:);
                    J_mean = mean(J, 2);
                    freq_interp = interp1(f, J_mean, f_vector_sparse(ff), 'linear');
                    interp_FFT(:,ch) = freq_interp(:);  % ensure column
                end

                m = scores * interp_FFT';
                s = svd(m, 0);
                current_coh = s(1)^2 / sum(s.^2);
                max_coh = max(max_coh, current_coh);
            end

            block_results(kk) = max_coh;
        end

        null_coherence(block_start:block_start+current_block_size-1) = block_results;

        elapsed = toc(total_start_time);
        done_iter = block_start + current_block_size - 1;
        progress = done_iter / n_iter * 100;
        remaining = elapsed / done_iter * (n_iter - done_iter);
        fprintf('Progress: %.1f%%, Elapsed: %.1fs, Remaining: %.1fs\n', progress, elapsed, remaining);
    end

    % Calculate threshold
    coh_threshold = prctile(null_coherence, 100*(1 - alpha));
    fprintf('\nCoherence threshold (α=%.2f): %.4f\n', alpha, coh_threshold);

    delete(gcp('nocreate'));
end
