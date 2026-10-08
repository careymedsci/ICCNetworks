
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


function result_one_anaylsis(im_data)
% 作者：Thomas Broggini 和 刘晓；法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Thomas Broggini and Xiao Liu
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

% RESULT_ONE_ANAYLSIS performs SVD-based analysis on a 3D image stack.
% The process includes decomposition, interactive component selection,
% temporal and spatial mode visualization, reconstruction, and energy 
% validation.
%
% Input:
%   im_data - 3D matrix of size [height, width, num_frames]
%             representing the image sequence

%% 1. Data Preparation & SVD Decomposition
% Input data dimensions: [height, width, num_frames]
[height, width, num_frames] = size(im_data);

% Reshape to space-time matrix [n_frames × n_pixels]
space_time_data = double(reshape(im_data, [], num_frames))';

% Perform economy-size SVD
[U, S, V] = svd(space_time_data, 'econ');
singular_values = diag(S);  % Correctly extract singular values

%% 2. Interactive Component Selection
% Create dual-scale scree plot
figure('Color','w','Position',[100 100 1000 500]);

% Linear scale plot
subplot(1,2,1);
plot(singular_values.^2/sum(singular_values.^2)*100, 'b-o', ...
    'MarkerSize',6,'LineWidth',1.5);
xlabel('Component Index');
ylabel('Energy Contribution (%)');
title('Scree Plot (Linear Scale)');
grid on;
xlim([1 min(100, length(singular_values))]);

% Log scale plot
subplot(1,2,2);
semilogy(singular_values.^2/sum(singular_values.^2)*100, 'r-o', ...
    'MarkerSize',6,'LineWidth',1.5);
xlabel('Component Index');
ylabel('Energy Contribution (%)');
title('Scree Plot (Log Scale)');
grid on;
xlim([1 min(100, length(singular_values))]);

% Interactive selection
[selected_x, ~] = ginput(1);
k = round(selected_x(1));
fprintf('Selected first %d components\n', k);

%% 3. Temporal Pattern Analysis
% Plot temporal modes
figure('Color','w','Name','Temporal Modes');
for comp = 5:min(5,k)
    subplot(5,1,comp);
    plot(U(:,comp), 'LineWidth',1.5);
    xlabel('Time (frames)');
    ylabel('Amplitude');
    title(sprintf('Temporal Mode %d', comp));
    xlim([1 num_frames]);
end

%% 4. Spatial Pattern Visualization
% Visualize spatial modes
figure('Color','w','Name','Spatial Modes');
for comp = 1:min(9,k)
    subplot(3,3,comp);
    imagesc(reshape(V(:,comp), height, width));
    axis image off;
    title(sprintf('Spatial Mode %d', comp));
    colormap(jet);
    colorbar;
end

%% 5. Data Reconstruction & Validation
% Reconstruct using selected components
recon_data = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';

% Calculate reconstruction error
recon_error = norm(space_time_data - recon_data, 'fro') ...
    / norm(space_time_data, 'fro') * 100;
fprintf('Reconstruction error: %.2f%%\n', recon_error);

% Visualize sample frame
sample_frame = 50; % Arbitrary frame selection
figure('Color','w','Name','Reconstruction Quality');
subplot(1,2,1);
imagesc(reshape(space_time_data(sample_frame,:), height, width));
title('Original Frame');
axis image off;

subplot(1,2,2);
imagesc(reshape(recon_data(sample_frame,:), height, width));
title(sprintf('Reconstructed (k=%d)', k));
axis image off;

%% 6. Energy Capture Validation
energy_ratio = singular_values.^2 / sum(singular_values.^2);  
cum_energy = cumsum(energy_ratio);  

figure;
plot(1:k, cum_energy(1:k)*100, 'bo-', 'LineWidth', 2);
hold on;
plot(k, cum_energy(k)*100, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Mode index');
ylabel('Cumulative energy (%)');
title(sprintf('Cumulative energy contribution of first %d modes', k));
grid on;
xlim([1 k]);
text(k, cum_energy(k)*100, sprintf('%.1f%%', cum_energy(k)*100), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

%% Plot Un
figure;
imagesc(U(:,1:k)');  
colormap jet;
colorbar;
xlabel('Time frame');
ylabel('Mode index');
title('Temporal components (U_n)');
axis tight;

end
