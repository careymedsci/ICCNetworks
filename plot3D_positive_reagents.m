
 
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

function plot3D_positive_reagents(params, savePath)
% plotRandom3DTraces - Plot a 3D view of randomly selected cell traces
%
% INPUTS:
%   params   - Struct containing:
%                .allTraces       [time x cells] matrix of fluorescence or signal
%                (optional) .peakLocations [time x cells] binary matrix of peak positions
%   nums     - Number of traces to highlight (e.g., 150)
%   savePath - Path to save the resulting .fig file
%
% OUTPUT:
%   fig - Handle to the created figure

if nargin < 3
    savePath = pwd; % Default to current directory if not provided
end

num_cells = size(params.allTraces, 2);
nums = num_cells;
num_frames = size(params.allTraces, 1);

% Select traces to highlight
if isfield(params, 'peakLocations')
    peak_binary = params.peakLocations;
    num_peaks_per_trace = sum(peak_binary, 1);
    [~, sorted_idx] = sort(num_peaks_per_trace, 'descend');

    multi_peak_idx  = sorted_idx(num_peaks_per_trace(sorted_idx) > 1);
    single_peak_idx = sorted_idx(num_peaks_per_trace(sorted_idx) == 1);
    no_peak_idx     = sorted_idx(num_peaks_per_trace(sorted_idx) == 0);

    selected_traces = multi_peak_idx;
    if length(selected_traces) < nums
        needed = nums - length(selected_traces);
        selected_traces = [selected_traces, single_peak_idx(1:min(needed, length(single_peak_idx)))];
    end
    if length(selected_traces) < nums
        needed = nums - length(selected_traces);
        selected_traces = [selected_traces, no_peak_idx(1:min(needed, length(no_peak_idx)))];
    end

    highlighted_cells = selected_traces(1:min(nums, length(selected_traces)));
else
    highlighted_cells = 1:min(nums, num_cells);
end

% Create time vector
total_duration = 10; % 10 minutes
time_vector = linspace(0, total_duration, num_frames);

% Create figure
fig = figure('Position', [200 200 1000 700], ...
             'Color', 'w', 'Name', '3D Cell Traces');
hold on

rng(52); % Fix random seed for reproducible x positions
x_positions = rand(1, num_cells) * num_cells;

% Plot each trace
for trace = 1:num_cells
    z_data = params.allTraces(:, trace);
    y_data = time_vector';
    x_data = x_positions(trace) * ones(size(z_data));

    if ismember(trace, highlighted_cells)
        color = [0.64 0.08 0.18];  % Dark red for highlighted cells
        z_scale = 10;               % Scale z-axis for visibility
        line_width = 1.5;
    else
        if rand < 0.30
            color = [1, 1, 1];     % White, almost invisible background traces
            z_scale = 1;
            line_width = 0.01;
        else
            color = [1, 1, 1];
            z_scale = 1;
            line_width = 0.1;
        end
    end

    plot3(x_data, y_data, z_data * z_scale, 'LineWidth', line_width, 'Color', color);
end

% Set view and formatting
view(60, 50)
axis off
set(gca, 'Visible', 'off')
set(gca, 'Color', 'none')
hold off
drawnow

% Save figure
if isfield(params, 'output')
    outputFigure = sprintf('%s_3D_traces_treatment.fig', params.output);
else
    outputFigure = '3D_traces_treatment.fig';
end
