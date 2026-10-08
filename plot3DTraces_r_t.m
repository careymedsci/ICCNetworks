
 
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

function plot3DTraces_r_t(params)
    savePath = params.currentFolder;
    load('all_ver_0.mat'); 
    nums = 200;
    
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

    total_duration = 10; % total time = 10 min
    num_frames = size(params.allTraces, 1);
    time_vector = linspace(0, total_duration, num_frames);  

    break_frame = 100;
    break_offset = 1; % add 1 min break indicateing the waiting time after addiation of drugs

    time_vector_shifted = time_vector;
    if break_frame < num_frames
        time_vector_shifted(break_frame+1:end) = time_vector_shifted(break_frame+1:end) + break_offset;
    end


   
    fig = figure('Position', [200 200 800 600], ...
                 'Color', 'w', ...
                 'Name', '3D Cell Traces');
    hold on

    num_cells = size(params.allTraces, 2);

    rng(42);  
    x_positions = rand(1, num_cells) * num_cells;

    for trace = 1:num_cells
        z_data_full = params.allTraces(:, trace);
        z_data = z_data_full;
        y_data = time_vector_shifted';
        x_data = x_positions(trace) * ones(size(z_data));

        if ismember(trace, highlighted_cells)
            color = [0.64 0.08 0.18];  
            z_scale = 35;
            line_width = 1.5;
        else
            if rand < 0.10
                color = [1, 1, 1];  
                z_scale = 1;
                line_width = 0.01;
            else
                color = [0.59 0.44 0.25];  
                z_scale = 1;
                line_width = 0.1;
            end
        end

        z_data_scaled = z_data * z_scale;

        plot3(x_data, y_data, z_data_scaled, ...
              'LineWidth', line_width, ...
              'Color', color);
    end

    plot3([0 num_cells], ...
          [time_vector(break_frame) + break_offset/2, time_vector(break_frame) + break_offset/2], ...
          [0 0], '--', 'Color', [0.6 0.6 0.6]);

    view(60, 50)
    axis off
    set(gca, 'Visible', 'off')
    set(gca, 'Color', 'none')

    hold off
    drawnow


    outputFigure = sprintf('%s_3D_traces_treatment.fig', params.output);
    savefig(fullfile(savePath, outputFigure));
end
