function plot3DTraces(params)

    savePath = params.currentFolder;
    load('all_ver_0.mat'); 
    nums = 150;

    peak_binary = params.peakLocations;
    num_peaks_per_trace = sum(peak_binary, 1);
    [~, sorted_idx] = sort(num_peaks_per_trace, 'descend');

    multi_peak_idx  = sorted_idx(num_peaks_per_trace(sorted_idx) > 1);
    single_peak_idx = sorted_idx(num_peaks_per_trace(sorted_idx) == 1);
    no_peak_idx     = sorted_idx(num_peaks_per_trace(sorted_idx) == 0);

    selected_traces = [multi_peak_idx];
    if length(selected_traces) < nums
        needed = nums - length(selected_traces);
        selected_traces = [selected_traces, single_peak_idx(1:min(needed, length(single_peak_idx)))];
    end
    if length(selected_traces) < nums
        needed = nums - length(selected_traces);
        selected_traces = [selected_traces, no_peak_idx(1:min(needed, length(no_peak_idx)))];
    end

    highlighted_cells = selected_traces(1:min(nums, length(selected_traces)));

    total_duration = 10;
    num_frames = size(params.allTraces,1);
    time_vector = linspace(0, total_duration, num_frames);


    fig = figure('Position', [200 200 800 600], ...
                 'Color', 'w', ...
                 'Name', '3D Cell Traces');
    hold on

    num_cells = size(params.allTraces, 2);

    for trace = 1:num_cells
        z_data_full = params.allTraces(:, trace);
        z_data = z_data_full;
        y_data = time_vector';
        x_data = trace * ones(size(z_data));

        if ismember(trace, highlighted_cells)
            color = [0.64 0.08 0.18];
            z_scale = 10.2;
            line_width = 1.5;
        else
            if rand < 0.30
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

    view(60, 70)

    
    axis off
    set(gca, 'Visible', 'off')
    set(gca, 'Color', 'none')


    hold off
    drawnow

    
    outputFigure = sprintf('%s_3D_traces.fig', params.output);
    savefig(fullfile(savePath, outputFigure));
end
