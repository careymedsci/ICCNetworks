function calculate_dfof_and_save(im_data, output_folder, baseline_frame_count, baseline_scale)
% CALCULATE_DFOF_AND_SAVE applies dF/F processing and saves the result as TIFF.


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



% Inputs:
%   im_data             - 3D image stack [height x width x time]
%   output_folder       - Path to save the result
%   baseline_frame_count - Number of frames to use as baseline (default: 50)
%   baseline_scale      - Scale factor for baseline (default: 0.8)
%
% Example:
%   calculate_dfof_and_save(im_data, 'C:\output\path')

    if nargin < 3
        baseline_frame_count = 50;
    end
    if nargin < 4
        baseline_scale = 0.8;
    end

    [height, width, num_frames] = size(im_data);
    
    % Reshape to 2D: [Time x Pixels]
    im_data_rs = reshape(im_data, [], num_frames)';
    waver = double(im_data_rs);

    % Calculate baseline using last N frames
    baseline = mean(waver(end - baseline_frame_count + 1:end, :), 1) * baseline_scale;

    % dF/F calculation
    wave = bsxfun(@rdivide, bsxfun(@minus, waver, baseline), baseline);
    rwave = reshape(wave', height, width, []);

    % Prepare output path
    [~, filename, ~] = fileparts(output_folder);
    dfof_folder = fullfile(output_folder, ['dfof_' filename]);
    if ~exist(dfof_folder, 'dir')
        mkdir(dfof_folder);
    end
    output_path = fullfile(dfof_folder, [filename '_dfof.tif']);

    % Save processed stack
    for n = 1:num_frames
        slice = rwave(:,:,n);
        if n == 1
            imwrite(slice, output_path);
        else
            imwrite(slice, output_path, 'WriteMode', 'append');
        end
    end

    disp(['dF/F stack saved to: ' output_path]);
end
