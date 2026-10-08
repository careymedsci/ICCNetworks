function [im_data, filename, output_folder, info] = load_tiff_stack()
% LOAD_TIFF_STACK prompts user to select a TIFF file and loads all its frames.

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

% Returns:
%   im_data        - 3D array of image data (height x width x frames)
%   filename       - name of the selected file
%   output_folder  - folder containing the selected file
%   info           - metadata from imfinfo

    msg = 'Please choose the .tiff data you saved to process.';
    h = msgbox(msg, 'Note:');
    uiwait(h);
    
    [filename, output_folder] = uigetfile('*.tif');
    if isequal(filename, 0)
        disp('User canceled file selection.');
        im_data = [];
        info = [];
        return;
    end

    % Get metadata
    file_path = fullfile(output_folder, filename);
    info = imfinfo(file_path);
    num_images = numel(info);

    % Get dimensions from first image
    first_frame = imread(file_path, 1);
    [height, width] = size(first_frame);

    % Initialize and load the full stack
    im_data = zeros(height, width, num_images, 'like', first_frame);
    for k = 1:num_images
        im_data(:,:,k) = imread(file_path, k);
    end

    disp(['Loaded ' num2str(num_images) ' frames from: ' filename]);
end
