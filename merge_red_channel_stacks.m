function merge_red_channel_stacks()
% MERGE_RED_CHANNEL_STACKS merges selected TIFF stacks' red channels into a single stack.

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


% Select input folder
input_folder = uigetdir('choose all the stacks');
if input_folder == 0
    disp('No input folder selected.');
    return;
end

% Get all .tif files in the folder
files = dir(fullfile(input_folder, '*.tif'));
if isempty(files)
    disp('No TIFF stack files found in the selected folder.');
    return;
end

file_names = {files.name};

% Let user choose which files to merge
[selected_indices, ok] = listdlg('ListString', file_names, ...
    'SelectionMode', 'multiple', ...
    'PromptString', 'Select TIFF stack files to merge:', ...
    'ListSize', [300, 300]);

if ~ok || isempty(selected_indices)
    disp('No files selected for merging.');
    return;
end

% Merge red channels
im_data = [];
count = 0;

for i = 1:numel(selected_indices)
    file_path = fullfile(input_folder, file_names{selected_indices(i)});
    stack_info = imfinfo(file_path);
    num_images = numel(stack_info);

    for j = 1:num_images
        count = count + 1;
        original_image = imread(file_path, j, 'Info', stack_info);
        im_data = cat(3, im_data, original_image(:,:,1)); % Red channel
    end
end

% Select output folder
output_folder = uigetdir('Choose a folder to save the merged stack');
if output_folder == 0
    disp('No output folder selected.');
    return;
end

% Construct output file path
[~, input_folder_name, ~] = fileparts(input_folder);
output_filename = fullfile(output_folder, [input_folder_name '_merged_red.tif']);

% Write merged stack
for n = 1:size(im_data, 3)
    if n == 1
        imwrite(im_data(:,:,n), output_filename);
    else
        imwrite(im_data(:,:,n), output_filename, 'WriteMode', 'append');
    end
end

disp(['Merged stack saved to: ' output_filename]);
disp(['Total number of images in the merged stack: ' num2str(count)]);
end
