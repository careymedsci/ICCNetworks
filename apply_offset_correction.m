function im_data = apply_offset_correction(im_data)
% APPLY_OFFSET_CORRECTION adjusts the first 100 frames based on intensity offset
% between two user-selected frames.

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


% Input:
%   im_data - 3D image data matrix (Height x Width x Time)
%
% Output:
%   im_data - Offset-corrected 3D image data

    % Ask user for reference frames
    answer1 = inputdlg('Which is the frame number before treatment?', 'Offset Correction');
    if isempty(answer1)
        disp('User canceled input.');
        return;
    end
    Alama1 = str2double(answer1{1});
    
    answer2 = inputdlg('Which is the frame number after treatment?', 'Offset Correction');
    if isempty(answer2)
        disp('User canceled input.');
        return;
    end
    Alama2 = str2double(answer2{1});
    
    % Validate frame numbers
    total_frames = size(im_data, 3);
    if any([Alama1, Alama2] < 1) || any([Alama1, Alama2] > total_frames)
        error('Frame numbers must be within the range of 1 to %d.', total_frames);
    end

    % Compute intensity difference between two frames
    mean_frame_1 = mean(im_data(:,:,Alama1), 'all');
    mean_frame_2 = mean(im_data(:,:,Alama2), 'all');
    mean_difference = mean_frame_2 - mean_frame_1;

    % Apply correction to first 100 frames
    num_corrected_frames = min(100, total_frames);
    im_data(:,:,1:num_corrected_frames) = im_data(:,:,1:num_corrected_frames) + mean_difference;

    disp(['Applied offset correction of ' num2str(mean_difference) ...
        ' to the first ' num2str(num_corrected_frames) ' frames.']);
end
