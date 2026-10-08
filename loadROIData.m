
%% 
% 作者：Thomas Broggini and 刘晓； 法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Thomas Broggini and Xiao Liu
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany
% Xiangyang No.1 people's Hospital, China
% Python version also exists
% Xiao.Liu@stud.uni-frankfurt.de
% 2024.10.09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = loadROIData(params)
%LOADROIDATA Load ROI intensity and coordinate data with robust format handling
%
% Input:
%   params - structure containing analysis parameters; must include 'mode'
%
% Output:
%   params - updated structure with loaded data and paths

    % --- Step 0: Confirm parameters ---
    msg = 'Parameters all set? Then press OK to begin';
    h = msgbox(msg, 'Confirmation');
    uiwait(h);

    % --- Step 1: Select ROI intensity file ---
    [filename, folder] = uigetfile({'*.*', 'All Files (*.*)'}, ...
        ['Hi! Please choose the file that contains all the ROI intensity information.', ...
         ' Typically _mean_values.csv or _mean_values.txt']);

    if isequal(filename,0) || isequal(folder,0)
        error('File selection cancelled by user.');
    end

    % Set working directory
    params.currentFolder = folder;
    cd(folder);

    % --- Step 2: Detect file type and read intensity data ---
    [~, baseName, ext] = fileparts(filename);
    ext = lower(ext);

    if ~ismember(ext, {'.csv', '.txt', '.tsv'})
        error('Unsupported file format: %s', ext);
    end

    % Read ROI intensity data
    try
        opts = detectImportOptions(filename);
        allTraces = readmatrix(filename, opts);
    catch ME
        error('Failed to read intensity file: %s\n%s', filename, ME.message);
    end

    % Store intensity data (skip header row/column)
    params.allTraces = allTraces(2:end,2:end);
    params.dataExtension = ext;

    % --- Step 3: Determine output name based on mode ---
    outputname = erase(baseName, "_mean_values"); % safer than string slicing
    switch params.mode
        case 0
            params.output = outputname;
        case 1
            params.output = outputname + "_circular_shift";
        case 2
            params.output = outputname + "_scrambled";
        case 3
            params.output = outputname + "_linear_shift";
        case 4
            params.output = outputname + "_distant_cells";
        otherwise
            warning('Unknown mode %d. Using default output name.', params.mode);
            params.output = outputname;
    end
    params.outputname = outputname;

    % --- Step 4: Load ROI coordinate file (auto-detect extension) ---
    coord_csv = outputname + "_coordintes.csv";
    coord_txt = outputname + "_coordintes.txt";

    if isfile(coord_csv)
        coord_file = coord_csv;
    elseif isfile(coord_txt)
        coord_file = coord_txt;
    else
        error('Coordinate file not found: expected %s or %s', coord_csv, coord_txt);
    end

    try
        opts = detectImportOptions(coord_file);
        coordinatesData = readmatrix(coord_file, opts);
        params.coordinatesData = coordinatesData(2,2:end);
    catch ME
        error('Failed to read coordinate file: %s\n%s', coord_file, ME.message);
    end

    % --- Step 5: Load max intensity projection image ---
    maxProjFile = outputname + "_max_z_proj.tif";
    if isfile(maxProjFile)
        params.maxZImages = imread(maxProjFile);
    else
        warning('Max Z projection image not found: %s', maxProjFile);
        params.maxZImages = [];
    end

    % --- Step 6: Basic parameters ---
    params.numFrames = size(params.allTraces,1);
    params.numRois   = size(params.allTraces,2);

    % --- Step 7: Adjust plotting option based on mode ---
    if params.mode > 0
        params.plot_networkx = 0;
    end

    % --- Step 8: Confirm data loaded successfully ---
    fprintf('ROI data loaded successfully.\n');
    fprintf('Current folder: %s\n', params.currentFolder);
    fprintf('Number of frames: %d\n', params.numFrames);
    fprintf('Number of ROIs: %d\n', params.numRois);
    fprintf('Data extension: %s\n', params.dataExtension);

end


%%

% 
% function params = loadROIData(params)
%     % Load ROI data and related images from files
%     % Input: params - structure containing analysis parameters; must include 'mode'
%     % Output: params - updated structure with loaded data
% 
%     % Display confirmation dialog for parameter setup
%     msg = 'Parameters all set? Then Press to begin';
%     h = msgbox(msg, 'Warning');
%     uiwait(h);
% 
%     % Select and read the intensity file
%     [filename, folder] = uigetfile('*.*', ['Hi! ################################# ', ...
%         'Please choose the file that contains all the ROIs intensities information ', ...
%         'the folder: * _roi_centre/ *  _mean_values.csv    ############################']);
% 
%     % Check if file selection was canceled
%     if isequal(filename,0) || isequal(folder,0)
%         disp('cancel selection');
%         error('File selection cancelled');
%     else
%         currentFolder = fileparts(fullfile(folder, filename));
%         clear folder
%     end
% 
%     % Store the current folder path in params
%     params.currentFolder = currentFolder;
%     cd(currentFolder);
% 
%     % Set plot_networkx based on mode
%     if params.mode > 0 
%         params.plot_networkx = 0;
%     end
% 
%     % Read and process intensity data
%     allTraces = readmatrix(filename); 
%     params.allTraces = allTraces(2:end,2:end);
% 
%     % Process output filename based on mode
%     outputname = filename(1:length(filename)-16);
%     if params.mode == 0
%         output = outputname;
%     elseif params.mode == 1
%         output = outputname + "_circular_shift";
%     elseif params.mode == 2
%         output = outputname + "_scrambled";
%     elseif params.mode == 3
%         output = outputname + "_linear_shift";
%     elseif params.mode == 4
%         output = outputname + "_distant_cells";
%     end
% 
%     % Store output-related information in params
%     params.output = output;
%     params.outputname = outputname;
% 
% 
%     % Read ROI coordinate data
%     coordinatesData = readmatrix(outputname + "_coordintes.csv");
%     params.coordinatesData = coordinatesData(2,2:end);
% 
%     % Read max intensity projection image
%     params.maxZImages = imread(outputname + "_max_z_proj.tif");
% 
%     % Clear temporary variables
%     clear h msg 
% 
%     % Compute and store basic parameters
%     params.numFrames = size(params.allTraces,1);
%     params.numRois = size(params.allTraces,2);
% 
%     % Verify that critical data were loaded correctly
%     if ~isfield(params, 'currentFolder') || isempty(params.currentFolder)
%         error('Failed to store current folder path');
%     end
%     if ~isfield(params, 'allTraces') || isempty(params.allTraces)
%         error('Failed to load trace data');
%     end
% 
%     % Display confirmation messages
%     fprintf('Data loaded successfully.\n');
%     fprintf('Current folder set to: %s\n', params.currentFolder);
% end
% 
% 
