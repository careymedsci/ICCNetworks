
% 作者：刘晓； 法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Xiao Liu
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany
% Xiangyang No.1 people's Hospital, China
% Also an python version coded 
% Xiao.Liu@stud.uni-frankfurt.de
% 2024.10.09


%%


function params = load200ROIData(params)
    % Load ROI data and related images from files
    % Input: params - structure containing analysis parameters; must include 'mode'
    % Output: params - updated structure with loaded data

    % Display confirmation dialog for parameter setup
    msg = 'Parameters all set? Then Press to begin';
    h = msgbox(msg, 'Warning');
    uiwait(h);

    % Select and read the intensity file
    [filename, folder] = uigetfile('*.*', ['Hi! ################################# ', ...
        'Please choose the file that contains all the ROIs intensities information ', ...
        'the folder: * _roi_centre/ *  _mean_values.csv   ############################']);

    % Check if file selection was canceled
    if isequal(filename,0) || isequal(folder,0)
        disp('cancel selection');
        error('File selection cancelled');
    else
        currentFolder = fileparts(fullfile(folder, filename));
        clear folder
    end

    % Store the current folder path in params
    params.currentFolder = currentFolder;
    cd(currentFolder);

    % Set plot_networkx based on mode
    if params.mode > 0 
        params.plot_networkx = 0;
    end

    % Read and process intensity data
    allTraces = readmatrix(filename); 
    params.allTraces = allTraces(102:end,2:end);

    % Process output filename based on mode
    outputname = filename(1:length(filename)-16);
    if params.mode == 0
        output = outputname;
    elseif params.mode == 1
        output = outputname + "_circular_shift";
    elseif params.mode == 2
        output = outputname + "_scrambled";
    elseif params.mode == 3
        output = outputname + "_linear_shift";
    elseif params.mode == 4
        output = outputname + "_distant_cells";
    end

    % Store output-related information in params
    params.output = output;
    params.outputname = outputname;

    % Read ROI coordinate data
    coordinatesData = readmatrix(outputname + "_coordintes.txt");
    params.coordinatesData = coordinatesData(102,2:end);

    % Read max intensity projection image
    params.maxZImages = imread(outputname + "_max_z_proj.tif");

    % Clear temporary variables
    clear h msg 

    % Compute and store basic parameters
    params.numFrames = size(params.allTraces,1);
    params.numRois = size(params.allTraces,2);

    % Verify that critical data were loaded correctly
    if ~isfield(params, 'currentFolder') || isempty(params.currentFolder)
        error('Failed to store current folder path');
    end
    if ~isfield(params, 'allTraces') || isempty(params.allTraces)
        error('Failed to load trace data');
    end

    % Display confirmation messages
    fprintf('Data loaded successfully.\n');
    fprintf('Current folder set to: %s\n', params.currentFolder);
end


