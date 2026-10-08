
% 作者：刘晓； 法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Xiao Liu
% Modified from DirkCHoffmann： 
% https://github.com/DirkCHoffmann/Hausmann-et-al.-2022-Nature/blob/main/RUN_getcorr.m
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany
% Xiangyang No.1 people's Hospital, China
% Also an python version coded 
% Xiao.Liu@stud.uni-frankfurt.de
% 2024.11.11


%%


function [params, analysis_summary] = findTriggerCells(params, analysis_summary)
% findTriggerCells - Identifies and visualizes trigger cells in the network
%
% Inputs:
%   params - Structure containing analysis parameters and data
%   analysis_summary - Structure containing analysis results
%
% Outputs:
%   params - Updated parameters structure
%   analysis_summary - Updated analysis summary
%   Also creates and saves visualization figures of trigger cells

%% SECTION 1: Initialize variables
% Extract variables from params structure
cacul_corr = params.cacul_corr;
find_triggercells = params.find_triggercells;
S = params.S;
C = params.C;
threshcorr = params.threshcorr;
define_trigger_cells = params.define_trigger_cells;
coordinatesData = params.coordinatesData;
scaling = params.scaling;
maxZImages = params.maxZImages;

numbering = params.numbering;
output = params.output;
mode = params.mode;
numCells = params.numRois;
x = params.x;
y = params.y;

if strcmp(params.testback, 'yes') 
 
    [selectionIndex, ok] = listdlg( ...
    'PromptString', 'Please select the frequency you want to analysis:', ...
    'SelectionMode', 'single', ...
    'ListString', {'1', '2', '3', '4'}, ...
    'Name', 'Select Image Group');

% If the user cancels, exit the script
if ~ok
    disp('User canceled selection.');
    return;
end

% Get the selected group index
selection = selectionIndex;
% Load the corresponding image files based on the selection
switch selection
    case 1
        mean_im_phase = imread( "Phase_9.tif");
        mean_im_maxp = imread( "Max_power_distribution_1.tif");
        mean_im_mag   = imread( "Magnitude_5.tif");
                file_list = dir('Freq_distribution_*.tif');
        if ~isempty(file_list)
            file_name = file_list(1).name;
            max_freq_global = imread(file_name);
        end
    case 2
        mean_im_phase = imread( "Phase_10.tif");
        mean_im_maxp = imread("Max_power_distribution_2.tif");
        mean_im_mag   = imread( "Magnitude_6.tif");
                file_list = dir('Freq_distribution_*.tif');
        if ~isempty(file_list)
            file_name = file_list(1).name;
            max_freq_global = imread(file_name);
        end
    case 3
        mean_im_phase = imread( "Phase_11.tif");
        mean_im_maxp = imread( "Max_power_distribution_3.tif");
        mean_im_mag   = imread( "Magnitude_7.tif");
                file_list = dir('Freq_distribution_*.tif');
        if ~isempty(file_list)
            file_name = file_list(1).name;
            max_freq_global = imread(file_name);
        end
    case 4
        mean_im_phase = imread( "Phase_12.tif");
        mean_im_maxp = imread( "Max_power_distribution_4.tif");
        mean_im_mag   = imread( "Magnitude_8.tif");
                file_list = dir('Freq_distribution_*.tif');
        if ~isempty(file_list)
            file_name = file_list(1).name;
            max_freq_global = imread(file_name);
        end
    otherwise
        error('Invalid selection.');
end

% Display selected group number in the command window
disp(['You selected picture group ', num2str(selection), '.']);

end



%% SECTION 2: Find trigger cells
if cacul_corr > 0
    if find_triggercells > 0
        % Initialize arrays for counting shifts and correlations
        posShifts = zeros(1,size(S,1));
        correlations = zeros(1,size(S,1));

        % Count positive shifts and correlations
        for i = 1:size(S,1)
            shifts_of_roi = S(:,i);
            for ii = 1:size(S,1)
                if C(ii,i) >= threshcorr
                    correlations(i) = correlations(i)+1;
                    if shifts_of_roi(ii)>0
                        posShifts(i) = posShifts(i)+1;
                    end
                end
            end
        end

        % Identify trigger cells
        n = 0;
        for i = 1:size(posShifts,2)
            if posShifts(i) >= define_trigger_cells
                n = n+1;
                trigger(i) = 1;
                xt(n) = coordinatesData(1,(i*2)-1)/scaling;
                yt(n) = coordinatesData(1,(i*2))/scaling;
            else
                trigger(i) = 0;
            end
        end

        %% SECTION 3: Visualize trigger cells on original image
        figure('Name','Visulaization of the Trigger Cells')
        imshow(maxZImages);
        hold on
        b = plot(x,y,'c.');
        set(b,'MarkerSize',8);
        if exist('xt', 'var')
            b = plot(xt,yt,'r.');
            set(b,'MarkerSize',15);
        end
        if numbering > 0
            for i = 1:numCells
                if any(x(i) == xt & y(i) == yt)
                    tx = 5; ty = -5;
                    t = text(x(i)+tx, y(i)+ty, num2str(i));
                    set(t, 'Color', 'blue', 'FontSize', 8);
                end
            end
        end
        title('Visulaization of the Trigger Cells','color', 'b','FontSize', 16)
        outputValue = sprintf('%s_triggercells_with_number%s.fig',output,num2str(mode));
        savefig(outputValue);

           if  strcmp(params.testback, 'yes') 
        %% SECTION 4: Visualize trigger cells on magnitude image
        figure('Name','Visulaization of the Trigger Cells')
        imshow(mean_im_mag);
        hold on
        b = plot(x,y,'c.');
        set(b,'MarkerSize',0.1);
        if exist('xt', 'var')
            b = plot(xt,yt,'w.');
            set(b,'MarkerSize',8);
        end
        for i = 1:numCells
            if any(x(i) == xt & y(i) == yt)
                tx = 5; ty = -5;
                t = text(x(i)+tx, y(i)+ty, num2str(i));
                set(t, 'Color', 'blue', 'FontSize', 8);
            end
        end
        title('Visulaization of the Trigger Cells on Magnatitude picture','color', 'b','FontSize', 16)
        outputValue = sprintf('%s_triggercells_without_number_Magnatitude%s.fig',output,num2str(mode));
        savefig(outputValue);

        %% SECTION 5: Visualize trigger cells on phase image
        figure('Name','Visulaization of the Trigger Cells on phase')
        imshow(mean_im_phase);
        hold on
        b = plot(x,y,'c.');
        set(b,'MarkerSize',0.1);
        if exist('xt', 'var')
            b = plot(xt,yt,'w.');
            set(b,'MarkerSize',8);
        end
        for i = 1:numCells
            if any(x(i) == xt & y(i) == yt)
                tx = 5; ty = -5;
                t = text(x(i)+tx, y(i)+ty, num2str(i));
                set(t, 'Color', 'blue', 'FontSize', 8);
            end
        end
        title('Visulaization of the Trigger Cells on Phase picture','color', 'b','FontSize', 16)
        outputValue = sprintf('%s_triggercells_without_number_phase%s.fig',output,num2str(mode));
        savefig(outputValue);

        %% SECTION 6: Visualize trigger cells on maxpower image
        figure('Name','Visulaization of the Trigger Cells on maxpower')
        imshow(mean_im_maxp);
        hold on
        b = plot(x,y,'c.');
        set(b,'MarkerSize',0.1);
        if exist('xt', 'var')
            b = plot(xt,yt,'w.');
            set(b,'MarkerSize',8);
        end
        for i = 1:numCells
            if any(x(i) == xt & y(i) == yt)
                tx = 5; ty = -5;
                t = text(x(i)+tx, y(i)+ty, num2str(i));
                set(t, 'Color', 'blue', 'FontSize', 8);
            end
        end
        title('Visulaization of the Trigger Cells on maxpower Picture','color', 'b','FontSize', 16)
        outputValue = sprintf('%s_triggercells_without_number_frequencey%s.fig',output,num2str(mode));
        savefig(outputValue);
           end
    end
end

%% SECTION 7: Update and export results
% Store trigger cell information in params
params.trigger = trigger;
params.triggerX = xt;
params.triggerY = yt;
params.posShifts = posShifts;
params.correlations = correlations;

end