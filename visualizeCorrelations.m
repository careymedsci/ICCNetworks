
% This part was modified from https://zenodo.org/badge/latestdoi/556336211
% Please contanct them for persmission.


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
% 2025.03.09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [params, analysis_summary] = visualizeCorrelations(params, analysis_summary)
% visualizeCorrelations - Visualizes correlations between cells on different image types
%
% Inputs:
%   params - Structure containing analysis parameters and data
%   analysis_summary - Structure containing analysis results
%
% Outputs:
%   params - Updated parameters structure
%   analysis_summary - Updated analysis summary
%   Also creates and saves multiple visualization figures

%% SECTION 1: Initialize variables
% Extract variables from params structure
cacul_corr = params.cacul_corr;
change_corr = params.change_corr;
C = params.C;
S = params.S;
threshcorr = params.threshcorr;
do_lines_on_image = params.do_lines_on_image;
coordinatesData = params.coordinatesData;
scaling = params.scaling;
maxZImages = params.maxZImages;
numbering = params.numbering;

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
%% SECTION 2: Process sparse correlation data
if cacul_corr > 0
    if change_corr == 1
        % Load and convert sparse correlation data
        sparseCC = "_svddata_sparseCC.dat";
        sparseCC = load(sparseCC);
        sparseCC = spconvert(sparseCC);
        CC = full(sparseCC);
        
        % Initialize and fill correlation matrix
        CCm = zeros(size(C,1));
        CCm(1:size(CC,1),1:size(CC,2)) = CC;
        
        % Update correlation and speed matrices based on threshold
        for i = 1:size(C,1)
            for ii = 1:size(C,1)
                if C(i,ii) >= threshcorr
                    if CCm(i,ii) ~= 1 || CCm(ii,i) ~= 1
                        C(i,ii) = NaN;
                        S(i,ii) = NaN;
                    end
                end
            end
        end
    end
end

%% SECTION 3: Visualization preparation
if cacul_corr > 0
    if do_lines_on_image > 0
        % Calculate cell coordinates
        numCells = size(coordinatesData,2)/2;
        x = zeros(numCells,1);
        y = zeros(numCells,1);
        for i = 1:numCells
            x(i) = coordinatesData(1,(i*2)-1)/scaling;
            y(i) = coordinatesData(1,(i*2))/scaling;
        end

        % Create color map for correlation visualization
        ColorMap = jet(100);
        
        %% SECTION 4: Plot on original image
        figure('Name','Visulaization of the Correlations')
        imshow(maxZImages);
        hold on
        a = plot(x,y,'c.');
        set(a,'MarkerSize',10);
        drawnow
        for i = 1:numCells
            for ii = 1:numCells
                if C(ii,i) >= threshcorr
                    LineColor = ColorMap(fix(C(ii,i)*100), :);
                    b = line([x(i) x(ii)], [y(i) y(ii)], 'Color', LineColor);
                    set(b,'LineWidth',2);
                end
            end
            if numbering > 0
                tx = 5; ty = -5;
                t = text(x(i)+tx, y(i)+ty, num2str(i));
                set(t, 'Color',[1, 0 ,0], 'FontSize', 6)
            end
        end
        drawnow
        title('Visulaization of the Correlations','Color','r')
        outputFigure = sprintf('_lines on image%s%s.fig',num2str(numbering));
        savefig(outputFigure);

        %% SECTION 5: Plot without background image
        figure('Name','Visulaization of the Correlations without Original Picture')
        empty_image = ones(size(maxZImages,1), size(maxZImages,2),3);
        imshow(empty_image, 'InitialMagnification', 'fit');
        hold on
        a = plot(x,y,'c.');
        set(a,'MarkerSize',10);
        drawnow
        for i = 1:numCells
            for ii = 1:numCells
                if C(ii,i) >= threshcorr
                    LineColor = ColorMap(fix(C(ii,i)*100), :);
                    b = line([x(i) x(ii)], [y(i) y(ii)], 'Color', LineColor);
                    set(b,'LineWidth',2);
                end
            end
            if numbering > 0
                tx = 5; ty = -5;
                t = text(x(i)+tx, y(i)+ty, num2str(i));
                set(t, 'Color',[1, 0 ,0], 'FontSize', 6)
            end
        end
        drawnow
        title('Visulaization of the Correlations without Original Picture','Color','r')
        outputFigure = sprintf('_lines on image_witout_picture%s%s.fig',num2str(numbering));
        savefig(outputFigure);


           if  strcmp(params.testback, 'yes') 
        %% SECTION 6: Plot on phase image
        figure('Name','Visulaization of the Correlations on phase picture')
        imshow(mean_im_phase);
        hold on
        a = plot(x,y,'c.');
        set(a,'MarkerSize',1);
        drawnow
        for i = 1:numCells
            for ii = 1:numCells
                if C(ii,i) >= threshcorr
                    LineColor = ColorMap(fix(C(ii,i)*100), :);
                    b = line([x(i) x(ii)], [y(i) y(ii)], 'Color', LineColor);
                    set(b,'LineWidth',2);
                end
            end
            if numbering > 0
                tx = 5; ty = -5;
                t = text(x(i)+tx, y(i)+ty, num2str(i));
                set(t, 'Color',[1, 1 ,1], 'FontSize', 6)
            end
        end
        drawnow

        figure(gcf) % 
            set(gcf, 'Color', 'none', 'InvertHardcopy', 'off'); % 
            print('correlation_arrows.svg', '-dsvg', '-vector'); % 

        title('Visulaization of the Correlations on Phase picture','Color','r')
        outputFigure = sprintf('_lines on Phase image%s%s.fig',num2str(numbering));
        savefig(outputFigure);



        %% SECTION 7: Plot on magnitude image
        figure('Name','Visulaization of the Correlations on magnitude picture')
        imshow(mean_im_mag);
        hold on
        a = plot(x,y,'c.');
        set(a,'MarkerSize',1);
        drawnow
        for i = 1:numCells
            for ii = 1:numCells
                if C(ii,i) >= threshcorr
                    LineColor = ColorMap(fix(C(ii,i)*100), :);
                    b = line([x(i) x(ii)], [y(i) y(ii)], 'Color', LineColor);
                    set(b,'LineWidth',2);
                end
            end
            if numbering > 0
                tx = 5; ty = -5;
                t = text(x(i)+tx, y(i)+ty, num2str(i));
                set(t, 'Color',[1, 1 ,1], 'FontSize', 6)
            end
        end
        drawnow
        title('Visulaization of the Correlations on Magnitude picture','Color','r')
        outputFigure = sprintf('_lines on Magnitude image%s%s.fig',num2str(numbering));
        savefig(outputFigure);

        %% SECTION 8: Plot on maxpower with specific image
        figure('Name','Visulaization of the Correlations on maxpower picture')
        imshow(mean_im_maxp);
        hold on
        a = plot(x,y,'c.');
        set(a,'MarkerSize',1);
        drawnow
        for i = 1:numCells
            for ii = 1:numCells
                if C(ii,i) >= threshcorr
                    LineColor = ColorMap(fix(C(ii,i)*100), :);
                    b = line([x(i) x(ii)], [y(i) y(ii)], 'Color', LineColor);
                    set(b,'LineWidth',2);
                end
            end
            if numbering > 0
                tx = 5; ty = -5;
                t = text(x(i)+tx, y(i)+ty, num2str(i));
                set(t, 'Color',[1, 1 ,1], 'FontSize', 6)
            end
        end
        drawnow
        title('Visulaization of the Correlations on maxpower picture','Color','r')
        outputFigure = sprintf('_lines on Frequency image%s%s.fig',num2str(numbering));
        savefig(outputFigure);
           end

    end
end

%% SECTION 9: Update and export results
% Update params structure with modified values
params.C = C;
params.S = S;
params.x = x;
params.y = y;


end