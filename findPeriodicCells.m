
 
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

function [params, analysis_summary] = findPeriodicCells(params, analysis_summary)
% Extract required variables from params
cacul_corr = params.cacul_corr;
C = params.C;
threshcorr = params.threshcorr;
find_periodicities = params.find_periodicities;
find_triggercells = params.find_triggercells;
S = params.S;
Dir = params.Dir;
peakLocations = params.peakLocations;
secondsPerFrame = params.secondsPerFrame;
pk_std = params.pk_std;
pk_num_define = params.pk_num_define;
coordinatesData = params.coordinatesData;
scaling = params.scaling;
output = params.output;
mode = params.mode;
maxZImages = params.maxZImages;
tracesleft = params.tracesleft;
amplitudes = params.amplitudes;
trigger = params.trigger;
hubs = params.hubs;
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
% Note that if you redo this code many times,loading pitures may
% go wrong, just change the name of the piture what you wangt
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


%% %%%%%%%%%%%%%%%%%%%%%%%find_perioficity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%🚩
if find_periodicities > 0
        if cacul_corr == 0 || find_triggercells == 0
            trigger = 0;
            hubs = 0;
        end
        if cacul_corr == 0
            C=0;
            S=0;
            Dir=0;
        end


lengthseq = 10;    %transfer the length of the sequences to be correlated in minutes
steps = 5;  %The time by which the sequences will be shifted, default: 5
%Indicates the number of sequences to be analyzed: 
% (length of the video - length of the last sequence) / the time by which the sequences will be shifted
numofseq = ceil(((size(peakLocations,1)*secondsPerFrame/60)-lengthseq) / steps); 
if numofseq < 1
    numofseq = 1;
end

for a=1 : numofseq
    start = floor(((a-1)*(steps*60/secondsPerFrame)))+1;
    lim = ceil(start-1+(lengthseq*60/secondsPerFrame));
    if a == numofseq
        lim = size(peakLocations,1);
    end
    di=lim-start;
    seqPeakLocations(a,1:di+1,:) = peakLocations(start:lim,:);
end

meanperiod = NaN(1,size(peakLocations,2));
deviation = NaN(1,size(peakLocations,2));
meanperiodsteps = NaN(numofseq,size(peakLocations,2));
deviationsteps = NaN(numofseq,size(peakLocations,2));

cellsabove3peaks = 0;
for roi = 1 : size(peakLocations,2)
    
    if sum(peakLocations(:,roi)) >= pk_num_define
        
        cellsabove3peaks = cellsabove3peaks + 1;
        %This loop (step) repeats the correlations for each sequence (seq)
        % and selects the lowest std at the end
        for step=1 : numofseq 
            
            pkpos = find(seqPeakLocations(step,:,roi));
            if sum(seqPeakLocations(step,:,roi)) > 3
                
                period = zeros(1,size(pkpos,1)-1);
                for i = 1 : size(pkpos,2)-1
                    period(i) = pkpos(i+1) - pkpos(i);
                end
                
                meanperiodsteps(step,roi) = mean(period);
                deviationsteps(step,roi) = std(period);
                
            else
                meanperiodsteps(step,roi) = NaN;
                deviationsteps(step,roi) = NaN;
            end
            
            clear pkpos period
            
        end
        
        [deviation(roi),loc] = min(deviationsteps(:,roi));
        meanperiod(roi) = meanperiodsteps(loc,roi);
        
        
    else
        meanperiod(roi) = NaN;
        deviation(roi) = NaN;
        
    end
end

if cacul_corr > 0
    % Here the number of correlations of the periodic and non-periodic cells is counted
    correlationsAll = NaN(1,size(S,1));
    posShiftsall = NaN(1,size(S,1));
    posShiftsperio = NaN(1,size(S,1)); % Number of cells triggered (all periodic cells)
    correlationsPeriodic = NaN(1,size(S,1)); % Number of connected cells (all periodic cells)
    posShiftsNOTperio = NaN(1,size(S,1)); % Number of cells triggered (all NON-periodic cells)
    correlationsNonPeriodic = NaN(1,size(S,1)); % Number of connected cells (all NON-periodic cells)
    directionAll = NaN(size(S,1)); % Directionality of all periodic cells
    directionPeriodic = NaN(size(S,1)); % Directionality of all periodic cells
    directionNonPeriodic = NaN(size(S,1)); % Directionality of all NON-periodic cells
    
    for i = 1 : size(S,1)
        if sum(peakLocations(:,i)) > 3 
            
            shifts_of_roi = S(i,:);
            if deviation(i) <= pk_std/secondsPerFrame 
                if isnan(correlationsAll(i))
                    correlationsAll(i)=0;
                    posShiftsall(i) = 0;
                    correlationsPeriodic(i) = 0;
                    posShiftsperio(i) = 0;
                end
                
                for ii = 1 : size(S,1)
                    if C(ii,i) >= threshcorr
                        correlationsPeriodic(i) = correlationsPeriodic(i)+1;
                        correlationsAll(i) = correlationsAll(i)+1;
                        if shifts_of_roi(ii)>0
                            posShiftsperio(i) = posShiftsperio(i)+1;
                            posShiftsall(i) = posShiftsall(i)+1;
                        end
                    end
                end
                directionPeriodic(i,:) = Dir(i,:);
                directionAll(i,:) = Dir(i,:);
                
            else % if it is a NON-periodic cell:
                if isnan(correlationsAll(i))
                    correlationsAll(i)=0;
                    posShiftsall(i) = 0;
                    correlationsNonPeriodic(i) = 0;
                    posShiftsNOTperio(i) = 0;
                end
                for ii = 1 : size(S,1)
                    if C(ii,i) >= threshcorr
                        correlationsNonPeriodic(i) = correlationsNonPeriodic(i)+1;
                        correlationsAll(i) = correlationsAll(i)+1;
                        if shifts_of_roi(ii)>0
                            posShiftsNOTperio(i) = posShiftsNOTperio(i)+1;
                            posShiftsall(i) = posShiftsall(i)+1;
                        end
                    end
                end
                directionNonPeriodic(i,:) = Dir(i,:);
                directionAll(i,:) = Dir(i,:);
                
            end
        end
    end
    
    % directionPeriodic = (posShiftsperio + posShiftsperio - correlationsPeriodic)./correlationsPeriodic;
    % directionNonPeriodic = (posShiftsNOTperio + posShiftsNOTperio - correlationsNonPeriodic)./correlationsNonPeriodic;
    
    
    analysis_summary(6,1) = "Mean number of co-active cells per periodic cell";
    analysis_summary(6,2) = mean(correlationsPeriodic,'omitnan');
    analysis_summary(7,1) = "Mean number of co-active cells per NOT-periodic cell (＞4 peaks）";
    analysis_summary(7,2) = mean(correlationsNonPeriodic,'omitnan');
    analysis_summary(8,1) = "Mean number of triggered cells per periodic cell";
    analysis_summary(8,2) = mean(posShiftsperio,'omitnan');
    analysis_summary(9,1) = "Mean number of triggered cells per NOT-periodic cell (＞4 peaks）";
    analysis_summary(9,2) = mean(posShiftsNOTperio,'omitnan');

    directionAll = directionAll(:);
    directionAll = directionAll(~isnan(directionAll));
    directionPeriodic = directionPeriodic(:);
    directionPeriodic = directionPeriodic(~isnan(directionPeriodic));
    directionNonPeriodic = directionNonPeriodic(:);
    directionNonPeriodic = directionNonPeriodic(~isnan(directionNonPeriodic));
    
    correlationsAll = correlationsAll(~isnan(correlationsAll));
    posShiftsall = posShiftsall(~isnan(posShiftsall));
    correlationsPeriodic = correlationsPeriodic(~isnan(correlationsPeriodic));
    correlationsNonPeriodic = correlationsNonPeriodic(~isnan(correlationsNonPeriodic));
    posShiftsperio = posShiftsperio(~isnan(posShiftsperio));
    posShiftsNOTperio = posShiftsNOTperio(~isnan(posShiftsNOTperio));
    networkAll = (2*posShiftsall)-correlationsAll;
    networkPeriodic = (2*posShiftsperio)-correlationsPeriodic;
    networkNonPeriodic = (2*posShiftsNOTperio)-correlationsNonPeriodic;
    
    if size(directionAll,1) > size(networkAll,2)
        mainMatrix = NaN(size(directionAll,1),10);
    else
        mainMatrix = NaN(size(networkAll,2),10);
    end
    mainMatrix(1:size(correlationsAll,2),1) = correlationsAll.';
    mainMatrix(1:size(correlationsPeriodic,2),2) = correlationsPeriodic.';
    mainMatrix(1:size(correlationsNonPeriodic,2),3) = correlationsNonPeriodic.';
    mainMatrix(1:size(networkAll,2),4) = networkAll.';
    mainMatrix(1:size(networkPeriodic,2),5) = networkPeriodic.';
    mainMatrix(1:size(networkNonPeriodic,2),6) = networkNonPeriodic.';
    mainMatrix(1:size(directionAll),8) = directionAll;
    mainMatrix(1:size(directionPeriodic),9) = directionPeriodic;
    mainMatrix(1:size(directionNonPeriodic),10) = directionNonPeriodic;
    
    writematrix (mainMatrix, output + "_ccc_nnn_ddd.txt");
end

numCells = size(coordinatesData,2)/2;
periodicity = zeros(1,numCells);
n=0;
m=0;
amplitudes_perio = NaN;
amplitudes_non_perio = NaN;

for i=1 : numCells
    if deviation(i) <= pk_std/secondsPerFrame
        n=n+1;
        periodicity(i) = 1;
        periodsofperiodiccells(n) = meanperiod(i)*secondsPerFrame;
        x(n) = coordinatesData(1,(i*2)-1)/scaling;
        y(n) = coordinatesData(1,(i*2))/scaling;
        p(n) = meanperiod(i)*secondsPerFrame;
        
        activityperiodic(n) = sum(peakLocations(:,i));
        amplitudes_perio(n) = mean(nonzeros(amplitudes(:,i)));
        
    elseif deviation(i) > pk_std/secondsPerFrame
        m=m+1;
        activitynotperiodic(m) = sum(peakLocations(:,i));
        
        amplitudes_non_perio(m) = mean(nonzeros(amplitudes(:,i)));

    end
end

writematrix (amplitudes_perio.', output + "_peakamplitudes_perio.txt");
writematrix (amplitudes_non_perio.', output + "_peakamplitudes_non-perio.txt");

if exist('periodicity','var')
else
    periodicity = 0;
end

x1=zeros(numCells,1);
y1=zeros(numCells,1);
for i=1 : numCells
    x1(i)= coordinatesData(1,(i*2)-1)/scaling;
    y1(i)= coordinatesData(1,(i*2))/scaling;
end
% ***plotting***
figure('Name','Periodic cells'); %.................................🚨🚩❗🗝️
imshow(maxZImages); 
hold on 
b=plot(x1,y1,'c.');
set(b,'MarkerSize',13);
if exist('x', 'var')
    b=plot(x,y,'g.');
    set(b,'MarkerSize',20);
    tx = 10; ty = -10;  % displacement so the text does not overlay the data points
    if isfield(params, 'numbering') && params.numbering == 1
    for i=1 : size(p,2)
        t = text(x(i)+tx, y(i)+ty, num2str(round(p(i))));
        set(t, 'Color',[1, 0 ,0], 'FontSize', 20)
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if  strcmp(params.testback, 'yes') 
% ***2. phase图***
figure('Name','Periodic cells'); %.................................🚨🚩❗🗝️
imshow(mean_im_phase); 
hold on 
b=plot(x1,y1,'c.');
set(b,'MarkerSize',2);
if exist('x', 'var')
    b=plot(x+3,y+3,'w.');
    set(b,'MarkerSize',18);
    
end
drawnow
outputValue = sprintf('%s_periodicity_with_number_Phase%s.fig',output,num2str(mode));
title ('Visulaization of Cells with Periodic activity on Phase Picture','color', 'b','FontSize', 16)
savefig(outputValue);
hold off


% ***3. Magtitude图***
figure('Name','Periodic cells'); %.................................🚨🚩❗🗝️
imshow(mean_im_mag); 
hold on 
b=plot(x1,y1,'c.');
set(b,'MarkerSize',2);
if exist('x', 'var')
    b=plot(x+3,y+3,'w.');
    set(b,'MarkerSize',18);
end
drawnow
outputValue = sprintf('%s_periodicity_with_number_Magtitude%s.fig',output,num2str(mode));
title ('Visulaization of Cells with Periodic activity on Magtitude Picture','color', 'b','FontSize', 16)
savefig(outputValue);
hold off


% ***4. max_power_global图***
figure('Name','Periodic cells'); %.................................🚨🚩❗🗝️
imshow(mean_im_maxp); 
hold on 
b=plot(x1,y1,'c.');
set(b,'MarkerSize',2);
if exist('x', 'var')
    b=plot(x+3,y+3,'w.');
    set(b,'MarkerSize',18);
end
drawnow
outputValue = sprintf('%s_periodicity_with_number_maxpower%s.fig',output,num2str(mode));
title ('Visulaization of Cells with Periodic activity on maxpower Picture','color', 'b','FontSize', 16)
savefig(outputValue);
hold off

% ***4. max_freq_global图***
figure('Name','max_power_global'); %.................................🚨🚩❗🗝️
imshow(max_freq_global); 
hold on 
b=plot(x1,y1,'c.');
set(b,'MarkerSize',2);
if exist('x', 'var')
    b=plot(x+3,y+3,'w.');
    set(b,'MarkerSize',18);
end
drawnow
outputValue = sprintf('%s_periodicity_with_number_global_freq%s.fig',output,num2str(mode));
title ('Visulaization of Cells with Periodic activity on global_freq Picture','color', 'b','FontSize', 16)
savefig(outputValue);
hold off


           end



analysis_summary(10,1) = "cells with at least 4 peaks （parameters: pk_num_define）";
analysis_summary(10,2) = cellsabove3peaks;
analysis_summary(11,1) = "periodic cells";
analysis_summary(11,2) = sum(periodicity(periodicity==1));
analysis_summary(6,3) = "mean of all periodic periods";
analysis_summary(7,3) = "std of all periodic periods";
if exist('periodsofperiodiccells','var')
    analysis_summary(6,4) = mean(periodsofperiodiccells);
    analysis_summary(7,4) = std(periodsofperiodiccells);
else
    analysis_summary(6,4) = "NaN";
    analysis_summary(7,4) = "NaN";
    periodsofperiodiccells = NaN;
end
writematrix (periodsofperiodiccells.', output + "_periodsofperiodiccells.txt");


analysis_summary(8,3) = "proportion of periodic cells in network";
if cellsabove3peaks == 0
    analysis_summary(8,4) = "NaN";
else
    analysis_summary(8,4) = sum(periodicity(periodicity==1)) / cellsabove3peaks;
end

sparsePeriodicity = sparse(periodicity);
sparse_meanperiod = sparse(meanperiod(~isnan(meanperiod)));
sparse_deviation = sparse(deviation(~isnan(deviation)));
combined = sparse_meanperiod*secondsPerFrame;
combined(2,:) = sparse_deviation*secondsPerFrame;
TF = isempty(combined);
if TF == 1
    combined = 0;
end
% writematrix (sparsePeriodicity, output + "_periodicity.txt");
% writematrix (combined, output + "_mean_sd_periods.txt");

if find_triggercells > 0 && cacul_corr > 0
    trigger = full(trigger);
    hubs = full(hubs);
    % Count cells that are both trigger and periodic
    overlaptri = sum(periodicity(trigger==1)); 
    % Count cells that are both hubs and periodic
    overlaphub = sum(periodicity(hubs==1)); 
    activityperiodiccell = size(activityperiodic,2);
    analysis_summary(9,3) = "periodic cells in trigger cells";
    analysis_summary(9,4) = overlaptri / sum(trigger(trigger==1));
    analysis_summary(10,3) = "trigger cells in periodic cells";
    analysis_summary(10,4) = overlaptri / sum(periodicity(periodicity==1));
    analysis_summary(12,3) = "number of periodic cells in hub cells";
    analysis_summary(12,4) = overlaphub;
    analysis_summary(11,3) = "Percentage of periodic cells in hub cells";
    analysis_summary(11,4) = overlaphub / sum(hubs(hubs==1));
    analysis_summary(17,3) = "Percentage of hub cells in  periodic cells";
    analysis_summary(17,4) = overlaphub / sum(periodicity(periodicity==1));

end
% Calculate mean number of peaks for periodic and non-periodic cells
if exist('activityperiodic','var')
    meanactperio = mean(activityperiodic(activityperiodic>3)); 
else
    meanactperio = NaN;
    activityperiodic = NaN;
end
% The average number of peaks of all non-periodic cells with at least 4 peaks.
if exist('activitynotperiodic','var')
    meanactnotperio = mean(activitynotperiodic(activitynotperiodic>3)); 
else
    meanactnotperio = NaN;
    activitynotperiodic = NaN;
end

analysis_summary(13,3) = "Mean number of peaks of periodic cells";
analysis_summary(13,4) = meanactperio;
analysis_summary(14,3) = "Mean number of peaks of NOT-periodic cells that have at least 4 peaks";
analysis_summary(14,4) = meanactnotperio;

perioactivity = NaN(size(peakLocations,2),2);
perioactivity(1:size(activityperiodic,2),1) = activityperiodic.';
perioactivity(1:size(activitynotperiodic,2),2) = activitynotperiodic.';
writematrix (perioactivity, output + "_perioactivity-perio-notperio.txt");


% Trigger coefficiant and Cluster coefficiant
if find_triggercells > 0 && cacul_corr > 0
    
    triggerCoefficient = NaN(size(peakLocations,2),1);
    clusterCoefficient = NaN(size(peakLocations,2),1);
    clear length
    allcorrelations = length(find(C>=threshcorr))/2;
    alltriggers = allcorrelations/2;
    
    mean_triggering = alltriggers / cellsabove3peaks;
    mean_clustering = allcorrelations / cellsabove3peaks;
    
    for i=1 : size(peakLocations,2)
        if sum(peakLocations(:,i))>=pk_num_define
            triggerCoefficient(i,1) = length(find(C(i,S(i,:)>0)>=threshcorr)) / (2*mean_triggering);
            clusterCoefficient(i,1) = length(find(C(i,:)>=threshcorr)) / (2*mean_clustering);
        end
    end
    
    triggerCoefficientPeriodic = triggerCoefficient(find(deviation <= pk_std/secondsPerFrame),1);
    triggerCoefficientNonPeriodic = triggerCoefficient(find(deviation > pk_std/secondsPerFrame),1);
    clusterCoefficientPeriodic = clusterCoefficient(find(deviation <= pk_std/secondsPerFrame),1);
    clusterCoefficientNonPeriodic = clusterCoefficient(find(deviation > pk_std/secondsPerFrame),1);
    
    periodicCoefficient = NaN(size(peakLocations,2),4);
    periodicCoefficient(1:size(clusterCoefficientPeriodic,1),1) = clusterCoefficientPeriodic;
    periodicCoefficient(1:size(clusterCoefficientNonPeriodic,1),2) = clusterCoefficientNonPeriodic;
    periodicCoefficient(1:size(triggerCoefficientPeriodic,1),3) = triggerCoefficientPeriodic;
    periodicCoefficient(1:size(triggerCoefficientNonPeriodic,1),4) = triggerCoefficientNonPeriodic;
    
    
    writematrix (periodicCoefficient, output + "_perio_coeff.txt");
 end
end

% proportion of periodic cells in non-hub cells
if exist('hubs','var') && ~isempty(hubs)

    isPeriodic = (sum(peakLocations) > 3) & (deviation <= pk_std/secondsPerFrame);
    isPeriodic = isPeriodic(:); 
    isHub = logical(hubs(:));  
    nonHubPeriodic = ~isHub & isPeriodic; % including peaks ≥ 4
    numNonHubPeriodic = sum(nonHubPeriodic);
    numNonHubs = tracesleft-sum(hubs(hubs==1)); 
    
    if numNonHubs > 0
        percentNonHubPeriodic = numNonHubPeriodic / numNonHubs * 100;
    else
        percentNonHubPeriodic = NaN;
    end
% std of at least 4 peaks cells
        min_peaks = params.pk_num_define;
        peakLocations = params.peakLocations;
        numCells = size(peakLocations, 2);
        deviation4peaks = NaN(1, numCells);  
        
        for i = 1:numCells
            if sum(peakLocations(:, i)) >= min_peaks
                deviation4peaks(i) = deviation(i); 
            else
                deviation4peaks(i) = NaN;  
            end
        end
        
        valid_deviation4peaks = deviation4peaks(~isnan(deviation4peaks));
        
        [min_std_val, idx_in_valid] = min(valid_deviation4peaks);
        
        
        original_idx = find(~isnan(deviation4peaks));
        cell_id = original_idx(idx_in_valid);
        
        
        fprintf('smallest std %.2f second，from cell %d \n', min_std_val, cell_id);


        % caculate the dircetion

        filesaaa = dir('*_ccc_nnn_ddd.txt');
        for k = 1:length(filesaaa)
            filename = filesaaa(k).name;
            data = readmatrix(filename);  % 
            %
        end
        direction_all = data(:, 8);
        dirperio = data(:, 9);
        dirnotperio = data(:, 10);
        pertt = str2num(analysis_summary{11, 2});
        direction_allpercell = sum(direction_all, 'omitnan')/params.tracesleft;
        dirextion_periopercell = sum(dirperio, 'omitnan')/pertt;
        direction_notperiopercell = sum(dirnotperio, 'omitnan')/(params.tracesleft-pertt);
end
%% Store modified/new variables back to params
params.periodsofperiodiccells = periodsofperiodiccells;
params.deviation4peaks=deviation4peaks;
params.periodicity = periodicity;
params.meanperiod = meanperiod;
params.deviation = deviation;
params.cellsabove3peaks = cellsabove3peaks;
params.correlationsAll = correlationsAll;
params.correlationsPeriodic = correlationsPeriodic;
params.correlationsNonPeriodic = correlationsNonPeriodic;
params.posShiftsperio = posShiftsperio;
params.posShiftsNOTperio = posShiftsNOTperio;
params.activityperiodic = activityperiodic;
params.activitynotperiodic = activitynotperiodic;
params.triggerCoefficient = triggerCoefficient;
params.clusterCoefficient = clusterCoefficient;
params.periodicCoefficient = periodicCoefficient;
analysis_summary(15,3) = "Number of periodic cells among non-hub cells";
analysis_summary(15,4) = numNonHubPeriodic;
analysis_summary(16,3) = "Percentage of periodic cells among non-hub cells";
analysis_summary(16,4) = percentNonHubPeriodic;
% Export important variables to workspace
assignin('base', 'periodicity', periodicity);
assignin('base', 'meanperiod', meanperiod);
assignin('base', 'deviation', deviation);
assignin('base', 'cellsabove3peaks', cellsabove3peaks);

analysis_summary(27,1) = "direction_allpercell";
analysis_summary(27,2) = direction_allpercell;
analysis_summary(28,1) = "direction_periopercell";
analysis_summary(28,2) = dirextion_periopercell;
analysis_summary(29,1) = "direction_notperiopercell";
analysis_summary(29,2) = direction_notperiopercell;
analysis_summary(30,1) = "percentage of periodic cells in network (above 3 peaks)";
analysis_summary(30,2) = sum(periodicity(periodicity==1))/tracesleft;
end


 filename = sprintf('all_for_3D_%d.mat', params.mode);
 save(filename);

end