

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


function [params, analysis_summary] = findHubCells(params, analysis_summary)
% findHubCells - Identifies hub cells and visualizes them with trigger cells
%
% Inputs:
%   params - Structure containing analysis parameters and data
%   analysis_summary - Structure containing analysis results
%
% Outputs:
%   params - Updated parameters structure
%   analysis_summary - Updated analysis summary
%   Also creates and saves visualization figures

%% SECTION 1: Initialize variables
% Extract variables from params structure
correlations = params.correlations;
define_hub_cells = params.define_hub_cells;
coordinatesData = params.coordinatesData;
scaling = params.scaling;
x = params.x;
y = params.y;
maxZImages = params.maxZImages;
numbering = params.numbering;
output = params.output;
mode = params.mode;
numCells =  params.numRois;
tracesleft = params.tracesleft;
% Get trigger cells data if exists
if isfield(params, 'triggerX')
    xt = params.triggerX;
    yt = params.triggerY;
    trigger = params.trigger;
end

%% SECTION 2: Find hub cells
n = 0;
for i = 1:size(correlations,2)
    if correlations(i) >= define_hub_cells
        n = n+1;
        hubs(i) = 1;
        xh(n) = coordinatesData(1,(i*2)-1)/scaling;
        yh(n) = coordinatesData(1,(i*2))/scaling;
    else
        hubs(i) = 0;
    end
end

%% SECTION 3: Visualize hub cells
figure('Name','Visulaization of the Hub Cells')
imshow(maxZImages);
hold on
b = plot(x,y,'c.');
set(b,'MarkerSize',8);

if exist('xh', 'var')
    b = plot(xh+2,yh+2,'b.');
    set(b,'MarkerSize',18);
end

if numbering > 0
    for i = 1:numCells
        % Check if current cell is in the list of interest
        if any(x(i) == xh & y(i) == yh)
            tx = 5; ty = -5;  % Displacement so the text does not overlay the data points
            t = text(x(i)+tx, y(i)+ty, num2str(i));
            set(t, 'Color', 'red', 'FontSize', 8);
        end
    end
end

drawnow
title('Visulaization of the Hub Cells','color', 'b','FontSize', 16)
outputValue = sprintf('%s_hubcells_withnumber%s.fig',output,num2str(mode));
savefig(outputValue);

%% SECTION 4: Visualize hub cells with trigger cells
figure('Name','Visulaization of the Hub Cells with Trigger Cells')
imshow(maxZImages);
hold on
b = plot(x,y,'c.');
set(b,'MarkerSize',10);
if exist('xt', 'var')
    l = plot(xt,yt,'r.');
    set(l,'MarkerSize',20);
end


if exist('xh', 'var')
    % Offset hub cells to prevent overlap with trigger cells
    b = plot(xh+1,yh+1,'b.');
    set(b,'MarkerSize',13);
end
drawnow
title('Visulaization of the Hub Cells with Trigger Cells','color', 'b','FontSize', 16)
outputValue = sprintf('%sTriggers_hubcells%s.fig',output,num2str(mode));
savefig(outputValue);
legend('Roi', 'Trigger Cells','Hub Cells', 'Location', 'northeastoutside', 'Box', 'off','FontSize', 12);

%% SECTION 5: Update analysis summary
analysis_summary(17,1) = "number of trigger cells";
analysis_summary(17,2) = sum(trigger(trigger==1));
analysis_summary(18,1) = "percentage of trigger cells in cells network (above 3 peaks)";
if tracesleft == 0
    analysis_summary(18,2) = "NaN";
else
    analysis_summary(18,2) = sum(trigger(trigger==1)) / tracesleft;
end

analysis_summary(19,1) = "hub cells";
analysis_summary(19,2) = sum(hubs(hubs==1));
analysis_summary(20,1) = "proportion of hub cells in network (above 3 peaks)";
if tracesleft == 0
    analysis_summary(20,2) = "NaN";
else
    analysis_summary(20,2) = sum(hubs(hubs==1)) / tracesleft;
end

%% SECTION 6: Update and export results
% Store hub cell information in params
params.hubs = hubs;
params.hubX = xh;
params.hubY = yh;


end