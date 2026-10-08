

 
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
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany
% Xiangyang No.1 people's Hospital, China
% Also an python version coded 
% Xiao.Liu@stud.uni-frankfurt.de
% 2024.10.09


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = visualizeCorrelationsWithArrows(params)
% visualizeCorrelationsWithArrows - Visualize cell correlations with optimized arrows
%
% Required input variables from params:
%   C - Correlation matrix
%   Dir - Direction matrix
%   scaling - Scaling factor for coordinates
%   coordinatesData - Cell coordinates
%   maxZImages - Maximum Z projection image
%   threshcorr - Correlation threshold
%   numbering - Flag for cell numbering

% Extract required variables from params
load('all_verVisulization_0.mat');

C = params.C;
Dir = params.Dir;
scaling = params.scaling;
coordinatesData = params.coordinatesData;
maxZImages = params.maxZImages;
threshcorr = params.threshcorr;
numbering = params.numbering;

%% SECTION 1: Process coordinates
numCells = size(coordinatesData,2)/2;
x = zeros(numCells,1);
y = zeros(numCells,1);
for i = 1:numCells
    x(i) = coordinatesData(1,(i*2)-1)/scaling;
    y(i) = coordinatesData(1,(i*2))/scaling;
end

%% SECTION 2: Create color map
ColorMap = jet(100);

%% SECTION 3: Create visualization
figure('Name','Visualization of the Correlations')
imshow(maxZImages);
hold on

% Plot cell positions
a = plot(x,y,'c.');
set(a,'MarkerSize',10);
drawnow

%% SECTION 4: Draw correlation arrows
for i = 1:numCells
    for ii = 1:numCells
        if C(ii,i) >= threshcorr
            % Get color (with boundary protection)
            colorIdx = max(1, min(100, fix(C(ii,i)*100)));
            LineColor = ColorMap(colorIdx, :);
            
            % Calculate arrow length (for head proportion adjustment)
            arrowLength = norm([x(ii)-x(i), y(ii)-y(i)]);
            
            % Determine direction based on Dir matrix
            if Dir(ii,i) > 0
                start = [x(i), y(i)];
                delta = [x(ii)-x(i), y(ii)-y(i)];
            else
                start = [x(ii), y(ii)];
                delta = [x(i)-x(ii), y(i)-y(ii)];
            end
            
            % Draw optimized arrow
            quiver(start(1), start(2), delta(1), delta(2), 0,... 
                'Color', LineColor,...
                'LineWidth', 1.2,...         % Original 2 → 1.2 (thinner line)
                'MaxHeadSize', 0.5,...       % Original 0.3 → 0.8 (larger head)
                'AutoScale', 'off',...
                'AlignVertexCenters', 'on');
        end
    end
    
    % Add cell numbering if enabled
    if numbering > 0
        tx = 5; ty = -5;
        t = text(x(i)+tx, y(i)+ty, num2str(i));
        set(t, 'Color',[1, 0 ,0], 'FontSize', 6)
    end
end

drawnow
hold off


figure(gcf) % 
    set(gcf, 'Color', 'none', 'InvertHardcopy', 'off'); % 
    print('correlation_arrows.svg', '-dsvg', '-vector'); % 

%% Store visualization data back to params
params.visualization = struct(...
    'x', x, ...
    'y', y, ...
    'ColorMap', ColorMap);