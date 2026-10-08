 
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


%%


function [params, analysis_summary] = visualize3DCellTraces(params, analysis_summary)
% visualize3DCellTraces - Create 3D visualization of cell activity traces
%
% Required input variables from params:
%   currentFolder - Save path for outputs
%   allTraces - Cell activity traces data
%   numRois - Number of ROIs/cells
%   hubs - Hub cells indicator array
%   trigger - Trigger cells indicator array
%   periodicity - Periodic cells indicator array
%   selected_ticks - Selected cell numbers for x-axis ticks

% Extract required variables from params
savePath = params.currentFolder;
hubs = params.hubs;
trigger = params.trigger;
periodicity = params.periodicity;
% 加载params结构体
loaded_data = load('all_ver_0.mat', 'params'); 
% 提取allTraces字段
allTraces = loaded_data.params.allTraces; 
load('all_for_3D.mat');
%% SECTION 1: Initialize parameters
total_duration = 10;      % Total duration in minutes
num_frames = size(allTraces,1);
cells = 1:params.numRois;
time_vector = linspace(0, total_duration, num_frames);
time_mask = time_vector >= 1.5 & time_vector <= 5.5;

%% SECTION 2: Create figure
fig = figure('Position', [200 200 800 600],...
            'Color', 'w',...
            'Name', '3D Cell Traces');
hold on

%% SECTION 3: Plot cell traces
for i = 1:length(cells)
    % Extract trace data
    trace = cells(i);
    z_data_full = allTraces(:,trace);
    z_data = z_data_full(time_mask);
    y_data = time_vector(time_mask)';
    x_data = trace * ones(size(z_data));
    
   % Set style parameters based on cell type
    if hubs(trace)
        color = [0.75, 0.75, 0.75];   % Blue - Hub cells
        z_scale = 0.5;
        line_width = 0.0001;
    elseif trigger(trace)
        color = [0.75, 0.75, 0.75];   % Green - Trigger cells
        z_scale = 0.5;
        line_width = 0.01;
    elseif periodicity(trace)
        color = [0.64 0.08 0.18];   % Red - Periodic cells
        z_scale = 1.2;
        line_width = 2;
    else
         if rand < 0.45
        color = [1, 1, 1];       % Background-like color (white)
        line_width = 0.1;      % Nearly invisible
        z_scale = 0.1;
          else
    color = [0.59 0.44 0.25]; % Muted earthy color - Other cells
    z_scale = 0.5;
    line_width = 0.1;
          end
    end
    
    % Plot scaled trace
    z_data_scaled = z_data * z_scale;
    plot3(x_data, y_data, z_data_scaled,...
         'LineWidth', line_width,...
         'Color', color,...
         'DisplayName', ['Cell ' num2str(trace)]);
end
view(60, 60)
%% SECTION 4: Configure view and axes
view(-37.5, 25)
ax = gca;
set(ax, 'XColor','none', 'YColor','none', 'ZColor','none')

% X-axis label system
ax.XAxis.TickLength = [0 0];
ax.XLabel.String = 'Cell Number';
ax.XLabel.FontSize = 12;
ax.XLabel.Color = [0.2 0.2 0.2];


%% SECTION 5: Create legend
leg_entries = gobjects(4,1);

% Plot legend example lines
leg_entries(1) = plot3([NaN],[NaN],[NaN],...
    'Color', [0.4 0.4 0.4], 'LineWidth',0.8, 'DisplayName','Non-Coactive Cells');
leg_entries(2) = plot3([NaN],[NaN],[NaN],...
    'Color', [0 1 0], 'LineWidth',1.1, 'DisplayName','Trigger Cells');
leg_entries(3) = plot3([NaN],[NaN],[NaN],...
    'Color', [0 0 1], 'LineWidth',1.1, 'DisplayName','Hub Cells');
leg_entries(4) = plot3([NaN],[NaN],[NaN],...
    'Color', [1 0 0], 'LineWidth',1.3, 'DisplayName','Periodic Cells');

% Configure legend
legend(leg_entries,...
    'Location', 'northeastoutside',...
    'Box', 'off',...
    'FontSize', 10,...
    'Color', [0.95 0.95 0.95],...
    'EdgeColor', [0.8 0.8 0.8]);
set(legend, 'Position', [0.78 0.65 0.15 0.25]);

% Final view settings
view(-37.5, 25)
box off
grid off

%% Store visualization data back to params
params.visualization3D = struct(...
    'time_vector', time_vector,...
    'time_mask', time_mask,...
    'figure_handle', fig);

% Export important variables to workspace
assignin('base', 'cell_traces_time_vector', time_vector);
assignin('base', 'cell_traces_time_mask', time_mask);

end