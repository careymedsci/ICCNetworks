
 
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

%%
% please mind that the C matrix was caculated under the rule that all cells
% has at least 4 peaks.

function [params, analysis_summary] = classifyCellsByNeighbors(params, analysis_summary)
% Extract required variables from params
cacul_corr = params.cacul_corr;
classify_neighbor_of_connections = params.classify_neighbor_of_connections;
find_periodicities = params.find_periodicities;
C = params.C;
threshcorr = params.threshcorr;
periodicity = params.periodicity;
output = params.output;
numRois = params.numRois;
tracesleft = params.tracesleft;
%% 
if cacul_corr > 0 && classify_neighbor_of_connections > 0
    if find_periodicities == 0
        periodicity = 0;
    end
    Cnum = zeros(size(C, 1), 1);
    numPeriodic = zeros(size(C, 1), 1);
    excludedCells = zeros(size(C, 1));

    for i = 1:size(C, 1)
        for ii = 1:size(C, 1)
            if C(i, ii) >= threshcorr && excludedCells(i, ii) == 0
                excludedCells(ii, i) = 1;
                Cnum(i, 1) = Cnum(i, 1) + 1;
                Cnum(ii, 1) = Cnum(ii, 1) + 1;
                if periodicity(i) == 1
                    numPeriodic(i, 1) = numPeriodic(i, 1) + 1;
                end
                if periodicity(ii) == 1
                    numPeriodic(ii, 1) = numPeriodic(ii, 1) + 1;
                end
            end
        end
    end

    % Probability distribution P(k): 
   probabilityDistribution = zeros(max(Cnum), 4);
for j = 1:max(Cnum)
    m = sum(Cnum == j);
    m_perio = sum(numPeriodic(Cnum == j) == j);  
    ratio = m / numRois;  
    probabilityDistribution(j, :) = [j, m, m_perio, ratio];
end


outfile = [char(output), '_probability_distribution.txt'];
writematrix(probabilityDistribution, outfile);

params.probabilityDistribution = probabilityDistribution;


figure('Name','Probability distribution');
h = histogram(Cnum, 'BinWidth', 1, 'BinLimits', [-0.5, max(Cnum)+0.5]);
grid on;
xlabel('Number of connections/cell');
ylabel('Cell Count');

counts = h.Values;
edges = h.BinEdges;
binCenters = edges(1:end-1) + diff(edges)/2;


for k = 1:length(counts)
    if counts(k) > 0
      
        ratioStr = sprintf('%.2f%%', (counts(k)/numRois)*100);
        text(binCenters(k), counts(k), {
            num2str(counts(k));   
            ratioStr              
            }, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 9);
    end

    
end

drawnow;
savefig([output, '_propablility distribution.fig']);

end
numCells = size(params.allTraces, 2);
excludedNodes = zeros(1, numCells); % initialize
for z=1 : numCells
    if sum(params.peakLocations(:,z)) < 1 
       excludedNodes(z) = 1;
    end
end

params.excludedNodes = excludedNodes; % store for future use
ca_cell_number = numRois - sum(excludedNodes);
ca_ratio = ca_cell_number/numRois *100; 
active_cells = tracesleft/numRois*100;
analysis_summary(15,1) = "number of Ca2+ transient cells";
analysis_summary(15,2) = ca_cell_number;
analysis_summary(12,1) = "percentage of Ca2+ transient cells";
analysis_summary(12,2) = ca_ratio;
analysis_summary(1,1) = "percentage of cells with at least 4 peaks(ca2+ active cells)";
analysis_summary(1,2) = active_cells;


% tracesleft = cells ≥ 4 peaks = active cells
% numRois - sum(Cnum == 0)： all co-active cells, that is cells (have at
% least 4 peaks) at least have coumunication with 1 other cells
% co_ratio: percentage of co-active cells in tracesleft
co_ratio = (numRois - sum(Cnum == 0) )/tracesleft *100; 

analysis_summary(8,3) = "percentage of co-active cells";
analysis_summary(8,4) = co_ratio;

filename = sprintf('all_for_analysis_%d.mat', params.mode);
save(filename);

end