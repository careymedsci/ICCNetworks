% 作者：刘晓；法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% Xiao Liu
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


function [params, analysis_summary] = analyze_peakwidths_of_p_and_active(params, analysis_summary)
% ANALYZE_PEAK_WIDTHS 计算周期细胞和峰值大于3的细胞的平均峰宽（单位：秒）
% ANALYZE_PEAK_WIDTHS computes the mean peak width (in seconds) for:
% - periodic cells
% - cells with more than 3 peaks

% 找到周期性细胞索引（周期性=1）
% Find indices of periodic cells (periodicity == 1)
periodic_cell_indices = find(params.periodicity == 1);

% 初始化变量
% Initialize width storage
periodic_widths = [];
all_cells_with_min3peaks = [];

peakLocations = params.peakLocations;

% 遍历每个 ROI（细胞）
% Loop through each ROI (cell)
for i = 1:params.numRois
    % 判断该细胞是否有超过3个峰值
    % Check if the cell has more than 2 peaks
    if sum(peakLocations(:, i)) > 2
        % 获取该细胞所有非零峰的位置
        % Get all non-zero peak positions
        cell_peaks = find(peakLocations(:, i));
        
        if ~isempty(cell_peaks)
            % 获取对应的峰宽度
            % Get corresponding peak widths
            cell_widths = params.widths(cell_peaks, i);
            
            % 添加到所有符合条件细胞列表中
            % Add to all-cells list
            all_cells_with_min3peaks = [all_cells_with_min3peaks; cell_widths];
            
            % 如果该细胞是周期性细胞，添加到周期细胞宽度列表
            % If cell is periodic, add to periodic list
            if ismember(i, periodic_cell_indices)
                periodic_widths = [periodic_widths; cell_widths];
            end
        end
    end
end

% 计算平均宽度（秒）
% Convert to seconds and compute mean widths
if ~isempty(periodic_widths)
    meanWidth_periodic = mean(periodic_widths) * params.secondsPerFrame;
else
    meanWidth_periodic = NaN;
end

if ~isempty(all_cells_with_min3peaks)
    meanWidth_min3peaks = mean(all_cells_with_min3peaks) * params.secondsPerFrame;
else
    meanWidth_min3peaks = NaN;
end

% 保存结果到输出和结构体中
% Store results to output and params struct
analysis_summary(end+1, 1) = "Mean peak width of periodic cells (at least 3 peaks; seconds)";
analysis_summary(end, 2) = meanWidth_periodic;

analysis_summary(end+1, 1) = "Mean peak width of cells (at least 3 peaks; seconds)";
analysis_summary(end, 2) = meanWidth_min3peaks;

params.meanWidth_periodic    = meanWidth_periodic;
params.meanWidth_min3peaks   = meanWidth_min3peaks;

save('_all_result.mat', 'params','analysis_summary' );
end
