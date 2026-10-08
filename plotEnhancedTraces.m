function plotEnhancedTraces(params)

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
% 2024.10.09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if params.plot_cells_traces <= 0
        return;
    end
    
    % 创建保存图像的文件夹
    savePath = fullfile(params.currentFolder, 'figures');
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    
    for trace = params.firstCellTrace : params.lastCellTrace
        fig = figure('Position', [100 100 1000 400], ...
                    'Color', [0.95 0.95 0.95], ...
                    'Name', sprintf('Enhanced Trace - Cell %d', trace));
        
        current_trace = params.allTraces(:,trace);
        
        % 绘制主信号
        plot(current_trace, ...
            'LineWidth', 1, ...
            'Color', [0.2 0.5 0.8], ...
            'DisplayName', 'Raw Signal');
        hold on
        
        % 绘制峰值点
        scatter(params.peaksInfo(trace).peakLocations, ...
                params.peaksInfo(trace).peakHeight, ...
                10, ...
                'Marker', 'o', ...
                'MarkerEdgeColor', [0.8 0.1 0.1], ...
                'MarkerFaceColor', 'none', ...
                'LineWidth', 1.2, ...
                'DisplayName', 'Detected Peaks');
        
        % 设置图形属性
        margin = 0.4 * (max(current_trace) - min(current_trace));
        ylim([min(current_trace)-margin, max(current_trace)+margin]);
        
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', round(xt*params.secondsPerFrame))
        xlabel(['Time (s), \Deltat = ' num2str(params.secondsPerFrame) 's/frame'])
        ylabel('Signal Intensity (a.u.)')
        
        tit = sprintf('Cell %d | %d Peaks Detected', ...
                     trace, length(params.peaksInfo(trace).peakLocations));
        title(tit, 'Color', 'b', 'FontSize', 12)
        legend('Location', 'best')
        grid on
        set(gca, 'GridAlpha', 0.3)
        
        hold off
        
        % 保存图像
        saveas(fig, fullfile(savePath, sprintf('Cell_%03d_enhanced.png', trace)));
    end
end