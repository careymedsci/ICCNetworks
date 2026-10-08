function plotSimpleTraces(params)
    
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



    required_fields = {'plot_cells_traces', 'firstCellTrace', 'lastCellTrace', ...
                      'allTraces', 'peaksInfo', 'secondsPerFrame', 'currentFolder'};
    for i = 1:length(required_fields)
        if ~isfield(params, required_fields{i})
            error('Missing required field: %s', required_fields{i});
        end
    end
    
    if params.plot_cells_traces <= 0
        return;
    end
    
    % 创建保存图像的文件夹
    savePath = fullfile(params.currentFolder, 'figures');
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    
    for trace = params.firstCellTrace : params.lastCellTrace
        fig = figure('Name', sprintf('Simple Trace - Cell %d', trace));
        plot(params.allTraces(:,trace))
        ylim([-100 100]);
        hold on
        % 使用params.peaksInfo而不是peaksInfo
        plot(params.peaksInfo(trace).peakLocations, ...
             params.peaksInfo(trace).peakHeight, 'o')
        hold off
        
        tit = ['Cell Number ' num2str(trace)];
        title(tit, 'Color', 'blue');
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', round(xt*params.secondsPerFrame))
        
        % 保存图像
        saveas(fig, fullfile(savePath, sprintf('Cell_%03d_simple.png', trace)));
    end
end