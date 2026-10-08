function plotChosenCellsonMaxZimage(params)


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
    
    numCells = size(params.coordinatesData,2)/2;
    x = zeros(numCells,1);
    y = zeros(numCells,1);
    

    for i = params.firstCellTrace : params.lastCellTrace
        x(i) = params.coordinatesData(1,(i*2)-1)/params.scaling;
        y(i) = params.coordinatesData(1,(i*2))/params.scaling;
        l(i) = i;
    end


    figure('Name','Chosen Cells')
    imshow(params.maxZImages);
    hold on
    b = plot(x,y,'g.');
    set(b,'MarkerSize',15);
    

    tx = 10; ty = -10;  
    for i = 1 : size(l,2)
        t = text(x(i)+tx, y(i)+ty, num2str(round(l(i))));
        set(t, 'Color',[1, 0 ,0], 'FontSize', 15)             
    end
    
    drawnow
    

    saveas(gcf, fullfile(params.currentFolder, [params.output '_chosen_cells.fig']));
end