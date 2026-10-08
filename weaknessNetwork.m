
 
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
% By Xiao Liu and Thomas Broggini;
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany
% Xiangyang No.1 people's Hospital, China
% Also an python version coded 
% Xiao.Liu@stud.uni-frankfurt.de
% 2024.10.09


function [params, analysis_summary] = weaknessNetwork(params, analysis_summary)

% Required input variables from params:
plot_networkx = params.plot_networkx;      % Flag for network plotting
cacul_corr = params.cacul_corr;           % Flag for correlation calculation
hubs = params.hubs;                       % Hub cells array
periodicity = params.periodicity;         % Periodicity array
CC = params.CC;                           % Connectivity matrix
find_periodicities = params.find_periodicities; % Flag for periodicity analysis
output = params.output;                   % Output file path
node = params.define_hub_cells;

%% Network visualization parameters
edge = 1.65;   % Edge width for visualization

%% Original code starts here with exact structure and order
if plot_networkx > 0 && cacul_corr > 0
    %% SECTION 1: Initialize hub and periodic cell indices
    % Get hub cell indices

    hh = NaN(1,size(hubs,2));
    for z=1 : size(hubs,2)
        if hubs(z) == 1
            hh(z) = z;
        end
    end
    
    % Get periodic cell indices

    pp = NaN(1,size(periodicity,2));
    for z=1 : size(periodicity,2)
        if periodicity(z) == 1
            pp(z) = z;
        end
    end

    %% SECTION 2: Process connectivity matrix
    % Remove unconnected cells from network

    stop = size(CC,1);
    CC_2=CC;
    kk=0;
    for k=1 : stop
        kk=kk+1;
        if sum(CC(:,k),'omitnan') == 0
            CC_2(kk,:)=[];
            CC_2(:,kk)=[];
            hh(1,kk:end)=hh(1,kk:end)-1;
            hh(:,kk)=[];
            pp(1,kk:end)=pp(1,kk:end)-1;
            pp(:,kk)=[];
            kk=kk-1;
        end
    end

    % Clean up indices and convert to sparse matrix

    pp=pp(~isnan(pp));
    pp_size = size(pp(1,:));
    hh=hh(~isnan(hh));
    CC_2 = sparse(CC_2);

    %% SECTION 3: Visualize original network with hub cells
    figure ('Name','Hub cells in Orignal network');
    G=graph(CC_2);
    colo=[0 0.5 0.6];
    g=plot(G, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge,'NodeColor',colo);
    layout(g, 'force','UseGravity',true)
    highlight(g, hh(1,:),'NodeColor','red','MarkerSize', node);
    drawnow
    title('Hub cells in Orignal network');
    outputFigure = sprintf('%s_network_hubs.fig', output);
    savefig(outputFigure);

  
    figure('Name','Periodic cells in Original Network');
    G=graph(CC_2);
    g=plot(G, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge,'NodeColor',colo);
    layout(g, 'force','UseGravity',true)
    highlight(g, pp(1,:),'NodeColor','#00ff00', 'MarkerSize', node*1.2);
    drawnow
    title('Periodic cells in Original Network');
    outputFigure = sprintf('%s_network_perio.fig', output);
    savefig(outputFigure);
    
    
    figure('Name','Full Founction Network Architecture ');
    G=graph(CC_2);
    g=plot(G, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge, 'NodeColor', colo);
    layout(g, 'force','UseGravity',true)
    highlight(g, pp(1,:),'NodeColor',[0, 0.95, 0], 'MarkerSize', node*2.2);
    hold
    %g2=plot(G, 'NodeLabel',{}, 'Marker','none', 'Edgecolor','none');
    g2=plot(G, 'NodeLabel',{}, 'NodeColor', [0 0.447 0.741],'Edgecolor',colo, 'MarkerSize', node);
    layout(g2, 'force','UseGravity',true)
    highlight(g2, hh(1,:),'NodeColor',[0.64, 0.08, 0.18],'MarkerSize', node*1);
    drawnow
    title('Full Founction Network Architecture ');
    outputFigure = sprintf('%s_network_hubsandperio.fig', output);
    savefig(outputFigure);



    %% SECTION 4: Analyze network without periodic cells
    if find_periodicities > 0
    periodicity2=periodicity;
    CC_notperio = full(CC_2);
    for k=1 : size(CC_notperio,1)
        for l=1 : size(CC_notperio,1)
            if ismember(k,pp) == 1
                periodicity2(k) = 1;
            else
                periodicity2(k) = 0;
            end
            if ismember(l,pp) == 1
                periodicity2(l) = 1;
            else
                periodicity2(l) = 0;
            end
            if periodicity2(k) == 1 || periodicity2(l) == 1
                CC_notperio(k,l) = 0;
                CC_notperio(l,k) = 0;
            end
        end
    end
    
  CC_notperio = sparse(CC_notperio); % *************************🚩❗🗝️
    figure('Name','After removal of Periodic cells');
    H=graph(CC_notperio);
    h=plot(H, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge);
layout(h, 'force', 'UseGravity', true);

    % layout(h, 'force', 'UseGravity', true);
    highlight(h, pp(1,:),'NodeColor','#646464');
    drawnow
    title('After removal of Periodic cells');
    outputFigure = sprintf('%s_network_notperio.fig', output);
    savefig(outputFigure);


    hubs2=hubs;
    CC_nothub = full(CC_2);
    for k=1 : size(CC_nothub,1)
        for l=1 : size(CC_nothub,1)
            if ismember(k,hh) == 1
                hubs2(k) = 1;
            else
                hubs2(k) = 0;
            end
            if ismember(l,hh) == 1
                hubs2(l) = 1;
            else
                hubs2(l) = 0;
            end
            if hubs2(k) == 1 || hubs2(l) == 1
                CC_nothub(k,l) = 0;
                CC_nothub(l,k) = 0;
            end
        end
    end
    
    CC_nothub = sparse(CC_nothub);
    figure('Name','After removal of Hub cells'); % ********************🚩❗🗝️
    H=graph(CC_nothub);
    h=plot(H, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge);
    layout(h, 'force','UseGravity',true)
    highlight(h, hh(1,:),'NodeColor','#646464');
    drawnow
    title('After removal of Hub cells');
    outputFigure = sprintf('%s_network_nothub.fig', output);
    savefig(outputFigure);

   
    rando=periodicity;
    CC_rando = full(CC_2);
    
    if size(CC_rando,2) >= sum(periodicity,'all','omitnan')
        r = randperm(size(CC_rando,2),pp_size(2));
        for k=1 : size(CC_2,1)
            for l=1 : size(CC_2,1)
                if ismember(k,r) == 1
                    rando(k) = 1;
                else
                    rando(k) = 0;
                end
                if ismember(l,r) == 1
                    rando(l) = 1;
                else
                    rando(l) = 0;
                end
                if rando(k) == 1 || rando(l) == 1
                    CC_rando(k,l) = 0;
                    CC_rando(l,k) = 0;
                end
            end
        end

        [~,binsizes_rando] = conncomp(graph(CC_rando));
        largestcluster_rand = max(binsizes_rando);
        
        CC_rando = sparse(CC_rando); %********************🗝️❗🚩
        % figure('Name','After removal of random cells');
        H=graph(CC_rando);
        h=plot(H, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge);
        layout(h, 'force','UseGravity',true);
        highlight(h, r,'NodeColor','#646464');

        drawnow
        title('After removal of random cells');
        outputFigure = sprintf('%s_network_notrand_%.0f.fig', output,largestcluster_rand);
        savefig(outputFigure);
    end
    
        %% SECTION 5: Calculate and store network metrics
    [~,binsizes] = conncomp(graph(CC_2));
    largestcluster = max(binsizes);
    [~,binsizes_notperio] = conncomp(graph(CC_notperio));
    largestcluster_perio = max(binsizes_notperio);
    [~,binsizes_rando] = conncomp(graph(CC_rando));
    largestcluster_rand = max(binsizes_rando);
    
    analysis_summary(3,3) = "Size of largest network";
    if sum(largestcluster,'all','omitnan') > 0
        analysis_summary(3,4) = largestcluster;
    else
        analysis_summary(3,4) = 0;
    end
    
    analysis_summary(4,3) = "Size of largest network after removing all periodic cells";
    if sum(largestcluster_perio,'all','omitnan') > 0
        analysis_summary(4,4) = largestcluster_perio;
    else
        analysis_summary(4,4) = 0;
    end
    
    analysis_summary(5,3) = "Size of largest network after removing random cells in the same number as periodic cells";
    if sum(largestcluster_rand,'all','omitnan') > 0
        analysis_summary(5,4) = largestcluster_rand;
    else
        analysis_summary(5,4) = 0;
    end
    

end
    
    writematrix (analysis_summary, output + "_analysis-summary.txt");
    diary off
end

%% Store modified/new variables back to params
params.CC_2 = CC_2;
params.CC_notperio = CC_notperio;
params.CC_nothub = CC_nothub;
params.CC_rando = CC_rando;
params.largestcluster = largestcluster;
params.largestcluster_perio = largestcluster_perio;
params.largestcluster_rand = largestcluster_rand;


end
