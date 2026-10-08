

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
% 2024.10.09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [params, analysis_summary] = network_theory_analysis(params, analysis_summary)
    % network_theory_analysis - Perform network theory analysis on correlation matrix
    % 网络理论分析函数，基于相关矩阵计算小世界网络指标
    %
    % 输入参数 params 是结构体，包含字段:
    %   do_networktheory - 是否执行网络分析 (0/1)
    %   cacul_corr       - 是否计算相关矩阵 (0/1)
    %   C                - 相关矩阵，方阵 (n×n)
    %   threshcorr       - 相关系数阈值，二值化用
    %   output           - 输出文件路径前缀（字符串）
   

    do_networktheory = params.do_networktheory;
    cacul_corr = params.cacul_corr;
    C = params.C;
    threshcorr = params.threshcorr;
    output = params.output;

    % 只有当计算相关矩阵且启用网络分析时才执行
    % Execute only if cacul_corr > 0 and do_networktheory > 0
    if cacul_corr > 0 && do_networktheory > 0
        n = size(C, 1);
        CC = zeros(n);

        % 构建无向二值连接矩阵CC
        % Build undirected binary adjacency matrix CC
        for i = 1:n
            for j = i+1:n
                if C(i,j) > threshcorr
                    CC(i,j) = 1;
                    CC(j,i) = 1;
                end
            end
        end
        
        % 计算至少有一个连接的节点数量
        % Count nodes with at least one connection
        conncells = sum(sum(CC,1) > 0);

        % 构建随机网络：Erdős–Rényi 随机图
        % Construct random network using Erdős–Rényi model (G(n,p))
        edge_count = sum(CC(:)) / 2; % edges witout direction
        p = edge_count / (conncells * (conncells - 1) / 2); % probroality of edges

        % 只考虑 connells 节点
        % only conncells knots
        randCC_small = rand(conncells) < p; 
        randCC_small = triu(randCC_small,1); % 去除对角线 / removal useless knots 
        randCC_small = randCC_small + randCC_small'; 
        randCC = zeros(n);
        randCC(1:conncells, 1:conncells) = randCC_small;

        % 计算原始网络的平均最短路径
        % Calculate mean shortest path in original network
        CCgraph = graph(CC);
        SP = distances(CCgraph); % two side
        SP(SP == Inf) = NaN;  % 不连通的路径设为NaN / unreachable set to NaN
        meanSP = mean(SP/2, 'all', 'omitnan');

        % 计算随机网络的平均最短路径
        % Calculate mean shortest path in random network
        randCCgraph = graph(randCC);
        randSP = distances(randCCgraph);
        randSP(randSP == Inf) = NaN;
        randmeanSP = mean(randSP, 'all', 'omitnan');

        % 保存平均最短路径结果 / save mean shortest path results
        analysis_summary(23,1) = "mean shortest path (original network)";
        analysis_summary(23,2) = meanSP;
        analysis_summary(24,1) = "mean shortest path (random network)";
        analysis_summary(24,2) = randmeanSP;

        % 计算原始网络的平均聚类系数
        % Calculate mean clustering coefficient of original network
        sparseCC = sparse(CC);
        Cl = clustercoeffs(sparseCC);
        meanCL = mean(Cl);

        % 计算随机网络的平均聚类系数
        % Calculate mean clustering coefficient of random network
        randsparseCC = sparse(randCC);
        randCl = clustercoeffs(randsparseCC);
        randmeanCL = mean(randCl);

        % 保存平均聚类系数结果 / save clustering coefficient results
        analysis_summary(25,1) = "mean clustering coefficient (original network)";
        analysis_summary(25,2) = meanCL;
        analysis_summary(26,1) = "mean clustering coefficient (random network)";
        analysis_summary(26,2) = randmeanCL;

        % 保存稀疏矩阵数据，格式：列号 行号 值
        % Save sparse adjacency matrix as [col row value]
        [i,j,val] = find(sparseCC);
        dlmwrite(output + "_sparseCC.dat", [j i val], 'delimiter', ' ', 'newline', 'pc');
    else
        warning('Parameters indicate no network analysis to perform.');
    end
    params.CC = CC;

    filename = sprintf('all_verVisulization_%d.mat', params.mode);
    save(filename);
end
