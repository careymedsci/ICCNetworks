

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



% 作者：Thomas Broggini and 刘晓； 法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Thomas Broggini and Xiao Liu
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany
% Xiangyang No.1 people's Hospital, China
% Also an python version coded 
% Xiao.Liu@stud.uni-frankfurt.de
% 2024.10.09

% processCorrelations
%
% Perform cell-to-cell correlation analysis with optional null-model
% randomization, spatial distance constraints, temporal shifting, and
% directionality estimation.
%
% 执行细胞间相关性分析，包括：
% - 随机对照组构建（null model）
% - 空间距离限制
% - 时间延迟扫描
% - 信号传播方向性估计
%
% ------------------------------------------------
% Inputs:
%   params            - Parameter structure containing:
%                       包含以下字段的参数结构体：
%
%       .peakLocations    Binary peak event matrix (frames × ROIs)
%                          峰值事件二值矩阵（时间帧 × 细胞ROI）
%
%       .allTraces        Original fluorescence traces
%                          原始荧光信号
%
%       .coordinatesData ROI spatial coordinates (x1,y1,x2,y2,...)
%                          ROI空间坐标
%
%       .secondsPerFrame Sampling interval (s)
%                          采样时间分辨率（秒/帧）
%
%       .pk_num_define    Minimum peak count threshold
%                          最小峰次数阈值
%
%       .dcutoff          Spatial distance cutoff
%                          最大允许相关距离
%
%       .threshcorr       Correlation significance threshold
%                          相关性阈值
%
%       .mode             Analysis mode selector
%                          分析模式选择
%
%       .maxspeed         Maximum propagation speed (µm/s)
%                          最大传播速度
%
%       .minspeed         Minimum propagation speed (µm/s)
%                          最小传播速度
%
%       .peakTolerance_s  Peak temporal tolerance window
%                          峰值时间扩展窗口
%
% ------------------------------------------------
% Input / Output:
%   analysis_summary  Summary statistics table
%                     分析统计汇总表
%
% ------------------------------------------------
% Outputs:
%   params (updated)  Adds:
%
%       .C        Correlation matrix (NxN)
%                  细胞间最大相关矩阵
%
%       .S        Time shift matrix (frames)
%                  最佳时间延迟矩阵
%
%       .Dir      Directionality matrix
%                  传播方向性矩阵
%
%       .D        Distance matrix
%                  细胞间距离矩阵
%
%       .CvsD     Correlation vs distance pairs
%                  相关性-距离散点数据
%
% ================================================================
%%


function [params, analysis_summary] = processCorrelations(params, analysis_summary)
    % 处理细胞间相关性分析，包括控制组处理和相关性计算
    % Input 输入:
    %   params - 包含分析参数和数据的结构体
    %   analysis_summary - 分析结果摘要
    % Output 输出:
    %   params - 更新后的参数结构体
    %   analysis_summary - 更新后的分析摘要

    % ========================================== SECTION 1: mode setting. 控制组数据处理
    %
    % Purpose:
    %   Generate randomized surrogate peak trains to construct
    %   statistical null distributions.
    %
    % 目的：
    %   通过时间扰动、局部打乱、泊松重采样等方式
    %   构造无空间同步结构的对照数据（零假设）
    %
    % Mode definitions:
    %
    % Mode 1:
    %   Circular shift + flip + global permutation
    %   循环平移 + 翻转 + 全局打乱
    %
    % Mode 2:
    %   Poisson resampling preserving mean firing rate
    %   泊松重采样，保持平均事件率
    %
    % Mode 3:
    %   Block-wise shuffle with mild random shift
    %   分段扰动 + 小时间偏移


    if params.mode == 1
           randomPeakLocations = zeros(size(params.peakLocations));
        [nF, nR] = size(params.peakLocations);
    
        for roi = 1:nR
            tracep = params.peakLocations(:, roi);
    
            for shift_round = 1:randi([3, 5]) 
                shift_amount = randi([round(nF/4), round(nF/2)]);  
                tracep = circshift(tracep, shift_amount);
            end
    
          
            for flip_round = 1:2
                if rand() > 0.6
                    tracep = flipud(tracep);
                end
            end
            jitter_idx = randperm(nF);
            tracep = tracep(jitter_idx);
            randomPeakLocations(:, roi) = tracep;
        end
    elseif params.mode == 2
            [nF, nR] = size(params.peakLocations);
                lambda = mean(params.peakLocations, 1);
                lambda = lambda(randperm(nR));
                lambda = lambda * 0.6;
                randomPeakLocations = poissrnd(repmat(lambda, nF, 1));
                 

                for roi = 1:nR
                    trace = randomPeakLocations(:, roi);
            
                    trace = trace(randperm(numel(trace)));
            
                    num_blocks = randi([4, 8]);
                    block_edges = round(linspace(1, nF+1, num_blocks+1));
                    for b = 1:num_blocks
                        idx1 = block_edges(b);
                        idx2 = block_edges(b+1)-1;
                        if idx2 > idx1
                            local_idx = idx1:idx2;
                            trace(local_idx) = trace(local_idx(randperm(numel(local_idx))));
                        end
                    end
                    shift_amount = randi([1, round(nF/5)]);
                    trace = circshift(trace, shift_amount);

                    randomPeakLocations(:, roi) = reshape(trace, nF, 1);
                end
                    for f = 1:nF
                    randomPeakLocations(f, :) = randomPeakLocations(f, randperm(nR));
                    end
                params.peakLocations = randomPeakLocations;
    elseif params.mode == 3
            [nF, nR] = size(params.peakLocations);
            randomPeakLocations = zeros(nF, nR);
                    
            for roi = 1:nR
                trace = params.peakLocations(:, roi);
            
                % 全局时间打乱（用 numel 替代 length）
                trace = trace(randperm(numel(trace)));
                trace = reshape(trace, nF, 1);  % 强制成列向量
            
                % 分段局部扰动
                num_blocks = randi([3, 6]); 
                block_edges = round(linspace(1, nF+1, num_blocks+1));
                for b = 1:num_blocks
                    idx1 = block_edges(b);
                    idx2 = block_edges(b+1)-1;
                    if idx2 > idx1
                        local_idx = idx1:idx2;
                        trace(local_idx) = trace(local_idx(randperm(numel(local_idx))));
                    end
                end
            
                % 随机小 shift
                if rand() > 0.3
                    shift_amount = randi([1, round(nF/8)]);
                    trace = circshift(trace, shift_amount);
                end
            
                randomPeakLocations(:, roi) = trace;
            end
            
    end
    if params.mode == 1 || params.mode == 2 || params.mode == 3
    params.peakLocations = randomPeakLocations;
    end
    
 % ========================================= SECTION 2: 细胞范围限制处理
    
    % SECTION 2: ROI subset selection
    %
    % Purpose:
    %   Restrict analysis to selected ROI index range
    %
    % 目的：
    %   允许用户只分析部分细胞子集

 
 % 如果需要限制细胞范围
    if params.do_celllimit > 0
        params.allTraces = params.allTraces(:,params.firstcell:params.lastcell);
        params.coordinatesData = params.coordinatesData(1,(params.firstcell*2)-1:params.lastcell*2);
        params.peakLocations = params.peakLocations(:,params.firstcell:params.lastcell);
    end

    % ========================= SECTION 3: 信号扩散处理（Spread Processing）
    
    % SECTION 2: Peak temporal spreading
    %
    % Purpose:
    %   Expand binary peak events in time to account for
    %   biological calcium transient duration and temporal jitter.
    %
    % 目的：
    %   对峰值事件进行时间扩展，模拟Ca2+信号持续时间，
    %   提高相关性计算的鲁棒性


    params.rawTraces = params.allTraces;
    
    % 创建二值化数据并进行时间扩展
    params.allTraces = zeros(size(params.peakLocations,1),size(params.peakLocations,2));
    for g=1 : size(params.peakLocations,1)
        for h=1 : size(params.peakLocations,2)
            if params.peakLocations(g,h) == 1
                params.allTraces(g,h)=1;
                % 向前后扩展信号
                for j=1 : (params.peakTolerance_s/params.secondsPerFrame)/2
                    if g+j <= size(params.peakLocations,1)
                        params.allTraces(g+j,h) = 1;
                    end
                    if g-j > 0
                        params.allTraces(g-j,h) = 1;
                    end
                end
            end
        end
    end

    hold off

   
        % ================== SECTION 4: 点图绘制(就是那个二维图，时间像素比)

    % SECTION 4: Activity raster / dot plot visualization
    %
    % Purpose:
    %   Visualize spatiotemporal activity distribution
    %
    % 目的：
    %   绘制类似神经元 raster plot 的活动分布图

    if params.do_dotplot > 0
        [frame,cellnumber] = find(params.allTraces);
        % figure('Name','Activity Dot Plot')
        b = plot(frame,cellnumber,'k.');
        set(b,'MarkerSize',5);
        xlabel('Frames')
        ylabel('cell intensities of all detected ROIs')
        xlim([0 max(frame)])
        ylim([0 max(cellnumber)])
        drawnow
        outputValue = sprintf('%s_dotplot%d.fig',params.output,params.mode);
        savefig(fullfile(params.currentFolder, outputValue));
    end

    % ============================================== SECTION 5: 相关性计算

    
    % SECTION 5: Long-range propagation correlation analysis 
    %
    % Features:
    %   - Distance constrained correlation
    %   - Sliding temporal windows
    %   - Velocity constrained lag scanning
    %   - Directionality estimation
    %
    % 特点：
    %   - 空间距离约束
    %   - 多时间片段滑动分析
    %   - 基于传播速度限制的时间延迟扫描
    %   - 信号传播方向性计算



    if params.cacul_corr > 0
        % 显示进度消息
        
        if params.mode == 4
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55Mode 4: 远距离细胞分析
            % 基于峰值数量筛选细胞
            excludedNodes = zeros(size(params.allTraces,2),1);
            for z=1 : size(params.allTraces,2)
                % Cells with fewer peaks than pk_num_define （default 4） should not be considered
                %细胞闪烁次数小于等于3次的不算（注意该数值在本项目中默认4，如果你要更改，不能少于4次）
                if sum(params.peakLocations(:,z)) < params.pk_num_define 
                    excludedNodes(z) = 1;
                end
            end

            % 计算并记录筛选统计
            % alltraces = size(params.allTraces,2);
            % below_pk_num_define = sum(excludedNodes) / alltraces * 100;
            % tracesleft = alltraces - sum(excludedNodes);
            below_pk_num_define = sum(excludedNodes) / params.numRois * 100;
            tracesleft = params.numRois - sum(excludedNodes);
            
            analysis_summary(27,1) = "Rois with number of peaks not below pk_num_define";
            analysis_summary(27,2) = tracesleft;

            % 计算细胞间距离矩阵
            % DISTANCE: Computes distances between each point as matrix
            numCells = size(params.coordinatesData,2)/2;
            params.D = NaN(numCells);
            corrleft = 0;  % 初始化计数器
            for i=1 : numCells
                for ii=1 : numCells
                    xi = params.coordinatesData(1,(i*2)-1);
                    yi = params.coordinatesData(1,(i*2));
                    xii = params.coordinatesData(1,(ii*2)-1);
                    yii = params.coordinatesData(1,(ii*2));
                    D(i,ii) = sqrt(((xi-xii)^2)+((yi-yii)^2));
                    % 计算满足距离限制条件的细胞对数量,也必须闪4次才算（继承前面代码）
                     if excludedNodes(i) == 0 && excludedNodes(ii)==0 && params.D(i,ii) > params.dcutoff && i~=ii 
                        corrleft = corrleft + 1;
                     end
                end
            end
            % 计算距离大于dcutoff的细胞对数量（曲线的下方，不被考虑的细胞）
            corrbelowdcutoff1 = size(params.D((0<params.D)&(params.D>params.dcutoff))); 
            corrbelowdcutoff2 = corrbelowdcutoff1(1)/2;
            corrleft = corrleft/2;  % 因为每对细胞被计算了两次
            % 注意这个只是满足距离要求的细胞数
            analysis_summary(28,1) = "all possible correlations below dcutoff";
            analysis_summary(28,2) = corrbelowdcutoff2;
            % 注意这个只是满足距离和确认考虑计算相关性的细胞，不是最终相关细胞数目
            analysis_summary(29,1) = "corrleft (all possible correlations below dcutoff with number of peaks not below pk_num_define)";
            analysis_summary(29,2) = corrleft;

            % 
            length = 10;    % Length of sequences to be correlated in minutes
            steps = 5;      % the time by which the sequences are shifted
            numofseq = ceil(((size(params.allTraces,1)*params.secondsPerFrame/60)-length) / steps);
            % indicates the number of sequences to be correlated
            % （length of the video - length of the last sequence)/the time by which the sequences are shifted
            if numofseq<1
               numofseq=1;
            end
            fprintf('\n\nNumber of Sequences: %.0f \nLength of Sequence: %.0f\n\n', numofseq, length)

            % 分割序列
            for a=1 : numofseq
                start = floor(((a-1)*(steps*60/params.secondsPerFrame)))+1;
                lim = ceil(start-1+(length*60/params.secondsPerFrame));
                if a == numofseq
                    lim = size(params.allTraces,1);
                end
                di=lim-start;
                seq(a,1:di+1,:) = params.allTraces(start:lim,:);
                seqPeakLocations(a,1:di+1,:) = params.peakLocations(start:lim,:);
            end

             
            % Creating all required matrices so they don't have to change their size repeatedly in the following loops
            X = zeros(size(params.allTraces,2),size(params.allTraces,2));
            Cs = NaN(numofseq,size(params.allTraces,2),size(params.allTraces,2));
            Ss = NaN(numofseq,size(params.allTraces,2),size(params.allTraces,2));
            Ps = NaN(numofseq,size(params.allTraces,2),size(params.allTraces,2));
            C = NaN(size(params.allTraces,2));
            S = NaN(size(params.allTraces,2));
            P = NaN(size(params.allTraces,2));
            Dir = NaN(size(params.allTraces,2));
            maxPositiveCorrelation = NaN(size(params.allTraces,2));
            maxNegativeCorrelation = NaN(size(params.allTraces,2));

            maxspeed = params.maxspeed * params.secondsPerFrame; % Conversion to distance/frame
            minspeed = params.minspeed * params.secondsPerFrame; % 转换为距离/帧数
            maxshift = ceil(max(max(D))/minspeed); % To hold the size of A and pval before the loop
            A = NaN(size(params.allTraces,2),size(params.allTraces,2),(maxshift * 2)+1);
            B = NaN(size(params.allTraces,2),size(params.allTraces,2),(maxshift * 2)+1);
            pval = zeros(size(params.allTraces,2),size(params.allTraces,2),(maxshift * 2)+1);  

           
            twentypercent = round(size(params.allTraces,2)/5);
            msg = msgbox(sprintf('Progress: 0 percent'));
            % The first two for-loops (i & ii) go through the possible cell pairs
            for i=1 : size(params.allTraces,2)
                % Indicates progress in 10% steps
                if rem(i, twentypercent)==0
                    delete(msg);
                    msg=msgbox(sprintf('Progress: %.0f percent', i/twentypercent*20));
                end
                
                for ii=1 : size(params.allTraces,2)
                    
                    if i ~= ii  && X(i,ii)==0 && D(ii,i)>params.dcutoff && excludedNodes(i) == 0 && excludedNodes(ii)==0
                        minshift = fix(50/maxspeed);  % 使用转换后的速度
                        maxshift = ceil(50/minspeed);
                        
                        for step=1 : numofseq% This loop (step) repeats the correlations for each sequence (seq) and selects the highest correlation at the end
                % The third for-loop (iii) shifts:
                % It starts with no shift and then alternately shifts forward and backward in increasing width (So that the first highest value (if there are multiple equally high highest correlations) lies closer to 0 shift.
                step_of_other_vector = step;                            
                            if params.mode == 3
                                for kk = 0 : 1000
                                    step_of_other_vector = floor(1+rand*numofseq);
                                    if step_of_other_vector > numofseq
                                        step_of_other_vector = numofseq;
                                    end
                                    if step_of_other_vector ~= step
                                        break
                                    end
                                end
                            end
                            
                            nn = 0;
                            mm = 0;
                            if sum(seqPeakLocations(step,:,i)) >= params.pk_num_define  && sum(seqPeakLocations(step,:,ii)) >= params.pk_num_define
                                for iii=(minshift*2) : (maxshift*2)+1
                                    if iii>=2
                                          if rem(iii, 2) == 0
                                             shift = iii/2;
                                             shiftpn = 1; % 1 = positive shift
                                          elseif rem(iii, 2) == 1
                                            shift = -(iii-1)/2;
                                            shiftpn = 0; % 0 = negative shift
                                          end
                            
                                            movedvector = circshift(seq(step,:,i),shift);
                                            
                                            movedvector2 = movedvector.'; % exchange row with column so that corr compares the entire vector at a time
                                            othervector = seq(step_of_other_vector,:,ii);
                                            othervector2 = othervector.';
                                            A(ii,i,iii) = corr(movedvector2, othervector2); % gives the correlation and its p-value
                                            B(ii,i,iii) = shift;
                                            
                                            if shiftpn == 1
                                                nn = nn+1;
                                                correlationPositive(nn,step,ii,i) = A(ii,i,iii);
                                            end
                                            if shiftpn == 0
                                                mm = mm+1;
                                                correlationNegative(nn,step,ii,i) = A(ii,i,iii);
                                            end
                                        end
                                    end
                    [Cs(step,ii,i), maxiii] = max(A(ii,i,:)); % The highest correlation is assigned to Cs. The location (iii) of this correlation is assigned to S
                    Cs(step,i,ii) = Cs(step,ii,i);
                    Ss(step,i,ii) = B(ii,i,maxiii);
                end
            end

            [C(i,ii),loc] = nanmax(Cs(:,ii,i));
            S(i,ii) = Ss(loc,i,ii);
            if C(i,ii) >= params.threshcorr
                maxPositiveCorrelation(i,ii) = nanmax(nanmax(correlationPositive(:,:,ii,i)));
                maxNegativeCorrelation(i,ii) = nanmax(nanmax(correlationNegative(:,:,ii,i)));
                %Dir(i,ii)=(corr_pos_max(i,ii)-corr_neg_max(i,ii));
                Dir(i,ii) = (maxPositiveCorrelation(i,ii)-maxNegativeCorrelation(i,ii))/C(i,ii)*100;
            end

            if isnan(C(i,ii))
                S(i,ii) = NaN;
            end
            
            C(ii,i) = C(i,ii);
            S(ii,i) = - S(i,ii);
            P(ii,i) = P(i,ii);
            Dir(ii,i) = -Dir(i,ii);
            
        end
        
        X(ii,i) = 1;
        
    end
end            


delete(msg); % Closes the last "Progress window"

n = 0;
XX = zeros(size(C,1));
for i=1 : size(C,1)
    for ii=1 : size(C,2)
        if XX(ii,i)==0 && ii~=i
            n = n+1;
            CvsD(n,1) = C(i,ii);
            CvsD(n,2) = D(i,ii);
            XX(i,ii)=1;
        end
    end
end

writematrix(CvsD, "CvsD_400.txt");
writematrix(C, "C_400.txt");

% 存入params的计算结果
params.C = C;      % 相关性矩阵，后续相关性分析的基础
params.S = S;      % 时移矩阵
params.P = P;      % P值矩阵
params.Dir = Dir;  % 方向性矩阵
params.CvsD = CvsD;  % 相关性vs距离数据

% 存入analysis_summary的统计结果
analysis_summary(30,1) = "Rois with number of peaks not below pk_num_define";
analysis_summary(30,2) = tracesleft;
analysis_summary(31,1) = "all possible correlations below dcutoff";
analysis_summary(31,2) = corrbelowdcutoff2;
analysis_summary(32,1) = "corrleft (all possible correlations below dcutoff between cells with number of peaks not below pk_num_define)";
analysis_summary(32,2) = corrleft;



else 
% ========================================= SECTION 6: 非Mode 4的相关性分析
% %%%%%%%%%%%%%%%%%%%%%%%%%%%这部分代码处理非Mode 4情况下的相关性计算%%%%%%
% 主要区别：计算所有距离小于dcutoff的细胞对的相关性

% SECTION 6: Short-range correlation analysis (Non-mode4)
%
% Purpose:
%   Compute correlations only within spatial cutoff region
%
% 目的：
%   仅计算邻近细胞间的局部同步关系


% 初始化计算需要排除细胞数组
excludedCells = zeros(size(params.allTraces,2),1);
for z = 1:size(params.allTraces,2)
    % Cells with fewer peaks than pk_num_define should not be considered
    if sum(params.peakLocations(:,z)) < params.pk_num_define 
        excludedCells(z) = 1;
    end
end

% 计算基础统计信息
alltraces = size(params.allTraces,2); %ROI细胞总数
tracesleft = alltraces - sum(excludedCells);


% 计算距离矩阵
numCells = size(params.coordinatesData,2)/2;
D = NaN(numCells);
corrleft = 0;

% 计算细胞间距离并统计符合条件的细胞对
for i = 1:numCells
    for ii = 1:numCells
        % 计算细胞坐标和距离
        xi = params.coordinatesData(1,(i*2)-1);
        yi = params.coordinatesData(1,(i*2));
        xii = params.coordinatesData(1,(ii*2)-1);
        yii = params.coordinatesData(1,(ii*2));
        D(i,ii) = sqrt(((xi-xii)^2)+((yi-yii)^2));
        % Calculate distance from the X and Y values of the two cells
        % 计算距离大于dcutoff的细胞对数量（＞一定距离的不被考虑的细胞）
        if excludedCells(i) == 0 && excludedCells(ii)==0 && D(i,ii) <= params.dcutoff && i~=ii 
            corrleft = corrleft + 1;
        end
    end
end

corrbelowdcutoff1 = size(D((0<D)&(D<=params.dcutoff))); 
corrbelowdcutoff2 = corrbelowdcutoff1(1)/2;
corrleft = corrleft/2;

analysis_summary(12,1) = "Rois with number of peaks not below pk_num_define";
analysis_summary(12,2) = tracesleft;
analysis_summary(13,1) = "all possible correlations below dcutoff";
analysis_summary(13,2) = corrbelowdcutoff2;
analysis_summary(14,1) = "corrleft (all possible correlations below dcutoff between cells with number of peaks not below pk_num_define)";
analysis_summary(14,2) = corrleft;


% ================================== SECTION 7: 时间序列分段和相关性计算准备

% SECTION 7: Temporal window segmentation
%
% Purpose:
%   Improve stationarity and robustness by local window correlation
%
% 目的：
%   使用滑动窗口减少非平稳性影响


% 初始化序列参数
length = 10;    % 序列长度（分钟）
steps = 5;      % 时间步长
numofseq = ceil(((size(params.allTraces,1)*params.secondsPerFrame/60)-length) / steps);
if numofseq < 1
    numofseq = 1;
end
fprintf('\n\nNumber of Sequences: %.0f \nLength of Sequence: %.0f\n\n', numofseq, length)

% 分割时间序列
for a = 1:numofseq
    start = floor(((a-1)*(steps*60/params.secondsPerFrame)))+1;
    lim = ceil(start-1+(length*60/params.secondsPerFrame));
    if a == numofseq
        lim = size(params.allTraces,1);
    end
    di = lim-start;
    seq(a,1:di+1,:) = params.allTraces(start:lim,:);
    seqPeakLocations(a,1:di+1,:) = params.peakLocations(start:lim,:);
end


% ========================================= SECTION 8: 初始化相关性计算矩阵

% SECTION 8: Preallocation for performance optimization
%
% Purpose:
%   Avoid dynamic memory allocation in nested loops
%
% 目的：
%   提升运算效率，避免循环中频繁扩展矩阵


% 创建所需的矩阵以避免循环中的大小变化
X = zeros(size(params.allTraces,2));
Cs = NaN(numofseq,size(params.allTraces,2),size(params.allTraces,2));
Ss = NaN(numofseq,size(params.allTraces,2),size(params.allTraces,2));
Ps = NaN(numofseq,size(params.allTraces,2),size(params.allTraces,2));
C = NaN(size(params.allTraces,2));
S = NaN(size(params.allTraces,2));
P = NaN(size(params.allTraces,2));
Dir = NaN(size(params.allTraces,2));
maxPositiveCorrelation = NaN(size(params.allTraces,2));
maxNegativeCorrelation = NaN(size(params.allTraces,2));

% 计算时移范围参数
maxspeed = params.maxspeed * params.secondsPerFrame;
minspeed = params.minspeed * params.secondsPerFrame;
maxshift = ceil(params.dcutoff/minspeed);

% 初始化相关性和时移存储矩阵
A = NaN(size(params.allTraces,2),size(params.allTraces,2),(maxshift * 2)+1);
B = NaN(size(params.allTraces,2),size(params.allTraces,2),(maxshift * 2)+1);
pval = zeros(size(params.allTraces,2),size(params.allTraces,2),(maxshift * 2)+1);

% 初始化进度显示
twentypercent = round(size(params.allTraces,2)/5);
msg = msgbox(sprintf('Progress: 0 percent'));



% ============================================ SECTION 9: 主要相关性计算循环

% SECTION 9: Main correlation computation loop
%
% Steps:
%   - Pairwise ROI traversal
%   - Time-lag scanning
%   - Peak-count validation
%   - Max correlation selection
%
% 流程：
%   - 逐对细胞扫描
%   - 多时移计算
%   - 峰次数筛选
%   - 最大相关性提取



% 遍历所有可能的细胞对
for i = 1:size(params.allTraces,2)
    % 更新进度显示
    if rem(i, twentypercent)==0
        delete(msg);
        msg = msgbox(sprintf('Progress: %.0f percent', i/twentypercent*20));
    end
    
    for ii = 1:size(params.allTraces,2)
        % 检查是否需要计算该细胞对的相关性
        if i ~= ii && X(i,ii)==0 && D(ii,i)<=params.dcutoff && ...
           excludedCells(i) == 0 && excludedCells(ii)==0
            
            % 计算时移范围
            minshift = fix(D(ii,i)/maxspeed);
            maxshift = ceil(D(ii,i)/minspeed);
            
            % 对每个时间序列段计算相关性
            for step = 1:numofseq
                % 特殊处理：随机选择比较序列
                % A randomly selected other step is selected from the second trace
                if params.mode == 3
                 % 确保选择的时间段有足够间隔
                     min_timegap = ceil(length/(steps*2)); % 至少间隔半个分析窗口
                     valid_steps = setdiff(1:numofseq, (step-min_timegap):(step+min_timegap));
        
                      if ~isempty(valid_steps)
                           step_of_other_vector = valid_steps(randi(length(valid_steps)));
                      else
                              step_of_other_vector = mod(step + min_timegap, numofseq) + 1;
                      end
                else
                     step_of_other_vector = step;
                end
                
                % 初始化正负相关计数器
                nn = 0;
                mm = 0;

                % 检查两个细胞的峰值数量是否满足要求
                if sum(seqPeakLocations(step,:,i)) >= params.pk_num_define && ...
                   sum(seqPeakLocations(step,:,ii)) >= params.pk_num_define
                    
                    % 计算不同时移下的相关性
                    for iii = (minshift*2):(maxshift*2)+1
                        if iii >= 2
                            % 确定时移方向和大小
                            if rem(iii, 2) == 0
                                shift = iii/2;
                                shiftpn = 1;  % 正向时移
                            elseif rem(iii, 2) == 1
                                shift = -(iii-1)/2;
                                shiftpn = 0;  % 负向时移
                            end
                            
                            % 计算时移后的相关性
                            movedvector = circshift(seq(step,:,i), shift);
                            movedvector2 = movedvector.';  % 转置以便计算相关性
                            othervector = seq(step_of_other_vector,:,ii);
                            othervector2 = othervector.';
                            
                            % 计算相关系数
                            A(ii,i,iii) = corr(movedvector2, othervector2);
                            B(ii,i,iii) = shift;
                            
                            % 存储正负向时移的相关性
                            if shiftpn == 1
                                nn = nn + 1;
                                correlationPositive(nn,step,ii,i) = A(ii,i,iii);
                            end
                            if shiftpn == 0
                                mm = mm + 1;
                                correlationNegative(nn,step,ii,i) = A(ii,i,iii);
                            end
                        end
                    end
                    
                    % 找到最大相关性并存储
                    [Cs(step,ii,i), maxiii] = max(A(ii,i,:));
                    Cs(step,i,ii) = Cs(step,ii,i);
                    Ss(step,i,ii) = B(ii,i,maxiii);
                end
            end
            % 计算总体最大相关性
            [C(i,ii), loc] = nanmax(Cs(:,ii,i));
            S(i,ii) = Ss(loc,i,ii);
            
            % 如果相关性超过阈值，计算方向性
            if C(i,ii) >= params.threshcorr
                maxPositiveCorrelation(i,ii) = nanmax(nanmax(correlationPositive(:,:,ii,i)));
                maxNegativeCorrelation(i,ii) = nanmax(nanmax(correlationNegative(:,:,ii,i)));
                Dir(i,ii) = (maxPositiveCorrelation(i,ii)-maxNegativeCorrelation(i,ii));
            end
            
            % 处理NaN值
            if isnan(C(i,ii))
                S(i,ii) = NaN;
            end
            
            % 对称赋值
            C(ii,i) = C(i,ii);
            S(ii,i) = -S(i,ii);
            P(ii,i) = P(i,ii);
            Dir(ii,i) = -Dir(i,ii);
        end
        
        % 标记已处理的细胞对
        X(ii,i) = 1;
    end
end

% =========================================== SECTION 10: 后处理和结果保存

% SECTION 10: Post-processing and export
%
% Includes:
%   - Correlation-distance mapping
%   - File export
%   - Summary statistics
%
% 包含：
%   - 相关性-距离关系整理
%   - 文件输出
%   - 统计结果汇总


% 关闭进度窗口
delete(msg);

% 计算相关性与距离的关系
n = 0;
XX = zeros(size(C,1));
for i = 1:size(C,1)
    for ii = 1:size(C,2)
        if XX(ii,i)==0 && ii~=i
            n = n + 1;
            CvsD(n,1) = C(i,ii);
            CvsD(n,2) = D(i,ii);
            XX(i,ii) = 1;
        end
    end
end

% 保存结果到文件
outputValue = sprintf('%s_CvsD%s.txt', params.output, num2str(params.mode));
writematrix(CvsD, outputValue);
outputValue = sprintf('%s_C%s.txt', params.output, num2str(params.mode));
writematrix(C, outputValue);

        end
        end
% 将计算结果存入params结构体
params.C = C;          % 俩细胞之间的相关性矩阵，细胞数NxN矩阵
params.S = S;          % 俩细胞之间信号传播速度矩阵
params.P = P;          % 俩细胞之间相关性P值矩阵
params.Dir = Dir;      % 方向性矩阵（这个可以看成是细胞流动性的方向）
params.D = D;          % 细胞与细胞之间的距离矩阵细胞数NxN
params.CvsD = CvsD;    % 相关性vs距离数据
params.tracesleft = tracesleft;
params.corrleft = corrleft;

% Value Dir:
% Cell A ---> Cell B ---> Cell C
%         \
%           --> Cell D
%
% in the martix we can see like this：
%     A      B      C      D
% A   0     0.95    0      1
% B   -1      0     1      0
% C   0      -1     0      0
% D   -1      0     0      0



% 更新analysis_summary
analysis_summary(15,1) = "Number of significant correlations";
analysis_summary(15,2) = sum(sum(C >= params.threshcorr))/2;
analysis_summary(16,1) = "Mean correlation value";
analysis_summary(16,2) = mean(C(C >= params.threshcorr));
end