function result = computePairCorrelation(i,ii,params,seq,seqPeakLocations,numofseq,D,maxspeed,minspeed)

    % 初始化
    Cs = NaN(numofseq,1);
    Ss = NaN(numofseq,1);
    correlationPositive = [];
    correlationNegative = [];

    for step = 1:numofseq
        % mode==3 时随机选另一个 step
        if params.mode == 3
            min_timegap = ceil(params.length/(params.steps*2));
            valid_steps = setdiff(1:numofseq, (step-min_timegap):(step+min_timegap));
            if ~isempty(valid_steps)
                step_of_other_vector = valid_steps(randi(length(valid_steps)));
            else
                step_of_other_vector = mod(step + min_timegap, numofseq) + 1;
            end
        else
            step_of_other_vector = step;
        end

        % 峰值数量过滤
        if sum(seqPeakLocations(step,:,i)) >= params.pk_num_define && ...
           sum(seqPeakLocations(step,:,ii)) >= params.pk_num_define

            minshift = fix(D(ii,i)/maxspeed);
            maxshift = ceil(D(ii,i)/minspeed);

            nn = 0; mm = 0;

            for iii = (minshift*2):(maxshift*2)+1
                if iii >= 2
                    if rem(iii,2)==0
                        shift = iii/2; shiftpn = 1; % 正向
                    else
                        shift = -(iii-1)/2; shiftpn = 0; % 负向
                    end

                    movedvector  = circshift(seq(step,:,i), shift);
                    othervector  = seq(step_of_other_vector,:,ii);

                    r = corr(movedvector.', othervector.');

                    if shiftpn==1
                        nn = nn+1;
                        correlationPositive(nn,step) = r;
                    else
                        mm = mm+1;
                        correlationNegative(mm,step) = r;
                    end
                end
            end

            % 找到该 step 下最大相关性
            [Cs(step), idx] = max(correlationPositive(:,step),[],'omitnan');
            if ~isempty(idx)
                Ss(step) = idx; % shift 对应位置
            end
        end
    end

    % 汇总结果
    [C, loc] = max(Cs,[],'omitnan');
    if isempty(loc), C = NaN; S = NaN;
    else, S = Ss(loc); end

    if C >= params.threshcorr
        maxPos = max(correlationPositive(:),[],'omitnan');
        maxNeg = max(correlationNegative(:),[],'omitnan');
        Dir    = maxPos - maxNeg;
    else
        maxPos = NaN; maxNeg = NaN; Dir = NaN;
    end

    % 输出结构
    result = struct('i',i,'ii',ii,'C',C,'S',S,...
                    'maxPos',maxPos,'maxNeg',maxNeg,'Dir',Dir);
end
