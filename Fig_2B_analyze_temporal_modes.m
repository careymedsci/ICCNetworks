function Fig_2B_analyze_temporal_modes()

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
% By Thomas Broggini and Xiao Liu
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


   clear; clc; close;
    % Loading Data

    str = 'n';
    Fs = 0.5;
    toplot.rate = Fs;

    [stackname, output_folder] = uigetfile('*.tif');
    if isequal(stackname, 0)
        disp('cancel');
        return;
    end

    info = imfinfo(fullfile(output_folder, stackname));
    num_images = numel(info);

    first_image = imread(fullfile(output_folder, stackname), 1);
    [height, width] = size(first_image);

    im_data_ori = zeros(height, width, num_images, 'uint8');

    for j = 1:num_images
        im_data_ori(:,:,j) = imread(fullfile(output_folder, stackname), j);
    end
    toplot.fname = stackname;
    cd(output_folder);

    clear stackname output_folder j info first_image num_images

    toplot.target = mean(im_data_ori, 3);  
    figure; hold on;
    imshow(toplot.target, []);  
    h = drawellipse();  
    disp('double left click to confirm');
    wait(h);  
    im_mask = createMask(h);
    close;

    figure
    A = toplot.target;
    Nmax = 1000;
    [ Avec, Ind ] = sort(A(:),1,'descend');
    max_values = Avec(1:Nmax);
    [ind_row, ind_col] = ind2sub(size(A),Ind(1:Nmax));
    A(ind_row, ind_col) = mean(A,'all');

    overlay = cat(3, normfunc(A).*255, 0.75.*255.*im_mask);
    overlay(:,:,3) = overlay(:,:,2);
    overlay(:,:,2) = overlay(:,:,1);
    overlay(:,:,1) = 0;
    imshow(uint8(overlay));
    daspect([1,1,1]);

    toplot.mask = im_mask;
    im_mask_ind = find(imresize(im_mask,1));
    toplot.mask_ind = im_mask_ind;
    num_skel_ind_unique = numel(toplot.mask_ind);
    skel_label_unique = 1 : num_skel_ind_unique;
    toplot.skel_label = skel_label_unique;
    toplot.mask_size = size(im_mask);
    im_mask_ind = find(imresize(im_mask,1));

    im_data_rs = reshape(im_data_ori, numel(im_data_ori(:, :, 1)),[]);
    wave = double(im_data_rs(im_mask_ind,:));
    wave = bsxfun(@minus, wave, mean(wave, 2));

    t = (0:size(wave,1)-1)/Fs;
    f = fit(t', mean(wave,2), 'exp1');
    wave = wave - f(t) + mean(f(t));

    [b, a] = butter(2, [0.0049 0.249]/(Fs/2), 'bandpass');
    space_time_data = filtfilt(b, a, wave);

    clear im_data_rs;

    [num_frame, num_pixel] = size(wave);
    [U, S, V] = svd(wave', "econ");
    for i = 1:size(S,2)
        lambda(i) = S(i,i)^2;
    end

    figure
    plot(log(lambda));
    xlim([0 100]);

    answer222 = inputdlg('Check lambda and input the sig-mode you want to plot?');
    Alama1 = str2double(answer222(1));
    close;

    k = Alama1;

    %%
    figure('Color','w','Position',[50 50 1200 800],'Name','Mode Power Spectra');

    rows = ceil(k/5);
    cols = min(k,5);

    gray_color = [0.5 0.5 0.5];
    highlight_color = [0.64, 0.08, 0.18];
    time_color = [0 0.45 0.74];

    t = tiledlayout(rows, cols, 'Padding','compact', 'TileSpacing','compact');

    for i = 1:k
        nexttile;

        temporal_mode = U(:,i);

        [pxx, f] = pwelch(temporal_mode, [], [], [], Fs);

        [peak_power, idx_peak] = max(pxx);
        peak_freq = f(idx_peak);

        if peak_freq < 0.005
            plot_color = gray_color;
        else
            plot_color = highlight_color;
        end

        plot(f, 10*log10(pxx), 'LineWidth', 1.5, 'Color', plot_color);
        hold on;

        time_mode_norm = (temporal_mode - min(temporal_mode)) / ...
                         (max(temporal_mode) - min(temporal_mode));
        time_interp = linspace(0, max(f), length(time_mode_norm));
        time_rescaled = interp1(time_interp, time_mode_norm, f, 'linear', 'extrap');
        time_rescaled = rescale(time_rescaled, min(10*log10(pxx)), max(10*log10(pxx)));

        plot(f, time_rescaled, 'LineWidth', 1.2, 'Color', time_color);

        xlim([0 0.2]);
        set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
        box on;

        text(peak_freq, 10*log10(peak_power), ...
             sprintf('%.3f Hz', peak_freq), ...
             'VerticalAlignment','bottom','HorizontalAlignment','left', ...
             'FontSize',8,'Color',plot_color);

        xlim_vals = xlim;
        ylim_vals = ylim;
        text(xlim_vals(2), ylim_vals(2), ...
             sprintf('Temporal Mode %d', i), ...
             'FontSize', 8, 'HorizontalAlignment','right', ...
             'VerticalAlignment','top', 'Color','k');

        if i == 1
            y_limits = ylim_vals;
        else
            ylim(y_limits);
        end
    end

    sgtitle(['Power Spectra with Temporal Evolution of First ', num2str(k), ' Modes'], ...
        'FontWeight','bold');
end
