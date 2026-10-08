
 
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

function toplot = calcium_fft_analysis_combine()

    % calcium_fft_analysis performs multitaper SVD-based frequency analysis 
    % on multi-frame TIF image data 
    %
    % Output:
    %   toplot — a structure containing analysis results, such as:
    %       toplot.rate     — sampling rate (input)
    %       toplot.mask     — selected ROI binary mask
    %       toplot.pwr      — dominant peak frequency per pixel
    %       toplot.coherence — coherence spectrum
    %       toplot.f_peak   — user-selected frequency peaks
    %       toplot.U        — spatial components per peak
    %       and more...
    
    
    % Modified by Xiao Liu 
    % Code By Spuervisor: Thomas Broggini
    % Department of Neurosurgery
    % Neuroscience Centre
    % University Hospital Frankfurt
    % Goethe University Frankfurt, Germany
    % Frankfurt Cancer Institute, Germany
    % Xiangyang No.1 people's Hospital, China
    % Modified from  UCSD
    % Xiao.Liu@stud.uni-frankfurt.de
    % 2025-04-19
    
    %% Loading Data
    str = 'n';   % note: for cell analysis. this code originally for nuron and vessel analysis
    Fs = 0.5; % the acquisition rate: 2s recorded 1 frame 
    toplot.rate = Fs;  % Store the acquisition rate into the 'toplot' structure
    
    % please check the size of the video
    height=1024;
    width  =1280;
    num_frames=300;
    im_data_ori = uint8(zeros(height, width, num_frames)); 
    [filename, folder] = uigetfile('*.*', 'Please choose the file that contains the ROIs intensities information!');
    if isequal(filename, 0) || isequal(folder, 0)
        disp('Canceled selection of ROI intensity file.');
        return;
    end
    currentfolder = fileparts(fullfile(folder, filename));
    cd(currentfolder);
    all_traces = readmatrix(filename); 
    all_traces = all_traces(2:end, 2:end); 
    [txtname, pathtxt] = uigetfile('*.*', 'Please choose the coordinates of target cells (_centre/_coordinates.txt)');
    if isequal(txtname, 0)
        disp('Canceled selection of ROI coordinates file.');
        return;
    end
    
    fulltxtpath = fullfile(pathtxt, txtname);
    disp(['Selected coordinates file: ', fulltxtpath]);
    cd(pathtxt);
    xyfile = readmatrix(txtname);
    xyfile = xyfile(2, 2:end); 
    
    
    roi_coords = reshape(xyfile, 2, []); 
    x_coords = roi_coords(1, :); 
    y_coords = roi_coords(2, :);
    
    
    num_frames = size(all_traces, 1); 
    num_rois = size(all_traces, 2);   
    
    for frame = 1:num_frames
        for roi = 1:num_rois
            pixel_value = all_traces(frame, roi);
            x = x_coords(roi);
            y = y_coords(roi);
            if x > 0 && y > 0 && x <= size(im_data_ori, 2) && y <= size(im_data_ori, 1)
                im_data_ori(y, x, frame) = uint8(pixel_value); 
            end
        end
    end
    clear x y
    
    im_mask_ind = sub2ind([height, width], xyfile(2:2:end), xyfile(1:2:end));
    
        im_mask=zeros(height, width);
        im_mask(im_mask_ind)=1;
        toplot.mask_ind = im_mask_ind;
        num_skel_ind_unique = numel(toplot.mask_ind);
    skel_label_unique = 1 : num_skel_ind_unique;
    toplot.skel_label = skel_label_unique;
    toplot.mask_size = size(im_mask);
    toplot.mask = im_mask;
    im_mask_ind = find(imresize(im_mask,1));
    im_data_rs = reshape(im_data_ori, numel(im_data_ori(:, :, 1)),[]);
    wave = double(im_data_rs(im_mask_ind,:));
    wave1 = wave';
        % 平滑处理
    for g = 1:num_rois
        wave1(:,g) = imgaussfilt(wave1(:,g), 7*Fs); 
    end
    wave= wave1';
    % Subtract mean from each pixel's time series
    % note that this is not ΔF/Fbase, which has already done before
    wave = bsxfun(@minus, wave, mean(wave, 2));  
    % wave = wave - repmat(mean(wave,2), 1, size(wave, 1); same
    
    
    % Apply bandpass filter (0.0049–0.249 Hz normalized to Nyquist)
    % Design Butterworth bandpass filter
    [b, a] = butter(2, [0.0049 0.249]/(Fs/2), 'bandpass');  
    wave = filtfilt(b, a, wave);  % Apply zero-phase digital filter
    
    clear im_data_rs;
    
    
    %% Dimensionality Reduction and Time-Frequency Resolution Preparation Using SVD and DPSS
    
    [num_frame, num_pixel] = size(wave);  % Get the number of frames (timepoints) and number of pixels (signals)
    [U, S, V] = svd(wave, "econ");  % Perform economical Singular Value Decomposition (SVD)
    
    for i = 1:size(S,2)
        lambda(i) = S(i,i)^2;  % Compute the squared singular values (eigenvalues of covariance matrix)
    end
    
    figure
    plot(log(lambda));  % Plot the logarithm of the eigenvalue spectrum (scree plot)
    xlim([0 100]);  % Limit x-axis to first 100 components
    
    % Prompt user to check the lambda plot and choose how many modes to retain
    answer222 = inputdlg('Check lambda and input the sig-mode you want to plot?');
    Alama1 = str2double(answer222(1));
    close;
    % Keep only significant modes according the Nyquist theorem and analysis before
    % 
    sig_modes = Alama1;  % Number of significant modes to retain
    % Spatial modes (left singular vectors) — each column of Un is a spatial mode
    % Note: Mind the shape of wave, U differs accordingly.
    Un = single(U(:, 1:sig_modes)); % U: left singular vectors (space modes)
    Sn = single(S(1:sig_modes, 1:sig_modes));  % Truncate singular values
    % Temporal modes (right singular vectors) — each column of Vn is a time course
    % Note: Mind the shape of wave, V differs accordingly.
    Vn = single(V(:, 1:sig_modes));   % V: right singular vectors (time modes)
    clear U S V  % Clear original full-size SVD components to save memory
    
    % Parameters for multitaper spectral analysis
    padding_ratio = 2;  % Time base will be zero-padded to 2× next power of 2 of frame count
    toplot.Delta_f = 0.008;  % Desired frequency resolution (half-bandwidth) in Hz
    num_pixel = size(Un,1);  % Number of pixels/signals
    num_frame = size(Vn,1);  % Number of timepoints (frames)
    num_frame_pad = (2 ^ ceil(log2(num_frame))) * padding_ratio;  % Padded number of timepoints for FFT
    toplot.pad = num_frame_pad;  % Store the padded frame length
    
    % Compute the time-bandwidth product p = Δf × T, rounded
    % Δf = frequency resolution, T = acquisition time = num_frame / rate
    p = round(num_frame * toplot.Delta_f / toplot.rate);
    toplot.Delta_f = p * toplot.rate / num_frame;  % Update Delta_f to the exact value based on integer p
    disp(['Bandwidth = ', num2str(toplot.Delta_f), ' Hz'])  % Display effective half-bandwidth
    toplot.num_tapers = 2 * p - 1;  % Number of DPSS tapers to use
    
    % Generate Slepian sequences (DPSS tapers)
    [slep, ~] = dpss(num_frame, p, toplot.num_tapers);  % Generate DPSS basis functions for multitaper analysis
    
    % Multitaper Spectral Estimation Using SVD Components
    
    % This section performs multitaper frequency analysis on the temporal modes (right singular vectors)
    % obtained from SVD, and reconstructs the total power spectrum projected onto spatial modes.
    
    if isunix
        addpath(genpath()) 
    else
        % If running on Windows, add your Chronux toolbox path for multitaper analysis
        addpath(genpath('C:\chronux_2_12'))  
    end
    
    % Format required by Chronux: [half_bandwidth, number_of_tapers]
    ntapers = [(toplot.num_tapers+1)/2, toplot.num_tapers];  
    nfft = toplot.pad;    % Zero-padding size for FFT
    
    % Generate tapers for multitaper spectral estimation
    % Returns DPSS tapers for specified parameters
    tapers = dpsschk(ntapers, size(Vn,1), Fs);  
    
    % Generate frequency grid for spectral estimation
    % Compute frequency vector and valid index range
    [f, findx] = getfgrid(Fs, nfft, [0, Fs/2]);  
    
    % Initialize matrix to store FFT of each tapered signal
    % Dimensions: [frequencies × tapers × temporal_modes]
    taperedFFT = complex(zeros(length(f), ntapers(2), size(Vn,2))); 
    
    % Loop through each temporal mode (right singular vector)
    for i = 1:size(Vn,2)
        J = mtfftc(Vn(:,i), tapers, nfft, Fs);  % Multitaper FFT for one temporal mode
        J = J(findx,:,:);                       % Restrict to valid frequency range
        taperedFFT(:,:,i) = J;                  % Store result
    end
    
    % Reconstruct weighted spatial modes (left singular vectors × singular values)
    scores = Un * Sn;  % [pixels × modes] × [modes × modes] = [pixels × modes]
    clear Un Sn i J
    
    % Initialize total power spectral matrix: [pixels × frequencies], 
    % giving the total power at each frequency for each spatial pixel.
    S_tot = zeros(size(scores,1), size(taperedFFT,1), 'single');
    
    % Loop over tapers to compute power spectrum projection per pixel
    for k = 1:ntapers(2)
        z = scores * squeeze(taperedFFT(:,k,:))';  % Project temporal FFT back onto spatial components
        S_tot = S_tot + conj(z) .* z;              % Accumulate power spectrum (magnitude squared)
    end
    
    % Peak Frequency Detection in Power Spectrum for Each Pixel
    % Note! please check the Start and end of the peak according the freq
    
    % This section visualizes the power spectrum and allows to manually 
    % define the frequency window to search for the peak in the spectrum.
    
    toplot.f = f;  % Store the frequency grid for later use
    % Create a full-screen figure for visualization
    figure('units','normalized','outerposition',[0 0 1 1]);  
    
    % Plot individual log power spectra for several frequency bins
    subplot(2,1,1)
    hold on
    for i = 1 :10: size(S_tot,2)  % Loop through S_tot, plot every 10th pixel's spectrum
        plot(f, log(S_tot(i,:)));  % Log of power spectrum for visualization
    end
    xlim([min(f) max(f)]);  % Set x-axis limits to match the frequency range
    xticks(min(f):0.005:max(f));  % Set x-tick intervals
    
    % Plot the average power spectrum across all pixels
    subplot(2,1,2)
    toplot.mpowr = mean(S_tot,1);  % Compute the mean power spectrum
    plot(log10(toplot.mpowr), 'k')  % Plot the mean spectrum in log scale
    prompt = 'Define Start and end of the peak search: ';  
    % User clicks to define start and end of frequency window
    [toplot.wind, ~] = ginput(2);  
    toplot.wind = fix(toplot.wind);  % Fix to integer values (for indices)
    toplot.pwr = [];  % Initialize power spectrum
    hold off  % Release the plot hold
    
    tic  % Start timing
    % Find the maximum power within the specified frequency window
    % Find peak power in the window
    [~, toplot.pwr] = max(S_tot(:, toplot.wind(1):toplot.wind(2)), [], 2);  
    % Convert index to corresponding frequency
    toplot.pwr = toplot.f(toplot.pwr + toplot.wind(1) - 1);  
    toc  % End timing
    disp(' Done');  % Display completion message
    close
    % Coherence Analysis and Peak Selection
    % choose 4 peaks of the coherence, which will be checked in terms of s.t.
    % This section analyzes the coherence across frequencies and allows to manually select peaks.
    
    % Number of frequency points for 4X oversampling
    plot_num_f_points = fix(4*toplot.rate/2 / toplot.Delta_f);  
    toplot.coherence = zeros(plot_num_f_points,1);  % Initialize coherence array
    % Generate frequency vector from 0 to Nyquist frequency
    toplot.f_vector = linspace(0, toplot.rate/2, plot_num_f_points);  
    
    tic  % Start timing the computation
    % Interpolate the FFT data to the desired frequency points (oversampled)
    interp_FFT = interp1(f, reshape(taperedFFT, size(taperedFFT,1), []), toplot.f_vector);
    % Reshape back to original dimensions
    interp_FFT = reshape(interp_FFT, [], size(taperedFFT,2), size(taperedFFT,3));  
    
    % Compute coherence for each frequency point
    for i = 1:length(toplot.f_vector)
        if str == 'v'  % Check if the user wants to compute with scoresC
            m = scoresC * squeeze(interp_FFT(i,:,ii))';  % Use scoresC for the calculation
        else
            m = scores * squeeze(interp_FFT(i,:,:))';  % Use scores for the calculation
        end
        s = svd(m, 0);  % Compute the singular value decomposition of the matrix
        % Coherence is the squared ratio of the first singular value to 
        % the sum of squared singular values
        toplot.coherence(i) = squeeze(s(1))^2 / sum(s.^2);  
    end
    toc  % End timing the computation
    disp('Done')  % Display a message indicating completion
    
    % Plot the coherence values
    figure('units','normalized','outerposition',[0 0 1 1]);
    % Plot coherence as a function of frequency
    plot(toplot.f_vector, toplot.coherence')  
    xlim([0.005 0.2]);  % Set x-axis limits (frequency range， according to the last ginput frq range)
    xticks(0:0.005:1);  % Set x-tick intervals
    % Allow to manually select 4 peaks on the coherence plot
    % Dont also forgot the target freq range [0.005~0.25hz]
    [toplot.f_peak, ~] = ginput(4);  
    hold on
    
    close
    clear i m s interp_FFT  % Clear variables used in this section
    
    % Singular Value Decomposition (SVD) for Selected Peaks
    % This section performs Singular Value Decomposition (SVD) on the 
    % interpolated FFT data for the selected frequency peaks.
    
    tic  % Start timing the computation
    close all  % Close any open figures
    
    f_global = toplot.f_peak;  % Selected frequency peaks for analysis
    
    % Interpolate the FFT data to the selected frequency points (f_global)
    interp_FFT = interp1(f, reshape(taperedFFT, size(taperedFFT, 1), []), f_global);
    % Reshape back to original dimensions
    interp_FFT = reshape(interp_FFT, [], size(taperedFFT, 2), size(taperedFFT, 3));  
    
    % Initialize the result matrix for SVD components
    toplot.U = zeros(length(toplot.skel_label), toplot.num_tapers, length(toplot.f_peak));  
    
    % Loop over each selected frequency peak
    for i = 1:length(toplot.f_peak)
        if str == 'v'  % Check if we are using scoresC
            m = scoresC * squeeze(interp_FFT(i, :, :))';  % Use scoresC for the computation
        else
            m = scores * squeeze(interp_FFT(i, :, :))';  % Use scores for the computation
        end
        [u, ~, ~] = svd(m, 0);  % Perform Singular Value Decomposition on matrix m
        % Store the left singular vectors (spatial modes) for each peak frequency
        toplot.U(:,:,i) = u(toplot.skel_label, :);  
    end
    % Rearrange dimensions of U for further analysis
    toplot.U = permute(toplot.U, [3, 1, 2]);  
    
    % Clear temporary variables
    clear i m s interp_FFT f_global  
    toc  % End timing the computation
    
    %% save mat
    clear S_tot wave wave1 z
      save([txtname,'.mat']); %'-append');
    
    
    
    %% parameters
    % load('dfof4_coordintes.txt.mat');
    im_mask = toplot.mask;  % Load mask from the 'toplot' structure
    tmp_mode = 1;  % Set the mode （≤ delta_f, ie. ）
    map = zeros(toplot.mask_size);  % Initialize an empty map
    im_size = size(toplot.mask);  % Get the size of the mask
    im_phase = zeros(size(im_mask));  % Initialize an array for phase values
    im_mag = zeros(size(im_mask));  % Initialize an array for magnitude values
    
    % Default parameter ranges for frames, phase, magnitude, power, and frequency
    tmp_max_frm = 10;  % Maximum frame value
    tmp_min_frm = -10;  % Minimum frame value
    tmp_max_phase = pi/2;  % Maximum phase value
    tmp_min_phase = -pi/2;  % Minimum phase value
    tmp_max_mag = 0.003;  % Maximum magnitude value
    tmp_min_mag = 0;  % Minimum magnitude value
    tmp_max_pwr = log10(130);  % Maximum power (log scale)
    tmp_min_pwr = log10(0.1);  % Minimum power (log scale)
    tmp_max_f0 = 0.13;  % Maximum frequency value
    tmp_min_f0 = 0.0;  % Minimum frequency value
    
    % Store the parameter values in the 'toplot' structure for later use
    toplot.prams.max_frm = tmp_max_frm;
    toplot.prams.min_frm = -tmp_min_frm;
    toplot.prams.max_phase = tmp_max_phase;
    toplot.prams.min_phase = tmp_min_phase;
    toplot.prams.max_mag = tmp_max_mag;
    toplot.prams.min_mag = tmp_min_mag;
    toplot.prams.max_pwr = tmp_max_pwr;
    toplot.prams.min_pwr = tmp_min_pwr;
    toplot.prams.binsize = linspace(-pi, pi, 360);  % Define the bin size for phase
    ext = 1;  % Extension parameter for further analysis
    
    % Initialize variables for frequency interpolation
    toplot.findx = findx;
    tapers_FT = fft(tapers, nfft) / Fs;  % Compute the Fourier Transform of the tapers
    tapers_FT = tapers_FT(1, :);  % Select the first row (frequency range)
    t_norm = sum(tapers_FT.^2);  % Normalize the tapers
    clear Un Sn i J
    
    % Define frequency peaks and handle cases with more than 3 peaks
    if size(toplot.f_peak, 1) > 3
        f0 = horzcat(toplot.f_peak(1), linspace(toplot.f_peak(2), toplot.f_peak(2) * fix(1 / toplot.f_peak(2)), fix(1 / toplot.f_peak(2))), linspace(toplot.f_peak(3), toplot.f_peak(3) * fix(1 / toplot.f_peak(3)), fix(1 / toplot.f_peak(3))), linspace(toplot.f_peak(4), toplot.f_peak(4) * fix(1 / toplot.f_peak(4)), fix(1 / toplot.f_peak(4))));
    elseif size(toplot.f_peak, 1) > 2
        f0 = horzcat(toplot.f_peak(1), linspace(toplot.f_peak(2), toplot.f_peak(2) * fix(1 / toplot.f_peak(2)), fix(1 / toplot.f_peak(2))), linspace(toplot.f_peak(3), toplot.f_peak(3) * fix(1 / toplot.f_peak(3)), fix(1 / toplot.f_peak(3))));
    else
        f0 = horzcat(toplot.f_peak(1), toplot.f_peak(2), toplot.f_peak(3));
    end
    
    % Interpolate FFT for each mode and calculate the power spectrum
    for mode = 1:size(scores, 2)
        for k = 1:ntapers(2)
            interp_FFT(:, k, mode) = interp1(f, taperedFFT(:, k, mode), f0);
        end
    end
    
    % Initialize array for storing coherence and power results
    A = zeros(length(f0), size(interp_FFT, 3), 'single');
    A = complex(A);
    
    % Sum the contributions from the tapers
    for k = 1:ntapers(2)
        A = A + tapers_FT(k) .* squeeze(interp_FFT(:, k, :));
    end
    A = A ./ t_norm;  % Normalize the result
    
    % Calculate the mean power for each mode
    mu = scores * A';  % Compute the power (could be sped up using a GPU)
    clear A k
    
    % Set frequency window for further analysis
    toplot.wind(2) = round(length(f) / 2);
    A = zeros(length(f(toplot.wind(1):toplot.wind(2))), size(taperedFFT, 3), 'single');
    A = complex(A);
    for k = 1:ntapers(2)
        A = A + tapers_FT(k) .* squeeze(taperedFFT(toplot.wind(1):toplot.wind(2), k, :));
    end
    A = A ./ t_norm;
    toplot.ampwind = scores * A';  % Store amplitude for the windowed frequencies
    clear A k
    
    % Compute frequency decomposition for each taper
    tic
    Fde = zeros(size(mu), 'single');
    for k = 1:ntapers(2)
        z = scores * squeeze(interp_FFT(:, k, :))' - mu * tapers_FT(k);  % Subtract the mean
        Fde = Fde + conj(z) .* z;  % Calculate the decomposition
        tmp_tic = tic;
        fprintf('Finish calculating %d/%d taper. Elapsed time is %f seconds\n', k, ntapers(2), toc(tmp_tic));  % Display progress
    end
    toc
    
    % Final frequency calculation using the decomposition
    F = (size(taperedFFT, 2) - 1) * (conj(mu) .* mu) * t_norm ./ Fde;
    clear Fde nsvd
    
    % Store results for analysis
    toplot.amps = mu;  % Store the amplitude spectrum
    toplot.fval = F;  % Store the frequency values
    rmpath(genpath('C:\chronux_2_12'));  % Remove the chronux toolbox path
    
    toplot.f0 = f0;  % Store the frequency peaks
    dfact = toplot.Delta_f * 2;  % Factor for scaling the frequency values
    toplot.maxf0 = zeros(1, length(toplot.amps));  % Initialize max frequencies
    
    % Find the maximum frequency for each mode
    for n = 1:length(toplot.amps)
        [val, idx] = max(abs(toplot.amps(n, :)).^2);
        toplot.maxf0(n) = toplot.f0(idx);  % Store the max frequency
    end
    
    % Initialize arrays for the power results
    toplot.maxf = zeros(1, length(toplot.ampwind));
    toplot.extpwr = zeros(1, length(toplot.ampwind));
    
    % Find the maximum power and corresponding frequencies
    for n = 1:length(toplot.ampwind)
        [val, idx] = max(abs(toplot.ampwind(n, :)).^2);
        toplot.maxf(n) = toplot.f(toplot.wind(1) + idx);  % Store the frequency at max power
        toplot.extpwr(n) = val / dfact;  % Store the normalized power
    end
    % Store the resonance（first peak）, puff(2nd peak), and harmonic powers(3rd peak)
    % Dont mind the name, for this project, they are just names.
    toplot.resonancepwr=abs(toplot.amps(:,1)).^2./dfact;
    toplot.puffpwr=abs(toplot.amps(:,2)).^2./dfact;
    toplot.puffharmpwr=abs(toplot.amps(:,3)).^2./dfact;
    % Adjust results based on skeletal label
    toplot.maxf=toplot.maxf(:,toplot.skel_label);
    toplot.maxf0=toplot.maxf0(:,toplot.skel_label);
    toplot.extpwr=toplot.extpwr(:,toplot.skel_label);
    toplot.pwrext=zeros(size(toplot.skel_label));
    toplot.pwrext(1,:)=toplot.resonancepwr(toplot.skel_label);
    toplot.pwrext(2,:)=toplot.puffpwr(toplot.skel_label);
    toplot.pwrext(3,:)=toplot.puffharmpwr(toplot.skel_label);
    % Handle additional peak frequency (5 peaks)
    if size(toplot.f_peak,2)==5
       audiopeak = 1+length(horzcat(toplot.f_peak(1),linspace(toplot.f_peak(2),toplot.f_peak(2)*fix(1/toplot.f_peak(2)),fix(1/toplot.f_peak(2))),linspace(toplot.f_peak(3),toplot.f_peak(3)*fix(1/toplot.f_peak(3)),fix(1/toplot.f_peak(3)))));
       toplot.audpwr=abs(toplot.amps(:,audiopeak)).^2./dfact;
       vispeak = 2+length(linspace(toplot.f_peak(2),toplot.f_peak(2)*fix(1/toplot.f_peak(2)),fix(1/toplot.f_peak(2))));
       toplot.vispwr=abs(toplot.amps(:,vispeak)).^2./dfact;
       toplot.pwrext(3,:)=toplot.vispwr(toplot.skel_label);
       toplot.pwrext(4,:)=toplot.audpwr(toplot.skel_label);
       toplot.pwrext(5,:)=toplot.puffharmpwr(toplot.skel_label);
    end
    if size(toplot.f_peak,2)==4 % (4 peaks)
       vispeak = 2+length(linspace(toplot.f_peak(2),toplot.f_peak(2)*fix(1/toplot.f_peak(2)),fix(1/toplot.f_peak(2))));
       toplot.vispwr=abs(toplot.amps(:,vispeak)).^2./dfact;
       toplot.pwrext(3,:)=toplot.vispwr(toplot.skel_label);
       toplot.pwrext(4,:)=toplot.puffharmpwr(toplot.skel_label);
    end
    % Final result for the frequency peaks
    toplot.f_peak = toplot.f_peak';
    %% create a blank pic
    fig=figure('units','inches','outerposition',[0 0 20 20]);
    ha = tight_subplot(size(toplot.f_peak,2),6,[.02 .0],[.02 .09],[.02 .01]);
    
    %% _________________Max Pwr distribution
    userdicision = true;
    prompt333 = {'For Max Pwr distribution, what is the parameter you want to put in?<eg. log10(130),you need to input 130>：'};
    
    while userdicision
        answer111 = inputdlg(prompt333);
        Valuemaxpwr = str2double(answer111{1});
        tmp_max_pwr = log10(Valuemaxpwr);
        tmp_min_pwr = -log10(Valuemaxpwr);
    
        disp(['Value you input：', num2str(Valuemaxpwr)]);
    
        if size(toplot.pwrext, 1) < size(toplot.f_peak, 2)
            toplot.pwrext(4, :) = toplot.puffharmpwr;
        end
    
        for ii = 1:size(toplot.f_peak, 2)
            axes(ha(ii));
            map = zeros(toplot.mask_size);
            im_size = size(toplot.mask);
    
    
            map(toplot.mask_ind) = log10(toplot.pwrext(ii, :));
            rescaled_pwr = nan(toplot.mask_size);
            tmp_pixel_value = map(toplot.mask_ind);
            tmp_pixel_value = min(tmp_max_pwr, max(tmp_min_pwr, tmp_pixel_value));
            rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;
    
    
            cmap = colormap('jet');
            int_to_cmap = linspace(tmp_min_pwr, tmp_max_pwr, size(cmap, 1));
            non_nan_ind = find(map);
            num_nonnan = numel(non_nan_ind);
            rbg_pwr_list = zeros(num_nonnan, 3);
            rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
            rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
            rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
            rgb_pwr = ones(3, prod(im_size)) * 0.65;
            rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
            rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
            toplot.rgb_pwr2 = permute(rgb_pwr, [2, 3, 1]);
    
    
            imagesc(toplot.rgb_pwr2);
            title(['@', num2str(toplot.f_peak(ii)), ' Hz']);
            axis off;
            box off;
            colormap jet;
            clim([tmp_min_pwr tmp_max_pwr]);
    
    
            hold on;
            [row, col] = find(toplot.mask);
            pixel_size = 10;
            for i = 1:length(row)
                rgb_value = toplot.rgb_pwr2(row(i), col(i), :);
                rectangle('Position', [col(i) - pixel_size/2, row(i) - pixel_size/2, pixel_size, pixel_size], ...
                          'FaceColor', squeeze(rgb_value)', 'EdgeColor', 'none');
            end
            hold off;
    
            daspect([1, 1, 1]);
    
            
    
            output_image = toplot.rgb_pwr2;
    
    
            for i = 1:length(row)
                rgb_value = squeeze(toplot.rgb_pwr2(row(i), col(i), :));
                y_start = max(1, row(i) - floor(pixel_size/2));
                y_end = min(size(output_image, 1), row(i) + floor(pixel_size/2));
                x_start = max(1, col(i) - floor(pixel_size/2));
                x_end = min(size(output_image, 2), col(i) + floor(pixel_size/2));
                for y = y_start:y_end
                    for x = x_start:x_end
                        output_image(y, x, :) = rgb_value;
                    end
                end
            end
    
            filename = sprintf('Max_power_distribution_%d.tif', ii);
            imwrite(output_image, filename, 'TIFF');
        end
    
    
        choice33333 = inputdlg('End and go to the next plotting? y for Yes; any other key for No');
    
        if strcmp(choice33333, 'y')
            userdicision = false;
        else
            userdicision = true;
        end
    end
    
    %% _____________ This is for the Magnitude
    
    
    userdicision = true;
    prompt2 = {'For the Magnitude, what is the parameter you want to put in <eg. 4e-1 = 0.003, you need to put in 0.003 or 4e-1>?'};
    
    while userdicision
    
        answer111 = inputdlg(prompt2);
        Valuemaxfredistri = str2double(answer111{1});
        disp(['Value you input：', num2str(Valuemaxfredistri)]);
        tmp_max_mag = Valuemaxfredistri;
    
    
        for ii = size(toplot.f_peak, 2) + 1 : 2 * size(toplot.f_peak, 2)
            axes(ha(ii));
    
    
            tmp_peak_idx = ii - size(toplot.f_peak, 2);
            tmp_mg_data = squeeze(abs(toplot.U(tmp_peak_idx, :, tmp_mode)));
            map = zeros(toplot.mask_size);
            map(toplot.mask_ind) = tmp_mg_data;
    
            rescaled_pwr = nan(toplot.mask_size);
            tmp_pixel_value = map(toplot.mask_ind);
            tmp_pixel_value = min(tmp_max_mag, max(tmp_min_mag, tmp_pixel_value));
            rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;
    
    
            cmap = colormap('jet');
            int_to_cmap = linspace(tmp_min_mag, tmp_max_mag, size(cmap, 1));
            non_nan_ind = find(map);
            num_nonnan = numel(non_nan_ind);
            rbg_pwr_list = zeros(num_nonnan, 3);
            rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
            rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
            rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
    
    
            rgb_pwr = ones(3, prod(im_size)) * 0.65;
            rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
            rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
            toplot.rgb_mag = permute(rgb_pwr, [2, 3, 1]);
    
     
            imagesc(toplot.rgb_mag);
            axis image; 
            box off;
            colormap jet;
            clim([tmp_min_mag tmp_max_mag]);
            daspect([1, 1, 1]);
            axis off;
    
    
            hold on;
            [row, col] = find(toplot.mask);
            pixel_size = 10;
            for i = 1:length(row)
                rgb_value = toplot.rgb_mag(row(i), col(i), :);
                rectangle('Position', [col(i) - pixel_size / 2, row(i) - pixel_size / 2, pixel_size, pixel_size], ...
                          'FaceColor', squeeze(rgb_value)', 'EdgeColor', 'none');
            end
            hold off;
    
    
            output_image = toplot.rgb_mag;
            for i = 1:length(row)
                rgb_value = squeeze(toplot.rgb_mag(row(i), col(i), :));
                y_start = max(1, row(i) - floor(pixel_size / 2));
                y_end = min(size(output_image, 1), row(i) + floor(pixel_size / 2));
                x_start = max(1, col(i) - floor(pixel_size / 2));
                x_end = min(size(output_image, 2), col(i) + floor(pixel_size / 2));
                for y = y_start:y_end
                    for x = x_start:x_end
                        output_image(y, x, :) = rgb_value;
                    end
                end
            end
            filename = sprintf('Magnitude_%d.tif', ii);
            imwrite(output_image, filename, 'TIFF');
        end
    
    
        choice3 = inputdlg('End and go to the next plotting? (y for Yes; any other key for No)');
        if strcmp(choice3, 'y')
            userdicision = false;
        else
            userdicision = true;
        end
    end
    
    %% _____________ The Phase
    
    for ii = 2*size(toplot.f_peak,2)+1 : 3*size(toplot.f_peak,2)
        axes(ha(ii));
    
    
        tmp_peak_idx = ii - (2 * size(toplot.f_peak,2));
        toplot.phase = squeeze(angle(toplot.U(tmp_peak_idx, :, tmp_mode))); % tmp_taper is last index
        toplot.tphase = angle(sum(toplot.U(tmp_peak_idx, :, tmp_mode)));
        toplot.ophase = squeeze(toplot.phase - toplot.tphase);
        toplot.ophase = mod(toplot.ophase + pi, 2*pi) - pi;
    
     
        map = zeros(toplot.mask_size);
        map(toplot.mask_ind) = toplot.ophase; 
        toplot.polarmap = map(toplot.mask == 1)';
        rescaled_pwr = nan(toplot.mask_size);
        tmp_pixel_value = map(toplot.mask_ind);
        tmp_pixel_value = min(tmp_max_phase, max(tmp_min_phase, tmp_pixel_value));
        rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;
    
        cmap = colormap('jet');
        int_to_cmap = linspace(tmp_min_phase, tmp_max_phase, size(cmap, 1));
        non_nan_ind = find(map);
        num_nonnan = numel(non_nan_ind);
        rbg_pwr_list = zeros(num_nonnan, 3);
        rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
        rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
        rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
    
      
        rgb_pwr = ones(3, prod(im_size)) * 0.65;
        rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
        rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
        toplot.rgb_phase = permute(rgb_pwr, [2, 3, 1]);
    
      
        clear i j pixcolor cmap int_to_cmap;
        imagesc(toplot.rgb_phase); 
        axis image; 
        box off;
        colormap jet; 
        clim([tmp_min_phase tmp_max_phase]);
        daspect([1,1,1]);
        axis off;
    
        hold on;
        [row, col] = find(toplot.mask);
        pixel_size = 10;
        for i = 1:length(row)
            rgb_value = toplot.rgb_phase(row(i), col(i), :);
            rectangle('Position', [col(i) - pixel_size/2, row(i) - pixel_size/2, pixel_size, pixel_size], ...
                      'FaceColor', squeeze(rgb_value)', 'EdgeColor', 'none');
        end
        hold off;
    
       
        if ii == 2*size(toplot.f_peak,2) + 1
            axis on;
            set(gca, 'YTickLabel', []);
            set(gca, 'XTickLabel', []);
            c = colorbar;
            c.LineWidth = 0.01;
            c.Label.String = '\bf Phase [rad]';
            c.Location = 'westoutside';
        end
    
       
        output_image = toplot.rgb_phase;
        for i = 1:length(row)
            rgb_value = squeeze(toplot.rgb_phase(row(i), col(i), :));
            y_start = max(1, row(i) - floor(pixel_size/2));
            y_end = min(size(output_image, 1), row(i) + floor(pixel_size/2));
            x_start = max(1, col(i) - floor(pixel_size/2));
            x_end = min(size(output_image, 2), col(i) + floor(pixel_size/2));
            for y = y_start:y_end
                for x = x_start:x_end
                    output_image(y, x, :) = rgb_value;
                end
            end
        end
    
       
        filename = sprintf('Phase_%d.tif', ii);
        imwrite(output_image, filename, 'TIFF');
    end
    
    %% _____________ plot global Freq distribution
    userdicision = true;
    prompt = {'For the Max Freq distribution, what phase parameter you want to put in<the highest frqency you choose,like 0.2?'};
    while userdicision
       answer = inputdlg(prompt);
       Valuemaxfredistri = str2num(answer{1}); 
    disp(['Value you input：', num2str(Valuemaxfredistri)]);
    prompt = {'For the Max Freq distribution, what phase parameter you want to put in<the lowest frqency you choose,like 0.04?'};
       answer = inputdlg(prompt);
       Valuemaxfredistri1 = str2num(answer{1}); 
    disp(['Value you input：', num2str(Valuemaxfredistri1)]);
    
    tmp_max_f0 = Valuemaxfredistri;
    tmp_min_f0 = Valuemaxfredistri1;
    axes(ha(4*size(toplot.f_peak,2)+1))
    map = zeros(toplot.mask_size);
    im_size = size(toplot.mask);
    map(toplot.mask_ind) = toplot.maxf; %Your Modification
    rescaled_pwr = nan(toplot.mask_size);
    tmp_pixel_value = map(toplot.mask_ind);
    tmp_pixel_value = min(tmp_max_f0, max(tmp_min_f0, tmp_pixel_value) );
    rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;
    cmap = colormap('jet');
    % int_to_cmap = linspace(0,1,size(cmap,1));
    int_to_cmap = linspace(tmp_min_f0, tmp_max_f0, size(cmap,1));
    non_nan_ind = find(map);
    num_nonnan = numel(non_nan_ind);
    rbg_pwr_list = zeros(num_nonnan, 3);
    rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
    % rgb_pwr = zeros(3, prod(im_size));
    rgb_pwr = ones(3, prod(im_size)) * 0.65;
    rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
    rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
    toplot.rgb_maxf0 = permute(rgb_pwr, [2, 3, 1]);
    clear i j pixcolor cmap int_to_cmap
    imagesc(toplot.rgb_maxf0); 
    axis image; 
    box off;
    colormap jet; 
    clim([tmp_min_f0 tmp_max_f0]);
    % caxis([min(toplot.pwr) max(toplot.pwr)]);
    daspect([1,1,1]);
    axis off
    c= colorbar;
        c.LineWidth=0.01;
    %    c.Label.String = ['\bf Phase [rad]'];
        c.Location = 'westoutside';
    c.AxisLocation ='out';
    title('Frequency at Max Power [Hz]')
     
    
    hold on;
            [row, col] = find(toplot.mask); 
            pixel_size = 10;
            for i = 1:length(row)
               
                rgb_value = toplot.rgb_maxf0(row(i), col(i), :);
               
                rectangle('Position', [col(i) - pixel_size/2, row(i) - pixel_size/2, pixel_size, pixel_size], ...
                          'FaceColor', squeeze(rgb_value)', 'EdgeColor', 'none');
            end
            hold off;
              
        output_image = toplot.rgb_phase;
        for i = 1:length(row)
            rgb_value = squeeze(toplot.rgb_phase(row(i), col(i), :));
            y_start = max(1, row(i) - floor(pixel_size/2));
            y_end = min(size(output_image, 1), row(i) + floor(pixel_size/2));
            x_start = max(1, col(i) - floor(pixel_size/2));
            x_end = min(size(output_image, 2), col(i) + floor(pixel_size/2));
            for y = y_start:y_end
                for x = x_start:x_end
                    output_image(y, x, :) = rgb_value;
                end
            end
        end
    
      
        filename = sprintf('Global_freq_distribution_%d.tif', ii);
        imwrite(output_image, filename, 'TIFF');
    
    choice= inputdlg('end and go to the next plotting？ y for Yes; anyother key for No');
    if strcmp(choice, 'y')
        userdicision = false;
      else
         userdicision = true; 
    end
    
    end
    % plot histgraph about frq distribution
        
       num_bins = 20; 
       bin_edges = linspace(tmp_min_f0, tmp_max_f0, num_bins + 1); 
       bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2; 
    
       
       frequency_values = toplot.maxf; 
       frequency_counts = histcounts(frequency_values, bin_edges); 
    
      
       figure;
       bar(bin_centers, frequency_counts, 'BarWidth', 0.9, 'FaceColor', 'r');
       xlabel('Frequency [Hz]');
       ylabel('Frequency Count');
       title('Frequency Distribution');
       grid off;set(gca, 'Box', 'off');
       
       for i = 1:length(frequency_counts)
          
           text(bin_centers(i), frequency_counts(i), num2str(frequency_counts(i)), ...
               'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
       end
         
       xtick_spacing = 0.001;
       xticks(tmp_min_f0:xtick_spacing:tmp_max_f0);
    
       saveas(gcf, ['global_Freq_distri_bar_plot', '.tif']); 
    
    
    im_size = size(toplot.mask); 
    rgb_pwr = ones(3, prod(im_size)) * 0.65;
    rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2)); 
    rgb_pwr = permute(rgb_pwr, [2, 3, 1]);
    
    % magnify the dot to display clearly
    pixel_size = 7; 
    
    
    output_stack = []; 
    for bin_idx = 1:num_bins
        
        current_min = bin_edges(bin_idx);
        current_max = bin_edges(bin_idx + 1);
        
        
        map = zeros(toplot.mask_size);
        map(toplot.mask_ind) = toplot.maxf; 
        mask_in_range = (map >= current_min & map < current_max); 
        
        
        figure('Visible', 'off'); 
        imshow(rgb_pwr); 
        hold on;
        [row, col] = find(mask_in_range); 
        for i = 1:length(row)
           
            rgb_value = toplot.rgb_maxf0(row(i), col(i), :); 
            rectangle('Position', [col(i) - pixel_size/2, row(i) - pixel_size/2, pixel_size, pixel_size], ...
                      'FaceColor', squeeze(rgb_value)', 'EdgeColor', 'none'); 
        end
        
        
        freq_label = sprintf('Freq: %.2f - %.2f Hz', current_min, current_max); 
        text(10, 10, freq_label, 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold', ...
             'BackgroundColor', 'black', 'Margin', 2); 
    
        hold off;
        
       
        frame = getframe(gca);
        img = frame.cdata; 
        output_stack = cat(4, output_stack, img); 
    end
    
    
    output_filename = 'global_freq_distri_stack.tif';
    for frame_idx = 1:size(output_stack, 4)
        imwrite(output_stack(:,:,:,frame_idx), output_filename, 'WriteMode', 'append');
    end
    implay(output_stack)
    
    
    
    %% _____________ Plot global Max  Power distribution
    axes(ha(4*size(toplot.f_peak,2)+2))
    map = zeros(toplot.mask_size);
    im_size = size(toplot.mask);
    map(toplot.mask_ind) = log10(toplot.extpwr); %Your Modification
    rescaled_pwr = nan(toplot.mask_size);
    tmp_pixel_value = map(toplot.mask_ind);
    tmp_pixel_value = min(tmp_max_pwr, max(tmp_min_pwr, tmp_pixel_value) );
    rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;
    cmap = colormap('jet');
    % int_to_cmap = linspace(0,1,size(cmap,1));
    int_to_cmap = linspace(tmp_min_pwr, tmp_max_pwr, size(cmap,1));
    non_nan_ind = find(map);
    num_nonnan = numel(non_nan_ind);
    rbg_pwr_list = zeros(num_nonnan, 3);
    rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
    % rgb_pwr = zeros(3, prod(im_size));
    rgb_pwr = ones(3, prod(im_size)) * 0.65;
    rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
    rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
    toplot.rgb_pwr2 = permute(rgb_pwr, [2, 3, 1]);
    clear i j pixcolor cmap int_to_cmap
    imagesc(toplot.rgb_pwr2); 
    axis image; 
    box off;
    colormap jet; 
    clim([tmp_min_pwr tmp_max_pwr]);
    % caxis([min(toplot.pwr) max(toplot.pwr)]);
    daspect([1,1,1]);
    axis off
    c= colorbar;
        c.LineWidth=0.01;
    %    c.Label.String = ['\bf Phase [rad]'];
        c.Location = 'westoutside';
    c.AxisLocation ='out';
    title('global_max Power distribution')
    
    
    hold on;
            [row, col] = find(toplot.mask); 
            pixel_size = 10; 
            for i = 1:length(row)
             
                rgb_value = toplot.rgb_pwr2(row(i), col(i), :);
                
                rectangle('Position', [col(i) - pixel_size/2, row(i) - pixel_size/2, pixel_size, pixel_size], ...
                          'FaceColor', squeeze(rgb_value)', 'EdgeColor', 'none');
            end
            hold off;

    
    %% _____________ Coherence
    figure; % Create a new figure window
    % Plot coherence over frequency vector
    plot(toplot.f_vector, toplot.coherence', 'Color', 'b');
    xlim([0.005 0.2]); % Limit x-axis from 0 to 1
    hax = gca; % Get current axis handle
    
    % Add vertical lines at peak frequencies with specified colors
    SP = toplot.f_peak(1);
    line([SP SP], get(hax, 'YLim'), 'Color', [1 0 0]); % Red line
    
    SP = toplot.f_peak(2);
    line([SP SP], get(hax, 'YLim'), 'Color', [0 1 0]); % Green line
    
    SP = toplot.f_peak(3);
    line([SP SP], get(hax, 'YLim'), 'Color', [0 1 1]); % Cyan line
    
    SP = toplot.f_peak(4);
    line([SP SP], get(hax, 'YLim'), 'Color', [0 0 1]); % Cyan line
    % Add legend using peak frequencies
    legend('coherence', ...
           [num2str(toplot.f_peak(1)), ' Hz'], ...
           [num2str(toplot.f_peak(2)), ' Hz'], ...
           [num2str(toplot.f_peak(3)), ' Hz'],...
            [num2str(toplot.f_peak(4)), ' Hz']);
    
    % Set plot title and y-axis label
    title('coherence');
    ylabel('coherence');
    
    % Create directory 'figures' if it doesn't exist
    if ~exist('figure', 'dir')
        mkdir('figure');
    end
    saveas(gcf, fullfile('figure', 'coherence.svg'));
    %% _____________ PowerSpect
    figure; % Create a new figure window
    % Add Chronux toolbox to path temporarily to get frequency grid
    addpath(genpath('C:\chronux_2_12'))
    [toplot.f, ~] = getfgrid(toplot.rate, toplot.pad, [0, toplot.rate/2]);
    rmpath(genpath('C:\chronux_2_12')); % Remove path after use
    % Compute log10 power (handling negative/zero values if needed)
    if mean(toplot.mpowr, 2) > 0
        tmp_pwr = log10(toplot.mpowr);
    else
        tmp_pwr = log10(exp(toplot.mpowr)); % Prevent log of negative
    end
    % Plot power spectrum
    if length(toplot.f) >= length(toplot.mpowr)
        plot(toplot.f(toplot.findx), tmp_pwr(1 : size(toplot.f(toplot.findx), 2)), 'k');
    else
        plot(toplot.f, tmp_pwr(1 : size(toplot.f, 2)), 'k');
    end
    
    % Labeling and axis configuration
    ylabel('log10 power [arb]');
    xlabel('frequency [Hz]');
    xlim([0.005 0.2]);
    title('mean power');
    
    % Add vertical lines at frequency peaks with specified colors
    hax = gca; % Get axis handle
    SP = toplot.f_peak(1);
    line([SP SP], get(hax, 'YLim'), 'Color', [1 0 0]); % Red
    
    SP = toplot.f_peak(2);
    line([SP SP], get(hax, 'YLim'), 'Color', [0 1 0]); % Green
    
    SP = toplot.f_peak(3);
    line([SP SP], get(hax, 'YLim'), 'Color', [0 1 1]); % Cyan
    
    SP = toplot.f_peak(4);
    line([SP SP], get(hax, 'YLim'), 'Color', [0 0 1]); % Cyan
    
    % Add legend with frequency labels
    legend('Power', ...
           [num2str(toplot.f_peak(1)),' Hz'], ...
           [num2str(toplot.f_peak(2)),' Hz'], ...
           [num2str(toplot.f_peak(3)),' Hz'], ...
           [num2str(toplot.f_peak(4)),' Hz']);
    
    % Create directory if it doesn't exist
    if ~exist('figure', 'dir')
        mkdir('figure');
    end
    print(gcf, fullfile('figure', 'mean power'), '-dsvg', '-r600');
   
  end