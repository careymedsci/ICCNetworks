function toplot = calcium_fft_analysis()
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors：Thomas Broggini
% Modified by Xiao Liu（刘晓）； 
% 法兰克福大学附属医院；湖北医药学院附属襄阳市第一人民医院
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany
% Xiangyang No.1 people's Hospital, China

% Xiao.Liu@stud.uni-frankfurt.de
% 2024-04-02

clear; clc; close;
%% Loading Data
str = 'n';  % note: for cell analysis. this code originally for nuron and vessel analysis
Fs = 0.5; 
% the acquisition rate: 2s recorded 1 frame 
toplot.rate = Fs;  
% Store the acquisition rate into the 'toplot' structure

% Prompt user to select a TIFF image file (image stack)
[stackname, output_folder] = uigetfile('*.tif');
if isequal(stackname, 0)
    disp('cancel');  
    return;  
end

% Get information about all frames in the selected TIFF file
info = imfinfo(fullfile(output_folder, stackname));
num_images = numel(info);  

% Read the first frame to get image dimensions
first_image = imread(fullfile(output_folder, stackname), 1);
[height, width] = size(first_image); 

% Preallocate array to store the full image stack
im_data_ori = zeros(height, width, num_images, 'uint8');

% Read each frame from the TIFF file and store in the array
for j = 1:num_images
    im_data_ori(:,:,j) = imread(fullfile(output_folder, stackname), j);
end
toplot.fname = stackname;  
cd(output_folder);  

% Clear temporary variables to clean workspace
clear stackname output_folder j info first_image num_images
%% taking ROI and Data preparation

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
% Use the average image for further processing
Nmax = 1000;  
% Number of maximum intensity pixels to consider
[ Avec, Ind ] = sort(A(:),1,'descend');  
% Sort all pixels in descending order of intensity
max_values = Avec(1:Nmax);  
% Extract top Nmax pixel values
[ind_row, ind_col] = ind2sub(size(A),Ind(1:Nmax));  
% Convert linear indices to row, column
A(ind_row, ind_col) = mean(A,'all');  
% Replace top Nmax pixels with the image mean

% Create an RGB overlay for visualizing the mask on modified image
overlay = cat(3, normfunc(A).*255, 0.75.*255.*im_mask);  
% Red channel = 0, Green and Blue = overlay of mask
overlay(:,:,3) = overlay(:,:,2);  
% Copy green to blue channel
overlay(:,:,2) = overlay(:,:,1);  
% Copy red to green channel (now all R=0, G=B=mask)
overlay(:,:,1) = 0;  
% Red channel remains 0
imshow(uint8(overlay));  
% Display the RGB overlay image
daspect([1,1,1]);  
% Set data aspect ratio to 1:1:1

toplot.mask = im_mask;  
% Store the binary mask
im_mask_ind = find(imresize(im_mask,1));  
% Get linear indices of selected ROI pixels
toplot.mask_ind = im_mask_ind;  
% Store mask indices
num_skel_ind_unique = numel(toplot.mask_ind);  
% Count number of selected pixels
skel_label_unique = 1 : num_skel_ind_unique;  
% Create label index
toplot.skel_label = skel_label_unique;  
% Store label indices
toplot.mask_size = size(im_mask);  
% Store size of the mask
im_mask_ind = find(imresize(im_mask,1));  
% Repeat: get ROI pixel indices again

% Reshape 3D image data into 2D (pixels × time)
im_data_rs = reshape(im_data_ori, numel(im_data_ori(:, :, 1)),[]);
% Extract signal (time series) from ROI pixels
wave = double(im_data_rs(im_mask_ind,:));

% Detrend via exponential fit
t = (0:size(wave,1)-1)/Fs; 
% Time vector for each pixel
f = fit(t', mean(wave,2), 'exp1');  
% Fit exponential to mean signal
wave = wave - f(t) + mean(f(t));  
% Detrend signal by subtracting exponential fit

% Subtract mean from each pixel's time series
% note that this is not ΔF/Fbase, which has already done before
wave = bsxfun(@minus, wave, mean(wave, 2));  
% wave = wave - repmat(mean(wave,2), 1, size(wave, 1); same


% Apply bandpass filter (0.0049–0.249 Hz normalized to Nyquist)
% Design Butterworth bandpass filter
[b, a] = butter(2, [0.0049 0.249]/(Fs/2), 'bandpass');  
% Apply zero-phase digital filter
wave = filtfilt(b, a, wave);  

clear im_data_rs;


%% Dimensionality Reduction and Time-Frequency Resolution Preparation Using
% SVD and DPSS

[num_frame, num_pixel] = size(wave);  
% Get the number of frames (timepoints) and number of pixels (signals)
[U, S, V] = svd(wave, "econ");  
% Perform economical Singular Value Decomposition (SVD)

for i = 1:size(S,2)
    lambda(i) = S(i,i)^2;  
    % Compute the squared singular values (eigenvalues of covariance matrix)
end

figure
plot(log(lambda));  
xlim([0 100]);  

% Prompt user to check the lambda plot and choose how many modes to retain
answer222 = inputdlg('Check lambda and input the sig-mode you want to plot?');
Alama1 = str2double(answer222(1));
close;
% significant modes according the Nyquist theorem and analysis before

sig_modes = Alama1;  
% Number of significant modes to retain
% Spatial modes (left singular vectors)
% each column of Un is a spatial mode
% Note: Mind the shape of wave, U differs accordingly.
Un = single(U(:, 3:sig_modes)); 
% U: left singular vectors (space modes)
Sn = single(S(3:sig_modes, 3:sig_modes));  
% Truncate singular values
% Temporal modes (right singular vectors) 
% each column of Vn is a time course
% Note: Mind the shape of wave, V differs accordingly.
Vn = single(V(:, 3:sig_modes));   
% V: right singular vectors (time modes)

clear U S V 

% Parameters for multitaper spectral analysis
padding_ratio = 2;  
% Time base will be zero-padded to 2× next power of 2 of frame count
% NOTE:   CHOOSE THE BEST DELTA F FOR YOUR DATA!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
toplot.Delta_f = 0.008;  
% Desired frequency resolution (half-bandwidth) in Hz
num_pixel = size(Un,1);  
% Number of pixels/signals
num_frame = size(Vn,1);  
% Number of timepoints (frames)
num_frame_pad = (2 ^ ceil(log2(num_frame))) * padding_ratio;  
% Padded number of timepoints for FFT
toplot.pad = num_frame_pad;  
% Store the padded frame length

% Compute the time-bandwidth product p = Δf × T, rounded
% Δf = frequency resolution, T = acquisition time = num_frame / rate
p = round(num_frame * toplot.Delta_f / toplot.rate);
toplot.Delta_f = p * toplot.rate / num_frame;  
% Update Delta_f to the exact value based on integer p
disp(['Bandwidth = ', num2str(toplot.Delta_f), ' Hz'])  
% Display effective half-bandwidth
toplot.num_tapers = 2 * p - 1;  
% Number of DPSS tapers to use

% Generate Slepian sequences (DPSS tapers)
[slep, ~] = dpss(num_frame, p, toplot.num_tapers);  
% Generate DPSS basis functions for multitaper analysis

%% Multitaper Spectral Estimation Using SVD Components

% This section performs multitaper frequency analysis on the temporal modes
% (right singular vectors)obtained from SVD, and reconstructs the total 
% power spectrum projected onto spatial modes.

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
    J = mtfftc(Vn(:,i), tapers, nfft, Fs);  
    % Multitaper FFT for one temporal mode
    J = J(findx,:,:);                       
    % Restrict to valid frequency range
    taperedFFT(:,:,i) = J;                  
    % Store result
end

% Reconstruct weighted spatial modes (left singular vectors × singular values)
scores = Un * Sn;  % [pixels × modes] × [modes × modes] = [pixels × modes]
clear Un Sn i J

% Initialize total power spectral matrix: [pixels × frequencies], 
% giving the total power at each frequency for each spatial pixel.
S_tot = zeros(size(scores,1), size(taperedFFT,1), 'single');

% Loop over tapers to compute power spectrum projection per pixel
for k = 1:ntapers(2)
    z = scores * squeeze(taperedFFT(:,k,:))';  
    % Project temporal FFT back onto spatial components
    S_tot = S_tot + conj(z) .* z;              
    % Accumulate power spectrum (magnitude squared)
end

%% Peak Frequency Detection in Power Spectrum for Each Pixel
% Note! please check the Start and end of the peak according the freq

% This section visualizes the power spectrum and allows to manually 
% define the frequency window to search for the peak in the spectrum.

toplot.f = f;  % Store the frequency grid for later use
% Create a full-screen figure for visualization
figure('units','normalized','outerposition',[0 0 1 1]);  

% Plot individual log power spectra for several frequency bins
subplot(2,1,1)
hold on
for i = 1 :10: size(S_tot,2)  
    % Loop through S_tot, plot every 10th pixel's spectrum
    plot(f, log(S_tot(i,:)));  
    % Log of power spectrum for visualization
end
xlim([min(f) max(f)]);  
% Set x-axis limits to match the frequency range
xticks(min(f):0.005:max(f));  
% Set x-tick intervals

% Plot the average power spectrum across all pixels
subplot(2,1,2)
toplot.mpowr = mean(S_tot,1);  
% Compute the mean power spectrum
plot(log10(toplot.mpowr), 'k')  
% Plot the mean spectrum in log scale
prompt = 'Define Start and end of the peak search: ';  
% User clicks to define start and end of frequency window
[toplot.wind, ~] = ginput(2);  
toplot.wind = fix(toplot.wind);  
% Fix to integer values (for indices)
toplot.pwr = []; 
% Initialize power spectrum
hold off  

tic  % Start timing
% Find the maximum power within the specified frequency window
% Find peak power in the window
[~, toplot.pwr] = max(S_tot(:, toplot.wind(1):toplot.wind(2)), [], 2);  
% Convert index to corresponding frequency
toplot.pwr = toplot.f(toplot.pwr + toplot.wind(1) - 1);  
toc  % End timing
disp(' Done');  
close
%% Coherence Analysis and Peak Selection
% choose 4 peaks of the coherence, which will be checked in terms of s.t.
% This section analyzes the coherence across frequencies and allows to 
% manually select peaks.

% Number of frequency points for 4X oversampling
plot_num_f_points = fix(4*toplot.rate/2 / toplot.Delta_f);  
toplot.coherence = zeros(plot_num_f_points,1);  
% Initialize coherence array
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
        m = scoresC * squeeze(interp_FFT(i,:,ii))';  
        % Use scoresC for the calculation
    else
        m = scores * squeeze(interp_FFT(i,:,:))';  
        % Use scores for the calculation
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
xlim([0.005 0.2]);  
% Set x-axis limits (frequency range，according to the last ginput frq range)
xticks(0:0.005:1);  
% Set x-tick intervals
% Allow to manually select 4 peaks on the coherence plot
% Dont also forgot the target freq range [0.005~0.25hz]
[toplot.f_peak, ~] = ginput(4);  
hold on

close
clear i m s interp_FFT  % Clear variables used in this section

%% Singular Value Decomposition (SVD) for Selected Peaks
% This section performs Singular Value Decomposition (SVD) on the 
% interpolated FFT data for the selected frequency peaks.

tic 
close all  

f_global = toplot.f_peak;  
% Selected frequency peaks for analysis

% Interpolate the FFT data to the selected frequency points (f_global)
interp_FFT = interp1(f, reshape(taperedFFT, size(taperedFFT, 1), []), f_global);
% Reshape back to original dimensions
interp_FFT = reshape(interp_FFT, [], size(taperedFFT, 2), size(taperedFFT, 3));  

% Initialize the result matrix for SVD components
toplot.U = zeros(length(toplot.skel_label), toplot.num_tapers, length(toplot.f_peak));  

% Loop over each selected frequency peak
for i = 1:length(toplot.f_peak)
    if str == 'v'  % Check if we are using scoresC
        m = scoresC * squeeze(interp_FFT(i, :, :))';  
        % Use scoresC for the computation
    else
        m = scores * squeeze(interp_FFT(i, :, :))';  
        % Use scores for the computation
    end
    [u, ~, ~] = svd(m, 0);  
    % Perform Singular Value Decomposition on matrix m
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
  save([toplot.fname,'.mat']); %'-append');
end