function coherence_power_analysis()
%% Coherence and Power Analysis for Peak Frequency Estimation

% 作者：Thomas Broggini 和 刘晓；法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% Authors：Thomas Broggini
% Modified by Xiao Liu（刘晓）；
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

% This section analyzes the coherence and power associated with specific frequency peaks.
% The code interpolates the FFT of the tapered data to estimate the power 
% and coherence at defined frequency peaks extracted from last part code.

% Initialize necessary parameters
% load([toplot.fname,'.mat']);
im_mask = toplot.mask;  % Load mask from the 'toplot' structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_mode = 1;  % Set the mode （≤ tappers, ie. ）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
tmp_max_f0 = 0.03;  % Maximum frequency value
tmp_min_f0 = 0.005;  % Minimum frequency value

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
ext = 1;  

% Initialize variables for frequency interpolation
toplot.findx = findx;
tapers_FT = fft(tapers, nfft) / Fs;  
% Compute the Fourier Transform of the tapers
tapers_FT = tapers_FT(1, :);  
% Select the first row (frequency range)
t_norm = sum(tapers_FT.^2);  
% Normalize the tapers
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
toplot.ampwind = scores * A';  
% Store amplitude for the windowed frequencies
clear A k

% Compute frequency decomposition for each taper
tic
Fde = zeros(size(mu), 'single');
for k = 1:ntapers(2)
    z = scores * squeeze(interp_FFT(:, k, :))' - mu * tapers_FT(k);  
    % Subtract the mean
    Fde = Fde + conj(z) .* z;  
    % Calculate the decomposition
    tmp_tic = tic;
    fprintf('Finish calculating %d/%d taper. Elapsed time is %f seconds\n', k, ntapers(2), toc(tmp_tic));  % Display progress
end
toc

% Final frequency calculation using the decomposition
F = (size(taperedFFT, 2) - 1) * (conj(mu) .* mu) * t_norm ./ Fde;
clear Fde nsvd

% Store results for analysis
toplot.amps = mu;  
% Store the amplitude spectrum
toplot.fval = F;  
% Store the frequency values
rmpath(genpath('C:\chronux_2_12'));  
% Remove the chronux toolbox path

toplot.f0 = f0;  
% Store the frequency peaks
dfact = toplot.Delta_f * 2;  
% Factor for scaling the frequency values
toplot.maxf0 = zeros(1, length(toplot.amps));  
% Initialize max frequencies

% Find the maximum frequency for each mode
for n = 1:length(toplot.amps)
    [val, idx] = max(abs(toplot.amps(n, :)).^2);
    toplot.maxf0(n) = toplot.f0(idx);  
    % Store the max frequency
end

% Initialize arrays for the power results
toplot.maxf = zeros(1, length(toplot.ampwind));
toplot.extpwr = zeros(1, length(toplot.ampwind));

% Find the maximum power and corresponding frequencies
for n = 1:length(toplot.ampwind)
    [val, idx] = max(abs(toplot.ampwind(n, :)).^2);
    toplot.maxf(n) = toplot.f(toplot.wind(1) + idx);  
    % Store the frequency at max power
    toplot.extpwr(n) = val / dfact;  
    % Store the normalized power
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
% %% create a blank pic
% fig=figure('units','inches','outerposition',[0 0 20 20]);
% ha = tight_subplot(size(toplot.f_peak,2),6,[.02 .0],[.02 .09],[.02 .01]);

%% _________________Max Pwr distribution_freq_related

% Create subfolder for output
if ~exist('figure', 'dir')
    mkdir('figure');
end

userdicision = true;
prompt333 = {'For Max Pwr distribution, what is the parameter you want to put in? <e.g., input 130 for log10(130)>：'};

while userdicision
    answer111 = inputdlg(prompt333);
    Valuemaxpwr = str2double(answer111{1});
    tmp_max_pwr = log10(Valuemaxpwr);
    tmp_min_pwr = -log10(Valuemaxpwr);
    disp(['Value you input: ', num2str(Valuemaxpwr)]);

    % Ensure pwrext has enough rows
    if size(toplot.pwrext,1) < size(toplot.f_peak,2)
        toplot.pwrext(4,:) = toplot.puffharmpwr;
    end

    % Prepare layout for subplots
    num_plots = size(toplot.f_peak,2);
    cols = ceil(sqrt(num_plots));
    rows = ceil(num_plots / cols);
    fig_all = figure;
    t = tiledlayout(rows, cols, 'Padding', 'compact', 'TileSpacing', 'compact');

    for ii = 1:num_plots
        % Prepare power map
        map = zeros(toplot.mask_size);
        im_size = size(toplot.mask);
        map(toplot.mask_ind) = log10(toplot.pwrext(ii,:));
        rescaled_pwr = nan(toplot.mask_size);
        tmp_pixel_value = map(toplot.mask_ind);
        tmp_pixel_value = min(tmp_max_pwr, max(tmp_min_pwr, tmp_pixel_value));
        rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;

        % Map to colormap
        cmap = colormap('jet');
        int_to_cmap = linspace(tmp_min_pwr, tmp_max_pwr, size(cmap,1));
        non_nan_ind = find(map);
        num_nonnan = numel(non_nan_ind);
        rbg_pwr_list = zeros(num_nonnan, 3);
        rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
        rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
        rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
        rgb_pwr = ones(3, prod(im_size)) * 0.5;
        rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
        rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
        toplot.rgb_pwr2 = permute(rgb_pwr, [2, 3, 1]);

        % Plot to multi-panel layout
        nexttile(t);
        imagesc(toplot.rgb_pwr2);
        title(['@', num2str(toplot.f_peak(ii), '%.3f'), ' Hz']);
        axis off; box off;
        colormap jet;
        caxis([tmp_min_pwr tmp_max_pwr]);
        daspect([1,1,1]);

        % Add colorbar only on first tile
        if ii == 1
            c = colorbar;
            c.LineWidth = 0.01;
            c.Label.String = '\bf Max Power distribution';
            c.Location = 'westoutside';
        end

        % -------- Save single image as SVG --------
        fig_single = figure('Visible','off');
        imagesc(toplot.rgb_pwr2);
        title(['@', num2str(toplot.f_peak(ii), '%.3f'), ' Hz']);
        axis off; box off;
        colormap jet;
        caxis([tmp_min_pwr tmp_max_pwr]);
        daspect([1,1,1]);
        c = colorbar;
        c.LineWidth = 0.01;
        c.Label.String = '\bf Max Power distribution';
        c.Location = 'eastoutside';

        filename_svg = sprintf('figure/MaxPwr_%.3fHz.svg', toplot.f_peak(ii));
        print(fig_single, filename_svg, '-dsvg', '-r600');
        close(fig_single);
    end

    % Ask if continue
    choice33333 = inputdlg('End and go to the next plotting? Input y for Yes; any other key for No');
    if strcmp(choice33333, 'y')
        userdicision = false;
    else
        userdicision = true;
    end
end
%% _____________ This is for the Magnitude_freq_related
userdicision = true;
prompt333 = {'For the Magnitude, what is the parameter you want to put in <e.g., 0.003 >:'};

while userdicision
    answer111 = inputdlg(prompt333);
    Valuemaxfredistri = str2double(answer111{1});
    tmp_max_mag = Valuemaxfredistri;
    tmp_min_mag = 0;  
    % Assuming zero as min for magnitude display
    disp(['Value you input: ', num2str(Valuemaxfredistri)]);

    num_freqs = size(toplot.f_peak, 2);
    im_size = size(toplot.mask);
    cols = ceil(sqrt(num_freqs));
    rows = ceil(num_freqs / cols);

    fig_all = figure;
    t = tiledlayout(rows, cols, 'Padding', 'compact', 'TileSpacing', 'compact');

    for ii = 1:num_freqs
        try
            tmp_peak_idx = ii;
            tmp_mg_data = squeeze(abs(toplot.U(tmp_peak_idx, :, tmp_mode)));

            % Build map
            map = zeros(toplot.mask_size);
            map(toplot.mask_ind) = tmp_mg_data;
            rescaled_pwr = nan(toplot.mask_size);
            tmp_pixel_value = map(toplot.mask_ind);
            tmp_pixel_value = min(tmp_max_mag, max(tmp_min_mag, tmp_pixel_value));
            rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;

            % RGB conversion
            cmap = colormap('jet');
            int_to_cmap = linspace(tmp_min_mag, tmp_max_mag, size(cmap,1));
            non_nan_ind = find(map);
            num_nonnan = numel(non_nan_ind);
            rbg_pwr_list = zeros(num_nonnan, 3);
            rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
            rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
            rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
            rgb_pwr = ones(3, prod(im_size)) * 0.5;
            rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
            rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
            toplot.rgb_mag = permute(rgb_pwr, [2, 3, 1]);

            % Multi-panel subplot
            nexttile(t);
            imagesc(toplot.rgb_mag);
            title(['@', num2str(toplot.f_peak(ii), '%.3f'), ' Hz']);
            axis off; box off;
            colormap jet;
            caxis([tmp_min_mag tmp_max_mag]);
            daspect([1,1,1]);

            if ii == 1
                c = colorbar;
                c.LineWidth = 0.01;
                c.Label.String = '\bf Magnitude';
                c.Location = 'westoutside';
            end

            % ---- Save single frequency figure ----
            fig_single = figure('Visible','off');
            imagesc(toplot.rgb_mag);
            title(['@', num2str(toplot.f_peak(ii), '%.3f'), ' Hz']);
            axis off; box off;
            colormap jet;
            caxis([tmp_min_mag tmp_max_mag]);
            daspect([1,1,1]);
            c = colorbar;
            c.LineWidth = 0.01;
            c.Label.String = '\bf Magnitude';
            c.Location = 'eastoutside';

            filename_svg = sprintf('figure/Mag_%.3fHz.svg', toplot.f_peak(ii));
            print(fig_single, filename_svg, '-dsvg', '-r600');
            close(fig_single);

        catch
            warning(['Skipping frequency index ', num2str(ii), ' due to error.']);
            continue;
        end
    end

    % Ask whether to continue
    choice33333 = inputdlg('End and go to the next plotting? Input y for Yes; any other key for No');
    if strcmp(choice33333, 'y')
        userdicision = false;
    else
        userdicision = true;
    end
end
%% _____________ The Phase_freq_related
% Define min/max phase limits 
tmp_max_phase = pi;
tmp_min_phase = -pi;

num_freqs = size(toplot.f_peak, 2);
im_size = size(toplot.mask);
cols = ceil(sqrt(num_freqs));
rows = ceil(num_freqs / cols);

% Create figure with tiled layout for all frequencies
fig_all = figure;
t = tiledlayout(rows, cols, 'Padding', 'compact', 'TileSpacing', 'compact');

for ii = 2*size(toplot.f_peak,2)+1 : 3*size(toplot.f_peak,2)
    try
        tmp_peak_idx = ii - 2 * size(toplot.f_peak,2);

        % Compute phase and normalize to [-pi, pi]
        toplot.phase = squeeze(angle(toplot.U(tmp_peak_idx, :, tmp_mode)));
        toplot.tphase = angle(sum(toplot.U(tmp_peak_idx, :, tmp_mode)));
        toplot.ophase = squeeze(toplot.phase - toplot.tphase);
        toplot.ophase = mod(toplot.ophase + pi, 2*pi) - pi;

        % Create phase map
        map = zeros(toplot.mask_size);
        map(toplot.mask_ind) = toplot.ophase;
        toplot.polarmap = map(toplot.mask == 1)';
        rescaled_pwr = nan(toplot.mask_size);
        tmp_pixel_value = map(toplot.mask_ind);
        tmp_pixel_value = min(tmp_max_phase, max(tmp_min_phase, tmp_pixel_value));
        rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;

        % Map phase to RGB using colormap
        cmap = colormap('jet');
        int_to_cmap = linspace(tmp_min_phase, tmp_max_phase, size(cmap,1));
        non_nan_ind = find(map);
        num_nonnan = numel(non_nan_ind);
        rbg_pwr_list = zeros(num_nonnan, 3);
        rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
        rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
        rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
        rgb_pwr = ones(3, prod(im_size)) * 0.5;
        rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
        rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
        toplot.rgb_phase = permute(rgb_pwr, [2, 3, 1]);

        % Display in tiled layout
        nexttile(t);
        imagesc(toplot.rgb_phase);
        title(['@', num2str(toplot.f_peak(tmp_peak_idx), '%.3f'), ' Hz']);
        axis image; box off;
        colormap jet;
        caxis([tmp_min_phase tmp_max_phase]);
        daspect([1,1,1]);
        axis off;

        if ii == 2*size(toplot.f_peak,2)+1
            axis on;
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            c = colorbar;
            c.LineWidth = 0.01;
            c.Label.String = '\bf Phase [rad]';
            c.Location = 'westoutside';
        end

        % Save each phase image as individual SVG file
        fig_single = figure('Visible','off');
        imagesc(toplot.rgb_phase);
        title(['@', num2str(toplot.f_peak(tmp_peak_idx), '%.3f'), ' Hz']);
        axis image; box off;
        colormap jet;
        caxis([tmp_min_phase tmp_max_phase]);
        daspect([1,1,1]);
        axis off;
        c = colorbar;
        c.LineWidth = 0.01;
        c.Label.String = '\bf Phase [rad]';
        c.Location = 'eastoutside';

        filename_svg = sprintf('figure/Phase_%.3fHz.svg', toplot.f_peak(tmp_peak_idx));
        print(fig_single, filename_svg, '-dsvg', '-r600');
        close(fig_single);

    catch
        warning(['Skipping index ', num2str(ii), ' due to error.']);
        continue;
    end
end

%% _____________ The Phase Plot_freq_related
% Create folder to save figures
if ~exist('figure', 'dir')
    mkdir('figure');
end

num_freqs = size(toplot.f_peak, 2);
cols = ceil(sqrt(num_freqs));
rows = ceil(num_freqs / cols);

% Create one big figure for all subplots
fig_all = figure;
t = tiledlayout(rows, cols, 'Padding', 'compact', 'TileSpacing', 'compact');

for ii = 3*num_freqs+1 : 4*num_freqs
    try
        tmp_peak_idx = ii - 3*num_freqs;

        % Calculate phase and centered phase
        toplot.phase = squeeze(angle(toplot.U(tmp_peak_idx, :, tmp_mode)));
        toplot.tphase = angle(sum(toplot.U(tmp_peak_idx, :, tmp_mode)));
        % toplot.ophase = squeeze(toplot.phase + tmp_angle);
        % toplot.ophase = mod(squeeze(toplot.phase - toplot.tphase) - pi, 2*pi) - pi;
        toplot.ophase = squeeze(toplot.phase - toplot.tphase);
        toplot.ophase = mod(toplot.ophase + pi, 2*pi) - pi;
       
        % Add subplot to main figure
        nexttile(t);
        h = polarhistogram(toplot.ophase, 'BinEdges', toplot.prams.binsize);
        hax = gca;
        hax.ThetaAxisUnits = 'radians';

        % Save individual polar histogram
        fig_single = figure('Visible', 'off');
        polarhistogram(toplot.ophase, 'BinEdges', toplot.prams.binsize);
        hax = gca;
        hax.ThetaAxisUnits = 'radians';

        % Save as SVG with 600 DPI
        filename_svg = sprintf('figure/Phase_plot_%.3fHz.svg', toplot.f_peak(tmp_peak_idx));
        print(fig_single, filename_svg, '-dsvg', '-r600');
        close(fig_single);

    catch
        continue
    end
end

%% _____________ Max Freq distribution
userdicision = true;
prompt = {'For the Max Freq distribution, what phase parameter you want to put in <the highest frequency you choose, e.g., 0.3>:'};
while userdicision
    answer = inputdlg(prompt);
    tmp_max_f0 = str2num(answer{1});
    disp(['Value you input: ', num2str(tmp_max_f0)]);

prompt = {'For the Min Freq distribution, what phase parameter you want to put in <the highest frequency you choose, e.g., 0.005>:'};
    answer = inputdlg(prompt);
    tmp_min_f0 = str2num(answer{1});
    disp(['Value you input: ', num2str(tmp_min_f0)]);   

    tmp_min_f0 = 0;  % Set the lower limit for frequency display
    im_size = size(toplot.mask);
    map = zeros(toplot.mask_size);
    map(toplot.mask_ind) = toplot.maxf;

    rescaled_pwr = nan(toplot.mask_size);
    tmp_pixel_value = map(toplot.mask_ind);
    tmp_pixel_value = min(tmp_max_f0, max(tmp_min_f0, tmp_pixel_value));
    rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;

    cmap = colormap('jet');
    int_to_cmap = linspace(tmp_min_f0, tmp_max_f0, size(cmap,1));
    non_nan_ind = find(map);
    num_nonnan = numel(non_nan_ind);
    rbg_pwr_list = zeros(num_nonnan, 3);
    rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
    rgb_pwr = ones(3, prod(im_size)) * 0.5;
    rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
    rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
    toplot.rgb_maxf0 = permute(rgb_pwr, [2, 3, 1]);

    % -------- Show figure in combined plot --------
    fig_all = figure;
    imagesc(toplot.rgb_maxf0);
    axis image;
    box off;
    colormap jet;
    caxis([tmp_min_f0 tmp_max_f0]);
    daspect([1,1,1]);
    axis off;
    c = colorbar;
    c.LineWidth = 0.01;
    c.Location = 'westoutside';
    c.AxisLocation = 'out';
    title('Frequency at Max Power [Hz]');

    % -------- Save single image --------
    fig_single = figure('Visible','off');
    imagesc(toplot.rgb_maxf0);
    axis image;
    box off;
    colormap jet;
    caxis([tmp_min_f0 tmp_max_f0]);
    daspect([1,1,1]);
    axis off;
    c = colorbar;
    c.LineWidth = 0.01;
    c.Location = 'eastoutside';
    c.AxisLocation = 'out';
    title('Frequency at Max Power [Hz]');
    filename_svg = sprintf('figure/MaxFreqAtPwr_%.3f.svg', tmp_max_f0);
    print(fig_single, filename_svg, '-dsvg', '-r600');
    close(fig_single);

    % Prompt for continuation
    choice33333 = inputdlg('end and go to the next plotting？ y for Yes; anyother key for No');
    if strcmp(choice33333, 'y')
        userdicision = false;
    else
        userdicision = true;
    end
end

%% _____________ Plot Max ext Power distribution
% Create new figure for displaying one map
fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.3, 0.4]);

% Define subplot layout: 1 row, 1 column, position 1
ha = tight_subplot(1, 1, 0.05, 0.05, 0.05); % [nRows, nCols, hgap, vgap, margins]
axes(ha(1)); % Select the only subplot

% Initialize map and get image size
map = zeros(toplot.mask_size);
im_size = size(toplot.mask);

% Apply log10 to external power data within the valid mask
map(toplot.mask_ind) = log10(toplot.extpwr); 

% Rescale values to color axis range
rescaled_pwr = nan(toplot.mask_size);
tmp_pixel_value = map(toplot.mask_ind);
tmp_pixel_value = min(tmp_max_pwr, max(tmp_min_pwr, tmp_pixel_value));
rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;

% Generate colormap
cmap = colormap('jet');
int_to_cmap = linspace(tmp_min_pwr, tmp_max_pwr, size(cmap, 1));

% Only use finite values to avoid interp1 error
non_nan_ind = find(map);
valid_ind = non_nan_ind(isfinite(rescaled_pwr(non_nan_ind)));

% Interpolate RGB from rescaled values
rbg_pwr_list = zeros(numel(valid_ind), 3);
rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(valid_ind));
rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(valid_ind));
rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(valid_ind));

% Assemble RGB image
rgb_pwr = ones(3, prod(im_size)) * 0.5;
rgb_pwr(:, valid_ind) = rbg_pwr_list.';
rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
toplot.rgb_pwr2 = permute(rgb_pwr, [2, 3, 1]);

% Display the RGB image
imagesc(toplot.rgb_pwr2);
axis image;
box off;
colormap jet;
caxis([tmp_min_pwr tmp_max_pwr]);
daspect([1, 1, 1]);
axis off;

% Add colorbar
c = colorbar;
c.LineWidth = 0.01;
c.Location = 'westoutside';
c.AxisLocation = 'out';

% Add title
title('Max ext Power distribution');

% Ensure figure directory exists
if ~exist('figure', 'dir')
    mkdir('figure');
end

% Save figure as high-resolution SVG
print(fig, 'figure/MaxPwr_ext.svg', '-dsvg', '-r600');%%
if size(toplot.f_peak,2)<4
    axes(ha(4*size(toplot.f_peak,2)+3))
    left_color = [0 0 1];
right_color = [0 0 0];
yyaxis left
hax=gca;
hold on
plot(toplot.f_vector,toplot.coherence','Color', 'b');
set(hax,'YColor',[left_color]);
xlim([0 2])
SP=toplot.f_peak(1); %your point goes here 
line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
SP=toplot.f_peak(2); %your point goes here 
line([SP SP],get(hax,'YLim'),'Color',[0 1 0])
SP=toplot.f_peak(3); %your point goes here 
line([SP SP],get(hax,'YLim'),'Color',[0 1 1])
legend('coherence',[num2str(toplot.f_peak(1)),' Hz'],[num2str(toplot.f_peak(2)),' Hz'],[num2str(toplot.f_peak(3)),' Hz']);
title('coherence / mean power')
%xlabel('frequency [Hz]');
ylabel('coherence');
yyaxis right
addpath(genpath('C:\chronux_2_12'))
[toplot.f,toplot.findx] = getfgrid(toplot.rate,toplot.pad,[0,toplot.rate/2]);
rmpath(genpath('C:\chronux_2_12'));
if mean(toplot.mpowr,2)>0
tmp_pwr = log10(toplot.mpowr);
else
tmp_pwr = log10(exp(toplot.mpowr));
end
if length(toplot.f) >= length(toplot.mpowr)==1;
    plot(toplot.f(toplot.findx), tmp_pwr(1 : size(toplot.f(toplot.findx),2)), 'k');
else
    plot(toplot.f, tmp_pwr(1 : size(toplot.f,2)), 'k');
end
ylabel('log10 power [arb]');
xlim([0 1])
ylim([-3.5 -2.5])
set(gca,'YColor',[right_color]);
hax.OuterPosition = [0.63   .01    0.35    0.16];
hax.XTickLabelMode ='auto';
else
end
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
