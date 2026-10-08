
 
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


function plot_csv_waveforms()
    % Manually select a CSV file
    [fileName, filePath] = uigetfile('*.csv', 'Select a CSV file with waveform data');
    if isequal(fileName, 0)
        disp('User canceled file selection.');
        return;
    end
    fullFileName = fullfile(filePath, fileName);

    % Read CSV data, starting from the second row and second column
    rawData = readmatrix(fullFileName);  % Assumes headers in row 1
    if size(rawData, 2) < 2
        error('Insufficient data columns. Waveform data should start from column 2.');
    end

    waveformData = rawData(:, 2:end);  % Skip the first column
    [numPoints, numTraces] = size(waveformData);

    % Create a time vector (frame-based index)
    timeVector = (1:numPoints)';

    % Create figure and plot waveforms
    figure('Color', 'w', 'Name', 'CSV Waveform Plot');
    hold on;
    colors = lines(numTraces);  % Generate distinguishable colors

    for i = 1:numTraces
        plot(timeVector, waveformData(:, i), 'LineWidth', 1.2, 'Color', colors(i, :));
    end

    xlabel('Time (frames)');
    ylabel('Signal amplitude');
    title('Waveforms from CSV File');
    legend(arrayfun(@(x) sprintf('Trace %d', x), 1:numTraces, 'UniformOutput', false));
    grid on;
    hold off;
end
