
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

function exportMathFitToExcel()
if ~evalin('base', 'exist(''analysis_summary'', ''var'')')
    errordlg('Variable analysis_summary not found in workspace', 'Variable Error');
    return;
end
analysis_summary = evalin('base', 'analysis_summary');

if ~evalin('base', 'exist(''params'', ''var'')')
    errordlg('Variable params not found in workspace', 'Variable Error');
    return;
end
params = evalin('base', 'params');

% select Excel to write
[filename, pathname] = uigetfile('*.xlsx', 'Select the Excel file to write to');
if isequal(filename, 0)
    disp('User cancelled operation');
    return;
end
excelFile = fullfile(pathname, filename);

% The data name and sequence
answer = inputdlg('Input data name:', 'Data Name', [1 50], {'mydata'});
if isempty(answer), return; end
dataname = answer{1};

% ---------- get data from mathlab workspace ----------
getValue = @(name) getValueFromSummary(analysis_summary, name);
% ---------- Scale-free-test ----------
sheet = 'Scale-free-test';
[~, ~, raw] = xlsread(excelFile, sheet);
col = size(raw,1) + 1; % next col 
coltext = sprintf('data%d:%s', col-1, dataname);
toWrite = {coltext, getValue('k_mean'), getValue('k_std')};
xlswrite(excelFile, toWrite, sheet, sprintf('A%d', col));

% ---------- small-word-test ----------
sheet = 'small-word-test';
[~, ~, raw] = xlsread(excelFile, sheet);
col = size(raw,1) + 1;
coltext = sprintf('Data%d:%s', col-2, dataname);
toWrite = {coltext, ...
    getValue('mean shortest path (original network)'), ...
    getValue('mean shortest path (random network)'), ...
    getValue('mean clustering coefficient (original network)'), ...
    getValue('mean clustering coefficient (random network)')};
xlswrite(excelFile, toWrite, sheet, sprintf('A%d', col));

% ---------- direction ----------
sheet = 'direction';
[~, ~, raw] = xlsread(excelFile, sheet);
col = size(raw,1) + 1;
coltext = sprintf('Data%d:%s', col-1, dataname);
toWrite = {coltext, ...
    getValue('direction_allpercell')*1000, ...
    getValue('direction_periopercell')*1000, ...
    getValue('direction_notperiopercell')*1000};
xlswrite(excelFile, toWrite, sheet, sprintf('A%d', col));

% ---------- periodic rate ----------
sheet = 'periodic rate';
[~, ~, raw] = xlsread(excelFile, sheet);
col = size(raw,1) + 1;
coltext = sprintf('Data%d:%s', col-1, dataname);
toWrite = {coltext, ...
    getValue('Percentage of periodic cells among non-hub cells'), ...
    getValue('Percentage of periodic cells in hub cells')*100, ...
    getValue('percentage of periodic cells in network (above 3 peaks)')*100};
xlswrite(excelFile, toWrite, sheet, sprintf('A%d', col));

% ---------- freq_distri_cross_correlat ----------
sheet = 'freq_distri_cross_correlat';
[~, ~, raw] = xlsread(excelFile, sheet);
col = size(raw,1) + 1;
coltext = sprintf('Data%d:%s', col-2, dataname);
toWrite = {coltext, ...
    getValue('periodic_min_frq')*1000, ...
    getValue('periodic_max_frq')*1000, ...
    getValue('general_min_frq')*1000, ...
    getValue('general_max_frq')*1000};
xlswrite(excelFile, toWrite, sheet, sprintf('A%d', col));

% ---------- layer ----------
sheet = 'layer';
[~, ~, raw] = xlsread(excelFile, sheet);
col = size(raw,1) + 1;
loglogname= col-1; % this is for loglog sheet
coltext = sprintf('Data%d:%s', col-1, dataname);
toWrite = {coltext, ...
    getValue('percentage of Ca2+ transient cells'), ...
    getValue('percentage of cells with at least 4 peaks(ca2+ active cells)'), ...
    getValue('percentage of co-active cells'), ...
    getValue('percentage of periodic cells in network (above 3 peaks)')*100, ...
    getValue('proportion of hub cells in network (above 3 peaks)')*100};
xlswrite(excelFile, toWrite, sheet, sprintf('A%d', col));

% ---------- loglog_dfof4 ----------
sheet = 'loglog_dfof4';
[~, ~, raw] = xlsread(excelFile, sheet);
col = size(raw,1) + 2; 
coltext = sprintf('Data%d:%s', loglogname, dataname);

% text for header title
rowData = {coltext, 'freq', 'Probability', 'Probability', 'Probability'};
xlswrite(excelFile, rowData, sheet, sprintf('A%d', col));

% title for the table
headers = {'Coactive cells / cell (k)', 'number of coactive', ...
    'Emperical Data', 'Powerlaw', 'Poisson'};
xlswrite(excelFile, headers, sheet, sprintf('A%d', col+1));

% get data
col1 = params.probabilityDistribution(:,1);
col2 = params.probabilityDistribution(:,2);
col3 = params.probabilityDistribution(:,4);
col4 = params.powerlaw_fit(:,1);
col5 = params.poisson_fit(:,1);

% get the size of the max rows
nRows = max([length(col1), length(col2), length(col3), length(col4), length(col5)]);

col1(end+1:nRows) = 0;
col2(end+1:nRows) = 0;
col3(end+1:nRows) = 0;
col4(end+1:nRows) = 0;
col5(end+1:nRows) = 0;

% preset
dataMat = zeros(nRows + sum(col2==0), 5);  % preset that col with 0
dataMat(1:nRows,1) = col1;
dataMat(1:nRows,2) = col2;
dataMat(1:nRows,3) = col3;
% if col 2 or 3 has 0, do offset to col 4 and 45
currRow = 1;
for i = 1:nRows
    if col2(i) == 0
    currRow = currRow + 1;
    dataMat(currRow,4:5) = [col4(i), col5(i)]; % offset to next row
    dataMat(currRow-1,4:5) = 0;               % set col 4-5 to 0
    else
        dataMat(currRow,4:5) = [col4(i), col5(i)];
    end
    currRow = currRow + 1;
end
while all(dataMat(end,:) == 0)
    dataMat(end,:) = [];
end
xlswrite(excelFile, dataMat, sheet, sprintf('A%d', col+2));

disp('Data successfully exported to Excel.');

% ===== support functions =====
function val = getValueFromSummary(analysis_summary, fieldName)
    %  analysis_summary's fieldName 
    
    if isstruct(analysis_summary)
        % ===  1: struct ===
        if isfield(analysis_summary, fieldName)
            val = analysis_summary.(fieldName);
        else
            error('Field "%s" not found in analysis_summary struct', fieldName);
        end

    elseif iscell(analysis_summary) || isstring(analysis_summary)
        % ===  2: cell  string  ===
        if isstring(analysis_summary)
            analysis_summary = cellstr(analysis_summary); % transform cell 
        end

        % find the correct fieldName
        [rowIdx, colIdx] = find(strcmpi(strtrim(analysis_summary), strtrim(fieldName)));

        if isempty(rowIdx)
            error('Field "%s" not found in analysis_summary cell/string array', fieldName);
        end

        % get the data to fieldName
        if colIdx(1) < size(analysis_summary, 2)
            rawVal = analysis_summary{rowIdx(1), colIdx(1) + 1};
            if ischar(rawVal) || isstring(rawVal)
                val = str2double(rawVal);
                if isnan(val), val = rawVal; end
            else
                val = rawVal;
            end
        else
            error('No value column found for field "%s"', fieldName);
        end

    else
        error('Unsupported analysis_summary type: %s', class(analysis_summary));
    end
end

end