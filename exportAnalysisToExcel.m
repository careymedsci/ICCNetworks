
 
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

function exportAnalysisToExcel()
    % Check if analysis_summary exists in workspace
    if ~evalin('base', 'exist(''analysis_summary'', ''var'')')
        errordlg('Variable analysis_summary not found in workspace', 'Variable Error');
        return;
    end
    analysis_summary = evalin('base', 'analysis_summary');

    % Let user choose an Excel file
    [filename, pathname] = uigetfile('*.xlsx', 'Select the Excel file to write to');
    if isequal(filename, 0)
        disp('User cancelled operation');
        return;
    end
    excelFile = fullfile(pathname, filename);

    % Ask whether to create a new sheet
    createSheetChoice = questdlg('Do you want to create a new sheet by copying?', ...
                                 'Create New Sheet', ...
                                 'Yes', 'No', 'Yes');
    if strcmp(createSheetChoice, 'Yes')
        % ========= Choose sheet to copy =========
        sheet_options = {'cbx', 'gap26','B16-F10.N3'};
        [selection, ok] = listdlg('PromptString', 'Select the sheet to copy:', ...
                                'SelectionMode', 'single', ...
                                'ListString', sheet_options);
        if ~ok
            disp('User cancelled operation');
            return;
        end
        copySheetName = sheet_options{selection};

        % ========= Enter new sheet name =========
        newSheetName = inputdlg('Enter the name for the new sheet:', 'New Sheet', 1);
        if isempty(newSheetName)
            disp('User cancelled operation');
            return;
        end
        newSheetName = newSheetName{1};

        % ========= Copy the entire sheet using COM =========
        excel = actxserver('Excel.Application');
        excel.Visible = false;
        wb = excel.Workbooks.Open(excelFile);

        try
            % Find the sheet to copy
            copySheet = wb.Sheets.Item(copySheetName);

            % Insert new sheet at the end
            newSheet = wb.Sheets.Add([], wb.Sheets.Item(wb.Sheets.Count));
            copySheet.Copy([], newSheet); % Full copy (including format)

            % Delete the empty sheet that was added
            newSheet.Delete;

            % Rename the copied sheet
            wb.Sheets.Item(wb.Sheets.Count).Name = newSheetName;

            % Ask whether to clear rows
            clearRowsChoice = questdlg('Do you want to clear specific rows?', ...
                                       'Clear Rows', ...
                                       'Yes', 'No', 'Yes');
            if strcmp(clearRowsChoice, 'Yes')
                rowsToClear = [4,7,10,13,16,19,22,25,26,29,30,33,34,37,38,41,42,45,48,51];
                targetSheet = wb.Sheets.Item(newSheetName);
                for r = rowsToClear
                    targetSheet.Rows.Item(r).ClearContents; % Clear content only, keep format
                end
            end

            % Save and close
            wb.Save;
            wb.Close;
            excel.Quit;
            delete(excel);

            msgbox(sprintf('Sheet copied successfully\nFile: %s\nNew Sheet: %s', filename, newSheetName), 'Completed');

        catch e
            wb.Close(false);
            excel.Quit;
            delete(excel);
            errordlg(['Operation failed: ' e.message], 'Error');
        end
    else
        % User chooses existing sheet
        [~, sheets] = xlsfinfo(excelFile);
        [selection, ok] = listdlg('PromptString', 'Select the target sheet:', ...
                                'SelectionMode', 'single', ...
                                'ListString', sheets);
        if ~ok
            disp('User cancelled operation');
            return;
        end
        newSheetName = sheets{selection};
    end

    % ========= Continue: Ask for column letter =========
    colLetter = inputdlg('Enter the column letter to write to (e.g., a, b, etc.):', 'Select Column', 1);
    if isempty(colLetter)
        disp('User cancelled operation');
        return;
    end
    colLetter = lower(colLetter{1}); % Convert to lowercase
    if ~isletter(colLetter) || length(colLetter) ~= 1
        errordlg('Please enter a single letter (a-z)', 'Input Error');
        return;
    end
    colNum = colLetter - 'a' + 1;

    % Mapping of fields to row numbers
    fieldMap = {
        'peaks / 100 cells in 600s', 4, @(x) x/10;
        'Mean peak width (seconds)', 7, @(x) x;
        'percentage of cells with at least 4 peaks(ca2+ active cells)', 10, @(x) x;
        'percentage of co-active cells', 13, @(x) x;
        'proportion of hub cells in network (above 3 peaks)', 16, @(x) x*100;
        'percentage of periodic cells in network (above 3 peaks)', 19, @(x) x*100;
        'Mean correlation value', 22, @(x) x;
        'periodic_min_frq', 25, @(x) x*1000;
        'periodic_max_frq', 26, @(x) x*1000;
        'general_min_frq', 29, @(x) x*1000;
        'general_max_frq', 30, @(x) x*1000;
        'mean shortest path (original network)', 33, @(x) x;
        'mean shortest path (random network)', 34, @(x) x;
        'mean clustering coefficient (original network)', 37, @(x) x;
        'mean clustering coefficient (random network)', 38, @(x) x;
        'k_mean', 41, @(x) x;
        'k_std', 42, @(x) x;
        'Mean peak width of periodic cells (at least 3 peaks; seconds)', 45, @(x) x;
        'Mean peak width of cells (at least 3 peaks; seconds)', 48, @(x) x;
        'connectivity (all correlations above threshcorr / corrleft)', 51, @(x) x*100;
    };

    % Prepare data to write
    dataToWrite = cell(51, 1);
    missingFields = {};

    % Check type
    if isstruct(analysis_summary)
        for i = 1:size(fieldMap, 1)
            fieldName = fieldMap{i, 1};
            if isfield(analysis_summary, fieldName)
                val = analysis_summary.(fieldName);
                if ~isempty(val)
                    dataToWrite{fieldMap{i, 2}} = fieldMap{i, 3}(val);
                else
                    missingFields{end+1} = [fieldName ' (empty value)'];
                end
            else
                missingFields{end+1} = fieldName;
            end
        end

    elseif isstring(analysis_summary) || ischar(analysis_summary)
        if ischar(analysis_summary)
            analysis_summary = string(analysis_summary);
        end
        for i = 1:size(fieldMap, 1)
            fieldName = fieldMap{i, 1};
            [rowIdx, colIdx] = find(strcmpi(strtrim(analysis_summary), strtrim(fieldName)));
            if ~isempty(rowIdx)
                if colIdx(1) < size(analysis_summary, 2)
                    valStr = analysis_summary(rowIdx(1), colIdx(1) + 1);
                    valNum = str2double(valStr);
                    if isnan(valNum)
                        val = valStr;
                    else
                        val = valNum;
                    end
                    dataToWrite{fieldMap{i, 2}} = fieldMap{i, 3}(val);
                else
                    missingFields{end+1} = [fieldName ' (no value column)'];
                end
            else
                missingFields{end+1} = fieldName;
            end
        end
    else
        errordlg('analysis_summary is neither a struct nor a string array', 'Type Error');
        return;
    end

    if ~isempty(missingFields)
        warndlg(['The following fields are missing or empty:' sprintf('\n') strjoin(missingFields, '\n')], 'Data Warning');
    end

    % ========= Write data into Excel =========
    try
        [~, ~, raw] = xlsread(excelFile, newSheetName);
        if size(raw, 1) < 51
            raw{51, 1} = [];
        end
        for row = 1:51
            if ~isempty(dataToWrite{row})
                raw{row, colNum} = dataToWrite{row};
            end
        end
        xlswrite(excelFile, raw, newSheetName);
        msgbox(sprintf('Data written successfully:\nFile: %s\nSheet: %s\nColumn: %s', filename, newSheetName, colLetter), 'Completed');
    catch e
        errordlg(['Write failed: ' e.message], 'Error');
    end

end
