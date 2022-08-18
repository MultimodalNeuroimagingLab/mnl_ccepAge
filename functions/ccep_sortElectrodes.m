%
%   Reorders an electrode table or file according to the row order of channel names in a channels.tsv table or file
%
%   This is relevant because the channels.tsv channel names match the electrophysiological data channels, but the rows in
%   electrodes.tsv are often arbitrarily sorted. Therefore the electrodes.tsv rows must be sorted in order to plot electrode
%   positions that correspond correctly to row indices of electrophysiological data.
%   
%   elecsOut = ccep_sortElectrodes(electrodes, channels, saveFile)
%
%       electrodes      = the electrodes to sort. Either as a table or as a filepath to electrodes.tsv file
%       channels        = the channels that the electrodes will be sorted to (by name). Either as a table or as a filepath to electrodes.tsv file
%       saveFile        = Whether to save the sorted electrodes to file (default = true). Will only
%                         store if the electrodes argument is a file.
%
%   Returns:
%       elecsOut        = Table with the rows of electrodes.tsv sorted in order of channel names in channels.tsv
%                         Channel names not found in electrodes.tsv (i.e. ends of leads without coordinates) will be
%                         NaN rows in elecsOut, in order to align overall table structure with channels
%                         NaN rows are saved to file as 'n/a' to match BIDS formatting
%
%
%   Harvey Huang, Multimodal Neuroimaging Lab, Mayo Clinic, 2021
%
function elecsOut = ccep_sortElectrodes(electrodes, channels, saveFile)

    if nargin < 3, saveFile = true; end
    
    if ~istable(channels)
        channels = readtable(channels, 'FileType', 'text', 'Delimiter', '\t'); % keep hyphens for filename
    end
    if istable(electrodes)
        saveFile = false;
        elecs = electrodes;
    else
        elecs = readtable(electrodes, 'FileType', 'text', 'Delimiter', '\t');
    end
    
    % delete hyphens from elecs and channels to ensure they match
    elecNames       = erase(strip(elecs.name), '-'); 
    channelNames    = erase(strip(channels.name), '-');
    
    assert(length(unique(elecNames)) == length(elecNames), 'Repeated names in electrodes file');
    assert(length(unique(channelNames)) == length(channelNames), 'Repeated names in channels file');

    [~, locb] = ismember(channelNames, elecNames);
    
    % varTypes in input cell, removing 
    varTypes = cellfun(@class, table2cell(elecs(1, :)), 'UniformOutput', false);
    varTypes(strcmp('char', varTypes)) = {'string'}; % so that Matlab doesn't scream at me for preallocating with 'char'
    varTypes{1} = 'name'; % ignore first col for indexing purposes
    
    elecsSave = table('Size', [length(channelNames), length(elecs.Properties.VariableNames)], ...
                            'VariableNames', elecs.Properties.VariableNames, ...
                            'VariableTypes', repmat({'string'}, [1, length(elecs.Properties.VariableNames)]));
    elecsSave.name = channels.name; % save with original channel names (whether or not they had hyphens)
    elecsSave(logical(locb), ~strcmp(varTypes, 'double')) = elecs(locb(locb > 0), ~strcmp(varTypes, 'double')); % directly put in non-double rows
    
    % For double variable types, save explicitly with 8 digits of precision
    numValues = table2cell(elecs(locb(locb > 0), strcmp(varTypes, 'double')));
    numValuesStr = cellfun(@(x) num2str(x, 8), numValues, 'UniformOutput', false);
    numValuesStr(strcmp(numValuesStr, 'NaN')) = {'n/a'}; % convert NaNs to 'n/a'
    elecsSave(logical(locb), strcmp(varTypes, 'double')) = numValuesStr;
    
    %nanRow = cell(1, length(varTypes)); % what to put into rows without input electrode info
    %nanRow(strcmpi(varTypes, 'double')) = {nan};
    %nanRow(~strcmpi(varTypes, 'double')) = {'n/a'};
    nanRow = repmat({'n/a'}, [1, length(elecs.Properties.VariableNames)]); % row of all 'n/a's
    elecsSave(~logical(locb), 2:end) = repmat(nanRow(2:end), sum(~logical(locb)), 1);
    
    % convoluted way of replacing <"missing"> values (previously NaN) with 'n/a' because fillmissing does not work with strings
%     C = table2cell(elecsOut);
%     C(ismissing(elecsOut)) = {'n/a'};
%     elecsOut = cell2table(C, 'VariableNames', elecsOut.Properties.VariableNames);
        
    elecsOut = elecsSave; % table to be returned has double columns maintained as double
    doubleCols = strcmp(varTypes, 'double');
    for ii = 1:length(varTypes)
        if doubleCols(ii), elecsOut.(elecs.Properties.VariableNames{ii}) = double(elecsOut.(elecs.Properties.VariableNames{ii})); end
    end

    % 
    if saveFile
        [saveDir, name] = fileparts(electrodes);
        outPath = fullfile(saveDir, sprintf('%s_sorted.tsv', name));
        if exist(outPath, 'file'), warning('Overwriting existing electrodes_sorted.tsv'); end
        writetable(elecsSave, outPath, 'FileType', 'text', 'Delimiter', '\t');
    end
    
end
 