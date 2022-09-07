%   
%   Determine the hemisphere for each electrode using the iEEGPlacementScheme field in the (first run's) JSON file
%
%   hemi = ccep_retrieveElecsHemisphere(jsonFilepathPattern, electrodeTable)
%
%       jsonFilePattern   = the filepath pattern used to search for JSON files. Each run has a JSON file
%                           (e.g. /<inputDir>/ccepAge/sub-<subj>/ses-1/ieeg/sub-<subj>_ses-1_task-SPESclin*_ieeg.json)
%       electrodeTable    = The electrodes table (loaded from a .tsv file) that corresponds to the JSON file. 
%                           For each electrode in the table the hemisphere will be determined.
% 
%
%   Returns: 
%       hemi              = cell array. Each cell represents one electrode and indicates whether
%                           it is either on the left 'L' or right 'R' hemisphere
%   
%
%   Dora Hermes, Max van den Boom, Dorien van Blooijs, 2022
%

function hemi = ccep_retrieveElecsHemisphere(jsonFilepathPattern, electrodeTable)

    % retrieve JSON files based on the given pattern (JSON files of the runs)
    jsonFiles = dir(jsonFilepathPattern);
    
    % load the JSON file that belongs to the first run
    jsonRunData = jsonread(fullfile(jsonFiles(1).folder, jsonFiles(1).name));
    
    % store then number of electrodes
    nElec = size(electrodeTable, 1);
    
    % determine the hemisphere placement based on the iEEGPlacementScheme in the JSON
    if contains(jsonRunData.iEEGPlacementScheme, 'left') && contains(jsonRunData.iEEGPlacementScheme, 'right')
        % electrodes are both left and right
        
        hemi = cell(nElec, 1);
        [hemi{:}] = deal('n/a');

        schemesplit = strsplit(jsonRunData.iEEGPlacementScheme, ';'); 
        rightcell = find(contains(schemesplit, 'right'));
        leftcell = find(contains(schemesplit, 'left'));

        if rightcell < leftcell
            leftcells = extractAfter(jsonRunData.iEEGPlacementScheme, 'left');
            rightcells = extractBetween(jsonRunData.iEEGPlacementScheme, 'right', 'left');
            rightcells = rightcells{:};
        else
            rightcells = extractAfter(jsonRunData.iEEGPlacementScheme, 'right');
            leftcells = extractBetween(jsonRunData.iEEGPlacementScheme, 'left', 'right');
            leftcells = leftcells{:};
        end

        leftelec = strsplit(leftcells, ';');
        leftelec =  leftelec(~cellfun('isempty', leftelec));
        rightelec = strsplit(rightcells, ';');
        rightelec = rightelec(~cellfun('isempty', rightelec));

        % set L in variable hemi for electrodes in the left hemisphere
        for elec=1:size(leftelec, 2)
            C = strsplit(leftelec{elec}, {'[', ']'});
            elecInd = find(contains(electrodeTable.name, C{1}));
            [hemi{elecInd}] = deal('L');
        end

        % set R in variable hemi for electrodes in the right hemisphere
        for elec=1:size(rightelec, 2)
            C = strsplit(rightelec{elec}, {'[', ']'});
            elecInd = find(contains(electrodeTable.name, C{1}));
            [hemi{elecInd}] = deal('R');
        end

    elseif contains(jsonRunData.iEEGPlacementScheme, 'left')
        
        % only in the left hemisphere
        hemi = num2cell(repmat('L', nElec, 1));
        
    elseif contains(jsonRunData.iEEGPlacementScheme, 'right')
        
        % only in the right hemisphere
        hemi = num2cell(repmat('R', nElec, 1));

    end
    
end
