%
% Script that aggregates data into one output file, containing:
%      - participants details (name, session, runs, age)
%      - electrode information with destrieux labels
%      - detected N1 responses
%      - Freesurfer labels for stimulation and recording pair
%
% Note: the subjects and runs that will be included in the output file are determined
%       by the dir/file structure that is the result of the 'ccep01_averageCCEPs.m' script.
% 
%
% Dora Hermes, Dorien van Blooijs, Max van den Boom, 2022
%



%% 
%  Set paths
clc
clear
myDataPath = setLocalDataPath(1);


%%
%  Get a list of datasets (output of the 'ccep01_averageCCEPs.m' script)

subjects = ccep_getSubFilenameInfo(myDataPath);


%% 
%  Initialize ccepData in ccepData_init.mat and add subject name, session, runs and age

if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_init.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_init.mat'), 'ccepData')
else
    
    ccepData = [];
    
    % load participants.tsv
    subjectsTsv = readtable(fullfile(myDataPath.input, 'participants.tsv'), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'});
    
    for iSubj = 1:length(subjects)
        disp(['subj ' num2str(iSubj) ' of ' num2str(length(subjects)), ' (', subjects(iSubj).name, ')']);
        
        % retrieve subject name based on dir/file structure
        subjName = subjects(iSubj).name;
        
        % find the subject (index) in the tsv (at first session, if there are multiple)
        [subjTsvIndex] = find(ismember(subjectsTsv.participant_id, subjName), 1);
        assert(~isempty(subjTsvIndex));
        
        % Store the name, session and age
        ccepData(iSubj).id = subjName;
        ccepData(iSubj).ses = subjects(iSubj).ses;
        ccepData(iSubj).age = subjectsTsv.age(subjTsvIndex);
        
        % get number of runs and electrodes info
        ccepData(iSubj).nrRuns = length(subjects(iSubj).run);
        ccepData(iSubj).electrodes = readtable(fullfile(myDataPath.input, subjects(iSubj).name, subjects(iSubj).ses, 'ieeg', [subjects(iSubj).name, '_', subjects(iSubj).ses, '_electrodes.tsv']), ...
                                                 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);
        
        % get hemisphere for each electrode and store in the structure
        hemi = ccep_retrieveElecsHemisphere(fullfile(myDataPath.input, subjects(iSubj).name, subjects(iSubj).ses, 'ieeg', [subjects(iSubj).name, '_', subjects(iSubj).ses, '_task-SPESclin*_ieeg.json']), ...
                                            ccepData(iSubj).electrodes);
        ccepData(iSubj).electrodes.jsonHemi = hemi;
        
        % check/assert if hemi or hemisphere field in electrodes.tsv reflects jsonHemi
        % Note: only check the electrodes that have a value in the hemisphere field
        if any(~strcmp(ccepData(iSubj).electrodes.hemisphere, '') & ~strcmp(ccepData(iSubj).electrodes.hemisphere, hemi))
            warning(['The hemisphere column in the electrodes table and the out of the ccep_retrieveElecsHemisphere differ, the output of the latter will be leading']);
        end
        
        
    end
    
    % optional save the ccepData structure, add more fields later as neccesary
    s = input('Do you want to save the ccepData structure? [y/n]: ', 's');
    if strcmp(s, 'y')
        save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_init.mat'), 'ccepData')
    end
    
end


%%
%  Load all N1 data and add to initialized ccepData_init, then save as ccepData_V1

for iSubj = 1:length(subjects)
    disp(['subj ' num2str(iSubj) ' of ' num2str(length(subjects)), ' (', subjects(iSubj).name, ')']);
    
    for iRun = 1:length(subjects(iSubj).run)
        
        runData = load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', subjects(iSubj).name ,subjects(iSubj).ses, subjects(iSubj).run{iRun}));
        ccepData(iSubj).run(iRun).runName            = subjects(iSubj).run{iRun};
        ccepData(iSubj).run(iRun).n1_peak_sample     = runData.n1_peak_sample;
        ccepData(iSubj).run(iRun).channel_names      = runData.channel_names;
        ccepData(iSubj).run(iRun).stimpair_names     = runData.stimpair_names;
        ccepData(iSubj).run(iRun).good_channels      = runData.good_channels;
        ccepData(iSubj).run(iRun).stimpair_currents  = runData.stimpair_currents;
        % loading all average cceps here makes it very heavy on the memory, do this later
        %ccepData(kk).run(ll).average_ccep            = runData.average_ccep;
        ccepData(iSubj).run(iRun).tt                 = runData.tt;
        clear runData
        
    end
end


%% 
%   Add Freesurfer labels for stimulation and recording pair to output struct

% loop over subjects
for iSubj = 1:length(ccepData)
   
    % loop over runs
    for iRun = 1:length(ccepData(iSubj).run)
        
        % pre-allocation: Destrieux labels and numbers for average CCEP stimulated pairs
        ccepData(iSubj).run(iRun).stimpair_DestrieuxLabel       = cell(size(ccepData(iSubj).run(iRun).stimpair_names, 1), 2);
        ccepData(iSubj).run(iRun).stimpair_DestrieuxNr          = cell(size(ccepData(iSubj).run(iRun).stimpair_names, 1), 2);
        
        % pre-allocation: Destrieux labels and numbers for measured channels
        ccepData(iSubj).run(iRun).channel_DestrieuxLabel         = cell(size(ccepData(iSubj).run(iRun).channel_names));
        ccepData(iSubj).run(iRun).channel_DestrieuxNr            = cell(size(ccepData(iSubj).run(iRun).channel_names));
        
        % loop through CCEP stimulated pairs
        for chPair = 1:length(ccepData(iSubj).run(iRun).stimpair_names)
            
            % get stimulated channels
            stimpchans = strsplit(ccepData(iSubj).run(iRun).stimpair_names{chPair}, '-');
            
            for ch = 1:2
                
                % get first stimulated channel number in_electrodes.tsv
                stim_el_nr = find(strcmpi(ccepData(iSubj).electrodes.name, stimpchans{ch}) == 1);
                
                % check whether the stimulated channel is found in the electrode tsv
                if isempty(stim_el_nr)
                    error(['Could not find stim channel ', stimpchans{ch}, ' in the eletrodes.tsv'])
                end
                
                % transfer the Destrieux codes and labels for the stimulated electrodes
                ccepData(iSubj).run(iRun).stimpair_DestrieuxLabel{chPair,ch} = ccepData(iSubj).electrodes.Destrieux_label_text{stim_el_nr};
                if isnumeric(ccepData(iSubj).electrodes.Destrieux_label)
                    ccepData(iSubj).run(iRun).stimpair_DestrieuxNr{chPair,ch} = int2str(ccepData(iSubj).electrodes.Destrieux_label(stim_el_nr));
                else
                    ccepData(iSubj).run(iRun).stimpair_DestrieuxNr{chPair,ch} = ccepData(iSubj).electrodes.Destrieux_label{stim_el_nr};
                end
                
            end
            clear stim_el_nr stimpchans
        end
        
        % loop through the channels
        for iChan = 1:length(ccepData(iSubj).run(iRun).channel_names)      
            
            % get channel number in_electrodes.tsv
            resp_ElecIndex = find(strcmpi(ccepData(iSubj).electrodes.name, ccepData(iSubj).run(iRun).channel_names{iChan}) == 1);
            
            if ~isempty(resp_ElecIndex)
                
                % transfer the Destrieux code and label for the response electrode
                ccepData(iSubj).run(iRun).channel_DestrieuxLabel{iChan} = ccepData(iSubj).electrodes.Destrieux_label_text{resp_ElecIndex};
                if isnumeric(ccepData(iSubj).electrodes.Destrieux_label)
                    ccepData(iSubj).run(iRun).channel_DestrieuxNr{iChan} = int2str(ccepData(iSubj).electrodes.Destrieux_label(resp_ElecIndex));
                else
                    ccepData(iSubj).run(iRun).channel_DestrieuxNr{iChan} = ccepData(iSubj).electrodes.Destrieux_label{resp_ElecIndex};
                end
                
                clear el1_nr
            else
                
                % warn
                warning(['Could not find response channel ', ccepData(iSubj).run(iRun).channel_names{iChan}, ' in the electrodes.tsv'])
                
                % leave empty
                ccepData(iSubj).run(iRun).channel_DestrieuxLabel{iChan} = NaN;
                ccepData(iSubj).run(iRun).channel_DestrieuxNr{iChan} = NaN;
            end 
            
        end
    end    
end


% optional save the ccepData structure, add more fields later as neccesary
s = input('Do you want to save the ccepData structure? [y/n]: ', 's');
if strcmp(s, 'y')
    save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'ccepData')
end
