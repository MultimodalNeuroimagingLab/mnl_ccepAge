%
% This script collects the metadata, extracts the average CCEPs and detects N1s for all the subjects
%
% Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, Max van den Boom, 2022
%


%% 
%  Configure

clc
clear
warning('on');
warning('backtrace', 'off')
myDataPath = setLocalDataPath(1);

% check fieldtrip availability and setup
if exist('ft_read_data') ~= 2 || exist('ft_read_header') ~= 2
   error('Could not find FieldTrip functions. FieldTrip is an external dependency, make sure it installed and added as a MATLAB path.');
end

% input whether to store output figures
outputFigures = 0;



%%
%  Collect all electrodes, channels and event metadata for all subjects
%  Extract the average CCEPS and detect the N1s

% inventorize and loop over the subjects
rootFiles = dir(fullfile(myDataPath.input, 'sub-ccepAge*'));
for iFile = 1:size(rootFiles, 1)
    bids_sub = rootFiles(iFile).name;
    
    % inventorize and loop over the sessions
    subjFiles = dir(fullfile(myDataPath.input, rootFiles(iFile).name, 'ses*'));
    for iSess = 1:size(subjFiles, 1)
        bids_ses = subjFiles(iSess).name;
        
        % check if there is a ieeg modality
        ieegSessPath = fullfile(myDataPath.input, rootFiles(iFile).name, subjFiles(iSess).name, 'ieeg');
        if ~exist(ieegSessPath, 'dir')
            error(['Could not find IEEG for participant ', rootFiles(iFile).name])
        end

        % check if the electrodes file is missing
        if ~exist(fullfile(ieegSessPath, [rootFiles(iFile).name, '_' subjFiles(iSess).name, '_electrodes.tsv']), 'file')
           warning(['Could not find electrodes file for subject ', bids_sub, ', session ', bids_ses, '. Skipping subject/session']);
           continue;
        end

        % load electrodes metadata
        electrodes = readtable(fullfile(ieegSessPath, [rootFiles(iFile).name, '_' subjFiles(iSess).name, '_electrodes.tsv']), ...
                                   'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);

        % check to make sure there are at least strips and grids
        if ~any(contains(electrodes.group, 'strip')) && ~any(contains(electrodes.group, 'grid'))
           warning(['No strip nor grid electrodes for subject ', bids_sub, ', session ', bids_ses, '. Skipping subject/session']);
           continue;
        end

        % check if there are no electrode positions available
        if ~any(electrodes.x ~= 0)
           warning(['No electrodes positions available for subject ', bids_sub, ', session ', bids_ses, '. Skipping subject/session']);
           continue;
        end

        % define the task we want to include
        bids_task = 'task-SPESclin';

        % inventorize and loop over the runs
        runFiles = dir(fullfile(myDataPath.input, bids_sub, bids_ses, 'ieeg', [bids_sub '_' bids_ses '_' bids_task '_*'  '_events.tsv']));
        for iRun = 1:size(runFiles, 1)

            % extract the run name
            bids_run = runFiles(iRun).name(strfind(runFiles(iRun).name, 'run-'):strfind(runFiles(iRun).name, 'run-') + 9);
            fprintf('- Run file %s\n', replace(runFiles(iRun).name, '_events.tsv', ''))


            %% 
            %  Load events and channels metadata, and extract essential information
            params = struct();

            % load the events.tsv
            events_name = fullfile(runFiles(iRun).folder,runFiles(iRun).name);
            ccep_events = readtable(events_name, 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);

            % extract the conditions from the events
            events_include = ismember(ccep_events.trial_type, {'electrical_stimulation'}) & ismember(ccep_events.sub_type, {'SPES', 'SPESclin'});
            params.mergeAmp = 1;
            params.mergePlusMin = 1;
            [stim_pair_nr, stim_pair_name, stim_pair_current] = ccep_bidsEvents2conditions(ccep_events, events_include, params);
            
            % detect errors in the events, events closer than 3s
            minOnsetDiff = 3;
            included_events = ccep_events(events_include, :);
            included_diffs = diff(included_events.onset);
            if any(included_diffs < minOnsetDiff)
                for iDiff = 1:length(included_diffs)
                    if included_diffs(iDiff) < minOnsetDiff
                        warning(['Events are less than ', num2str(minOnsetDiff), 's apart. Event onsets: ', num2str(included_events(iDiff, :).onset) ' - ' num2str(included_events(iDiff + 1, :).onset), ' (dist: ', num2str(included_diffs(iDiff)), ')' ]);
                    end 
                end
            end
            
            % check and remove stim-pairs with less than 5 trials
            for iCond = max(stim_pair_nr):-1:1
                trialIndices = find(stim_pair_nr == iCond);
                if length(trialIndices) < 5
                    warning(['Less than 5 trials for stim-pair (', num2str(length(trialIndices)), ' trial(s) - ', stim_pair_name{trialIndices(1)}, '), ignoring stim-pair']);
                    
                    % remove condition
                    stim_pair_nr(trialIndices) = [];
                    stim_pair_name(trialIndices) = [];
                    stim_pair_current(trialIndices) = [];
                    
                    % move other condition
                    lowerIndices = stim_pair_nr > iCond;
                    stim_pair_nr(lowerIndices) = stim_pair_nr(lowerIndices) - 1;
                    
                end
            end
            
            
            if isnan(max(stim_pair_nr))
                warning(['No more trials left in this run, consider removing run!']);
                continue;
            end
            
            
            
            % read the channel 
            channels_tsv_name = replace(fullfile(runFiles(iRun).folder, runFiles(iRun).name), 'events.tsv', 'channels.tsv');
            channels_table = readtable(channels_tsv_name, 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);

            % list of channel names
            channel_names = channels_table.name;

            % find good ECoG channels
            good_channels = find(ismember(upper(channels_table.type), {'ECOG'}) & ismember(channels_table.status, 'good'));



            %%
            %  Load the data (BrainVision BIDS format)

            ieeg_name = replace(fullfile(runFiles(iRun).folder, runFiles(iRun).name), 'events.tsv', 'ieeg.eeg');
            data = ft_read_data(ieeg_name, 'dataformat', 'brainvision_eeg');
            data_hdr = ft_read_header(ieeg_name, 'dataformat', 'brainvision_eeg');
            srate = data_hdr.Fs;



            %%
            %  Epoch and average each condition (stim pair)

            params.epoch_length         = 5;    % total epoch length in sec, default = 5
            params.epoch_prestim_length = 2;    %: prestimulus epoch length in sec, default = 2
            params.baseline_subtract    = 1;    % subtract median baseline from each trial

            [average_ccep, stimpair_names, stimpair_currents, ccep, tt] = ccep_averageConditions(data, srate, ccep_events, channel_names, stim_pair_nr, stim_pair_name, params, 0);



            %% 
            %  Detect N1 in the average CCEPs

            params.amplitude_thresh = 3.4;
            params.n1_peak_range    = 100;
            params.srate            = srate;

            % detect the N1s
            % Note: passing the good channels will result in the non-good channels (i.e. indices in average_ccep that 
            %       are not in the variable good_channels) to be excluded from N1 detection and NaN'ed in the output.
            [n1_peak_sample, n1_peak_amplitude] = ccep_detect_n1peak_ECoG(average_ccep, good_channels, params);
            
            % clear the N1s where either of the stim-pair channels is marked as bad
            good_channel_names = upper(channel_names(good_channels));
            for iStimPair = 1:length(stimpair_names)
                stimPiarElecs = split(stimpair_names{iStimPair}, '-');
                
                % check if either of the stimulation electrodes is not included as good
                if isempty(find(ismember(good_channel_names, upper(stimPiarElecs{1})), 1)) || isempty(find(ismember(good_channel_names, upper(stimPiarElecs{2})), 1))
                    
                    % NaN the N1s for this stim-pair
                    n1_peak_sample(:, iStimPair)    = nan;
                    n1_peak_amplitude(:, iStimPair) = nan;
                    
                end
                
            end
            
            % save files the average CCEPs and N1 detection output (create folder if needed
            saveName = fullfile(myDataPath.output,'derivatives','av_ccep', bids_sub, bids_ses, replace(runFiles(iRun).name, 'events.tsv', 'averageCCEPs.mat'));
            if ~exist(fullfile(myDataPath.output,'derivatives', 'av_ccep', bids_sub,bids_ses), 'dir')
                mkdir(fullfile(myDataPath.output,'derivatives', 'av_ccep', bids_sub,bids_ses));
                disp(['making dir: ', fullfile(myDataPath.output, 'derivatives', 'av_ccep', bids_sub, bids_ses)]);
            end
            save(saveName, 'average_ccep', 'stimpair_names', 'stimpair_currents', 'tt', 'channel_names', 'good_channels', 'n1_peak_sample', 'n1_peak_amplitude')



            %% 
            %  Visually inspect detected N1 in each averaged signal

            % [n1_peak_amplitude_check, n1_peak_sample_check ] = ccep_visualcheck_n1peak_ECoG(average_ccep, ccep,average_ccep_names,channel_names,tt,n1_peak_amplitude,n1_peak_sample);



            %% 
            %  Plot and save averages per channel (optionally)

            if outputFigures == 1
                params.save_fig = 1;
                ccep_plot_av(average_ccep, tt, n1_peak_sample, n1_peak_amplitude, stimpair_names, ...
                             channel_names, good_channels, myDataPath, bids_sub, bids_ses, bids_task, bids_run, params)
            end
            
        end     % end run loop
    end     % end sessions loop
end     % end subjects loop
