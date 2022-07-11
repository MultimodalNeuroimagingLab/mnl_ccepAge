%
% This script can be used as workflow script to create average CCEPs
% for the CCEP data in the RESPect database.
%
% Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, 2019
%

%% Set paths
clc
clear
myDataPath = setLocalDataPath(1);

%% Metadata: run all participants and all runs in a dataset
% if ECoG and if electrode positions are determined.

files = dir(fullfile(myDataPath.input, 'sub-ccepAge*'));

for n = 1:size(files,1)

    filessub = dir(fullfile(myDataPath.input, files(n).name, 'ses*'));
    bids_sub = files(n).name;

    for m = 1:size(filessub, 1)
        
        if exist(fullfile(myDataPath.input, files(n).name, filessub(m).name, 'ieeg'), 'dir')
            filesses = dir(fullfile(myDataPath.input, files(n).name, filessub(m).name, 'ieeg'));
            bids_ses = filessub(m).name;

            if exist(fullfile(filesses(1).folder, [files(n).name, '_' filessub(m).name, '_electrodes.tsv']), 'file')
                electrodes_tsv = read_tsv(fullfile(filesses(1).folder, [files(n).name, '_' filessub(m).name, '_electrodes.tsv']));

                if any(contains(electrodes_tsv.group, 'strip')) || any(contains(electrodes_tsv.group, 'grid'))

                    if any(~isnan(str2double(electrodes_tsv.x)))
                        electrodes_tsv.x = str2double(electrodes_tsv.x);
                    end

                    if any(electrodes_tsv.x ~= 0) % only run the rest if electrode positions are available

                        bids_task = 'task-SPESclin';

                        filesrun = dir(fullfile(myDataPath.input, bids_sub, bids_ses, 'ieeg', ...
                            [bids_sub '_' bids_ses '_' bids_task '_*'  '_events.tsv']));

                        for ii = 1:size(filesrun, 1)

                            bids_runs = filesrun(ii).name(strfind(filesrun(ii).name, 'run-'):strfind(filesrun(ii).name, 'run-') + 9);
                            fprintf('Run file %s!\n', replace(filesrun(ii).name, '_events.tsv', ''))

                            %% load events

                            % load the events.tsv
                            events_name = fullfile(filesrun(ii).folder,filesrun(ii).name);
                            ccep_events = readtable(events_name,'FileType', 'text', 'Delimiter', '\t');

                            % generate vector for averaging across trials
                            events_include = ismember(ccep_events.sub_type, {'SPES', 'SPESclin'});
                            params.mergeAmp = 1;
                            params.mergePlusMin = 1;

                            [stim_pair_nr, stim_pair_name] = ccep_bidsEvents2conditions(ccep_events, events_include, params);

                            %% load data and channels

                            % load the data, as BrainVision BIDS format
                            ieeg_name = replace(fullfile(filesrun(ii).folder, filesrun(ii).name),'events.tsv', 'ieeg.eeg');
                            data = ft_read_data(ieeg_name, 'dataformat', 'brainvision_eeg');
                            data_hdr = ft_read_header(ieeg_name, 'dataformat', 'brainvision_eeg');

                            channels_tsv_name = replace(fullfile(filesrun(ii).folder, filesrun(ii).name), 'events.tsv', 'channels.tsv');
                            channels_table = readtable(channels_tsv_name,'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'});

                            %% get necessary parameters from data

                            % sampling frequency
                            srate = data_hdr.Fs;

                            % list of channel names
                            channel_names = channels_table.name;

                            % find good ECoG channels
                            good_channels = find(ismember(channels_table.type, {'ECOG'}) & ismember(channels_table.status, 'good'));

                            %% load data for each condition and save averages

                            params.epoch_length = 5;          % total epoch length in sec, default = 5
                            params.epoch_prestim_length = 2;  %: prestimulus epoch length in sec, default = 2
                            params.baseline_subtract = 1;     % subtract median baseline from each trial

                            [average_ccep, average_ccep_names, ccep, tt] = ccep_averageConditions(data, srate, ccep_events, channel_names, stim_pair_nr, stim_pair_name, params);

                            %% detect N1 in each averaged signal

                            params.amplitude_thresh = 3.4;
                            params.n1_peak_range = 100;
                            params.srate = srate;

                            % including no detection possible in bad channels
                            [n1_peak_sample, n1_peak_amplitude] = ccep_detect_n1peak_ECoG(average_ccep, good_channels,params);

                            % save files
                            saveName = fullfile(myDataPath.output,'derivatives','av_ccep', bids_sub, bids_ses, ...
                                                replace(filesrun(ii).name, 'events.tsv', 'averageCCEPs.mat'));

                            if ~exist(fullfile(myDataPath.output,'derivatives', 'av_ccep',bids_sub,bids_ses), 'dir')
                                mkdir(fullfile(myDataPath.output,'derivatives', 'av_ccep',bids_sub,bids_ses))
                                sprintf(['making dir:\n', fullfile(myDataPath.output, 'derivatives', 'av_ccep', bids_sub, bids_ses)])
                            end

                            save(saveName, 'average_ccep', 'average_ccep_names', 'tt', 'channel_names', 'good_channels', ...
                                 'n1_peak_sample', 'n1_peak_amplitude')

                            %% check detected N1 in each averaged signal

                            % [n1_peak_amplitude_check, n1_peak_sample_check ] = ccep_visualcheck_n1peak_ECoG(average_ccep, ccep,average_ccep_names,channel_names,tt,n1_peak_amplitude,n1_peak_sample);

                            %% plot and save averages per channel
                            params.save_fig = 0;%str2double(input('Do you want to save the figures? [yes = 1, no = 0]: ','s'));

                            % plotting with N1 peak detection:
                            if params.save_fig == 1 % only save figure if also plotting
                                ccep_plot_av(average_ccep, tt, n1_peak_sample, n1_peak_amplitude, average_ccep_names, ...
                                             channel_names, good_channels, myDataPath, bids_sub, bids_ses, bids_task, bids_runs, params)
                            end

                            fprintf('File %s has run!\n', replace(filesrun(ii).name, '_events.tsv', ''))
                        end
                    end
                end
            end
        end
    end
end