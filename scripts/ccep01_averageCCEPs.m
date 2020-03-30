%
% This script can be used as workflow script to create average CCEPs 
% for the CCEP data in the RESPect database.
%
% Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, 2019 
%

%% Set paths

myDataPath = setLocalDataPath(1);


%% Metadata: fill in yourself

% add subject(s) information
bids_sub = 'RESP0276';
bids_ses = '1';
bids_task = 'SPESclin';
bids_runs = {'021448'};


%% load events

run_nr = 1;

% load the events.tsv
events_name = fullfile(myDataPath,['sub-' bids_sub], ...
    ['ses-' bids_ses],'ieeg',...
    ['sub-' bids_sub '_ses-' bids_ses ...
    '_task-' bids_task '_run-' bids_runs{run_nr} '_events.tsv']);
ccep_events = readtable(events_name,'FileType','text','Delimiter','\t');

% generate vector for averaging across trials
% TODO: add an exclusion of trials with noise yet !
events_include = ismember(ccep_events.sub_type,'SPES');
[stim_pair_nr,stim_pair_name] = ccep_bidsEvents2conditions(ccep_events,events_include);


%% load data and channels

% load the data, as BrainVision BIDS format
ieeg_name = fullfile(myDataPath,['sub-' bids_sub], ...
    ['ses-' bids_ses],'ieeg',...
    ['sub-' bids_sub '_ses-' bids_ses ...
    '_task-' bids_task '_run-' bids_runs{run_nr} '_ieeg.eeg']);
data = ft_read_data(ieeg_name,'dataformat','brainvision_eeg');
data_hdr = ft_read_header(ieeg_name,'dataformat','brainvision_eeg');

channels_tsv_name = fullfile(myDataPath,['sub-' bids_sub], ...
    ['ses-' bids_ses],'ieeg',...
    ['sub-' bids_sub '_ses-' bids_ses ...
    '_task-' bids_task '_run-' bids_runs{run_nr} '_channels.tsv']);
channels_table = readtable(channels_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

%% get necessary parameters from data

% sampling frequency
srate = data_hdr.Fs;

% list of channel names
channel_names = channels_table.name;

% find good sEEG/ECoG channels
good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));

%% load data for each condition and save averages

params.epoch_length = 5; % total epoch length in sec, default = 5
params.epoch_prestim_length = 2;%: prestimulus epoch length in sec, default = 2
params.baseline_subtract = 0; % subtract median baseline from each trial

[average_ccep,average_ccep_names] = ccep_averageConditions(data,srate,ccep_events,stim_pair_nr,stim_pair_name,params);

