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

%% Metadata: fill in yourself

% add subject(s) information
bids_sub = ['sub-' input('Patient number: sub- (RESPXXXX): ','s')];
bids_ses = ['ses-' input('Session number: ses- (X): ','s')];
bids_task = 'task-SPESclin';

% choose between available runs
files = dir(fullfile(myDataPath.input,bids_sub, bids_ses,'ieeg',...
    [bids_sub '_' bids_ses '_' bids_task '_*'  '_events.tsv']));
names = {files.name};
strings = cellfun(@(x) x(strfind(names{1},'run-'):strfind(names{1},'run-')+9), names, 'UniformOutput', false);
stringsz = [repmat('%s, ',1,size(strings,2)-1),'%s'];

bids_runs = input(sprintf(['Choose one of these runs: \n' stringsz '\n'],strings{:}),'s');

clear files names strings stringsz

%% load events

% load the events.tsv
events_name = fullfile(myDataPath.input,bids_sub, bids_ses,'ieeg',...
    [bids_sub '_' bids_ses '_' bids_task '_' bids_runs '_events.tsv']);
ccep_events = readtable(events_name,'FileType','text','Delimiter','\t');

% generate vector for averaging across trials
% TODO: add an exclusion of trials with noise yet !
% TODO: add option to separate based on stimulation currents?
events_include = ismember(ccep_events.sub_type,'SPES');
params.mergePlusMin = 1;
[stim_pair_nr,stim_pair_name] = ccep_bidsEvents2conditions(ccep_events,events_include,params);

%% load data and channels

% load the data, as BrainVision BIDS format
ieeg_name = fullfile(myDataPath.input, bids_sub, bids_ses,'ieeg',...
    [ bids_sub '_' bids_ses '_' bids_task '_' bids_runs '_ieeg.eeg']);
data = ft_read_data(ieeg_name,'dataformat','brainvision_eeg');
data_hdr = ft_read_header(ieeg_name,'dataformat','brainvision_eeg');

channels_tsv_name = fullfile(myDataPath.input,bids_sub, bids_ses,'ieeg',...
    [ bids_sub '_' bids_ses ...
    '_' bids_task '_' bids_runs '_channels.tsv']);
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

% TODO: remove bad channels
[average_ccep,average_ccep_names,tt] = ccep_averageConditions(data,srate,ccep_events,channel_names,stim_pair_nr,stim_pair_name,params);

saveName = fullfile(myDataPath.output,'derivatives','av_ccep',bids_sub,bids_ses,...
    [ bids_sub '_' bids_ses '_' bids_task '_' bids_runs '_averageCCEPs.mat']);

if ~exist(fullfile(myDataPath.output,'derivatives','av_ccep',bids_sub,bids_ses),'dir')
    mkdir(fullfile(myDataPath.output,'derivatives','av_ccep',bids_sub,bids_ses))
    sprintf(['making dir:\n',...
        fullfile(myDataPath.output,'derivatives','av_ccep',bids_sub,bids_ses)])
end

save(saveName,'average_ccep','average_ccep_names','tt','channel_names','good_channels')

%% detect N1 in each averaged signal

params.amplitude_thresh = 2.6;
params.n1_peak_range = 100;
params.srate = srate;

[n1_peak_sample,n1_peak_amplitude] = ccep_detect_n1peak_ECoG(average_ccep,params);


%% plot and save averages per channel

ccep_plot_av(average_ccep,tt,n1_peak_sample, n1_peak_amplitude,average_ccep_names,channel_names,good_channels,myDataPath,bids_sub,bids_ses,bids_task,bids_runs)
