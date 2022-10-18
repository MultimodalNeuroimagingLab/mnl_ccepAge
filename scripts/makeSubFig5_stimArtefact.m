%
% This script loads CCEPs of one subject and plots responses to stimulation
% Dorien van Blooijs, Dora Hermes 2022
%

%% Set paths
clc
clear
myDataPath = setLocalDataPath(1);

%% Metadata: fill in yourself: run one run of one selected participant

% add subject(s) information
bids_sub = 'sub-ccepAgeUMCU50';
bids_ses = 'ses-1';
bids_task = 'task-SPESclin';
bids_runs = 'run-021222';

files = dir(fullfile(myDataPath.input,bids_sub, bids_ses,'ieeg',...
    [bids_sub '_' bids_ses '_' bids_task '_*'  '_events.tsv']));

electrodes_tsv = read_tsv(fullfile(files(1).folder,[bids_sub, '_' bids_ses,'_electrodes.tsv']));

% load events
events_name = fullfile(files(1).folder,[bids_sub ,'_',...
    bids_ses,'_',bids_task,'_',bids_runs,'_events.tsv']);
ccep_events = readtable(events_name,'FileType','text','Delimiter','\t');

% generate vector for averaging across trials
events_include = ismember(ccep_events.sub_type,{'SPES','SPESclin'});
params.mergeAmp = 1;
params.mergePlusMin = 0; % separate the anodal and cathodal stimuli

[stim_pair_nr,stim_pair_name] = ccep_bidsEvents2conditions(ccep_events,events_include,params);

% load the data, as BrainVision BIDS format
ieeg_name = fullfile(files(1).folder,[bids_sub ,'_',...
    bids_ses,'_',bids_task,'_',bids_runs,'_ieeg.eeg']);
data = ft_read_data(ieeg_name,'dataformat','brainvision_eeg');
data_hdr = ft_read_header(ieeg_name,'dataformat','brainvision_eeg');

channels_tsv_name = replace(fullfile(ieeg_name),'ieeg.eeg', 'channels.tsv');
channels_table = readtable(channels_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

% sampling frequency
srate = data_hdr.Fs;

% list of channel names
channel_names = channels_table.name;

% find good ECoG channels
good_channels = find(ismember(channels_table.type,{'ECOG'}) & ismember(channels_table.status,'good'));

% load data for each condition and save averages

params.epoch_length = 5; % total epoch length in sec, default = 5
params.epoch_prestim_length = 2;%: prestimulus epoch length in sec, default = 2
params.baseline_subtract = 1; % subtract median baseline from each trial

[average_ccep,average_ccep_names,ccep,tt] = ccep_averageConditions(data,srate,ccep_events,channel_names,stim_pair_nr,stim_pair_name,params);

%% detect N1s

params.amplitude_thresh = 3.4;
params.n1_peak_range = 100;
params.srate = srate;

% including no detection possible in bad channels
[n1_peak_sample,n1_peak_amplitude] = ccep_detect_n1peak_ECoG(average_ccep,good_channels,params);

%% plot

respChan = [1:64];

figure,
for stimTrial = [1:2:96]
    for chan = respChan

        if ~isnan(n1_peak_sample(chan,stimTrial)) || ~isnan(n1_peak_sample(chan,stimTrial+1)) 
            plot(tt,squeeze(average_ccep(chan,stimTrial,:)))
            hold on
            plot(tt,squeeze(average_ccep(chan,stimTrial+1,:)))

            hold off
            title(sprintf('Stimpair %s, response Channel %s',average_ccep_names{stimTrial},channel_names{chan}))
            xlim([-0.1 0.5])
            pause
        end
    end
end

%% make figure of selected stimulus pair and response channel

stimTrial = [17,18];
respChan = [15, 20];

maxVal = 1.2 * max(max(max(squeeze(average_ccep(respChan,stimTrial,:)))));
minVal = 1.2 * min(min(min(squeeze(average_ccep(respChan,stimTrial,:)))));

figure(1),
h3 = fill([0 19*1/2048 19/2048 0], ...
    [maxVal maxVal, minVal, minVal],[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
hold on
h1 = plot(tt,squeeze(average_ccep(respChan(1),stimTrial(1),:)),'b');
h2 = plot(tt,squeeze(average_ccep(respChan(1),stimTrial(2),:)),'r');
h4 = plot(tt,mean(squeeze(average_ccep(respChan(1),stimTrial(:),:))),'k','LineWidth',2);
hold off
xlim([-0.1 0.5])
ylim([minVal maxVal])
xlabel('Time (s)')
title(sprintf('response Channel %s',channel_names{respChan(1)}))
legend([h3, h1, h2, h4],'stimulus artefact', ...
    average_ccep_names{stimTrial(1)}, average_ccep_names{stimTrial(2)}, ...
    'mean signal')

% save figure
figureName = sprintf('%sderivatives/av_ccep_figures/supfig5_stimArtefact_%s',...
    myDataPath.output,channel_names{respChan(1)});

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

figure(2),
h3 = fill([0 19*1/2048 19/2048 0], ...
    [maxVal maxVal, minVal, minVal],[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
hold on
h1 = plot(tt,squeeze(average_ccep(respChan(2),stimTrial(1),:)),'b');
h2 = plot(tt,squeeze(average_ccep(respChan(2),stimTrial(2),:)),'r');
h4 = plot(tt,mean(squeeze(average_ccep(respChan(2),stimTrial(:),:))),'k','LineWidth',2);
hold off
xlim([-0.1 0.5])
ylim([minVal maxVal])
xlabel('Time (s)')
title(sprintf('response Channel %s',channel_names{respChan(2)}))
legend([h3,h1, h2, h4],'stimulus artefact', ...
    average_ccep_names{stimTrial(1)}, average_ccep_names{stimTrial(2)}, ...
    'mean signal')

% save figure
figureName = sprintf('%sderivatives/av_ccep_figures/supfig5_stimArtefact_%s',...
    myDataPath.output,channel_names{respChan(2)});

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)


