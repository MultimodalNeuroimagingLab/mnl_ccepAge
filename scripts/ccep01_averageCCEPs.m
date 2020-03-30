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

[average_ccep,average_ccep_names,tt] = ccep_averageConditions(data,srate,ccep_events,channel_names,stim_pair_nr,stim_pair_name,params);

%%
saveName = fullfile(myDataPath,'derivatives','av_ccep',['sub-' bids_sub],...
    ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_runs{run_nr} '_averageCCEPs.mat']);

if ~exist(fullfile(myDataPath,'derivatives','av_ccep',['sub-' bids_sub]),'dir')
    mkdir(fullfile(myDataPath,'derivatives','av_ccep',['sub-' bids_sub]))
    sprintf(['making dir:\n',...
        fullfile(myDataPath,'derivatives','av_ccep',['sub-' bids_sub])])
end

save(saveName,'average_ccep','average_ccep_names','tt','channel_names','good_channels')

%%

elnrs_plot = good_channels;

for ll = 20%:length(elnrs_plot)
    el_plot = elnrs_plot(ll);
    figure('Position',[0 0 700 700]),hold on
    for kk = 1:length(average_ccep_names)        
        this_ccep_plot = squeeze(average_ccep(el_plot,kk,:));
%         this_ccep_plot(tt>-0.010 & tt<0.010) = NaN;
        
        plot(tt,kk*500+zeros(size(tt)),'Color',[.8 .8 .8])
        plot(tt,kk*500+this_ccep_plot)
    end
    xlim([-.2 1])
    set(gca,'YTick',500*[1:length(average_ccep_names)],'YTickLabel',average_ccep_names)
    title([channel_names{el_plot}])
    
    ylabel('stimulated electrodes')
    xlabel('time(s)')
    
    % add amplitude bar
    plot([0.9 0.9],[1000 1500],'k','LineWidth',2)
    text(0.91,1250,['500 ' native2unicode(181,'latin1') 'V'])

    % filename
    figureName = fullfile(myDataPath,'derivatives','av_ccep_figures',['sub-' bids_sub],...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_runs{run_nr} '_incomingCCEP_el' channel_names{el_plot}]);
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',figureName)
    print('-depsc','-r300',figureName)
    close all
end