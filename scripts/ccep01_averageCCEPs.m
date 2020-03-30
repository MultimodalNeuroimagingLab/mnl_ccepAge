%
% This script can be used as workflow script to create average CCEPs 
% for the CCEP data in the RESPect database.


% Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, 2019 

%% Set paths
% set rootpath
dataRootPath = fullfile('~/Documents/data/ccep/');

% add folder with functions and scripts
addpath(fullfile('~/Documents','git','ccep'))
% add ecogBasicCode to path
addpath(fullfile('~/Documents','git','ecogBasicCode')) 
% be sure the JSONio toolbox is added to the path
addpath(fullfile('~/Documents','git','JSONio'))


%% Metadata: fill in yourself

% add subject(s) information
subjects = {'RESP0574'};
sessions = {'1'};
tasks = {'SPESclin'};
runs = {'021358'};

% do I need this?
hemi_smalls = {'l'};
hemi_caps = {'L'};
hemi_small = hemi_smalls{s};
hemi_cap = hemi_caps{s};


%% load data

s = 1;
ccepSet.subj = subjects{s};
ccepSet.ses = sessions{s};
ccepSet.task = tasks{s};
ccepSet.run = runs{s};

% load the data, as BrainVision BIDS format
ieeg_name = fullfile(dataRootPath,['sub-' ccepSet.subj], ...
    ['ses-' ccepSet.ses],'ieeg',...
    ['sub-' ccepSet.subj '_ses-' ccepSet.ses ...
    '_task-' ccepSet.task '_run-' ccepSet.run '_ieeg.eeg']);
data = ft_read_data(ieeg_name,'dataformat','brainvision_eeg');
data_hdr = ft_read_header(ieeg_name,'dataformat','brainvision_eeg');

% load the events.tsv
events_name = fullfile(dataRootPath,['sub-' ccepSet.subj], ...
    ['ses-' ccepSet.ses],'ieeg',...
    ['sub-' ccepSet.subj '_ses-' ccepSet.ses ...
    '_task-' ccepSet.task '_run-' ccepSet.run '_events.tsv']);
ccep_events = readtable(events_name,'FileType','text','Delimiter','\t');

% % load the electrodes.tsv
% electrodes_name = fullfile(dataRootPath,['sub-' ccepSet.subj], ...
%     ['ses-' ccepSet.ses],'ieeg',...
%     ['sub-' ccepSet.subj '_ses-' ccepSet.ses ...
%     '_electrode_positions_fouratlases.tsv']);
% electrodes_table = readtable(electrodes_name,'Filetype','text','Delimiter','\t', 'ReadVariableNames', true);

%% epoch

% set epoch parameters
epoch_prestim_length = round((3*data_hdr.Fs))-1; % 3 seconds pre stimulation to samples
epoch_poststim_length = round((5*data_hdr.Fs)); % 5 seconds post stimulation

% count how many stimulations there are
total_stim_count = sum(strcmp(ccep_events.trial_type,'electrical_stimulation'));

% define time vector for all epochs (for plotting and knowing when tt==0 for onset stimulation)
tt = [-epoch_prestim_length : epoch_poststim_length]./data_hdr.Fs;

% create ccep_events_onlystims, which makes table of only the
% stimulations and not e.g. artifact data
ll_counter = 1;
ccep_events_onlystims = table();
for ll = 1:length(ccep_events.onset)
    % only include electrical stimulation trials:
    if strcmp(ccep_events.trial_type(ll),'electrical_stimulation')
        % remove stimulations with less than prestim length before and 
        % stimulations with less than poststim length after
        if ccep_events.sample_start(ll)-epoch_prestim_length>0 && ...
                ccep_events.sample_start(ll)+epoch_poststim_length<size(data,2)
            ccep_events_onlystims(ll_counter,:) = ccep_events(ll,:);
            ll_counter = ll_counter + 1;
        end
    end
end

% define the output structure
data_epoch = NaN(ll_counter,size(data,1),epoch_poststim_length+epoch_prestim_length+1);

% loop through all stimulations and add to the output structure
for ll = 1:length(ccep_events_onlystims.onset) % for all epochs
    data_epoch(ll,:,:) = ...
        data(:,ccep_events_onlystims.sample_start(ll)-epoch_prestim_length:...
        ccep_events_onlystims.sample_start(ll)+epoch_poststim_length);
    disp(['data epoched for epoch ' int2str(ll) ' of ' int2str(length(ccep_events_onlystims.onset))])
end

%% matrix with stimulation sites:

% loop through stimulation events
for ee = 1:size(ccep_events_onlystims,1)
    %extract the first stimulation site and add as column
    ccep_events_onlystims.electrical_stimulation_site_num_1(ee) = str2double(extractBefore(ccep_events_onlystims.electrical_stimulation_site_num(ee),'  '));
    % extract the second stimulation site and add as column
    ccep_events_onlystims.electrical_stimulation_site_num_2(ee) = str2double(extractAfter(ccep_events_onlystims.electrical_stimulation_site_num(ee),'  '));
end

%% plot ERP

% stimulated electrode #1
stim_el = 29;

% measured electrode
meas_elec = 31;

% find epochs with stim_el
plot_epoch_nrs = find(ccep_events_onlystims.electrical_stimulation_site_num_1==stim_el);

figure
plot(tt,squeeze(data_epoch(plot_epoch_nrs,meas_elec,:)))

% plot ERSP

addpath(genpath('/Users/dora/Documents/m-files/Chronux/'))

% Make a nice colormap
cm1 = [repmat([0 0 0],100,1)];
cm1(1:40,1) = [0.7]';
cm1(1:40,2) = [0.7:-0.6/39:0.1]';
cm1(1:40,3) = [0.7:-0.7/39:0]';
cm1(40:100,1) = [0.7:(1-0.7)/60:1]';
cm1(40:100,2) = [0.1:.9/60:1]';
cm2 = [repmat([0 0 0],100,1)];
cm2(1:30,3) = [0.7]';
cm2(1:30,1) = [0.7:-0.7/29:0]';
cm2(1:30,2) = [0.7:-0.7/29:0]';
cm2(30:100,3) = [0.7:(1-0.7)/70:1]';
cm2(30:100,2) = [0:1/70:1]';
cm = [cm2(end:-1:1,:); cm1];

params.pad = -1; % default-1
movingwin = [.500 .05]; % default [.200 .05]
params.tapers = [3 5]; % default [3 5]
params.fpass = [0 200];
params.Fs = data_hdr.Fs;
params.trialave = 0;

data_temp = squeeze(data_epoch(plot_epoch_nrs,meas_elec,:));

% calculate baseline
params.trialave = 1;
[S1b,t_tf_b,f_b] = mtspecgramc(data_temp(:,tt>-2 & tt<-0.5)',movingwin,params);
S1b = mean(S1b,1);

% calculate spectgram trial
params.trialave = 0;
[S1,t_tf,f] = mtspecgramc(data_temp',movingwin,params);
t_tf = t_tf+tt(1);
% normalize wrt baseline
S1 = nanmean(S1,3);
S1 = S1./repmat(S1b,[size(S1,1),1]);

figure('Color',[1 1 1],'Position',[0 0 600 300])
subplot(4,1,1)
plot(tt,squeeze(data_epoch(plot_epoch_nrs,meas_elec,:)))
xlim([t_tf(1) t_tf(end)]),ylim([-800 800])
colorbar
title(['stim el 1 = ' int2str(stim_el) ' measure el = ' int2str(meas_elec)])

subplot(4,1,2:4)
imagesc(t_tf,f,log10(S1)',[-1.5 1.5])
axis xy
colormap(cm)

colorbar

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',fullfile(dataRootPath,'derivatives','figures',...
%     ['ERSP_sub-' ccepSet.subj '_stimel' int2str(stim_el) 'measureel' int2str(meas_elec)]))
% print('-depsc','-r300',fullfile(dataRootPath,'derivatives','figures',...
%     ['ERSP_sub-' ccepSet.subj '_stimel' int2str(stim_el) 'measureel' int2str(meas_elec)]))
