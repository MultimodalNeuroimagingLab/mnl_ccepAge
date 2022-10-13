%
% This script loads CCEPs of one subject and plots a distribution of
% latencies to distances from the stimulus pair
% 
% Dorien van Blooijs 2022
%

%% Set paths
clc
clear
myDataPath = setLocalDataPath(1);

%% add subject information for five different patients
close all

% bids_sub = 'sub-ccepAgeUMCU50'; % 7y, M
% bids_ses = 'ses-1';
% bids_task = 'task-SPESclin';
% bids_runs = 'run-021222';

% bids_sub = 'sub-ccepAgeUMCU55'; % 14y, M
% bids_ses = 'ses-1';
% bids_task = 'task-SPESclin';
% bids_runs = 'run-031717';

bids_sub = 'sub-ccepAgeUMCU19';  % 25y, F
bids_ses = 'ses-1';
bids_task = 'task-SPESclin';
bids_runs = 'run-040955';

% bids_sub = 'sub-ccepAgeUMCU70'; % 38y, F
% bids_ses = 'ses-1';
% bids_task = 'task-SPESclin';
% bids_runs = 'run-021404';

% bids_sub = 'sub-ccepAgeUMCU59'; % 50y, F
% bids_ses = 'ses-1';
% bids_task = 'task-SPESclin';
% bids_runs = 'run-041501';

load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', bids_sub ,bids_ses, ...
    [bids_sub, '_', bids_ses,'_',bids_task,'_',bids_runs, '_averageCCEPs.mat']));

tb_elec = readtable(fullfile(myDataPath.input,bids_sub,bids_ses,'ieeg',...
    [bids_sub, '_', bids_ses,'_','electrodes.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);

%% calculate distance between stimpair and response electrode

distance = NaN(size(n1_peak_sample));
elecPos = [tb_elec.x tb_elec.y tb_elec.z];

for iStimp = 1:size(n1_peak_sample,2)
    for iChan = 1:size(n1_peak_sample,1)
        if ~isnan(n1_peak_sample(iChan,iStimp))
            stimp1 = extractBefore(stimpair_names{iStimp},'-');
            stimp2 = extractAfter(stimpair_names{iStimp},'-');
            stimpChan1 = find(strcmp(stimp1, tb_elec.name) == 1);
            stimpChan2 = find(strcmp(stimp2, tb_elec.name) == 1);
            
            midStimp = (elecPos(stimpChan1,:) + elecPos(stimpChan2,:))/2;
            distance(iChan,iStimp) = sqrt(sum((midStimp - elecPos(iChan,:)).^2));
        end
    end
end

%% figure of all latencies of 1 subject

all_n1_peak_sample = n1_peak_sample(:);
all_n1_peak_sample(isnan(all_n1_peak_sample)) = [];

figure(1),
histogram(tt(all_n1_peak_sample)*1000,'BinWidth',1)
xlabel('N1 latency (ms)')
ylabel('Number of responses')
title('Distribution of N1 latencies')
xlim([0 100])

% save the image
figureName = fullfile(myDataPath.output,'derivatives','age',...
            [bids_sub,'_histogram_N1latencies']);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc','-r300',figureName)


%% figure of all latencies for different distances from stimulus pair

cmap = jet(size(n1_peak_sample,2));
stimnum = repmat(1:size(n1_peak_sample,2),size(n1_peak_sample,1),1);

all_distance = distance(:);
all_stimnum = stimnum(:);
all_stimnum(isnan(all_distance)) = [];
all_distance(isnan(all_distance)) = [];

figure(2),
hold on,
for iCh = 1:size(all_distance,1)
    plot(all_distance(iCh),tt(all_n1_peak_sample(iCh))*1000,'.','Color',cmap(all_stimnum(iCh),:))
end
hold off
xlabel('Distance (mm)')
ylabel('N1 latency (ms)')
title('Distribution of N1 latencies')
xlim([0 125])
ylim([0 100])

% save the image
figureName = fullfile(myDataPath.output,'derivatives','age',...
            [bids_sub,'_distance_N1latencies']);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc','-r300',figureName)

%% figure

counts =[];
for iStim = 1:size(n1_peak_sample,2)
    [n,edges] = histcounts(tt(all_n1_peak_sample(all_stimnum == iStim))*1000,[5:10:105]);
    
    counts(iStim,:) = n;
end

sz = 1:5:100;

figure(3),
hold on,
for iStim = 1:size(n1_peak_sample,2)
    plot([iStim iStim],[0 100],'Color',[0.9 0.9 0.9])
for iLat = 1:10
    if counts(iStim,iLat) >0
    scatter(iStim,iLat*10,sz(counts(iStim,iLat)),cmap(iStim,:),"filled")
    end
end
end

xlim([0 size(n1_peak_sample,2)+1])
ylim([0 100])
xlabel('Stimulus pair')
ylabel('N1 latency (ms)')

%% figure

minVal =-0.45;
maxVal = 0.45;
randVal = (maxVal - minVal).*rand(size(all_n1_peak_sample,1),1) + minVal;

figure,
plot(all_stimnum+randVal,(tt(all_n1_peak_sample)*1000),'.')

