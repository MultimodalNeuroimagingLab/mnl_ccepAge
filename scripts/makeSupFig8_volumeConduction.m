%
% This script loads CCEPs of one subject and plots a distribution of
% latencies to distances from the stimulus pair
% 
% 
% Dora Hermes, Max van den Boom, Dorien van Blooijs 2022
%


clear
close all
warning('on');
warning('backtrace', 'off')

myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData')
else
    disp('Run scripts ccep02_loadN1.m and ccep03_addtracts.m first')
end

stimStimElec_excludeDist = 18;     % the distance between the stimulated electrodes (in mm) above which N1s are excluded, 0 = not excluding
respStimElec_excludeDist = 13;     % the distance between a stimulated and response electrode (in mm) within which N1s are excluded, 0 = not excluding


%%
%  pull distances for all ROIs

stim_roi = [];
resp_roi = [];

% load tracts and their corresponding end-point ROIs
rois = ccep_categorizeAnatomicalRegions();

% loop over the (sub-)tracts and directions
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract) 
        % direction along tract
        for iDir = [false true] 
            
            stim_roi = [stim_roi, rois(iTr).sub_tract(iSubTr).(['roi', num2str(iDir + 1)])];
            resp_roi = [resp_roi, rois(iTr).sub_tract(iSubTr).(['roi', num2str(~iDir + 1)])];
             
        end
    end
end

% extract the latencies and number of N1s/CCEPs between the end-point ROIs for a specific (sub-)tract and direction
out = ccep_N1sBetweenRegions(ccepData, stim_roi, resp_roi, stimStimElec_excludeDist, respStimElec_excludeDist);

%%

figure('Position',[0 0 1600 600])

subplot(5,1,1:4),hold on

all_varLatencies = [];

for iSubj = 1:length(ccepData) 

    stim_pairs = unique(out(iSubj).StimPairNr);
    for kk = 1:length(stim_pairs)
        this_plot = find(out(iSubj).StimPairNr==stim_pairs(kk));
        if length(this_plot)>5 % several stim/record pairs
            plot(std(1000*out(iSubj).latencies(this_plot)),iSubj,'k.')
            all_varLatencies = [all_varLatencies; std(1000*out(iSubj).latencies(this_plot))];
        end
    end
end
plot([0 0],[0 75],'k')
xlim([-1 50])
set(gca,'XTick', [0:10:50], 'XTickLabel', [])
ylim([0 75])
ylabel('Subject number')

subplot(5, 1, 5), hold on
histogram(all_varLatencies, [0:50], 'FaceColor', [1 1 1]);
% plot([0 0], [0 200], 'k')
xlim([-1 50])
xlabel('Standard deviation (ms)')
ylabel('Number of stimulation pairs')

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS8_volumeConductionCheck');
set(gcf, 'PaperPositionMode', 'auto')
print('-dpng', '-r300', '-painters', figureName)
print('-depsc2', '-r300', '-painters', figureName)



