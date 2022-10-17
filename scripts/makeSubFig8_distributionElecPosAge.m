% distributionElecPosAge
% This is an extra analysis to analyse the distribution of electrode
% positions in different age groups.

%   Dorien van Blooijs, UMCU 2022

%% load all N1 latencies
clear
close all

myDataPath = setLocalDataPath(1);

% get a list of datasets
theseSubs = ccep_getSubFilenameInfo(myDataPath);

if exist(fullfile(myDataPath.output,'derivatives','av_ccep','ccepData_V1.mat'),'file')

    % if the ccepData_V1.mat was saved after running ccep02_loadN1, load
    % the ccepData structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','ccepData_V1.mat'))
else
    disp('Run first ccep02_loadN1.mat')
end


%% lobe specification

rois(1).name = 'frontal';
rois(1).regions = [12, 14, 15, 52, 53];
rois(2).name = 'parietal';
rois(2).regions = [25, 26, 27, 56];
rois(3).name = 'temporal';
rois(3).regions = [33, 34, 36, 37, 38, 60, 72, 73, 74];
rois(4).name = 'central';
rois(4).regions = [3, 4, 28, 29, 45];

%% run electrode position distribution for all subjects

% pre-allocation
% overview = struct([]);
agegroups = struct([]);

minAge = 1:10:51;
maxAge = 10:10:60;

for n = 1:size(minAge,2)
    agegroups(n).minAge = minAge(n);
    agegroups(n).maxAge = maxAge(n);
    agegroups(n).elecF = 0;
    agegroups(n).elecP = 0;
    agegroups(n).elecT = 0;
    agegroups(n).elecC = 0;
end

for iSubj = 1:size(ccepData,2)

    % pre-allocation
    good_channels = cell(size(ccepData(iSubj).run));
    countF = 0; countP = 0; countT = 0; countC = 0;

    for iRun = 1:size(ccepData(iSubj).run,2)
        good_channels{iRun} = ccepData(iSubj).run(iRun).channel_names(ccepData(iSubj).run(iRun).good_channels);
    end

    unique_channels = unique(lower(vertcat(good_channels{:})));

    for iElec = 1:size(ccepData(iSubj).electrodes,1)
        if any(strcmpi(ccepData(iSubj).electrodes.name(iElec),unique_channels))

            idx = strcmpi(ccepData(iSubj).electrodes.name(iElec),unique_channels);
            if ismember(ccepData(iSubj).electrodes.Destrieux_label(idx),rois(1).regions)
                countF = countF +1;
            elseif ismember(ccepData(iSubj).electrodes.Destrieux_label(idx),rois(2).regions)
                countP = countP +1;
            elseif ismember(ccepData(iSubj).electrodes.Destrieux_label(idx),rois(3).regions)
                countT = countT +1;
            elseif ismember(ccepData(iSubj).electrodes.Destrieux_label(idx),rois(4).regions)
                countC = countC +1;
            end
        end
    end

%     overview(iSubj).age = ccepData(iSubj).age;
%     overview(iSubj).elecF = countF;
%     overview(iSubj).elecP = countP;
%     overview(iSubj).elecT = countT;
%     overview(iSubj).elecC = countC;

    idxAge = find(ccepData(iSubj).age >= minAge & ccepData(iSubj).age <= maxAge) ;
    agegroups(idxAge).elecF = agegroups(idxAge).elecF + countF;
    agegroups(idxAge).elecP = agegroups(idxAge).elecP + countP;
    agegroups(idxAge).elecT = agegroups(idxAge).elecT + countT;
    agegroups(idxAge).elecC = agegroups(idxAge).elecC + countC;  

end

%% figures

figure, 
bar(mean([agegroups(:).minAge; agegroups(:).maxAge]),[agegroups(:).elecF])
xlabel('Age (years)')
ylabel('Number of electrodes')
title('Frontal lobe')
ylim([0 350])

% save the image
figureName = fullfile(myDataPath.output,'derivatives','age',...
            'ageDistribution_frontal');
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc','-r300',figureName)

figure, 
bar(mean([agegroups(:).minAge; agegroups(:).maxAge]),[agegroups(:).elecT])
xlabel('Age (years)')
ylabel('Number of electrodes')
title('Temporal lobe')
ylim([0 350])

% save the image
figureName = fullfile(myDataPath.output,'derivatives','age',...
            'ageDistribution_temporal');
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc','-r300',figureName)

figure, 
bar(mean([agegroups(:).minAge; agegroups(:).maxAge]),[agegroups(:).elecP])
xlabel('Age (years)')
ylabel('Number of electrodes')
title('Parietal lobe')
ylim([0 350])

% save the image
figureName = fullfile(myDataPath.output,'derivatives','age',...
            'ageDistribution_parietal');
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc','-r300',figureName)

figure, 
bar(mean([agegroups(:).minAge; agegroups(:).maxAge]),[agegroups(:).elecC])
xlabel('Age (years)')
ylabel('Number of electrodes')
title('Central lobe')
ylim([0 350])

% save the image
figureName = fullfile(myDataPath.output,'derivatives','age',...
            'ageDistribution_central');
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc','-r300',figureName)
