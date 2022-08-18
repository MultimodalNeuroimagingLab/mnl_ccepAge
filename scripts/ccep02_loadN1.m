%
% Script that creates one output file that includes:
%      - participants details (name, session, runs, age)
%      - electrode information with destrieux labels
%      - detected N1 responses
%      - Freesurfer labels for stimulation and recording pair
%
% Note: the subjects and runs that will be included in the output file are 
%       determined by the dir/file structure that is the result of the 'ccep01_averageCCEPs.m' script.
% 
%
% Dora Hermes, Dorien van Blooijs, Max van den Boom, 2022
%



%% 
%  Set paths
clc
clear
myDataPath = setLocalDataPath(1);


%%
%  Get a list of datasets (output of the 'ccep01_averageCCEPs.m' script)

subjects = ccep_getSubFilenameInfo(myDataPath);


%% 
%  Initialize ccepData in ccepData_init.mat and add subject name, session, runs and age

if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_init.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_init.mat'), 'ccepData')
else
    
    ccepData = [];
    
    % load participants.tsv
    subjectsTsv = readtable(fullfile(myDataPath.input, 'participants.tsv'), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'});
    
    for iSubj = 1:length(subjects)
        disp(['subj ' num2str(iSubj) ' of ' num2str(length(subjects)), ' (', subjects(iSubj).name, ')']);
        
        % retrieve subject name based on dir/file structure
        subjName = subjects(iSubj).name;
        
        % find the corresponding details in the tsv (at first session, if there are multiple)
        [subjTsvIndex] = find(ismember(subjectsTsv.participant_id, subjName), 1);
        assert(~isempty(subjTsvIndex));
        
        % Store the name, session and age
        ccepData(iSubj).id = subjName;
        ccepData(iSubj).ses = subjects(iSubj).ses;
        ccepData(iSubj).age = subjectsTsv.age(subjTsvIndex);
        
        % get number of runs and electrodes info
        ccepData(iSubj).nrRuns = length(subjects(iSubj).run);
        ccepData(iSubj).elecs = readtable(fullfile(myDataPath.input, subjects(iSubj).name, subjects(iSubj).ses, 'ieeg', [subjects(iSubj).name, '_', subjects(iSubj).ses, '_electrodes.tsv']), ...
                                                 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);
        
        % get hemisphere for each electrode and store in the structure
        hemi = ccep_retrieveElecsHemisphere(fullfile(myDataPath.input, subjects(iSubj).name, subjects(iSubj).ses, 'ieeg', [subjects(iSubj).name, '_', subjects(iSubj).ses, '_task-SPESclin*_ieeg.json']), ...
                                            ccepData(iSubj).elecs);
        ccepData(iSubj).elecs.jsonHemi = hemi;
        
        % TODO: check/assert if hemi or hemisphere field in electrodes.tsv reflects jsonHemi
        
        
    end
    
    % optional save the ccepData structure, add more fields later as neccesary
    s = input('Do you want to save the ccepData structure? [y/n]: ', 's');
    if strcmp(s, 'y')
        save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_init.mat'), 'ccepData')
    end
    
end


%%
%  Load all N1 data and add to initialized ccepData_init, then save as ccepData_V1

for iSubj = 1:length(subjects)
    disp(['subj ' num2str(iSubj) ' of ' num2str(length(subjects)), ' (', subjects(iSubj).name, ')']);
    
    for iRun = 1:length(subjects(iSubj).run)
        
        runData = load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', subjects(iSubj).name ,subjects(iSubj).ses, subjects(iSubj).run{iRun}));
        ccepData(iSubj).run(iRun).runName            = subjects(iSubj).run{iRun};
        ccepData(iSubj).run(iRun).allLatencies       = runData.tt(runData.n1_peak_sample(~isnan(runData.n1_peak_sample)));
        ccepData(iSubj).run(iRun).n1_peak_sample     = runData.n1_peak_sample;
        ccepData(iSubj).run(iRun).channel_names      = runData.channel_names;
        ccepData(iSubj).run(iRun).average_ccep_names = runData.average_ccep_names;
        ccepData(iSubj).run(iRun).good_channels      = runData.good_channels;
        % loading all average cceps here makes it very heavy on the memory, do this later
        %ccepData(kk).run(ll).average_ccep            = runData.average_ccep;
        ccepData(iSubj).run(iRun).tt                 = runData.tt;
        clear runData
        
    end
end


%% 
%   Add Freesurfer labels for stimulation and recording pair to output struct

for iSubj = 1:length(ccepData) % loop subjects
   
    for iRun = 1:length(ccepData(iSubj).run) % loop runs
        
        % pre-allocation: Destrieux labels and numbers for average CCEP stimulated pairs
        ccepData(iSubj).run(iRun).average_ccep_DestrieuxLabel    = cell(size(ccepData(iSubj).run(iRun).average_ccep_names,1),2);
        ccepData(iSubj).run(iRun).average_ccep_DestrieuxNr       = cell(size(ccepData(iSubj).run(iRun).average_ccep_names,1),2);
        
        % pre-allocation: Destrieux labels and numbers for measured channels
        ccepData(iSubj).run(iRun).channel_DestrieuxLabel         = cell(size(ccepData(iSubj).run(iRun).channel_names));
        ccepData(iSubj).run(iRun).channel_DestrieuxNr            = cell(size(ccepData(iSubj).run(iRun).channel_names));
        
        % loop through CCEP stimulated pairs
        for chPair = 1:length(ccepData(iSubj).run(iRun).average_ccep_names)
            
            % get stimulated channels
            stimpchans = strsplit(ccepData(iSubj).run(iRun).average_ccep_names{chPair}, '-');
            
            for ch = 1:2
                
                % get first stimulated channel number in_electrodes.tsv
                stim_el_nr = find(strcmpi(ccepData(iSubj).elecs.name, stimpchans{ch}) == 1);
                
                % sometimes the stim pair is called TP1 and the channel name is
                % TP01, we need to check for this
                if isempty(stim_el_nr)
                    
                    % insert a zero and check
                    newName = insertBefore(stimpchans{1},length(stimpchans{1}), '0');
                    stim_el_nr = find(strcmpi(ccepData(iSubj).elecs.name, newName) == 1);
                    if isempty(stim_el_nr)
                        disp(['no match for ' stimpchans{1}])
                    end
                    
                end
                
                ccepData(iSubj).run(iRun).average_ccep_DestrieuxLabel{chPair,ch} = ccepData(iSubj).elecs.Destrieux_label_text{stim_el_nr};
                if isnumeric(ccepData(iSubj).elecs.Destrieux_label)
                    ccepData(iSubj).run(iRun).average_ccep_DestrieuxNr{chPair,ch} = int2str(ccepData(iSubj).elecs.Destrieux_label(stim_el_nr));
                else
                    ccepData(iSubj).run(iRun).average_ccep_DestrieuxNr{chPair,ch} = ccepData(iSubj).elecs.Destrieux_label{stim_el_nr};
                    
                end
            end
            clear stim_el_nr stimpchans % housekeeping
        end
        
        % loop through channels
        for chSig = 1:length(ccepData(iSubj).run(iRun).channel_names)      
            
            % get channel number in_electrodes.tsv
            el1_nr = find(strcmpi(ccepData(iSubj).elecs.name,ccepData(iSubj).run(iRun).channel_names{chSig}) == 1);
            if ~isempty(el1_nr)
                ccepData(iSubj).run(iRun).channel_DestrieuxLabel{chSig} = ccepData(iSubj).elecs.Destrieux_label_text{el1_nr};
                
                if isnumeric(ccepData(iSubj).elecs.Destrieux_label)
                    ccepData(iSubj).run(iRun).channel_DestrieuxNr{chSig} = int2str(ccepData(iSubj).elecs.Destrieux_label(el1_nr));
                else
                    ccepData(iSubj).run(iRun).channel_DestrieuxNr{chSig} = ccepData(iSubj).elecs.Destrieux_label{el1_nr};
                end
                
                clear el1_nr
            else
                ccepData(iSubj).run(iRun).channel_DestrieuxLabel{chSig} = NaN;
                ccepData(iSubj).run(iRun).channel_DestrieuxNr{chSig} = NaN;
            end 
            
        end
    end    
end


% optional save the ccepData structure, add more fields later as neccesary
s = input('Do you want to save the ccepData structure? [y/n]: ', 's');
if strcmp(s, 'y')
    save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'ccepData')
end




%% 
%  Some plots to check:
%
%  plot means, var etc

% initialize output: age, mean and variance in latency per subject
my_output = NaN(length(ccepData) - 1, 3);

% get variable per subject
for iSubj = 1:length(ccepData)
    my_output(iSubj,1) = ccepData(iSubj).age;
    allLatencies = [];
    for iRun = 1:length(ccepData(iSubj).run)
        allLatencies = [allLatencies ccepData(iSubj).run(iRun).allLatencies];
    end
    my_output(iSubj, 2) = mean(allLatencies);
    my_output(iSubj, 3) = var(allLatencies);
    clear allLatencies
end

figure
subplot(2, 1, 1),
plot(my_output(:, 1), 1000 * my_output(:, 2), '.')
xlabel('age (years)'), ylabel('mean latency (ms)')
[r, p] = corr(my_output(:, 1), my_output(:, 2), 'Type', 'Pearson');
title(['r=' num2str(r, 3) ' p=' num2str(p, 3)])

[P, S] = polyfit(my_output(:, 1), 1000 * my_output(:, 2), 1);
[y_fit, ~] = polyval(P, my_output(:, 1), S);
            
hold on
% Plot polyfit throught data points
plot(my_output(:,1), y_fit, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
hold off

subplot(2, 1, 2),
plot(my_output(:, 1),my_output(:, 3), '.')
xlabel('age (years)'), ylabel('variance in latency')
[r, p] = corr(my_output(:, 1), my_output(:, 3), 'Type', 'Pearson');
title(['r=' num2str(r, 3) ' p=' num2str(p, 3)])

[P, S] = polyfit(my_output(:, 1),my_output(:, 3), 1);
[y_fit, ~] = polyval(P, my_output(:, 1), S);
            
hold on
% Plot polyfit throught data points
plot(my_output(:, 1), y_fit, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2)
hold off
sgtitle('Pearson correlation between age and N1-latency')

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'corrAgeVsN1latency');

set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)


%%
%  Plot all under 40

figure
subplot(2, 1, 1),
plot(my_output(my_output(:, 1) < 40, 1),1000 * my_output(my_output(:, 1) < 40, 2),'.')
xlabel('age (years)'), ylabel('mean latency (ms)')
[r, p] = corr(my_output(my_output(:, 1) < 40, 1), my_output(my_output(:, 1) < 40, 2), 'Type', 'Pearson');
title(['r=' num2str(r, 3) ' p=' num2str(p, 3)])

[P, S] = polyfit(my_output(my_output(:, 1) < 40, 1), 1000 * my_output(my_output(:,1) < 40, 2), 1);
[y_fit, ~] = polyval(P, my_output(my_output(:, 1) < 40, 1), S);
            
hold on
% Plot polyfit throught data points
plot(my_output(my_output(:, 1) < 40, 1), y_fit, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2)
hold off

subplot(2,1,2),
plot(my_output(my_output(:, 1) < 40, 1),my_output(my_output(:,1) < 40,3),'.')
xlabel('age (years)'),ylabel('variance in latency')
[r,p] = corr(my_output(my_output(:, 1) < 40, 1), my_output(my_output(:, 1) < 40, 3),'Type','Pearson');
title(['r=' num2str(r, 3) ' p=' num2str(p, 3)])

[P,S] = polyfit(my_output(my_output(:, 1) < 40, 1),my_output(my_output(:, 1) < 40, 3), 1);
[y_fit, ~] = polyval(P, my_output(my_output(:, 1) < 40, 1), S);

hold on
% Plot polyfit throught data points
plot(my_output(my_output(:, 1) < 40, 1), y_fit, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2)
hold on

sgtitle('Pearson correlation between age(<40 years) and N1-latency')

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'corrAgeVsN1latency_40yrs');

set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)



