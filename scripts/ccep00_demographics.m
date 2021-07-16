%
% Script to show a demographic overview of the included patients in this
% study. This is based on the patients from whom the detected N1s are saved
% in a folder called derivatives/av_ccep/. This is performed in script
% ccep01_averageCCEPs.m
% 
% Dora Hermes, Dorien van Blooijs, 2020
%

%% Set paths
clc
clear
myDataPath = setLocalDataPath(1);

%% get a list of datasets

% list with all subject names, sessions and runs
theseSubs = [];

% get all subjects
tempSubs = dir(fullfile(myDataPath.output,'derivatives','av_ccep'));
sub_counter = 0;
for kk = 1:length(tempSubs)
    if strcmp(extractBefore(tempSubs(kk).name,'-'),'sub')
        sub_counter = sub_counter+1;
        theseSubs(sub_counter).name = tempSubs(kk).name; %#ok<SAGROW>
    end
end
clear tempSubs sub_counterjj

% load the participants.tsv
participants_info = readtable(fullfile(myDataPath.input,'participants.tsv'),'FileType','text','Delimiter','\t');

% copy the information in participants.tsv if the patient is included in
% this study
Race = repmat({'n/a'},length(theseSubs),1);
Ethnicity = repmat({'n/a'},length(theseSubs),1);
Gender = repmat({'n/a'},length(theseSubs),1);
Age_Unit = repmat({'years'},length(theseSubs),1);
Age = NaN(length(theseSubs),1);
RESP_nr = repmat({''},length(theseSubs),1);

for kk = 1:length(theseSubs)
    RESP_nr{kk,1} = theseSubs(kk).name;
    
    [~,thisInd] = ismember(RESP_nr{kk,1},participants_info.participant_id);
    Age(kk,1) = participants_info.age(thisInd);
    Gender(kk,1) = participants_info.sex(thisInd);
end

demographics_table = table(Race,Ethnicity,Gender,Age,Age_Unit,RESP_nr);

% print some information of the included patients in the command window
fprintf('Mean age + SD = %2.1f +- %2.1f \n',...
    mean(Age), std(Age))

fprintf('Median age (min - max) = %2.1f (%2.1f - %2.1f) \n',...
    median(Age), min(Age), max(Age))

fprintf('Gender: %1.0f Male, %1.0f Female \n',...
    sum(strcmpi(Gender,'male')), sum(strcmpi(Gender,'female')))

% writetable(demographics_table,'./demographics_table74.tsv','FileType','text','Delimiter','\t')

                                
