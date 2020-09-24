%
% Script to load detected N1 responses, combine all in 1 file and assign
% the Destrieux labels
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
        theseSubs(sub_counter).name = tempSubs(kk).name;
    end
end
clear tempSubs sub_counterjj

% load the participants.tsv
participants_info = readtable(fullfile(myDataPath.input,'participants.tsv'),'FileType','text','Delimiter','\t');

Race = repmat({'n/a'},length(theseSubs),1);
Ethnicity = repmat({'n/a'},length(theseSubs),1);
Gender = repmat({'n/a'},length(theseSubs),1);
Age_Unit = repmat({'years'},length(theseSubs),1);
Age = NaN(length(theseSubs),1);
RESP_nr = repmat({''},length(theseSubs),1);

for kk = 1:length(theseSubs)
    RESP_nr{kk,1} = extractAfter(theseSubs(kk).name,'-');
    
    [~,thisInd] = ismember(RESP_nr{kk,1},participants_info.participant_id);
    Age(kk,1) = participants_info.age(thisInd);
end

demographics_table = table(Race,Ethnicity,Gender,Age,Age_Unit,RESP_nr);

% writetable(demographics_table,'./demographics_table74.tsv','FileType','text','Delimiter','\t')

                                
