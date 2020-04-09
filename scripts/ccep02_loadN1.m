%
% Script to load detecte N1 responses.
% 
% Dora Hermes, Dorien van Blooijs, 2020
%

%% Set paths
clc
clear
myDataPath = setLocalDataPath(1);

%% get a list of datasets

% list with all subject names
theseSubs = [];
tempSubs = dir(fullfile(myDataPath.output,'derivatives','av_ccep'));
sub_counter = 0;
for kk = 1:length(tempSubs)
    if strcmp(extractBefore(tempSubs(kk).name,'-'),'sub')
        sub_counter = sub_counter+1;
        theseSubs(sub_counter).name = tempSubs(kk).name;
    end
end
clear tempSubs sub_counter
%%
% for each subject get N1s
for kk = 1:length(theseSubs)
    % get the first session
    tempSes = dir(fullfile(myDataPath.output,'derivatives','av_ccep',theseSubs(kk).name));
    for ll = 1:length(tempSes)
        if strcmp(extractBefore(tempSes(ll).name,'-'),'ses')
            theseSubs(kk).ses = tempSes(ll).name;
            break
        end
    end
    clear tempSes
    
    % get all runs
    tempRun = dir(fullfile(myDataPath.output,'derivatives','av_ccep',theseSubs(kk).name,theseSubs(kk).ses));
    run_counter = 0;
    for ll = 1:length(tempRun)
        if strcmp(extractBefore(tempRun(ll).name,'-RESP'),'sub')
            run_counter = run_counter+1;
            theseSubs(kk).run{run_counter} = tempRun(ll).name;
        end
    end
    clear tempSes

    
end

