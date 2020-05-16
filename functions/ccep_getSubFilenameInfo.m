function theseSubs = ccep_getSubFilenameInfo(myDataPath)
% input: bids datapath where derivatives from av_ccep are saved
% returns list with all subject names, sessions and runs
%
% dhermes 2020

theseSubs = [];

% get all subjects
tempSubs = dir(fullfile(myDataPath.output,'derivatives','av_ccep','sub-*'));
for kk = 1:length(tempSubs)
    theseSubs(kk).name = tempSubs(kk).name;
end
clear tempSubs

% for each subject get N1s
for kk = 1:length(theseSubs)
    % get the session (there is only one)
    tempSes = dir(fullfile(myDataPath.output,'derivatives','av_ccep',theseSubs(kk).name,'ses-*'));
    theseSubs(kk).ses = tempSes.name;
    clear tempSes
    
    % get all runs
    tempRun = dir(fullfile(myDataPath.output,'derivatives','av_ccep',theseSubs(kk).name,theseSubs(kk).ses,'sub-*'));
    for ll = 1:length(tempRun)
        theseSubs(kk).run{ll} = tempRun(ll).name;
    end
    clear tempRun
end
