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
clear tempSubs sub_counter

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

%% initialize N1latencies and get subject label and age

n1Latencies = [];

% load participants.tsv
sub_info = readtable(fullfile(myDataPath.input,'participants.tsv'),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

for kk = 1:length(theseSubs) 
    disp(['subj ' int2str(kk) ' of ' int2str(length(theseSubs))])
    
    % add subject age
    thisSubName = extractAfter(theseSubs(kk).name,'sub-');
    [thisSubInd] = find(ismember(sub_info.name,thisSubName),1); % first session age
    n1Latencies(kk).id = thisSubName;
    n1Latencies(kk).ses = theseSubs(kk).ses;
    n1Latencies(kk).age = sub_info.age(thisSubInd);
    
    % get number of runs and electrodes info
    n1Latencies(kk).nrRuns = length(theseSubs(kk).run);
    n1Latencies(kk).elecs_tsv = read_tsv(fullfile(myDataPath.input,theseSubs(kk).name,theseSubs(kk).ses,'ieeg',...
        [theseSubs(kk).name,'_',theseSubs(kk).ses,'_electrodes.tsv']));
    
end

% optional save the n1Latencies structure, add more fields later as neccesary
% save(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_init.mat'),'n1Latencies')

%% or skip the previous box and load

load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_init.mat'),'n1Latencies')


%% load all N1 data
for kk = 1:length(theseSubs)
    disp(['subj ' int2str(kk) ' of ' int2str(length(theseSubs))])
        
    for ll = 1:length(theseSubs(kk).run)
        
        clear thisData
        thisRun = fullfile(myDataPath.output,'derivatives','av_ccep',theseSubs(kk).name,theseSubs(kk).ses,...
            theseSubs(kk).run{ll});
        thisData = load(thisRun);
        n1Latencies(kk).run(ll).allLatencies = thisData.tt(thisData.n1_peak_sample(~isnan(thisData.n1_peak_sample)));
        n1Latencies(kk).run(ll).n1_peak_sample = thisData.n1_peak_sample;
        n1Latencies(kk).run(ll).channel_names = thisData.channel_names;
        n1Latencies(kk).run(ll).average_ccep_names = thisData.average_ccep_names;
        n1Latencies(kk).run(ll).good_channels = thisData.good_channels;
        % loading all average cceps here makes it very heavy on the memory,
        %   do this later
        % n1Latencies(kk).run(ll).average_ccep = thisData.average_ccep;
        n1Latencies(kk).run(ll).tt = thisData.tt;
    end
    
end

%% get Freesurfer labels for stimulation and recording pair

for kk = 1:length(n1Latencies) % loop subjects
   
    for ll = 1:length(n1Latencies(kk).run) % loop runs
        
        % pre-allocation: Destrieux labels and numbers for average CCEP stimulated pairs
        n1Latencies(kk).run(ll).average_ccep_DestrieuxLabel = cell(size(n1Latencies(kk).run(ll).average_ccep_names,1),2);
        n1Latencies(kk).run(ll).average_ccep_DestrieuxNr = cell(size(n1Latencies(kk).run(ll).average_ccep_names,1),2);
        
        % pre-allocation: Destrieux labels and numbers for measured channels
        n1Latencies(kk).run(ll).channel_DestrieuxLabel = cell(size(n1Latencies(kk).run(ll).channel_names));
        n1Latencies(kk).run(ll).channel_DestrieuxNr = cell(size(n1Latencies(kk).run(ll).channel_names));
        
        % loop through CCEP stimulated pairs
        for chPair = 1:length(n1Latencies(kk).run(ll).average_ccep_names)
            % get stimulated channels
            stimpchans = strsplit(n1Latencies(kk).run(ll).average_ccep_names{chPair},'-');
            
            for ch = 1:2
                % get first stimulated channel number in_electrodes.tsv
                stim_el_nr = find(strcmpi(n1Latencies(kk).elecs_tsv.name,stimpchans{ch})==1);
                
                % sometimes the stim pair is called TP1 and the channel name is
                % TP01, we need to check for this
                if isempty(stim_el_nr)
                    % insert a zero and check
                    newName = insertBefore(stimpchans{1},length(stimpchans{1}),'0');
                    stim_el_nr = find(strcmpi(n1Latencies(kk).elecs_tsv.name,newName)==1);
                    if isempty(stim_el_nr)
                        disp(['no match for ' stimpchans{1}])
                    end
                end
                
                n1Latencies(kk).run(ll).average_ccep_DestrieuxLabel{chPair,ch} = ...
                    n1Latencies(kk).elecs_tsv.Destrieux_label_text{stim_el_nr};
                if isnumeric(n1Latencies(kk).elecs_tsv.Destrieux_label)
                    n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{chPair,ch} = ...
                        int2str(n1Latencies(kk).elecs_tsv.Destrieux_label(stim_el_nr));
                else
                    n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{chPair,ch} = ...
                        n1Latencies(kk).elecs_tsv.Destrieux_label{stim_el_nr};
                    
                end
            end
            clear stim_el_nr stimpchans % housekeeping
        end
        
        % loop through channels
        for chSig = 1:length(n1Latencies(kk).run(ll).channel_names)            
            % get channel number in_electrodes.tsv
            el1_nr = find(strcmpi(n1Latencies(kk).elecs_tsv.name,n1Latencies(kk).run(ll).channel_names{chSig})==1);
            if ~isempty(el1_nr)
                n1Latencies(kk).run(ll).channel_DestrieuxLabel{chSig} = ...
                    n1Latencies(kk).elecs_tsv.Destrieux_label_text{el1_nr};
                if isnumeric(n1Latencies(kk).elecs_tsv.Destrieux_label)
                    n1Latencies(kk).run(ll).channel_DestrieuxNr{chSig} = ...
                        int2str(n1Latencies(kk).elecs_tsv.Destrieux_label(el1_nr));
                else
                    n1Latencies(kk).run(ll).channel_DestrieuxNr{chSig} = ...
                        n1Latencies(kk).elecs_tsv.Destrieux_label{el1_nr};
                end
                clear el1_nr
            else
                n1Latencies(kk).run(ll).channel_DestrieuxLabel{chSig} = NaN;
                n1Latencies(kk).run(ll).channel_DestrieuxNr{chSig} = NaN;
            end            
        end
    end    
end


%% some plots to check:
%%
%% plot means, var etc

% initialize output: age, mean and variance in latency per subject
my_output = NaN(length(n1Latencies)-1,3);

% get variable per subject
for kk = 1:length(n1Latencies)
    my_output(kk,1) = n1Latencies(kk).age;
    allLatencies = [];
    for ll = 1:length(n1Latencies(kk).run)
        allLatencies = [allLatencies n1Latencies(kk).run(ll).allLatencies]; %#ok<AGROW>
    end
    my_output(kk,2) = mean(allLatencies);
    my_output(kk,3) = var(allLatencies);
    clear allLatencies
end

figure
subplot(2,1,1),hold on
plot(my_output(:,1),1000*my_output(:,2),'.')
xlabel('age (years)'),ylabel('mean latency (ms)')
[r,p] = corr(my_output(:,1),my_output(:,2),'Type','Pearson');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

subplot(2,1,2),hold on
plot(my_output(:,1),my_output(:,3),'.')
xlabel('age (years)'),ylabel('variance in latency')
[r,p] = corr(my_output(:,1),my_output(:,3),'Type','Pearson');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

%% plot all under 40

figure
subplot(2,1,1),hold on
plot(my_output(my_output(:,1)<40,1),1000*my_output(my_output(:,1)<40,2),'.')
xlabel('age (years)'),ylabel('mean latency (ms)')
[r,p] = corr(my_output(my_output(:,1)<40,1),my_output(my_output(:,1)<40,2),'Type','Pearson');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

subplot(2,1,2),hold on
plot(my_output(my_output(:,1)<40,1),my_output(my_output(:,1)<40,3),'.')
xlabel('age (years)'),ylabel('variance in latency')
[r,p] = corr(my_output(my_output(:,1)<40,1),my_output(my_output(:,1)<40,3),'Type','Pearson');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

%%

% optional save the n1Latencies structure, add more fields later as neccesary
save(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
