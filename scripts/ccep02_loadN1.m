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

%% run through all N1 data

n1Latencies = [];

for kk = 1%:length(theseSubs) 
    disp(['subj ' int2str(kk) ' of ' int2str(length(theseSubs))])
    
    n1Latencies(kk).nrRuns = length(theseSubs(kk).run);
    n1Latencies(kk).elecs_tsv = read_tsv(fullfile(myDataPath.input,theseSubs(kk).name,theseSubs(kk).ses,'ieeg',...
        [theseSubs(kk).name,'_',theseSubs(kk).ses,'_electrodes.tsv']));
        
    for ll = 1%:length(theseSubs(kk).run)
        
        clear thisData
        thisRun = fullfile(myDataPath.output,'derivatives','av_ccep',theseSubs(kk).name,theseSubs(kk).ses,...
            theseSubs(kk).run{ll});
        thisData = load(thisRun);
        n1Latencies(kk).run(ll).allLatencies = thisData.tt(thisData.n1_peak_sample(~isnan(thisData.n1_peak_sample)));
        n1Latencies(kk).run(ll).n1_peak_sample = thisData.n1_peak_sample;
        n1Latencies(kk).run(ll).channel_names = thisData.channel_names;
        n1Latencies(kk).run(ll).average_ccep_names = thisData.average_ccep_names;
        n1Latencies(kk).run(ll).good_channels = thisData.good_channels;
        % loading all average cceps makes it too heavy on the memory
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
                n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{chPair,ch} = ...
                    n1Latencies(kk).elecs_tsv.Destrieux_label{stim_el_nr};
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
                n1Latencies(kk).run(ll).channel_DestrieuxNr{chSig} = ...
                    n1Latencies(kk).elecs_tsv.Destrieux_label{el1_nr};
                clear el1_nr
            else
                n1Latencies(kk).run(ll).channel_DestrieuxLabel{chSig} = NaN;
                n1Latencies(kk).run(ll).channel_DestrieuxNr{chSig} = NaN;
            end            
        end
    end    
end


%% get ages

% load participants.tsv
sub_info = readtable(fullfile(myDataPath.input,'participants.tsv'),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    
for kk = 1:length(theseSubs)
    thisSubName = extractAfter(theseSubs(kk).name,'sub-');
    [thisSubInd] = find(ismember(sub_info.name,thisSubName),1); % first session age
    n1Latencies(kk).age = sub_info.age(thisSubInd);
end


%% here, we should probably save the N1 latency structure

save(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies.mat'),'n1Latencies')

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
% temporal areas:
% G_temporal_inf, G_temporal_middle
% roi = {'37','38'};
% G_temporal_inf, G_temporal_middle, G_temp_sup-Lateral,
% G_oc-temp_med-Parahip, G_oc-temp_lat-fusifor
roi{1} = {'37','38','34','23','21'};
roi_name{1} = 'temporal';
% frontal areas:
% G_front_inf-Triangul, G_front_middle, G_front_inf-Opercular
roi{2} = {'14';'15';'12'}; % maybe add 16: G_front_sup
roi_name{2} = 'frontal';% parietal areas:
% G_pariet_inf-Angular, G_pariet_inf-Supramar, G_parietal_sup
roi{3} = {'25','26','27'};
roi_name{3} = 'parietal';% occipital areas:
% G_occipital_middle G_oc-temp_med-Lingual Pole_occipital
roi{4} = {'19','22','42'};
roi_name{4} = 'occipital';% sensorimotor:
% G_postcentral G_precentral S_central
roi{5} = {'28','29','46'};
roi_name{5} = 'central';

average_ccep_age = cell(max([n1Latencies.age]),1);
for kk = 1: size(n1Latencies,2)
    
   age = n1Latencies(kk).age;
   average_ccep_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2),5*2048); % [roi_start, roi_end,run,tt]
   
   for ll = 1:size(n1Latencies(kk).run,2)
       
       if size(n1Latencies(kk).run(ll).average_ccep,3) == 5*2048 % ERROR WHEN fs is not 2048Hz....
           
           for rr = 1:size(roi,2)
               % find stimulation pair within specific region
               chanPair = find(sum(contains(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr,roi{rr}),2)>0);
               for rr2 = 1:size(roi,2)
                   % find response electrode within specific region
                   chanSig = find(ismember(str2double(n1Latencies(kk).run(ll).channel_DestrieuxNr),str2double(roi{rr2}))>0);
                   
                   % collect all signals with stimulation pair and response
                   % electrode within specific region
                   
                   average_ccep_select = NaN(size(chanPair,1),size(chanSig,2),5*2048);
                   for cp = 1:size(chanPair,1)
                       for cs = 1:size(chanSig,1)
                           if ~isnan(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)))
                                average_ccep_select(cp,cs,:) = n1Latencies(kk).run(ll).average_ccep(chanSig(cs),chanPair(cp),:);
                           end
                       end
                   end
                               
                   average_ccep_run(rr,rr2,ll,:) = squeeze(mean(mean(average_ccep_select,'omitnan'),'omitnan'))';
               end
           end
       end
   end
   
   average_ccep_pat = squeeze(mean(average_ccep_run,3,'omitnan'));
   
   if ~isempty(average_ccep_age{age})
       n = size(average_ccep_age{age},3);
       average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
   else
       n = 0;
       average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
   end
end
      

average_ccep_age_mean = cell(size(average_ccep_age))   ; 
for age = 1: max([n1Latencies.age])
       average_ccep_age_mean{age} = mean(average_ccep_age{age},3,'omitnan');
end
    
%% make figure with all ccep signals
% define rois
rr = 1;
rr2 = 5;

tt = n1Latencies(74).run(1).tt;
ttmin = -0.02;
ttmax = 0.1;
amp = 500;
ymin = (min([n1Latencies.age])-1)*amp;
ymax = (max([n1Latencies.age])+1)*amp;

figure(1),
hold on
for age = 1:max([n1Latencies.age])
    plot(tt(tt>ttmin & tt< ttmax),zeros(size(tt(tt>ttmin & tt<ttmax)))+amp*age,'Color',[.8 .8 .8])
    if ~isempty(average_ccep_age_mean{age})
        plot(tt(tt>ttmin & tt<ttmax),squeeze(average_ccep_age_mean{age}(rr,rr2,:,tt>ttmin & tt<ttmax))+amp*age)
    end
end
hold off

xlabel('Time (s)')
ylabel('Age (years)')
ylim([ymin,ymax])
ax = gca;
ax.YTick = (1:max([n1Latencies.age]))*amp;
ax.YTickLabel = num2cell(1:max([n1Latencies.age]));
title(sprintf('N1 from %s to %s in increasing age',roi_name{rr},roi_name{rr2}))

