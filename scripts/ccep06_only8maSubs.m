% ccep06_only8mASubs
% this is an extra analysis to make sure that the variable pulse current
% (4mA/8mA or not known) does not have a great impact on the results. 
% We only include data of subjects of which we are certain that 8mA was applied.

%% load all N1 latencies

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
    
    % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
else
    disp('Run first ccep02_loadN1.mat')
end

%% include only the runs in which we are certain that 8mA was applied.

n1Latencies8ma = struct;

CountSub = 1;
for n=1:size(theseSubs,2)
    CountRun = 1;
    
    for m = 1:size(theseSubs(n).run,2)
        
        % load events.tsv
        events_tsv = read_tsv(fullfile(myDataPath.input, theseSubs(n).name, theseSubs(n).ses,'ieeg',...
            replace(theseSubs(n).run{m},'_averageCCEPs.mat','_events.tsv')));
                
        % find events of stimulation
        idx =  ismember(events_tsv.sub_type,{'SPES','SPESclin'}) & ismember(events_tsv.trial_type,{'electrical_stimulation'});
        
        if sum(idx) == 0
           warning('%s does not have any stimulation events',replace(theseSubs(n).run{m},'_averageCCEPs.mat','_events.tsv'))
        end
        
        if iscell(events_tsv.electrical_stimulation_current)
            stimcur = str2double(events_tsv.electrical_stimulation_current(idx));
        else
            stimcur = events_tsv.electrical_stimulation_current(idx);
        end
        
        % include only the runs with 8mA that are certain
        if all(~contains(events_tsv.notes(idx),'Stimulation intensity is suggested to be 0.008 A but may differ when applied in eloquent tissue')) && ...
                all(stimcur == 0.008)
            
            n1Latencies8ma(CountSub).id        = n1Latencies(n).id;
            n1Latencies8ma(CountSub).ses       = n1Latencies(n).ses;
            n1Latencies8ma(CountSub).age       = n1Latencies(n).age;
            n1Latencies8ma(CountSub).elecs_tsv = n1Latencies(n).elecs_tsv;
            n1Latencies8ma(CountSub).run(CountRun) = n1Latencies(n).run(m);  
            
            CountRun = CountRun +1;
            if m == size(theseSubs(n).run,2)
                    CountSub = CountSub +1;
            end
        end
        
    end
end

%% change n1Latencies in only the subjects with 8mA to run the rest of the code

n1Latenciesall = n1Latencies;
n1Latencies = n1Latencies8ma;

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

% plot all under 40
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

%% run other scripts with only patients with 8mA stimulation

ccep03_plotN1_meanacrossage

