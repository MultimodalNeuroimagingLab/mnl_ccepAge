% ccep06_only8mASubs
% this is an extra analysis to make sure that the variable pulse current
% (1-8mA or not known) does not have a great impact on the results. 
% We only include data of subjects of which we are certain that 8mA was applied.

%   Dorien van Blooijs, UMCU 2021

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

%% include only the runs in which we are certain that 8mA was applied.

ccepData8ma = struct;

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
            
            ccepData8ma(CountSub).id        = ccepData(n).id;
            ccepData8ma(CountSub).ses       = ccepData(n).ses;
            ccepData8ma(CountSub).age       = ccepData(n).age;
            ccepData8ma(CountSub).electrodes = ccepData(n).electrodes;
            ccepData8ma(CountSub).run(CountRun) = ccepData(n).run(m);  
            
            CountRun = CountRun +1;
            if m == size(theseSubs(n).run,2)
                    CountSub = CountSub +1;
            end
        end
        
    end
end

%% plot means for all subjects and the ones with only 8mA

% initialize output: age, mean and variance in latency per subject
my_output_all = NaN(length(ccepData)-1,3);
my_output_8ma = NaN(length(ccepData8ma)-1,3);

% get variable per subject --> all subjects
for kk = 1:length(ccepData)
    my_output_all(kk,1) = ccepData(kk).age;
    allLatenciesall = [];
    for ll = 1:length(ccepData(kk).run)
        allLatenciesall = [allLatenciesall ccepData(kk).run(ll).allLatencies]; %#ok<AGROW>
    end
    my_output_all(kk,2) = mean(allLatenciesall);
    my_output_all(kk,3) = var(allLatenciesall);
    clear allLatenciesall
end

% get variable per subject --> only 8mA
for kk = 1:length(ccepData8ma)
    my_output_8ma(kk,1) = ccepData8ma(kk).age;
    allLatencies8ma = [];
    for ll = 1:length(ccepData8ma(kk).run)
        allLatencies8ma = [allLatencies8ma ccepData8ma(kk).run(ll).allLatencies]; %#ok<AGROW>
    end
    my_output_8ma(kk,2) = mean(allLatencies8ma);
    my_output_8ma(kk,3) = var(allLatencies8ma);
    clear allLatencies8ma
end

figure
subplot(1,2,1),hold on
plot(my_output_all(:,1),1000*my_output_all(:,2),'.k')
xlabel('age (years)'),ylabel('mean latency (ms)')
[r,p] = corr(my_output_all(:,1),my_output_all(:,2),'Type','Spearman');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])
[P,S] = polyfit(my_output_all(:,1),1000*my_output_all(:,2),1);
[y_fit, ~] = polyval(P,my_output_all(:,1),S);
hold on
% Plot polyfit throught data points
plot(my_output_all(:,1),y_fit,'Color','r','LineWidth',1)
hold off
xlim([0 50]), ylim([10 50])

subplot(1,2,2),hold on
plot(my_output_8ma(:,1),1000*my_output_8ma(:,2),'.k')
xlabel('age (years)'),ylabel('mean latency (ms)')
[r,p] = corr(my_output_8ma(:,1),my_output_8ma(:,2),'Type','Spearman');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])
[P,S] = polyfit(my_output_8ma(:,1),1000*my_output_8ma(:,2),1);
[y_fit, ~] = polyval(P,my_output_8ma(:,1),S);
hold on
% Plot polyfit throught data points
plot(my_output_8ma(:,1),y_fit,'Color','r','LineWidth',1)
hold off
xlim([0 50]), ylim([10 50])

if ~exist(fullfile(myDataPath.output,'derivatives','age'),'dir')
    mkdir(fullfile(myDataPath.output,'derivatives','age'));
end

figureName = fullfile(myDataPath.output,'derivatives','age',...
    'corrAgeVsN1latency_allSub_8mA');

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-depsc','-r300',figureName)

%% plot means, var etc

% initialize output: age, mean and variance in latency per subject
my_output = NaN(length(ccepData8ma)-1,3);

% get variable per subject
for kk = 1:length(ccepData8ma)
    my_output(kk,1) = ccepData8ma(kk).age;
    allLatencies = [];
    for ll = 1:length(ccepData8ma(kk).run)
        allLatencies = [allLatencies ccepData8ma(kk).run(ll).allLatencies]; %#ok<AGROW>
    end
    my_output(kk,2) = mean(allLatencies);
    my_output(kk,3) = var(allLatencies);
    clear allLatencies
end

figure
subplot(1,2,1),hold on
plot(my_output(:,1),1000*my_output(:,2),'.')
xlabel('age (years)'),ylabel('mean latency (ms)')
[r,p] = corr(my_output(:,1),my_output(:,2),'Type','Spearman');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])
[P,S] = polyfit(my_output(:,1),1000*my_output(:,2),1);
[y_fit, ~] = polyval(P,my_output(:,1),S);
hold on
% Plot polyfit throught data points
plot(my_output(:,1),y_fit,'Color',[0.7,0.7,0.7],'LineWidth',2)
hold off

subplot(1,2,2),hold on
plot(my_output(:,1),my_output(:,3),'.')
xlabel('age (years)'),ylabel('variance in latency')
[r,p] = corr(my_output(:,1),my_output(:,3),'Type','Spearman');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])
[P,S] = polyfit(my_output(:,1),my_output(:,3),1);
[y_fit, ~] = polyval(P,my_output(:,1),S);
            
hold on
% Plot polyfit throught data points
plot(my_output(:,1),y_fit,'Color',[0.7,0.7,0.7],'LineWidth',2)
hold off
sgtitle('Pearson correlation between age and N1-latency')

% plot all under 40
figure
subplot(2,1,1),hold on
plot(my_output(my_output(:,1)<40,1),1000*my_output(my_output(:,1)<40,2),'.')
xlabel('age (years)'),ylabel('mean latency (ms)')
[r,p] = corr(my_output(my_output(:,1)<40,1),my_output(my_output(:,1)<40,2),'Type','Spearman');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

subplot(2,1,2),hold on
plot(my_output(my_output(:,1)<40,1),my_output(my_output(:,1)<40,3),'.')
xlabel('age (years)'),ylabel('variance in latency')
[r,p] = corr(my_output(my_output(:,1)<40,1),my_output(my_output(:,1)<40,3),'Type','Spearman');
title(['r=' num2str(r,3) ' p=' num2str(p,3)]) 


%% run other scripts with only patients with 8mA stimulation

% optional save the n1Latencies structure, add more fields later as neccesary
s = input('Do you want to save the n1Latencies structure with subjects in whom only 8mA stimulation is applied? [y/n]: ','s');
if strcmp(s,'y')
    save(fullfile(myDataPath.output,'derivatives','av_ccep','ccepData_V1_8ma.mat'),'ccepData8ma')
end

% run makeFig2_plotResponsesAge
% run makeFig3_plotN1_meanacrossage
