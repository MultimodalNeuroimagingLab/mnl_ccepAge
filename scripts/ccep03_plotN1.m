%% load the n1Latencies from the derivatives

% if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')

%% connections from one region to another

region_start = input('Choose roi where connections start [temporal, frontal, parietal, occipital]: ','s');
region_end = input('Choose roi where connections end [temporal, frontal, parietal, occipital]: ','s');

% region_start = 'occipital';
% region_end = 'occipital';

out = ccep_connectRegions(n1Latencies,region_start,region_end);

% initialize output: age, mean and variance in latency per subject
my_output = NaN(length(out.sub),3);

% get variable per subject
for kk = 1:length(out.sub)
    my_output(kk,1) = out.sub(kk).age;
    my_output(kk,2) = nanmean(out.sub(kk).latencies);
    my_output(kk,3) = nanvar(out.sub(kk).latencies);
end

figure
plot(my_output(:,1),1000*my_output(:,2),'k.','MarkerSize',10)
xlabel('age (years)'),ylabel('mean dT (ms)')
[r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),2),'Type','Pearson');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

%% figure latency connections from one region to another
clear out
regions_connect = {'temporal','central','parietal','frontal'};
conn_matrix = [1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];

for outInd = 1:size(conn_matrix,1)
    region_start = regions_connect{conn_matrix(outInd,1)};
    region_end = regions_connect{conn_matrix(outInd,2)};
    temp = ccep_connectRegions(n1Latencies,region_start,region_end);
    out{outInd} = temp;
    out{outInd}.name = {region_start,region_end};
end

%%
figure('position',[0 0 700 700])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean and variance in latency per subject
    my_output = NaN(length(out{outInd}.sub),3);

    % get variable per subject
    for kk = 1:length(out{outInd}.sub)
        my_output(kk,1) = out{outInd}.sub(kk).age;
        my_output(kk,2) = nanmean(out{outInd}.sub(kk).latencies);
        my_output(kk,3) = nanvar(out{outInd}.sub(kk).latencies);
    end
    
    subplot(4,4,outInd),hold on
    plot(my_output(:,1),1000*my_output(:,2),'k.','MarkerSize',10)
    xlabel('age (years)'),ylabel('mean dT (ms)')
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),2),'Type','Pearson');
%     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    title(['r=' num2str(r,3) ' p=' num2str(p,3)])
    xlim([0 60]), ylim([0 100])
    
    % Yeatman et al., fit a second order polynomial:
    % y  = w1* age^2 * w2*age + w3
    [P,S] = polyfit(my_output(~isnan(my_output(:,2)),1),1000*my_output(~isnan(my_output(:,2)),2),2);
    x_age = [1:1:50];
    y_fit = P(1)*x_age.^2 + P(2)*x_age + P(3);
    plot(x_age,y_fit,'r')

    % Let's fit a first order polynomial:
    % y  =  w1*age + w2
    [P,S] = polyfit(my_output(~isnan(my_output(:,2)),1),1000*my_output(~isnan(my_output(:,2)),2),1);
    x_age = [1:1:50];
    y_fit = P(1)*x_age + P(2);
    plot(x_age,y_fit,'r')
end

if ~exist(fullfile(myDataPath.output,'derivatives','age'),'dir')
    mkdir(fullfile(myDataPath.output,'derivatives','age'));    
end

% figureName = fullfile(myDataPath.output,'derivatives','age','AgeVsLatency_N1');

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)


%% overview of number of connections from one region to another

DestAmatall = zeros(size(n1Latencies,2),75,75);
for kk = 1:length(n1Latencies) % loop subjects    
    % initialize connections as empty
    DestAmat = zeros(75,75); % roi_start --> roi_end
    
    for ll = 1:length(n1Latencies(kk).run) % loop runs
        % run through all first stimulated channels (nr 1 in pair)
        for chPair = 1:size(n1Latencies(kk).run(ll).n1_peak_sample,2)
            
            % find region label of chPair
            chPairDest = str2double(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr(chPair,:));
            
            % find N1 responders
            chSig = find(~isnan(n1Latencies(kk).run(ll).n1_peak_sample(:,chPair))==1);
            % find region label of N1 responders
            chSigDest = str2double(n1Latencies(kk).run(ll).channel_DestrieuxNr(chSig));
            
            
            DestAmat(chSigDest(~isnan(chSigDest) & chSigDest ~= 0),chPairDest(~isnan(chPairDest) & chPairDest ~= 0)) = ...
                DestAmat(chSigDest(~isnan(chSigDest) & chSigDest ~= 0),chPairDest(~isnan(chPairDest) & chPairDest ~= 0))+1;
        end
    end
    
    DestAmatall(kk,:,:) = DestAmat; 
end

DestAmatsum = squeeze(sum(DestAmatall,1));

%% visualize all existing connections per patient and in total

subj = 1;

figure,
subplot(1,2,1),
imagesc(squeeze(DestAmatall(subj,:,:)))
title(sprintf('Connections between regions in subject %1.0f',subj))
xlabel('Responding Destrieux region')
ylabel('Stimulated Destrieux region')

subplot(1,2,2),
imagesc(DestAmatsum)
title('Connections between regions in all subjects')
xlabel('Responding Destrieux region')
ylabel('Stimulated Destrieux region')

%% 
