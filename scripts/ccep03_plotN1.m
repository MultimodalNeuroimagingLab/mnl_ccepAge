%% add: load the n1Latencies from the derivatives


%% connections from one region to another

%%% need to wrap this in a function where we can enter 2 areas of interest
region_start = input('Choose roi where connections start [temporal, frontal, parietal, occipital]: ','s');
region_end = input('Choose roi where connections end [temporal, frontal, parietal, occipital]: ','s');

ccep_connectRegions(n1Latencies,region_start,region_end)

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

