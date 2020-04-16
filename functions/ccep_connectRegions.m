function ccep_connectRegions(n1Latencies,region_start,region_end)

if strcmp(region_start,'temporal')
    % temporal areas:
    % G_temporal_inf, G_temporal_middle
    roi_start = {'37','38'};
elseif strcmp(region_start,'frontal')
    % frontal areas:
    % G_front_inf-Triangul, G_front_middle, G_front_inf-Opercular
    roi_start = {'14';'15';'12'}; % maybe add 16: G_front_sup
elseif strcmp(region_start,'parietal')
    % parietal areas:
    % G_pariet_inf-Angular, G_pariet_inf-Supramar, G_parietal_sup
    roi_start = {'25','26','27'};
elseif strcmp(region_start,'occipital')
    % occipital areas:
    % G_occipital_middle G_oc-temp_med-Lingual Pole_occipital
    roi_start = {'19','22','42'};
end

if strcmp(region_end,'temporal')
    roi_end = {'37','38'};
    
elseif strcmp(region_end,'frontal')
    roi_end =  {'14';'15';'12'}; % maybe add 16: G_front_sup
    
elseif strcmp(region_end,'parietal')
    roi_end = {'25','26','27'};
    
elseif strcmp(region_end,'occipital')
    roi_end = {'19','22','42'};
end

out = [];

for kk = 1:length(n1Latencies) % loop subjects
    out.sub(kk).age = n1Latencies(kk).age;
    
    % initialize connections as empty
    out.sub(kk).samples = []; % roi_start --> roi_end
    
    out.sub(kk).latencies = []; % roi_start --> roi_end
    
    for ll = 1:length(n1Latencies(kk).run) % loop runs
        % run through all first stimulated channels (nr 1 in pair)
        for chPair = 1:size(n1Latencies(kk).run(ll).n1_peak_sample,2)
            
            % if stimulated channel 1 or 2 is on roi_start
            if ismember(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{chPair,1},roi_start) || ...
                    ismember(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{chPair,2},roi_start)
                
                % find N1 responders
                chSig = find(~isnan(n1Latencies(kk).run(ll).n1_peak_sample(:,chPair))==1);
                % find region label of N1 responders
                chDestrieux = str2double(n1Latencies(kk).run(ll).channel_DestrieuxNr(chSig));
                
                % if N1 responders are on roi_end
                if any(ismember(chDestrieux,str2double(roi_end)))
                    
                    thisSample = n1Latencies(kk).run(ll).n1_peak_sample(chSig(ismember(chDestrieux,str2double(roi_end))),chPair);
                    
                    out.sub(kk).samples = [out.sub(kk).samples; thisSample];
                    out.sub(kk).latencies = [out.sub(kk).latencies n1Latencies(kk).run(ll).tt(thisSample(~isnan(thisSample)))];
                end
            end
        end
    end
end

%% figure

outInd = 1;

% initialize output: age, mean and variance in latency per subject
my_output = NaN(length(out(outInd).sub),3);

% get variable per subject
for kk = 1:length(out(outInd).sub)
    my_output(kk,1) = out(outInd).sub(kk).age;
    my_output(kk,2) = nanmean(out(outInd).sub(kk).latencies);
    my_output(kk,3) = nanvar(out(outInd).sub(kk).latencies);
end

figure
subplot(2,1,1),hold on
plot(my_output(:,1),1000*my_output(:,2),'b.')
xlabel('age (years)'),ylabel('mean latency (ms)')
[r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),2),'Type','Pearson');
title([region_start ' to ' region_end   ', r=' num2str(r,3) ' p=' num2str(p,3)])

subplot(2,1,2),hold on
plot(my_output(:,1),my_output(:,3),'r.')
xlabel('age (years)'),ylabel('variance in latency')
[r,p] = corr(my_output(~isnan(my_output(:,3)),1),my_output(~isnan(my_output(:,3)),3),'Type','Pearson');
title([region_start ' to ' region_end   ', r=' num2str(r,3) ' p=' num2str(p,3)])

end