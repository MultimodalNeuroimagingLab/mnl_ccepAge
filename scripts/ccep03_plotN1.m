% temporal areas:
% G_temporal_inf, G_temporal_middle
roi(1).areas = {'37','38'};

% frontal areas:
% G_front_inf-Triangul, G_front_middle, G_front_inf-Opercular
roi(2).areas = {'14';'15';'12'}; % maybe add 16: G_front_sup

% parietal areas:
% G_pariet_inf-Angular, G_pariet_inf-Supramar, G_parietal_sup
roi(3).areas = {'25','26','27'};

% occipital areas:
% G_occipital_middle G_oc-temp_med-Lingual Pole_occipital
roi(4).areas = {'19','22','42'};

%% add: load the n1Latencies from the derivatives


%% frontal to temporal connections

%%% need to wrap this in a function where we can enter 2 areas of interest

roi_start = roi(1).areas;
roi_end = roi(3).areas;

out = [];

for kk = 1:length(n1Latencies) % loop subjects
    out.sub(kk).age = n1Latencies(kk).age;
    
    % initialize connections as empty
    out.sub(kk).samples = []; % roi_start --> roi_end
    out.sub(kk).latencies = []; % roi_start --> roi_end
    
    for ll = 1:length(n1Latencies(kk).run) % loop runs       
        % run through all first stimulated channels (nr 1 in pair)
        for chPair = 1:size(n1Latencies(kk).run(ll).n1_peak_sample,2)
            
            % temporal --> frontal
            % if stimulated channel 1 or 2 is on roi_start
            if ismember(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{chPair,1},roi_start) || ...
                    ismember(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{chPair,2},roi_start)
                % run through all measured channels
                for chSig = 1:size(n1Latencies(kk).run(ll).n1_peak_sample,1)
                    thisChanDestrieuxNr = n1Latencies(kk).run(ll).channel_DestrieuxNr{chSig};
                    if isnan(thisChanDestrieuxNr)
                        thisChanDestrieuxNr = '';
                    end
                    if ~isempty(thisChanDestrieuxNr) 
                        % if roi_end channel is measured
                        if ismember(thisChanDestrieuxNr,roi_end)
                            % we have a roi_start stim, roi_end measured pair
                            thisSample = n1Latencies(kk).run(ll).n1_peak_sample(chSig,chPair);
                            
                            %%% this is still in samples, need to put this
                            %%% in time
                            out.sub(kk).samples = [out.sub(kk).samples thisSample];
                            if ~isnan(thisSample)
                                out.sub(kk).latencies = [out.sub(kk).latencies n1Latencies(kk).run(ll).tt(thisSample)];
                            else
                                out.sub(kk).latencies = [out.sub(kk).latencies NaN];
                            end
                        end
                    end
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
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

subplot(2,1,2),hold on
plot(my_output(:,1),my_output(:,3),'r.')
xlabel('age (years)'),ylabel('variance in latency')
[r,p] = corr(my_output(~isnan(my_output(:,3)),1),my_output(~isnan(my_output(:,3)),3),'Type','Pearson');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])



