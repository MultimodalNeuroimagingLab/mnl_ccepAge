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

%% frontal to temporal connections

%%% need to wrap this in a function where we can enter 2 areas of interest

out = [];

for kk = 1:length(n1Latencies) % loop subjects
    out(1).sub(kk).age = n1Latencies(kk).age;
    out(2).sub(kk).age = n1Latencies(kk).age;
    
    % initialize connections as empty
    out(1).sub(kk).latencies = []; % temporal --> frontal
    out(2).sub(kk).latencies = []; % frontal --> temporal
    
    for ll = 1:length(n1Latencies(kk).run) % loop runs       
        % run through all first stimulated channels (nr 1 in pair)
        for ch = 1:size(n1Latencies(kk).run(ll).n1_peak_sample,2)
            
            % temporal --> frontal
            % if the stimulated channel 1 is temporal
            if ismember(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{ch,1},roi(1).areas)
                % run through all measured channels
                for chSig = 1:size(n1Latencies(kk).run(ll).n1_peak_sample,1)
                    thisChanDestrieuxNr = n1Latencies(kk).run(ll).channel_DestrieuxNr{chSig};
                    if ~isempty(thisChanDestrieuxNr)
                        % if a frontal channel is measured
                        if ismember(thisChanDestrieuxNr,roi(2).areas)
                            % we have a temporal stim, frontal measured pair
                            thisSample = n1Latencies(kk).run(ll).n1_peak_sample(chSig,ch);
                            
                            %%% this is still in samples, need to put this
                            %%% in time
                            out(1).sub(kk).latencies = [out(1).sub(kk).latencies thisSample];
                        end
                    end
                end
            end
            
            %%%% add here also if stimulated channel 2 is temporal!
            
            % frontal --> temporal
            % if the stimulated channel 1 is frontal
            if ismember(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr{ch,1},roi(2).areas)
                % run through all measured channels
                for chSig = 1:size(n1Latencies(kk).run(ll).n1_peak_sample,1)
                    thisChanDestrieuxNr = n1Latencies(kk).run(ll).channel_DestrieuxNr{chSig};
                    if ~isempty(thisChanDestrieuxNr)
                        % if a frontal channel is measured
                        if ismember(thisChanDestrieuxNr,roi(1).areas)
                            % we have a temporal stim, frontal measured pair
                            thisSample = n1Latencies(kk).run(ll).n1_peak_sample(chSig,ch);
                            
                            %%% this is still in samples, need to put this
                            %%% in time
                            out(2).sub(kk).latencies = [out(2).sub(kk).latencies thisSample];
                        end
                    end
                end
            end
            
            %%%% add here also if stimulated channel 2 is frontal!
        end
        
    end
end

%% figure

outInd = 2;

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
plot(my_output(:,1),my_output(:,2),'.')
xlabel('age (years)'),ylabel('mean latency (ms)')
[r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),2),'Type','Spearman');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

subplot(2,1,2),hold on
plot(my_output(:,1),my_output(:,3),'.')
xlabel('age (years)'),ylabel('variance in latency')
[r,p] = corr(my_output(~isnan(my_output(:,3)),1),my_output(~isnan(my_output(:,3)),3),'Type','Spearman');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])




%%

%                    allLatencies: [1×411 double]
%                  n1_peak_sample: [133×46 double]
%                   channel_names: {133×1 cell}
%              average_ccep_names: {46×1 cell}
%                   good_channels: [104×1 double]
%     average_ccep_DestrieuxLabel: {46×2 cell}
%        average_ccep_DestrieuxNr: {46×2 cell}
%          channel_DestrieuxLabel: {133×1 cell}
%             channel_DestrieuxNr: {133×1 cell}


