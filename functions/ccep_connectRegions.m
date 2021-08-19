function [out] = ccep_connectRegions(n1Latencies,roi_start, roi_end)
% determine for each region (start) to another region (end), the latency of
% CCEPs, the number of CCEPs and the relative number of CCEPs (corrected
% for total number of electrodes on region (end)). 

out = [];

for kk = 1:length(n1Latencies) % loop subjects
    out.sub(kk).age = n1Latencies(kk).age;
    
    out.sub(kk).el_roi_end = sum(ismember(str2double(n1Latencies(kk).run(1).channel_DestrieuxNr(:)),str2double(roi_end(:))));
    
    % initialize connections as empty
    out.sub(kk).samples = []; % roi_start --> roi_end
    out.sub(kk).latencies = []; % roi_start --> roi_end
    out.sub(kk).numCCEPs = [];
    out.sub(kk).relCCEPs = [];
    
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
                    
                    % number of CCEPs per stimulus                   
                    out.sub(kk).numCCEPs = [out.sub(kk).numCCEPs, size(thisSample,1)];
                    
                    % total number of electrodes on roi_end minus stimulated
                    % electrodes on roi_end, because in stimulated
                    % electrodes, no CCEP can be detected.
                    totCh_roi_end = sum(ismember(str2double(n1Latencies(kk).run(ll).channel_DestrieuxNr),str2double(roi_end(:)))) - ...
                        sum(ismember(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr(chPair,:),roi_end));
                    
                    % relative number of CCEPs per stimulus
                    out.sub(kk).relCCEPs = [out.sub(kk).relCCEPs, size(thisSample,1)/totCh_roi_end];
                end
            end
        end
    end    
end


end