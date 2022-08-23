%   
%   Extract the latencies and number of N1s/CCEPs between one stimulated ROI and a response ROI
%   [out] = ccep_connectRegions(ccepData, roiStim, roiResp)
%
%       ccepData        = ...
%       roiStim         = An array of Destrieux codes that defines the ROI in which stimulation channels/electrodes will be included
%       roiResp         = An array of Destrieux codes that defines the ROI in which response channels/electrodes will be included
% 
%
%   Returns: 
%       out             = ...
%   
%
%   Dora Hermes, Max van den Boom, Dorien van Blooijs, 2022
%
function [out] = ccep_connectRegions(ccepData, roiStim, roiResp)
    out = struct;
    
    % loop over subjects
    for iSubj = 1:length(ccepData) 
        out.sub(iSubj).age = ccepData(iSubj).age;
        
        % count the number of response electodes that are in the response ROI
        out.sub(iSubj).numElecRespROI = sum(ismember(str2double(ccepData(iSubj).run(1).channel_DestrieuxNr), roiResp));

        % initialize connections as empty
        out.sub(iSubj).samples      = [];
        out.sub(iSubj).latencies    = [];
        out.sub(iSubj).numCCEPs     = [];
        out.sub(iSubj).relCCEPs     = [];

        % loop over runs
        for iRun = 1:length(ccepData(iSubj).run)
            
            % loop over all stim-pairs
            for iStimpair = 1:size(ccepData(iSubj).run(iRun).n1_peak_sample, 2)
                
                % check if either of the electrodes in the stim-pair is on the stim ROI
                if ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr{iStimpair, 1}), roiStim) || ...
                   ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr{iStimpair, 2}), roiStim)
                    
                    % retrieve the channels that have a N1 response for this stim-pair
                    respChans = find(~isnan(ccepData(iSubj).run(iRun).n1_peak_sample(:, iStimpair)));
                    
                    % retrieve the region label for the channels with a N1 response
                    respChanDestrieux = str2double(ccepData(iSubj).run(iRun).channel_DestrieuxNr(respChans));

                    % check if any of the electrodes with a N1 response is on the reponse ROI
                    if any(ismember(respChanDestrieux, roiResp))
                        
                        % pick the N1 sample-indices of the channels that have an N1 response which are on the response ROI
                        n1SampleIndices = ccepData(iSubj).run(iRun).n1_peak_sample(respChans(ismember(respChanDestrieux, roiResp)), iStimpair);

                        % store the N1 sample-indices in the output structure, as index and and as latency (in ms)
                        out.sub(iSubj).samples = [out.sub(iSubj).samples; n1SampleIndices];
                        out.sub(iSubj).latencies = [out.sub(iSubj).latencies ccepData(iSubj).run(iRun).tt(n1SampleIndices(~isnan(n1SampleIndices)))];

                        % store the number of N1s/CCEPs per stimulus                   
                        out.sub(iSubj).numCCEPs = [out.sub(iSubj).numCCEPs, size(n1SampleIndices, 1)];

                        % TODO: check if metric is used...
                        % total number of electrodes on roi_end minus stimulated
                        % electrodes on roi_end, because in stimulated
                        % electrodes, no CCEP can be detected.
                        % <number of response channels that are on the response ROI> / <number of the current stim-pair channels are on the response ROI>
                        totCh_roi_end = sum(ismember(str2double(ccepData(iSubj).run(iRun).channel_DestrieuxNr), roiResp)) - ...
                                        sum(ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr(iStimpair, :)), roiResp));

                        % relative number of CCEPs per stimulus
                        out.sub(iSubj).relCCEPs = [out.sub(iSubj).relCCEPs, size(n1SampleIndices, 1) / totCh_roi_end];
                    end
                end
            end
            
        end     % end run loop   
    end     % end subject loop

end