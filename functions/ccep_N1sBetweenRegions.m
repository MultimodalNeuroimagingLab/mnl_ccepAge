%   
%   Extract the latencies and (relative) number of N1s between one stimulated ROI and a response ROI
%   [out] = ccep_connectRegions(ccepData, roiStim, roiResp)
%
%       ccepData            = ...
%       roiStim             = An array of Destrieux codes that defines the ROI in which stimulation channels/electrodes will be included
%       roiResp             = An array of Destrieux codes that defines the ROI in which response channels/electrodes will be included
% 
%
%   Returns: 
%       out                 = A struct with N1 metrics for all subjects, each row represents one subject
%
%       out(x).samples      = A vector that holds all the sample-indices of N1s that occured in the response electrodes (found to be on the response ROI). 
%                             This vector concatenates values from all of the subject's runs and stim-pairs (of which one electrode is on the stimulus ROI)
%
%       out(x).latencies    = A vector that holds all the latencies of N1s that occured in the response electrodes (found to be on the response ROI). 
%                             This vector concatenates values from all of the subject's runs and stim-pairs (of which one electrode is on the stimulus ROI)
%
%       out(x).numN1s       = A vector that - for each stim-pair - holds the number of N1s that occur in the response electrodes (found to be on the response ROI). 
%                             This vector concatenates values from all of the subject's runs
%   
%       out(x).relNumN1s     = A vector that - for each stim-pair and within the set of response electrodes that are on the response ROI - holds the ratio
%                             between the number of N1s and response electrodes. This vector concatenates values from all of the subject's runs
%   

% divide the number of N1s for this stim-pair by the number of response channels that are on the response ROI
%
%   Dora Hermes, Max van den Boom, Dorien van Blooijs, 2022
%
function [out] = ccep_connectRegions(ccepData, roiStim, roiResp)
    out = struct;
    
    % loop over subjects
    for iSubj = 1:length(ccepData) 
        out(iSubj).age = ccepData(iSubj).age;
        
        % count the number of response electodes that are in the response ROI
        out(iSubj).numElecRespROI = sum(ismember(str2double(ccepData(iSubj).run(1).channel_DestrieuxNr), roiResp));

        % initialize connections as empty
        out(iSubj).samples      = [];
        out(iSubj).latencies    = [];
        out(iSubj).numN1s     = [];
        out(iSubj).relNumN1s     = [];

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

                    % check if any of the electrodes with a N1 response is on the response ROI
                    if any(ismember(respChanDestrieux, roiResp))
                        % at least one electrode of the current stim-pair is on the stim-ROI and
                        % at least one electrodes with a N1 response is on the response ROI
                        
                        % pick the N1 sample-indices of the channels that have an N1 response which are on the response ROI
                        n1SampleIndices = ccepData(iSubj).run(iRun).n1_peak_sample(respChans(ismember(respChanDestrieux, roiResp)), iStimpair);

                        % store the N1 sample-indices and latency (in ms) in the output structure
                        out(iSubj).samples = [out(iSubj).samples; n1SampleIndices];
                        out(iSubj).latencies = [out(iSubj).latencies ccepData(iSubj).run(iRun).tt(n1SampleIndices(~isnan(n1SampleIndices)))];

                        % store the number of N1s for this stim-pair
                        out(iSubj).numN1s = [out(iSubj).numN1s, size(n1SampleIndices, 1)];

                        
                        %
                        % relative number of N1s for this stim-pair
                        %
                        
                        % get number of electrodes in the response ROI minus the number of stimulated electrodes 
                        % on the response ROI (because in stimulated electrodes no N1 can be detected)
                        numRespCh = sum(ismember(str2double(ccepData(iSubj).run(iRun).channel_DestrieuxNr), roiResp)) - ...
                                        sum(ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr(iStimpair, :)), roiResp));
                        
                        % divide the number of N1s for this stim-pair by the number of response channels that are on the response ROI
                        out(iSubj).relNumN1s = [out(iSubj).relNumN1s, size(n1SampleIndices, 1) / numRespCh];
                        
                    end
                end
            end
            
        end     % end run loop   
    end     % end subject loop

end