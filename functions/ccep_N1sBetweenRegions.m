%   
%   Extract the latencies and (relative) number of N1s between one stimulated ROI and a response ROI
%   [out] = ccep_N1sBetweenRegions(ccepData, roiStim, roiResp)
%
%       ccepData                    = ...
%       roiStim                     = An array of Destrieux codes that defines the ROI in which stimulation channels/electrodes will be included
%       roiResp                     = An array of Destrieux codes that defines the ROI in which response channels/electrodes will be included
%       stimStimElec_excludeDist    = The distance between the stimulated electrodes (in mm) above which N1s are excluded, 0 = not excluding
%       respStimElec_excludeDist    = The distance between a stimulated and response electrode (in mm) within which N1s are excluded, 0 = not excluding
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
%       out(x).ratioN1s     = A vector that - for each stim-pair and within the set of response electrodes that are on the response ROI - holds the ratio
%                             between the number of N1s and response electrodes. This vector concatenates values from all of the subject's runs
%  
%
%   Dora Hermes, Max van den Boom, Dorien van Blooijs, 2022
%
function [out] = ccep_N1sBetweenRegions(ccepData, roiStim, roiResp, stimStimElec_excludeDist, respStimElec_excludeDist)
    out = struct;
    
    % loop over subjects
    for iSubj = 1:length(ccepData) 
        out(iSubj).age = ccepData(iSubj).age;
        
        % count the number of response electodes that are in the response ROI
        out(iSubj).numElecRespROI = sum(ismember(str2double(ccepData(iSubj).run(1).channel_DestrieuxNr), roiResp));

        % initialize connections as empty
        out(iSubj).samples      = [];
        out(iSubj).latencies    = [];
        out(iSubj).numN1s       = [];
        out(iSubj).ratioN1s     = [];
        out(iSubj).distRespStim = [];
        out(iSubj).StimPairNr   = [];

        % loop over runs
        for iRun = 1:length(ccepData(iSubj).run)
            
            % loop over all stim-pairs
            for iStimPair = 1:size(ccepData(iSubj).run(iRun).n1_peak_sample, 2)

                % check if native electrode distances are available
                if isfield(ccepData(iSubj), 'nativeElecDistances')

                    % retrieve the stimulated electrodes and their respective indices in the electrodes table
                    stimPiarElecs = split(ccepData(iSubj).run(iRun).stimpair_names{iStimPair}, '-');
                    stim1_elecIndex = find(ismember(upper(ccepData(iSubj).electrodes.name), upper(stimPiarElecs{1})));
                    stim2_elecIndex = find(ismember(upper(ccepData(iSubj).electrodes.name), upper(stimPiarElecs{2})));
                    if isempty(stim1_elecIndex)
                        error(['Stim electrode 1 name (' , stimPiarElecs{1}, ') does not match any of the electrode names in ', ccepData(iSubj).id, ' - run ', num2str(iRun), ', trying case-insensitive.']); 
                    end
                    if isempty(stim2_elecIndex)
                        error(['Stim electrode 2 name (' , stimPiarElecs{2}, ') does not match any of the electrode names in ', ccepData(iSubj).id, ' - run ', num2str(iRun), ', trying case-insensitive.']); 
                    end

                    % retrieve the distance between stimulated electrodes
                    stim_stim_dist = ccepData(iSubj).nativeElecDistances(stim1_elecIndex, stim2_elecIndex);

                    % check if the distance between the stimulated electrodes is larger than exclusion threshold
                    if stimStimElec_excludeDist ~= 0 && stim_stim_dist > stimStimElec_excludeDist
                        %warning(['Distance between two stimulated electrodes (', ccepData(iSubj).run(iRun).stimpair_names{iStimPair}, ') is larger than ', num2str(stimStimElec_excludeDist), ' (', num2str(stim_stim_dist), ')']);

                        % skip this stim-pair
                        continue;

                    end

                end
                
                % check if either of the electrodes in the stim-pair is on the stim ROI
                if ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr{iStimPair, 1}), roiStim) || ...
                   ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr{iStimPair, 2}), roiStim)
                    
                    % retrieve the channels that have a N1 response for this stim-pair
                    respChans = find(~isnan(ccepData(iSubj).run(iRun).n1_peak_sample(:, iStimPair)));
                    
                    % retrieve the region label for the channels with a N1 response
                    respChanDestrieux = str2double(ccepData(iSubj).run(iRun).channel_DestrieuxNr(respChans));

                    % check if native electrode distances are available and an exclusion distance is set
                    if isfield(ccepData(iSubj), 'nativeElecDistances') && respStimElec_excludeDist ~= 0
                        
                        % find response electrodes that are to close to the stimulated electrodes
                        exclChans = [];
                        all_respStimDist = [];

                        for iRespChan = 1:length(respChans)

                            % retrieve the response electrode and it's respective index in the electrodes table
                            respElec = ccepData(iSubj).run(iRun).channel_names{respChans(iRespChan)};
                            resp_elecIndex = find(ismember(upper(ccepData(iSubj).electrodes.name), upper(respElec)));
                            if isempty(resp_elecIndex)
                                error(['Response electrode channel name (' , respElec, ') does not match any of the electrode names in ', ccepData(iSubj).id, ' - run ', num2str(iRun), ', trying case-insensitive.']); 
                            end

                            % retrieve the distances between the stimulated and response electrodes
                            resp_stim1_dist = ccepData(iSubj).nativeElecDistances(stim1_elecIndex, resp_elecIndex);
                            resp_stim2_dist = ccepData(iSubj).nativeElecDistances(stim2_elecIndex, resp_elecIndex);
                            all_respStimDist(end+1,:) = [resp_stim1_dist resp_stim2_dist];

                            % check whether either of the electrodes of the stimulus pair is within x mm of the response channel/electrode, skip if so
                            if respStimElec_excludeDist ~= 0 && resp_stim1_dist < respStimElec_excludeDist
                                %warning(['Distance between stim1 electrode (', stimPiarElecs{1}, ') and response electrode (', respElec, ') is smaller than ', num2str(respStimElec_excludeDist), ' (', num2str(resp_stim1_dist), '), skipping']);
                                exclChans(end + 1) = iRespChan;
                            end
                            if respStimElec_excludeDist ~= 0 && resp_stim2_dist < respStimElec_excludeDist
                                %warning(['Distance between stim2 electrode (', stimPiarElecs{2}, ') and response electrode (', respElec, ') is smaller than ', num2str(respStimElec_excludeDist), ' (', num2str(resp_stim2_dist), '), skipping']);
                                exclChans(end + 1) = iRespChan;
                            end
                        end
                        
                        % exclude response electrodes
                        respChans(exclChans) = [];
                        respChanDestrieux(exclChans) = [];
                        all_respStimDist(exclChans,:) = [];
                        
                        clear exclChans respElec resp_elecIndex iRespChan respElec resp_stim1_dist resp_stim2_dist;
                        
                    end

                    % check if any of the electrodes with a N1 response is on the response ROI
                    if any(ismember(respChanDestrieux, roiResp))
                        % at least one electrode of the current stim-pair is on the stim-ROI and
                        % at least one electrodes with a N1 response is on the response ROI
                        
                        % pick the N1 sample-indices of the channels that have an N1 response which are on the response ROI
                        n1SampleIndices = ccepData(iSubj).run(iRun).n1_peak_sample(respChans(ismember(respChanDestrieux, roiResp)), iStimPair);

                        % store the N1 sample-indices and latency (in ms) in the output structure
                        out(iSubj).samples = [out(iSubj).samples; n1SampleIndices];
                        out(iSubj).latencies = [out(iSubj).latencies ccepData(iSubj).run(iRun).tt(n1SampleIndices(~isnan(n1SampleIndices)))];

                        % store the number of N1s for this stim-pair
                        out(iSubj).numN1s = [out(iSubj).numN1s, size(n1SampleIndices, 1)];

                        % store the distances
                        out(iSubj).distRespStim = [out(iSubj).distRespStim; all_respStimDist(ismember(respChanDestrieux, roiResp),:)];
                                    
                        % store the stim pair nr, add 100 to the second
                        % run, because there are not >100 stim pairs
                        out(iSubj).StimPairNr = [out(iSubj).StimPairNr; (iRun-1)*100+iStimPair*ones(size(n1SampleIndices))];

                        %
                        % relative number of N1s for this stim-pair
                        %
                        
                        % get number of electrodes in the response ROI, minus the number of stimulated electrodes 
                        % on the response ROI (because in stimulated electrodes no N1 can be detected)
                        numRespCh = sum(ismember(str2double(ccepData(iSubj).run(iRun).channel_DestrieuxNr), roiResp)) - ...
                                        sum(ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr(iStimPair, :)), roiResp));
                        
                        % divide the number of N1s for this stim-pair by the number of response channels that are on the response ROI
                        out(iSubj).ratioN1s = [out(iSubj).ratioN1s, size(n1SampleIndices, 1) / numRespCh];
                        
                    end
                end
            end
            
        end     % end run loop   
    end     % end subject loop

end