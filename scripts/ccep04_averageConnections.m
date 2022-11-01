%
%  Scripts that averages the CCEP responses and N1 properties over the connections (from stimulated region to responding regions) for each run and each patient
% 
%  Note: the cells in the output variables ('subjectResponses', 'subjectResponses_nonnorm', 'subjectN1s', 'subjectN1sFWHM', 
%        'average_L_trk_length' and 'average_R_trk_length') represent the age (so cell 4 contains two values because it we have two subjects ages 4, etc).
%
%
%  Max van den Boom, Dorien van Blooijs, Dora Hermes. MultimodalNeuroimaging Lab (MNL), 2022
%

clear
close all
warning('on');
warning('backtrace', 'off')

myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

stimStimElec_excludeDist = 18;     % the distance between the stimulated electrodes (in mm) above which N1s are excluded, 0 = not excluding
respStimElec_excludeDist = 13;     % the distance between a stimulated and response electrode (in mm) within which N1s are excluded, 0 = not excluding


%%
%  

% load V2 of the data
if ~exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'file')
    error(['Could not find/load _V2 file (', fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), ')']);
end
load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData');

% load tracts and their corresponding end-point ROIs
rois = ccep_categorizeAnatomicalRegions();

% retrieve a time vector in ms (this patient has fs-2048)
tt = ccepData(2).run(1).tt;

% initialize variables to store ...
% (each cell represents a year of age, upto the highest age of all participants)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        subjectResponses{iTr}{iSubTr}           = cell(max([ccepData.age]), 1);
        subjectResponses_nonnorm{iTr}{iSubTr}   = cell(max([ccepData.age]), 1);
        subjectN1s{iTr}{iSubTr}                 = cell(max([ccepData.age]), 1);
        subjectN1FWHM{iTr}{iSubTr}              = cell(max([ccepData.age]), 1);
        average_L_trk_length{iTr}{iSubTr}       = cell(max([ccepData.age]), 1);
        average_R_trk_length{iTr}{iSubTr}       = cell(max([ccepData.age]), 1);
    end
end

% loop over the subjects
for iSubj = 1:size(ccepData, 2)

    %
    age = ccepData(iSubj).age;

    %
    disp('--------');
    disp(['  subj ', num2str(iSubj), ' of ', num2str(size(ccepData, 2))]);
    disp(['  id: ', ccepData(iSubj).id, ' - age: ', num2str(age), ' - ']);

    % create the cells and empty matrices to store CCEP traces and N1 quantifications in, and eventually contain the average
    for iTr = 1:length(rois)
        for iSubTr = 1:length(rois(iTr).sub_tract)
            runResponses{iTr}{iSubTr}           = NaN(2, size(ccepData(iSubj).run, 2), 5 * 2048);    % [roiDir, run, tt]
            runResponses_nonnorm{iTr}{iSubTr}   = NaN(2, size(ccepData(iSubj).run, 2), 5 * 2048);    % [roiDir, run, tt]
            runN1s{iTr}{iSubTr}                 = NaN(2, size(ccepData(iSubj).run, 2));              % [roiDir, run]
            runN1FWHM{iTr}{iSubTr}              = NaN(2, size(ccepData(iSubj).run, 2));              % [roiDir, run]
        end
    end

    % get the session directory name
    sesDir = dir(fullfile(myDataPath.output, 'derivatives', 'av_ccep', ccepData(iSubj).id, 'ses-*'));
    sesDir = sesDir.name;

    % loop over the runs
    for iRun = 1:size(ccepData(iSubj).run, 2)

        % load the run data
        runData = load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', ccepData(iSubj).id, sesDir, ccepData(iSubj).run(iRun).runName));

        % we will resample all cceps to 2048Hz, no need to preselect

        % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
        for iTr = 1:length(rois)
            for iSubTr = 1:length(rois(iTr).sub_tract)

                % use both directions (consider stim on roi1 and response on roi2, and vice versa)
                for iDir = [false true]
                    stimRoi = rois(iTr).sub_tract(iSubTr).(['roi', num2str(iDir + 1)]);
                    respRoi = rois(iTr).sub_tract(iSubTr).(['roi', num2str(~iDir + 1)]);

                    % find which stimulation pairs have at least one electrode on top of the stimulated ROI
                    stimPairs = find(sum(ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr), stimRoi), 2) > 0);

                    % find which response electrodes are on top of the response ROI
                    respChan = find(ismember(str2double(ccepData(iSubj).run(iRun).channel_DestrieuxNr), respRoi) > 0);


                    %
                    % collect all responses for the stimulation pair and response electrode combination that have a N1 response
                    % Note: that stim-pair and response channels that are bad will remain NaN in the matrix and skipped later
                    %

                    average_ccep_select     = NaN(size(stimPairs, 1), size(respChan, 1), size(runData.average_ccep, 3));       % stores evoked signal (in epoch)
                    average_N1_select       = NaN(size(stimPairs, 1), size(respChan, 1));                                      % stores the N1 peak latency
                    average_N1_FWHM_select  = NaN(size(stimPairs, 1), size(respChan, 1));                                      % stores the N1 width (FWHM in ms)

                    % loop over the stim-pairs
                    for iStimPair = 1:size(stimPairs, 1)

                        % check if native electrode distances are available
                        if isfield(ccepData(iSubj), 'nativeElecDistances')

                            % retrieve the stimulated electrodes and their respective indices in the electrodes table
                            stimPiarElecs = split(ccepData(iSubj).run(iRun).stimpair_names{stimPairs(iStimPair)}, '-');
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
                                warning(['Distance between two stimulated electrodes (', ccepData(iSubj).run(iRun).stimpair_names{stimPairs(iStimPair)}, ') is larger than ', num2str(stimStimElec_excludeDist), ' (', num2str(stim_stim_dist), ')']);

                                % skip this stim-pair
                                continue;

                            end

                        end

                        % loop over the response electrodes
                        for iRespChan = 1:size(respChan, 1)

                            % check if there was no N1 response between this stim-pair and the response electrode, goto next if no N1
                            % Note: stim-pair and response channels that are bad should have a NaN in the matrix and should be skipped here
                            if isnan(ccepData(iSubj).run(iRun).n1_peak_sample(respChan(iRespChan), stimPairs(iStimPair)))
                                continue;
                            end

                            % check if native electrode distances are available and an exclusion distance is set
                            if isfield(ccepData(iSubj), 'nativeElecDistances') && respStimElec_excludeDist ~= 0

                                % retrieve the response electrode and it's respective index in the electrodes table
                                respElec = ccepData(iSubj).run(iRun).channel_names{respChan(iRespChan)};
                                resp_elecIndex = find(ismember(upper(ccepData(iSubj).electrodes.name), upper(respElec)));
                                if isempty(resp_elecIndex)
                                    error(['Response electrode channel name (' , respElec, ') does not match any of the electrode names in ', ccepData(iSubj).id, ' - run ', num2str(iRun), ', trying case-insensitive.']); 
                                end

                                % retrieve the distances between the stimulated and response electrodes
                                resp_stim1_dist = ccepData(iSubj).nativeElecDistances(stim1_elecIndex, resp_elecIndex);
                                resp_stim2_dist = ccepData(iSubj).nativeElecDistances(stim2_elecIndex, resp_elecIndex);

                                % check whether either of the electrodes of the stimulus pair is within x mm of the response channel/electrode, skip if so
                                if respStimElec_excludeDist ~= 0 && resp_stim1_dist < respStimElec_excludeDist
                                    %warning(['Distance between stim1 electrode (', stimPiarElecs{1}, ') and response electrode (', respElec, ') is smaller than ', num2str(respStimElec_excludeDist), ' (', num2str(resp_stim1_dist), '), skipping']);
                                    continue;
                                end
                                if respStimElec_excludeDist ~= 0 && resp_stim2_dist < respStimElec_excludeDist
                                    %warning(['Distance between stim2 electrode (', stimPiarElecs{2}, ') and response electrode (', respElec, ') is smaller than ', num2str(respStimElec_excludeDist), ' (', num2str(resp_stim2_dist), '), skipping']);
                                    continue;
                                end

                            end

                            % extact and add response
                            response = squeeze(runData.average_ccep(respChan(iRespChan), stimPairs(iStimPair), :));
                            average_ccep_select(iStimPair, iRespChan, :) = response;

                            % extract and add N1 latency (in ms)
                            n1PeakSample = ccepData(iSubj).run(iRun).n1_peak_sample(respChan(iRespChan), stimPairs(iStimPair));
                            average_N1_select(iStimPair, iRespChan) = runData.tt(n1PeakSample);

                            % calculate and add N1 width
                            n1PeakAmpHalfWidth = .5 * response(n1PeakSample);
                            steps_back = find(response(n1PeakSample:-1:1) > n1PeakAmpHalfWidth, 1) - 1;     % walk back from min
                            steps_forw = find(response(n1PeakSample:1:end) > n1PeakAmpHalfWidth, 1) - 1;    % walk forward from min
                            if ~isempty(steps_back) && ~isempty(steps_forw)
                                % offset in response

                                if runData.tt(n1PeakSample - steps_back) > 0 && (runData.tt(n1PeakSample + steps_forw) - runData.tt(n1PeakSample - steps_back)) < .1 
                                    % can't onset <0, may occur if not at baseline at 9 ms
                                    % width must be smaller than 100ms, also offset issue

                                    average_N1_FWHM_select(iStimPair, iRespChan) = runData.tt(n1PeakSample + steps_forw) - runData.tt(n1PeakSample - steps_back);

                                end
                            end
                            clear response n1PeakSample n1PeakAmpHalfWidth steps_back steps_forw;

                        end
                        clear stimPiarElecs respElec resp_elecIndex stim1_elecIndex stim2_elecIndex;

                    end

                    % if responses are added to average_ccep_select
                    if any(~isnan(average_ccep_select(:)))

                        % average over both channels in the stim pair, and average over all response channels 
                        % is average response on a stim-pair that is one ROI, over all electrodes in the other ROI
                        runAverageResp      = squeeze(mean(mean(average_ccep_select,    1, 'omitnan'), 2, 'omitnan'));    
                        runAverageN1        = squeeze(mean(mean(average_N1_select,      1, 'omitnan'), 2, 'omitnan'));
                        runAverageN1FWHM    = squeeze(mean(mean(average_N1_FWHM_select, 1, 'omitnan'), 2, 'omitnan'));

                        % upsample if necessary
                        if length(runAverageResp) ~= 5 * 2048
                            runAverageResp = resample(runAverageResp, 5 * 2048, length(runAverageResp));
                        end

                        % store a non-normalized copy of the average response
                        runAverageRespNonNorm = runAverageResp; 

                        % 'normalize' 0.15-0.100 by unit length
                        runAverageResp = runAverageResp ./ sqrt(sum(runAverageResp(tt > .015 & tt < 0.100) .^ 2));

                        % add this run's response and N1
                        runResponses{iTr}{iSubTr}(iDir + 1, iRun, :)           = runAverageResp;           % [roiDir, run, samples]
                        runResponses_nonnorm{iTr}{iSubTr}(iDir + 1, iRun, :)   = runAverageRespNonNorm;    % [roiDir, run, samples]
                        runN1s{iTr}{iSubTr}(iDir + 1, iRun)                    = runAverageN1;             % [roiDir, run]
                        runN1FWHM{iTr}{iSubTr}(iDir + 1, iRun)                 = runAverageN1FWHM;         % [roiDir, run]

                        clear runAverageResp runAverageRespNonNorm runAverageN1 runAverageN1FWHM;
                    else

                        runResponses{iTr}{iSubTr}(iDir + 1, iRun, :)           = NaN(5 * 2048, 1);
                        runResponses_nonnorm{iTr}{iSubTr}(iDir + 1, iRun, :)   = NaN(5 * 2048, 1);
                        runN1s{iTr}{iSubTr}(iDir + 1, iRun)                    = NaN;
                        runN1FWHM{iTr}{iSubTr}(iDir + 1, iRun)                 = NaN;

                    end
                    clear average_ccep_select average_N1_select average_N1_FWHM_select;

                end

            end
        end

        clear runData averageN1;
    end


    %
    % Average over the runs and store in an 'age' specific cell
    %

    % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
    for iTr = 1:length(rois)
        for iSubTr = 1:length(rois(iTr).sub_tract)

            % mean response per subject (i.e. average over runs)
            subjectMeanResponse          = squeeze(mean(runResponses{iTr}{iSubTr}, 2, 'omitnan'));         %[roiDir, samples]
            subjectMeanResponse_nonnorm  = squeeze(mean(runResponses_nonnorm{iTr}{iSubTr}, 2, 'omitnan')); %[roiDir, samples]
            subjectMeanN1                = squeeze(mean(runN1s{iTr}{iSubTr}, 2, 'omitnan'));               %[roiDir]
            subjectMeanN1FWHM            = squeeze(mean(runN1FWHM{iTr}{iSubTr}, 2, 'omitnan'));            %[roiDir]

            % determine the index of the latest subject
            n = 0;
            if ~isempty(subjectResponses{iTr}{iSubTr}{age})
                n = size(subjectResponses{iTr}{iSubTr}{age}, 2);
            end

            % 
            subjectResponses{iTr}{iSubTr}{age}(:, n + 1, :)         = subjectMeanResponse;          % [roiDir, subjects, samples]
            subjectResponses_nonnorm{iTr}{iSubTr}{age}(:, n + 1, :) = subjectMeanResponse_nonnorm;  % [roiDir, subjects, samples]
            subjectN1s{iTr}{iSubTr}{age}(:, n + 1)                  = subjectMeanN1;                % [roiDir, subjects]
            subjectN1sFWHM{iTr}{iSubTr}{age}(:, n + 1)              = subjectMeanN1FWHM;            % [roiDir, subjects]

            % add the distance (if there are electrodes on that side
            %disp(['     tract: ', rois(iTr).tract_name, ' - ', rois(iTr).sub_tract(iSubTr).name]);
            if any(contains(ccepData(iSubj).electrodes.jsonHemi, 'L'))
                average_L_trk_length{iTr}{iSubTr}{age}(n + 1)           = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{1};
                %disp(['          left hemi - dist: ', num2str(ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{1})]);
            else
                average_L_trk_length{iTr}{iSubTr}{age}(n + 1)           = nan;
            end
            if any(contains(ccepData(iSubj).electrodes.jsonHemi, 'R'))
                average_R_trk_length{iTr}{iSubTr}{age}(n + 1)           = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{2};
                %disp(['          right hemi - dist: ', num2str(ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{2})]);
            else
                average_R_trk_length{iTr}{iSubTr}{age}(n + 1)           = nan;
            end

        end
    end
    clear runResponses runResponses_nonnorm runN1s
    clear subjectMeanResponse subjectMeanResponse_nonnorm subjectMeanN1 subjectMeanN1FWHM ROIsDist ROIsTrcs;
end

% save
s = input('Do you want to save the ccepAverages structure? [y/n]: ', 's');
if strcmp(s, 'y')
    save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepAverages.mat'), ...
         'subjectResponses', 'subjectResponses_nonnorm', 'subjectN1s', 'subjectN1sFWHM', 'tt', 'rois', 'average_L_trk_length', 'average_R_trk_length');
end

