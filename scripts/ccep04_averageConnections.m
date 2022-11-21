%
%  Scripts that averages the CCEP responses and N1 properties over the connections (from stimulated region to responding regions) for each run and each patient
% 
%  Note: the cells in the output variables ('subjectResponses', 'subjectResponses_nonnorm', 'subjectN1s', 'subjectN1Widths', 
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
assert(length(tt) == 10240, 'tt is expected to be 10240 (5*2048), all runs are upsampled to this length')

% initialize variables to store ...
% (each cell represents a year of age, upto the highest age of all participants)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)

        subjectResponses{iTr}{iSubTr}               = cell(max([ccepData.age]), 1);
        subjectResponses_nonnorm{iTr}{iSubTr}   = cell(max([ccepData.age]), 1);
        subjectN1s{iTr}{iSubTr}                 = cell(max([ccepData.age]), 1);
        subjectN1Widths{iTr}{iSubTr}              = cell(max([ccepData.age]), 1);
        average_L_trk_length{iTr}{iSubTr}       = cell(max([ccepData.age]), 1);
        average_R_trk_length{iTr}{iSubTr}       = cell(max([ccepData.age]), 1);
        
        % debug, variables to store the subject/run averages
        for iDir = [false true]
            allRunIds{iTr}{iSubTr}{iDir + 1}            = {};
            allRunElecConns{iTr}{iSubTr}{iDir + 1}      = {};
            allRunResponses{iTr}{iSubTr}{iDir + 1}      = [];
            allRunN1s{iTr}{iSubTr}{iDir + 1}      = [];
        end
        
    end
end

% loop over the subjects
for iSubj = 1:size(ccepData, 2)

    %
    age = ccepData(iSubj).age;

    %
    disp('--------');
    disp(['  subj ', num2str(iSubj), ' of ', num2str(size(ccepData, 2))]);
    disp(['  id: ', ccepData(iSubj).id, ' - age: ', num2str(age)]);

    % create the cells and empty matrices to store CCEP traces and N1 quantifications in, and eventually contain the average
    for iTr = 1:length(rois)
        for iSubTr = 1:length(rois(iTr).sub_tract)
            for iDir = [false true]
                allResponses{iTr}{iSubTr}{iDir + 1}    = [];
                allN1s{iTr}{iSubTr}{iDir + 1}          = [];
                allN1FWHM{iTr}{iSubTr}{iDir + 1}       = [];
            end
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
                    clear stimRoi respRoi;
                    

                    %
                    % collect all responses for the stimulation pair and response electrode combination that have a N1 response
                    % Note: that stim-pair and response channels that are bad will remain NaN in the matrix and skipped later
                    %

                    runResponses  = NaN(size(stimPairs, 1), size(respChan, 1), size(runData.average_ccep, 3));       % stores evoked signal (in epoch)
                    runN1s        = NaN(size(stimPairs, 1), size(respChan, 1));                                      % stores the N1 peak latency
                    runN1Widths  = NaN(size(stimPairs, 1), size(respChan, 1));                                       % stores the N1 width (FWHM in ms)
                    runElecConns  = cell(size(stimPairs, 1), size(respChan, 1));

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
                                %warning(['Distance between two stimulated electrodes (', ccepData(iSubj).run(iRun).stimpair_names{stimPairs(iStimPair)}, ') is larger than ', num2str(stimStimElec_excludeDist), ' (', num2str(stim_stim_dist), ')']);

                                % skip this stim-pair
                                continue;

                            end
                            clear stim_stim_dist;

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
                            runResponses(iStimPair, iRespChan, :) = response;

                            % extract and add N1 latency (in ms)
                            n1PeakSample = ccepData(iSubj).run(iRun).n1_peak_sample(respChan(iRespChan), stimPairs(iStimPair));
                            runN1s(iStimPair, iRespChan) = runData.tt(n1PeakSample);

                            % calculate and add N1 width
                            n1PeakAmpHalfWidth = .5 * response(n1PeakSample);
                            steps_back = find(response(n1PeakSample:-1:1) > n1PeakAmpHalfWidth, 1) - 1;     % walk back from min
                            steps_forw = find(response(n1PeakSample:1:end) > n1PeakAmpHalfWidth, 1) - 1;    % walk forward from min
                            if ~isempty(steps_back) && ~isempty(steps_forw)
                                % offset in response

                                if runData.tt(n1PeakSample - steps_back) > 0 && (runData.tt(n1PeakSample + steps_forw) - runData.tt(n1PeakSample - steps_back)) < .1 
                                    % can't onset <0, may occur if not at baseline at 9 ms
                                    % width must be smaller than 100ms, also offset issue

                                    runN1Widths(iStimPair, iRespChan) = runData.tt(n1PeakSample + steps_forw) - runData.tt(n1PeakSample - steps_back);

                                end
                            end
                            
                            runElecConns{iStimPair, iRespChan} = [stimPiarElecs{1} '-' stimPiarElecs{2}, ' --> ', respElec];
                            
                            
                            clear response n1PeakSample n1PeakAmpHalfWidth steps_back steps_forw resp_stim1_dist resp_stim2_dist;
                        end
                        
                        clear stimPiarElecs respElec resp_elecIndex stim1_elecIndex stim2_elecIndex;
                    end

                    % check if any N1 responses were found (added to average_ccep_select)
                    assert(any(~isnan(runResponses(:))) == any(~isnan(runN1s(:))), 'different');
                    if any(~isnan(runN1s(:)))
                        
                        % reshape all average CCEPs to a 2D matrix <all stim/pairs> x <time>
                        runResponses = reshape(runResponses, [], size(runResponses, 3));                        
                        runResponses(all(isnan(runResponses), 2), :) = [];
                        
                        % upsample if necessary
                        if size(runResponses, 2) ~= 5 * 2048
                            runResponses = resample(runResponses', 5 * 2048, size(runResponses, 2))';
                        end
                        
                        % retrieve the N1s and N1 widths
                        runN1s = runN1s(~isnan(runN1s));
                        runN1s = runN1s(:);
                        runN1Widths = runN1Widths(~isnan(runN1Widths));
                        runN1Widths = runN1Widths(:);

                        % collect/store all CCEP (properties)
                        if isempty(allN1s{iTr}{iSubTr}{iDir + 1})
                            allResponses{iTr}{iSubTr}{iDir + 1} = runResponses;
                            allN1s{iTr}{iSubTr}{iDir + 1} = runN1s;
                            allN1FWHM{iTr}{iSubTr}{iDir + 1} = runN1Widths;
                        else
                            allResponses{iTr}{iSubTr}{iDir + 1} = cat(1, allResponses{iTr}{iSubTr}{iDir + 1}, runResponses);
                            allN1s{iTr}{iSubTr}{iDir + 1} = [allN1s{iTr}{iSubTr}{iDir + 1}; runN1s];
                            allN1FWHM{iTr}{iSubTr}{iDir + 1} = [allN1s{iTr}{iSubTr}{iDir + 1}; runN1Widths];
                        end
                        
                        
                        %
                        % debug, average the run and store the subject/run
                        %
                        
                        % calculate run average
                        averageRunResponse = mean(runResponses, 1, 'omitnan');
                        if isempty(averageRunResponse)
                            averageRunResponse = NaN(5 * 2048, 1);
                        end
                        averageRunResponse = averageRunResponse ./ sqrt(sum(averageRunResponse(tt > .015 & tt < 0.100) .^ 2));
            
                        runId = ['(', num2str(ccepData(iSubj).age), ') ', ...
                                 strrep(strrep(strrep(strrep(ccepData(iSubj).run(iRun).runName, ...
                                                        '_ses-1_task-SPESclin_run-', '_run'), ...
                                                        '_ses-1b_task-SPESclin_run-', '_run'), ...
                                                        '_averageCCEPs.mat', ''), ...
                                                        'sub-ccepAgeUMCU', 's')];
                        
                        % store subject/run average
                        allRunIds{iTr}{iSubTr}{iDir + 1}            = [allRunIds{iTr}{iSubTr}{iDir + 1}; strrep(runId, '_', '\_')];
                        cellConn = runElecConns(~cellfun(@isempty, runElecConns));
                        allRunElecConns{iTr}{iSubTr}{iDir + 1}      = [allRunElecConns{iTr}{iSubTr}{iDir + 1}; strjoin(cellConn, '   ')];
                        allRunResponses{iTr}{iSubTr}{iDir + 1}      = [allRunResponses{iTr}{iSubTr}{iDir + 1}; averageRunResponse];
                        allRunN1s{iTr}{iSubTr}{iDir + 1}            = [allRunN1s{iTr}{iSubTr}{iDir + 1}; mean(runN1s)];

                        
                        %
                        % print
                        %
                        %{
                        figure('Position', [0 0 1600 800]);     hold on;
                        for iConn = 1:size(runResponses, 1)
                            this_ccep_plot = squeeze(runResponses(iConn, :));
                            this_ccep_plot(tt > -0.010 & tt <= 0.009) = NaN;
                            plot(tt, iConn * 500 + zeros(size(tt)), 'Color', [.8 .8 .8]);
                            plot(tt, iConn * 500 + this_ccep_plot);
                            plot(runN1s(iConn), iConn * 500, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 3)
                        end
                        plot(tt, (iConn + 1) * 500 + averageRunResponse * 1000, 'Color', [0 0 0], 'LineWidth', 2);
                        plot(mean(runN1s), (iConn + 1) * 500, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 3)
                        xlim([-1 1])
                        ylim([-500, (iConn + 2) * 500])
                        set(gca, 'YTick', 500 * (1:size(runResponses, 1)), 'YTickLabel', cellConn)
                        ylabel('stimulated electrodes')
                        xlabel('time(s)')
                        
                        subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
                        strTitle = [rois(iTr).tract_name, ' - ', subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];
                        title(strrep([strTitle, '   ', runId], '_', '\_'));

                        % save
                        if ~exist(fullfile(myDataPath.output,'derivatives', 'inspect'), 'dir')
                            mkdir(fullfile(myDataPath.output,'derivatives', 'inspect'));
                        end
                        figureName = fullfile(myDataPath.output,'derivatives', 'inspect', ['averConnSignals_', runId, '_', strrep(strrep(strTitle, ' -> ', '_'), ' ', '')]);
                        set(gcf, 'PaperPositionMode', 'auto')
                        print('-dpng', '-r300', figureName)
                        close(gcf)
                        %}
                        
                    end
                    
                    clear runResponses runN1s runN1Widths;
                end

            end
        end

        clear runData;
    end

    
    %
    % Average over the runs and store in an 'age' specific cell
    %

    % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
    for iTr = 1:length(rois)
        for iSubTr = 1:length(rois(iTr).sub_tract)
            
            % determine the index of the latest subject
            n = 0;
            if ~isempty(subjectResponses{iTr}{iSubTr}{age})
                n = size(subjectResponses{iTr}{iSubTr}{age}, 2);
            end
            
            % 
            for iDir = [false true]
                
                % average over all the runs, stim-pairs and response electrodes
                averageResp = mean(allResponses{iTr}{iSubTr}{iDir + 1}, 1, 'omitnan');
                if isempty(averageResp)
                    averageResp = NaN(5 * 2048, 1);
                end
                
                % store a non-normalized copy of the average responses
                averageRespNonNorm = averageResp; 
                
                % 'normalize' 0.15-0.100 by unit length
                averageResp = averageResp ./ sqrt(sum(averageResp(tt > .015 & tt < 0.100) .^ 2));

                % 
                subjectResponses{iTr}{iSubTr}{age}(iDir + 1, n + 1, :)          = averageResp;                              % [roiDir, subjects, samples]
                subjectResponses_nonnorm{iTr}{iSubTr}{age}(iDir + 1, n + 1, :)  = averageRespNonNorm;                       % [roiDir, subjects, samples]
                subjectN1s{iTr}{iSubTr}{age}(iDir + 1, n + 1)                   = mean(allN1s{iTr}{iSubTr}{iDir + 1});      % [roiDir, subjects]
                subjectN1Widths{iTr}{iSubTr}{age}(iDir + 1, n + 1)              = mean(allN1FWHM{iTr}{iSubTr}{iDir + 1});   % [roiDir, subjects]
                
            end
            
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
    clear allResponses allResponses_nonnorm allN1s
    clear subjectMeanResponse subjectMeanResponse_nonnorm subjectMeanN1 subjectMeanN1FWHM ROIsDist ROIsTrcs;
end

%{
%%
% Debug, visualize subject/runs to inspect for artifacts

ttmin = -0.3;
ttmax = .500;

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            strSubTitle = [subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];

            %
            figure('Position', [0 0 2400 1200])
            hold on;

            imagesc(1000 * tt(tt > ttmin & tt < ttmax), ...
                    1:length(allRunIds{iTr}{iSubTr}{iDir + 1}), ...
                    -allRunResponses{iTr}{iSubTr}{iDir + 1}(:, tt > ttmin & tt < ttmax), ...
                    [-0.1 0.1]);
            colormap(parula)
            
            axis tight
            set(gca, 'YTick', 1:length(allRunIds{iTr}{iSubTr}{iDir + 1}), 'YTickLabel', allRunIds{iTr}{iSubTr}{iDir + 1})
            for iElecConn = 1:length(allRunIds{iTr}{iSubTr}{iDir + 1})
               text(ttmin * 1000, iElecConn, allRunElecConns{iTr}{iSubTr}{iDir + 1}{iElecConn});
               plot(allRunN1s{iTr}{iSubTr}{iDir + 1}(iElecConn) * 1000, iElecConn, 'or');
            end
            title(strrep([rois(iTr).tract_name, ' - ', strSubTitle], '_', '\_'));

            %
            % save
            %
            if ~exist(fullfile(myDataPath.output, 'derivatives', 'inspect'), 'dir')
                mkdir(fullfile(myDataPath.output, 'derivatives', 'inspect'));
            end
            figureName = fullfile(myDataPath.output,'derivatives', 'inspect', ['averConn_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_')]);

            set(gcf, 'PaperPositionMode', 'auto')
            print('-dpng', '-r300', figureName)
            close(gcf)
            
        end
    end
end
%}



%%
%  Sort CCEPs and their properties according to age (for each connection)
%
%  Note: this does not average anything over multiple subjects at the same age.
%

sortedCCEPs = {};
numSubjectsWithROICoverage = {};
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]
            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp = [];
            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm = [];
            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age = [];
            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1 = [];
            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1FWHM = [];
            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.L_trkLength = [];
            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.R_trkLength = [];
            numSubjectsWithROICoverage{iTr}{iSubTr}{iDir + 1} = [];
        end
    end
end

for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        for age = 1:max([ccepData.age])
            if ~isempty(subjectResponses{iTr}{iSubTr}{age})
                
                % for both directions (consider stim on roi1 and response on roi2, and vice versa)
                for iDir = [false true]
                    
                    % 
                    ageResponse = squeeze(subjectResponses{iTr}{iSubTr}{age}(iDir + 1, :, :));                  % [subjects, samples]
                    ageResponseNonnorm = squeeze(subjectResponses_nonnorm{iTr}{iSubTr}{age}(iDir + 1, :, :));   % [subjects, samples]
                    ageN1 = squeeze(subjectN1s{iTr}{iSubTr}{age}(iDir + 1, :));                                 % [subjects]
                    ageN1Width = squeeze(subjectN1Widths{iTr}{iSubTr}{age}(iDir + 1, :));                       % [subjects]
                    ageLTrkLength = average_L_trk_length{iTr}{iSubTr}{age};                                     % [subjects]         (note, no direction in average_L_trk_length, since same length)
                    ageRTrkLength = average_R_trk_length{iTr}{iSubTr}{age};                                     % [subjects]         (note, no direction in average_R_trk_length, since same length)
                    
                    % if squeeze changed the orientation of the matrix, undo
                    if size(ageResponse, 1) > size(ageResponse, 2)
                        ageResponse = ageResponse';
                        ageResponseNonnorm = ageResponseNonnorm';
                    end
                    
                    % determine if there are nans (meaning no N1s - for subjects - and between a specific tract, sub-tract and direction)
                    excludeN1s = find(isnan(ageResponse(:, 1)));
                    ageResponse(excludeN1s, :) = [];
                    ageResponseNonnorm(excludeN1s, :) = [];
                    ageN1(excludeN1s) = [];
                    ageN1Width(excludeN1s) = [];
                    ageLTrkLength(excludeN1s) = [];
                    ageRTrkLength(excludeN1s) = [];
                    
                    if ~isempty(ageResponse) % there are subjects with electrodes on ROI
                        nr_subs = size(ageResponse, 1);
                        
                        numSubjectsWithROICoverage{iTr}{iSubTr}{iDir + 1} = [numSubjectsWithROICoverage{iTr}{iSubTr}{iDir + 1}, (zeros(1, nr_subs) + age)];

                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp; ageResponse];                             % [all subjects, samples]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm; ageResponseNonnorm];      % [all subjects, samples]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1 = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1, ageN1];                                       % [all subjects]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1FWHM = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1FWHM, ageN1Width];                          % [all subjects]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age, zeros(1, nr_subs) + age];                                 % [all subjects]
                        
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.L_trkLength = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.L_trkLength, ageLTrkLength];                           % [all subjects]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.R_trkLength = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.R_trkLength, ageRTrkLength];                           % [all subjects]
                        
                        clear addThis addThisNonnorm addN1 ageN1Width LTrkLength RTrkLength
                    end

                end
            end
        end
    end
end

%%
% save
s = input('Do you want to save the ccepAverages structure? [y/n]: ', 's');
if strcmp(s, 'y')
    save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepAverages.mat'), ...
         'subjectResponses', 'subjectResponses_nonnorm', 'subjectN1s', 'subjectN1Widths', 'tt', 'rois', 'average_L_trk_length', 'average_R_trk_length', 'sortedCCEPs');
end

