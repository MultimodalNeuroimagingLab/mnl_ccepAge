% 
%  Load the ccepData from the derivatives
%  This code is used to plot the normalized ccep of all patients in order of age. 
%  This figure is displayed as Figure 2 in the article.
%


%% 
%  1. Load the ccepData from the derivatives

clear
close all

selectPat = input('Would you like to include all patients, or only the ones for whom it is certain that 8mA was applied (supplemental material)? [all/8] ','s');

myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData')
    filename_averageCCEP = fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepAverageResponses.mat');
else
    disp('Run scripts ccep02_loadN1.m and ccep03_addtracts.m first')
end



%% 
%  2. Load CCEP responses and categorize into connections from stimulated region to responding regions
%
%  the N1s are averaged for each run, and then averaged N1s per patient are collected for all subjects. 
%  Multiple subjects with the same age are collected for each age (subjectResponses_nonnorm) and normalized (subjectResponses)

if exist(filename_averageCCEP, 'file')
    load(filename_averageCCEP);
    
else
    
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
                        
                        % find stimulation pairs within specific region
                        % Note: will include a stimulated electrodes also if only one of the pair is inside of the ROI. TODO: discuss with Dora
                        stimPairs = find(sum(ismember(str2double(ccepData(iSubj).run(iRun).stimpair_DestrieuxNr), stimRoi), 2) > 0);
                        
                        % find response electrodes within specific region
                        respChan = find(ismember(str2double(ccepData(iSubj).run(iRun).channel_DestrieuxNr), respRoi) > 0);
                        

                        
                        %
                        % collect all signals with stimulation pair and response electrode within specific region
                        %
                        
                        average_ccep_select  = NaN(size(stimPairs, 1), size(respChan, 1), size(runData.average_ccep, 3));
                        average_N1_select    = NaN(size(stimPairs, 1), size(respChan, 1)); % N1latency
                        for cp = 1:size(stimPairs, 1)
                            for cs = 1:size(respChan, 1)
                                if ~isnan(ccepData(iSubj).run(iRun).n1_peak_sample(respChan(cs),stimPairs(cp)))
                                    if runData.tt(ccepData(iSubj).run(iRun).n1_peak_sample(respChan(cs), stimPairs(cp))) == 0, disp('what???'), end
                                    tempResp = squeeze(runData.average_ccep(respChan(cs),stimPairs(cp), :));
                                    average_ccep_select(cp, cs, :) = tempResp;
                                    average_N1_select(cp, cs) = runData.tt(ccepData(iSubj).run(iRun).n1_peak_sample(respChan(cs), stimPairs(cp)));
                                    clear tempResp
                                end
                            end
                        end

                        % if responses are added to average_ccep_select
                        if any(~isnan(average_ccep_select(:)))
                            
                            % average over both channels in the stim pair, and average over all response channels 
                            % is average response on a stim-pair that is one ROI, over all electrodes in the other ROI
                            runAverageResp = squeeze(mean(mean(average_ccep_select, 1, 'omitnan'), 2, 'omitnan'));    
                            runAverageN1   = squeeze(mean(mean(average_N1_select,   1, 'omitnan'), 2, 'omitnan'));

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

                            clear runAverageResp runAverageRespNonNorm runAverageN1;
                        else

                            runResponses{iTr}{iSubTr}(iDir + 1, iRun, :)           = NaN(5 * 2048, 1);
                            runResponses_nonnorm{iTr}{iSubTr}(iDir + 1, iRun, :)   = NaN(5 * 2048, 1);
                            runN1s{iTr}{iSubTr}(iDir + 1, iRun)                    = NaN;

                        end
                        clear average_ccep_select average_N1_select;
                    
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

                % determine the index of the latest subject
                n = 0;
                if ~isempty(subjectResponses{iTr}{iSubTr}{age})
                    n = size(subjectResponses{iTr}{iSubTr}{age}, 2);
                end
                
                % 
                subjectResponses{iTr}{iSubTr}{age}(:, n + 1, :)         = subjectMeanResponse;          % [roiDir, subjects, samples]
                subjectResponses_nonnorm{iTr}{iSubTr}{age}(:, n + 1, :) = subjectMeanResponse_nonnorm;  % [roiDir, subjects, samples]
                subjectN1s{iTr}{iSubTr}{age}(:, n + 1)                  = subjectMeanN1;                % [roiDir, subjects]
                
                % add the distance (if there are electrodes on that side
                disp(['     tract: ', rois(iTr).tract_name, ' - ', rois(iTr).sub_tract(iSubTr).name]);
                if any(contains(ccepData(iSubj).electrodes.jsonHemi, 'L'))
                    average_L_trk_length{iTr}{iSubTr}{age}(n + 1)           = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{1};
                    disp(['          left hemi - dist: ', num2str(ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{1})]);
                else
                    average_L_trk_length{iTr}{iSubTr}{age}(n + 1)           = nan;
                end
                if any(contains(ccepData(iSubj).electrodes.jsonHemi, 'R'))
                    average_R_trk_length{iTr}{iSubTr}{age}(n + 1)           = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{2};
                    disp(['          right hemi - dist: ', num2str(ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{2})]);
                else
                    average_R_trk_length{iTr}{iSubTr}{age}(n + 1)           = nan;
                end
                
            end
        end
        clear runResponses runResponses_nonnorm runN1s
        clear subjectMeanResponse subjectMeanResponse_nonnorm subjectMeanN1 ROIsDist ROIsTrcs;
    end

    save(filename_averageCCEP, 'subjectResponses', 'subjectResponses_nonnorm', 'subjectN1s', 'tt', 'rois', 'average_L_trk_length', 'average_R_trk_length');
end



%%
%  3. Sort CCEPs and their properties according to age (for each connection)
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
                    ageResponse = squeeze(subjectResponses{iTr}{iSubTr}{age}(iDir + 1, :, :));                 % [subjects, samples]
                    ageResponseNonnorm = squeeze(subjectResponses_nonnorm{iTr}{iSubTr}{age}(iDir + 1, :, :));  % [subjects, samples]
                    ageN1 = squeeze(subjectN1s{iTr}{iSubTr}{age}(iDir + 1, :));                                % [subjects]
                    ageLTrkLength = average_L_trk_length{iTr}{iSubTr}{age};                                    % [subjects]         (note, no direction in average_L_trk_length, since same length)
                    ageRTrkLength = average_R_trk_length{iTr}{iSubTr}{age};                                    % [subjects]         (note, no direction in average_R_trk_length, since same length)
                    
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
                    ageLTrkLength(excludeN1s) = [];
                    ageRTrkLength(excludeN1s) = [];
                    
                    if ~isempty(ageResponse) % there are subjects with electrodes on ROI
                        nr_subs = size(ageResponse, 1);
                        
                        numSubjectsWithROICoverage{iTr}{iSubTr}{iDir + 1} = [numSubjectsWithROICoverage{iTr}{iSubTr}{iDir + 1}, (zeros(1, nr_subs) + age)];

                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp; ageResponse];                               % [all subjects, samples]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm; ageResponseNonnorm];        % [all subjects, samples]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1 = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1, ageN1];                                         % [all subjects]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age, zeros(1, nr_subs) + age];                                   % [all subjects]
                        
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.L_trkLength = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.L_trkLength, ageLTrkLength];                            % [all subjects]
                        sortedCCEPs{iTr}{iSubTr}{iDir + 1}.R_trkLength = [sortedCCEPs{iTr}{iSubTr}{iDir + 1}.R_trkLength, ageRTrkLength];                            % [all subjects]
                        
                        clear addThis addThisNonnorm addN1 LTrkLength RTrkLength
                    end

                end
            end
        end
    end
end


%% 
%  Prepare data for plotting and calculate some statistics

% count the number of tests
numTests = 0;
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        numTests = numTests + 1;
    end
end
numTests = numTests * 2;    % test both directions

%
% Prepare summary matrices
% 
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        all_pBonf{iTr}{iSubTr} = zeros(1, 2);
        all_p{iTr}{iSubTr} = zeros(1, 2);
        all_r{iTr}{iSubTr} = zeros(1, 2);
        
        for iDir = [false true]
            
            n1Latencies = 1000 * sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1;
            age = sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age;
            
            % check if there are at least two subjects to correlate age and latency
            if numel(age) < 2
                all_pBonf{iTr}{iSubTr}(iDir + 1) = nan;
                all_p{iTr}{iSubTr}(iDir + 1) = nan;
                all_r{iTr}{iSubTr}(iDir + 1) = nan;
            else
                [r, p] = corr(n1Latencies', age', 'type', 'Spearman');
                all_pBonf{iTr}{iSubTr}(iDir + 1) = p / numTests;
                all_p{iTr}{iSubTr}(iDir + 1) = p;
                all_r{iTr}{iSubTr}(iDir + 1) = r;
                
            end
            
            
        end
    end
end



%%
%  Age descriptive plots

for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]
            
            % construct sub-tract string
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            strSubTitle = [subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];
            
            figure('Position',[0 0 600 300]);
            s = histogram(numSubjectsWithROICoverage{iTr}{iSubTr}{iDir + 1});
            title([rois(iTr).tract_name, ' - ', strSubTitle, '  (', num2str(length(numSubjectsWithROICoverage{iTr}{iSubTr}{iDir + 1})), ' subjects)']);

            if ~exist(fullfile(myDataPath.output,'derivatives', 'age'), 'dir')
                mkdir(fullfile(myDataPath.output,'derivatives', 'age'));
            end
            figureName = fullfile(myDataPath.output,'derivatives', 'age', ['AllSortAge_descr_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_')]);

            set(gcf,'PaperPositionMode','auto')
            print('-dpng','-r300',figureName)
            close(gcf)
            
        end
    end
end



%%
%  Plot age vs latencies, and age vs speed

% loop over the tracts
for iTr = 1:length(rois)

    %
    % each tract is one figure
    %
    figure('Position',[0 0 length(rois(iTr).sub_tract) * 600 800])

    % loop over the subtracts and direction
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]
            
            % construct sub-tract string
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            strSubTitle = [subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];
            
            
            %
            % each sub-tract/direction is one column (3 plots per column)
            %
            
            plotIndex = ((iSubTr - 1) * 2 + 1) + iDir;
            
            
            % age vs mean latency
            subplot(3, length(rois(iTr).sub_tract) * 2, plotIndex);
            x = sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age;
            y = 1000 * sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1;
            plot(x, y, '.')
            if iSubTr == 1 && iDir == 0
                ylabel('latency (ms)');
            end
            [r, p] = corr(x', y', 'Type', 'Pearson');
            title({rois(iTr).tract_name, strSubTitle, ' ', ['r=' num2str(r, 3) ' p=' num2str(p, 3)]})
            
            [P, S] = polyfit(x, y, 1);
            [y_fit, ~] = polyval(P, x, S);
            hold on
            plot(x, y_fit, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
            hold off
            
            
            
            % age vs tract length
            subplot(3, length(rois(iTr).sub_tract) * 2, length(rois(iTr).sub_tract) * 2 + plotIndex);
            
            x = sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age;
            y = mean([sortedCCEPs{iTr}{iSubTr}{iDir + 1}.L_trkLength; sortedCCEPs{iTr}{iSubTr}{iDir + 1}.R_trkLength], 'omitnan');
            x(isnan(y)) = [];
            y(isnan(y)) = [];
            
            plot(x, y, '.')
            if iSubTr == 1 && iDir == 0
                ylabel('tract length (mm)');
            end
            [r, p] = corr(x', y', 'Type', 'Pearson');
            title(['r=' num2str(r, 3) ' p=' num2str(p, 3)])
            
            [P, S] = polyfit(x, y, 1);
            [y_fit, ~] = polyval(P, x, S);
            hold on
            plot(x, y_fit, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
            hold off
            
            
            % age vs speed
            subplot(3, length(rois(iTr).sub_tract) * 2, 2 * length(rois(iTr).sub_tract) * 2 + plotIndex);
            
            x = sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age;
            y = mean([sortedCCEPs{iTr}{iSubTr}{iDir + 1}.L_trkLength, sortedCCEPs{iTr}{iSubTr}{iDir + 1}.R_trkLength], 'omitnan') ./ ...
                (1000 * sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1);
            x(isnan(y)) = [];
            y(isnan(y)) = [];
            
            plot(x, y, '.')
            xlabel('age (years)');
            if iSubTr == 1 && iDir == 0
                ylabel('speed (mm per ms)');
            end
            [r, p] = corr(x', y', 'Type', 'Pearson');
            title(['r=' num2str(r, 3) ' p=' num2str(p, 3)])
            
            [P, S] = polyfit(x, y, 1);
            [y_fit, ~] = polyval(P, x, S);
            hold on
            plot(x, y_fit, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
            hold off
            
        end
    end
    
    %
    % save
    %
    if ~exist(fullfile(myDataPath.output,'derivatives', 'age'), 'dir')
        mkdir(fullfile(myDataPath.output,'derivatives', 'age'));
    end
    figureName = fullfile(myDataPath.output,'derivatives', 'age', ['correlations_', rois(iTr).tract_name]);
    
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',figureName)
    print('-depsc','-r300',figureName)
    close(gcf)
    
end



%% 
%  figure of subplots for each stimulated and responding region with normalized CCEPs + N1 sorted by age


ttmin = 0.010;
ttmax = .100;

for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]
            
            % construct sub-tract string
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            strSubTitle = [subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];
            
            if ~isempty(sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age)
                
                %
                figure('Position',[0 0 600 300])
                hold on;

                imagesc(1000 * tt(tt > ttmin & tt < ttmax), ...
                        1:length(sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age), ...
                        -sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp(:, tt > ttmin & tt < ttmax), ...
                        [-0.1 0.1]);

                plot(1000 * sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1, 1:length(sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age), 'k.')
                colormap(parula)
                hold on
                set(gca,'XTick',20:20:80,'YTick',[])
                axis tight

                strSign = [' (pbonf = ', num2str(all_pBonf{iTr}{iSubTr}(iDir + 1)), ')'];
                if all_pBonf{iTr}{iSubTr}(iDir + 1) < .05,  strSign = [strSign, ' *'];     end
                title([rois(iTr).tract_name, ' - ', strSubTitle, strSign]);

                %
                % save
                %
                if ~exist(fullfile(myDataPath.output,'derivatives', 'age'), 'dir')
                    mkdir(fullfile(myDataPath.output,'derivatives', 'age'));
                end
                figureName = fullfile(myDataPath.output,'derivatives', 'age', ['sortedAge_tmax' int2str(ttmax * 1000), '_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_')]);

                set(gcf,'PaperPositionMode','auto')
                print('-dpng','-r300',figureName)
                print('-depsc','-r300',figureName)
                close(gcf)
                
            end
        end
    end
end


