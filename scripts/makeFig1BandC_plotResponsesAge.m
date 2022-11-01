% 
% 
%  This code is used to plot Fig 1C - the normalized ccep of all patients in order of age, 
%  and Fig 1B - example CCEPs for a few older and younger subjects 
%


%% 
%  Load the ccepData and ccepAverages from the derivatives

clear
close all
warning('on');
warning('backtrace', 'off')

myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData')
else
    disp('Run scripts ccep02_aggregateToStruct.m and ccep03_addtracts.m first')
end

if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepAverages.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepAverages.mat'))
else
    disp('Run script ccep04_averageConnections.m first')
end



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

% Prepare summary matrices
all_Ps = [];
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        all_r{iTr}{iSubTr} = zeros(1, 2);
        all_p{iTr}{iSubTr} = zeros(1, 2);
        all_p_fdr{iTr}{iSubTr} = zeros(1, 2);
        
        for iDir = [false true]
            
            n1Latencies = 1000 * sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageN1;
            age = sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age;
            
            % check if there are at least two subjects to correlate age and latency
            if numel(age) < 2
                all_r{iTr}{iSubTr}(iDir + 1) = nan;
                all_p{iTr}{iSubTr}(iDir + 1) = nan;
            else
                [r, p] = corr(n1Latencies', age', 'type', 'Spearman');
                
                all_r{iTr}{iSubTr}(iDir + 1) = r;
                all_p{iTr}{iSubTr}(iDir + 1) = p;
                
            end
            all_p_fdr{iTr}{iSubTr}(iDir + 1) = nan;
            
            % store all P values for FDR correction later
            all_Ps(end + 1, :) = [iTr, iSubTr, iDir + 1, all_p{iTr}{iSubTr}(iDir + 1)];
            
        end
    end
end

% FDR corrected
[~, ~, ~, all_Ps(:, 5)]  = fdr_bh(all_Ps(:, 4), 0.05, 'pdep');
for iP = 1:size(all_Ps, 1)
    all_p_fdr{all_Ps(iP, 1)}{all_Ps(iP, 2)}(all_Ps(iP, 3)) = all_Ps(iP, 5);

end



%% 
% Figure 1C - Plots for each stimulated and responding region with normalized CCEPs + N1 sorted by age

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
                figure('Position', [0 0 600 300])
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

                strSign = [' (p\_fdr = ', num2str(all_p_fdr{iTr}{iSubTr}(iDir + 1)), ')'];
                if all_p_fdr{iTr}{iSubTr}(iDir + 1) < .05,  strSign = [strSign, ' *'];     end
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



%%
% Figure 1B - Plot example CCEPs for a few older and younger subjects 
%

ttmin = -0.300;
ttmax = .400;

for iTr = 3                 % SLF
    for iSubTr = 1          % parietal->frontal/frontal->parietal subtract
        for iDir = [true]   % parietal to frontal
            
            % construct sub-tract string
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            strSubTitle = [subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];
            
            if ~isempty(sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age)
                
                %
                figure('Position', [0 0 500 150])
                plot(1000 * tt(tt > ttmin & tt < ttmax), ...
                        zeros(size(tt(tt > ttmin & tt < ttmax))),'Color',[.5 .5 .5]);
                hold on
                
                % plot traces for a couple of young subjects
                sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age([1 4 5])
                for iSubj = [1 4 5]
                    plot(1000 * tt(tt > ttmin & tt < ttmax), ...
                            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm(iSubj, tt > ttmin & tt < ttmax), 'k', 'LineWidth', 1);
                end
                
                % plot traces for a couple of older subjects
                sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age([28 30:32])
                for iSubj = [28 30:32]
                    plot(1000 * tt(tt > ttmin & tt < ttmax), ...
                            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm(iSubj, tt > ttmin & tt < ttmax), 'b', 'LineWidth', 1);
                end

                hold on
                set(gca, 'XTick', 0:50:100, 'YTick', [-250 0 250])
                axis tight
                ylim([-350 150])

                fill([-10 10 10 -10], [500 500 -500 -500], [.5 .5 .5], 'EdgeColor', 'w', 'FaceAlpha', .8)

                strSign = [' (p\_fdr = ', num2str(all_p_fdr{iTr}{iSubTr}(iDir + 1)), ')'];
                if all_p_fdr{iTr}{iSubTr}(iDir + 1) < .05,  strSign = [strSign, ' *'];     end
                title([rois(iTr).tract_name, ' - ', strSubTitle, strSign]);

                %
                % save
                %
                if ~exist(fullfile(myDataPath.output,'derivatives', 'age'), 'dir')
                    mkdir(fullfile(myDataPath.output,'derivatives', 'age'));
                end
                figureName = fullfile(myDataPath.output,'derivatives', 'age', ['CCEPexamples_6subjects', '_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_')]);

                set(gcf,'renderer','Painters')
                set(gcf,'PaperPositionMode','auto')
                print('-dpng','-r300',figureName)
                print('-depsc','-r300',figureName)
                close(gcf)
                
            end
        end
    end
end

