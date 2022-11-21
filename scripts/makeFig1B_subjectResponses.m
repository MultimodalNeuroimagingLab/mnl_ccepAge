% 
% 
%  This code is used to create Figure 1B - example CCEPs for a few older and younger subjects 
%


%% 
%  Load the ccepData and ccepAverages from the derivatives

clear
close all
warning('on');
warning('backtrace', 'off')

myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepAverages.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepAverages.mat'))
else
    disp('Run script ccep04_averageConnections.m first')
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
                figure('Position', [0 0 1800 200])
                plot(1000 * tt(tt > ttmin & tt < ttmax), ...
                        zeros(size(tt(tt > ttmin & tt < ttmax))),'Color',[.5 .5 .5]);
                hold on
                
                % plot traces for a couple of young subjects
                disp(['Ages young: ', num2str(sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age([1 4 5]))]);
                for iSubj = [1 4 5]
                    plot(1000 * tt(tt > ttmin & tt < ttmax), ...
                            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm(iSubj, tt > ttmin & tt < ttmax), 'k', 'LineWidth', 1);
                end
                
                % plot traces for a couple of older subjects
                %disp(['Ages older: ', num2str(sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age([28 30 32]))]);
                disp(['Ages older: ', num2str(sortedCCEPs{iTr}{iSubTr}{iDir + 1}.age([28 29 30]))]);
                for iSubj = [28 30 32]
                    plot(1000 * tt(tt > ttmin & tt < ttmax), ...
                            sortedCCEPs{iTr}{iSubTr}{iDir + 1}.averageResp_nonnorm(iSubj, tt > ttmin & tt < ttmax), 'b', 'LineWidth', 1);
                end

                hold on
                set(gca, 'XTick', 0:50:100, 'YTick', [-250 0 250])
                axis tight
                ylim([-400 180])
                fill([-10 10 10 -10], [500 500 -500 -500], [.5 .5 .5], 'EdgeColor', 'w', 'FaceAlpha', .8)
                title([rois(iTr).tract_name, ' - ', strSubTitle]);

                
                %
                % save
                %
                if ~exist(fullfile(myDataPath.output,'derivatives', 'age'), 'dir')
                    mkdir(fullfile(myDataPath.output,'derivatives', 'age'));
                end
                figureName = fullfile(myDataPath.output,'derivatives', 'age', ['CCEPexamples_6subjects', '_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_')]);

                set(gcf, 'renderer', 'Painters')
                set(gcf, 'PaperPositionMode', 'auto')
                print('-dpng', '-r300', figureName)
                print('-depsc', '-r300', figureName)
                close(gcf)
                
            end
            
        end
    end
end

