% 
% 
%  This code is used to create Figure 1C - the normalized CCEPs of all patients in order of age, 
%  
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

                set(gcf, 'PaperPositionMode', 'auto')
                print('-dpng', '-r300', figureName)
                print('-depsc', '-r300', figureName)
                close(gcf)
                
            end
        end
    end
end



