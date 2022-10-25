%
% This script produces Figure 3 of the article
%


%%
%  Load the n1Latencies from the derivatives

clear
close all
warning('on');
warning('backtrace', 'off')

myDataPath = setLocalDataPath(1);

if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData')
else
    disp('Run scripts ccep02_loadN1.m and ccep03_addtracts.m first')
end

stimStimElec_excludeDist = 18;     % the distance between the stimulated electrodes (in mm) above which N1s are excluded, 0 = not excluding
respStimElec_excludeDist = 13;     % the distance between a stimulated and response electrode (in mm) within which N1s are excluded, 0 = not excluding


%%
%  extract for all (sub-)tract ROIs and directions the latencies and number of N1s/CCEPs for each of the subject

clear out
clc

% load tracts and their corresponding end-point ROIs
rois = ccep_categorizeAnatomicalRegions();

% out will get a cell array in tracts X subtracts X directions
out = {}; 

% loop over the (sub-)tracts and directions
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract) 
        % direction along tract
        for iDir = [false true] 
            
            % construct sub-tract string (that includes the stim-resp direction) and store
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            strSubTitle = [subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];
            out{iTr}{iSubTr}{iDir + 1}.name = [rois(iTr).tract_name, ' - ', strSubTitle];
            
            % extract the latencies and number of N1s/CCEPs between the end-point ROIs for a specific (sub-)tract and direction
            out{iTr}{iSubTr}{iDir + 1}.metrics = ccep_N1sBetweenRegions(ccepData, ...
                                                                        rois(iTr).sub_tract(iSubTr).(['roi', num2str(iDir + 1)]), ...
                                                                        rois(iTr).sub_tract(iSubTr).(['roi', num2str(~iDir + 1)]), ...
                                                                        stimStimElec_excludeDist, respStimElec_excludeDist);
            
            % add the native tract length 
            for iSubj = 1:length(out{iTr}{iSubTr}{iDir + 1}.metrics)
                if any(contains(ccepData(iSubj).electrodes.jsonHemi, 'L')) && any(contains(ccepData(iSubj).electrodes.jsonHemi, 'R'))
                    % left and right
                    out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).nativeTractDist = mean(cell2mat(ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances), 'omitnan');
                    
                elseif any(contains(ccepData(iSubj).electrodes.jsonHemi, 'L'))
                    % only left
                    out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).nativeTractDist = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{1};
                    
                else
                    % only right
                    out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).nativeTractDist = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{2};
                    
                end
            end
            
        end
    end
end



%%
%   Generate the images for figure 2 and 3


% loop over the (sub-)tracts and directions
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]

            % construct sub-tract string (that includes the stim-resp direction)
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            strSubTitle = [subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];

            % initialize output: age, mean, variance in latency, and number of connections per subject
            nsubs = length(out{iTr}{iSubTr}{iDir + 1}.metrics);
            subsValues = NaN(nsubs, 5);              % [subject, <age, mean latency, standard error, number of latency values, tract dist>]

            % retrieve the age, mean latency, standard error and number of latency values per subject
            for iSubj = 1:nsubs
                subsValues(iSubj, 1) = out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).age;
                subsValues(iSubj, 2) = mean(out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).latencies, 'omitnan');
                numLatencies = sum(~isnan(out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).latencies));
                subsValues(iSubj, 3) = std(out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).latencies, 'omitnan') ./ sqrt(numLatencies);
                subsValues(iSubj, 4) = numLatencies;
                
                % take the average distance of the tract in the left and right hemisphere
                subsValues(iSubj, 5) = out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).nativeTractDist;
                
                clear numLatencies nativeTractDist
            end

            % retrieve the unique ages and calculate mean across same age latencies
            ages = unique(sort(subsValues(~isnan(subsValues(:, 2)), 1)));
            n1LatencyMeans = zeros(size(ages));
            n1SpeedMeans = zeros(size(ages));
            for iAge = 1:length(ages)
                subjInclIndices = ismember(subsValues(:, 1), ages(iAge));
                
                % latency
                n1LatencyMeans(iAge) = 1000 * mean(subsValues(subjInclIndices, 2), 'omitnan');
                
                % speed - divide each subject by it's own length, then average across subjects of the same age in m/s
                n1SpeedMeans(iAge) = mean(.001 * subsValues(subjInclIndices, 5) ./ subsValues(subjInclIndices, 2), 'omitnan');
                
            end
            
            % exclude 
            exclMeans = isnan(n1LatencyMeans) & isnan(n1SpeedMeans);
            ages(exclMeans) = [];
            n1LatencyMeans(exclMeans) = [];
            n1SpeedMeans(exclMeans) = [];
            
            lat_cross_linear = NaN(length(n1LatencyMeans), 4);      % <age> x <size latency (ms), prediction (ms), p1 (slope), p2 (intercept) of left out>
            spd_cross_linear = NaN(length(n1LatencyMeans), 4);      % <age> x <size latency (ms), prediction (ms), p1 (slope), p2 (intercept) of left out>
            lat_cross_second = NaN(length(n1LatencyMeans), 5);      % <age> x <size latency (ms) X prediction (ms) X p1 (age^2) X p2 (age) X p3 (intercept) of left out>
            spd_cross_second = NaN(length(n1LatencyMeans), 5);      % <age> x <size latency (ms) X prediction (ms) X p1 (age^2) X p2 (age) X p3 (intercept) of left out>
            
            
            
            %
            % fit first order polynomial
            %
            
            % loop over the ages (for leave-one-out)
            age_counter = 0;
            for iAge = 1:length(n1LatencyMeans)
                age_counter = age_counter + 1;
                
                % determine training set of N1s (leaving one age out)
                subsTrain = ~ismember(1:length(n1LatencyMeans), iAge)';
                
                % fit linear on training set of N1s
                lat_polyfit = robustfit([ages(subsTrain)], n1LatencyMeans(subsTrain),'bisquare',4.685);% default
                lat_polyfit = lat_polyfit(end:-1:1); % reverse order to match previous polyfit
                spd_polyfit = robustfit([ages(subsTrain)], n1SpeedMeans(subsTrain),'bisquare',4.685);% default
                spd_polyfit = spd_polyfit(end:-1:1); % reverse order to match previous polyfit

                % 
                lat_cross_linear(age_counter, 3:4) = lat_polyfit;                                   % linear parameters
                lat_cross_linear(age_counter, 1) = n1LatencyMeans(iAge);                            % measured N1 for iAge
                lat_cross_linear(age_counter, 2) = lat_polyfit(1) * ages(iAge) + lat_polyfit(2);    % left out prediction for iAge
                
                % 
                spd_cross_linear(age_counter, 3:4) = spd_polyfit;                                   % linear parameters
                spd_cross_linear(age_counter, 1) = n1SpeedMeans(iAge);                              % measured N1 for iAge
                spd_cross_linear(age_counter, 2) = spd_polyfit(1) * ages(iAge) + spd_polyfit(2);    % left out prediction for iAge
                
            end
            
            % coefficient of determination between prediction and left out measurement
            out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1) = calccod(lat_cross_linear(:, 2), lat_cross_linear(:, 1), 1, [], 1);
            out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1) = calccod(spd_cross_linear(:, 2), spd_cross_linear(:, 1), 1, [], 1);
            
            %
            out{iTr}{iSubTr}{iDir + 1}.lat_sp_out(1) = corr(lat_cross_linear(:, 2), lat_cross_linear(:, 1), 'type', 'Pearson');
            out{iTr}{iSubTr}{iDir + 1}.spd_sp_out(1) = corr(spd_cross_linear(:, 2), spd_cross_linear(:, 1), 'type', 'Pearson');
            
            % average parameters for linear fit across ages
            out{iTr}{iSubTr}{iDir + 1}.lat_linear_avparams = mean(lat_cross_linear(:, 3:4));
            out{iTr}{iSubTr}{iDir + 1}.lat_linear_CIinter = [quantile(lat_cross_linear(:, 4), .025, 1) quantile(lat_cross_linear(:, 4), .0975, 1)];
            
            out{iTr}{iSubTr}{iDir + 1}.spd_linear_avparams = mean(spd_cross_linear(:, 3:4));
            out{iTr}{iSubTr}{iDir + 1}.spd_linear_CIinter = [quantile(spd_cross_linear(:, 4), .025, 1) quantile(spd_cross_linear(:, 4), .0975, 1)];
            
            % 
            if isnan(out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1)) || isnan(out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1))
                error('nan for cod');
            end
            
            
            
            %
            % fit second order polynomial
            %
            
            % Like Yeatman et al., for DTI fit a second order polynomial:
            
            % 
            age_counter = 0;
            for iAge = 1:length(n1LatencyMeans)
                age_counter = age_counter + 1;
                
                % determine training set of N1s (leaving one age out)
                subsTrain = ~ismember(1:length(n1LatencyMeans), iAge)';
                
                % fit second-order polynomial on training set of N1s
                lat_polyfit = robustfit([ages(subsTrain) ages(subsTrain).^2], n1LatencyMeans(subsTrain),'bisquare',4.685);% default
                lat_polyfit = lat_polyfit(end:-1:1); % reverse order to match previous polyfit
                spd_polyfit = robustfit([ages(subsTrain) ages(subsTrain).^2], n1SpeedMeans(subsTrain),'bisquare',4.685);% default
                spd_polyfit = spd_polyfit(end:-1:1); % reverse order to match previous polyfit
                
                %
                lat_cross_second(age_counter, 3:5) = lat_polyfit;
                lat_cross_second(age_counter, 1) = n1LatencyMeans(iAge);
                lat_cross_second(age_counter, 2) = lat_polyfit(1) * ages(iAge) .^2 + lat_polyfit(2) * ages(iAge) + lat_polyfit(3);
                
                %
                spd_cross_second(age_counter, 3:5) = spd_polyfit;
                spd_cross_second(age_counter, 1) = n1SpeedMeans(iAge);
                spd_cross_second(age_counter, 2) = spd_polyfit(1) * ages(iAge) .^2 + spd_polyfit(2) * ages(iAge) + spd_polyfit(3);
                
            end
            
            %
            out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(2) = calccod(lat_cross_second(:, 2), lat_cross_second(:, 1), 1, [], 1);
            out{iTr}{iSubTr}{iDir + 1}.lat_sp_out(2) = corr(lat_cross_second(:, 2), lat_cross_second(:, 1), 'type', 'Pearson');
            
            out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(2) = calccod(spd_cross_second(:, 2), spd_cross_second(:, 1), 1, [], 1);
            out{iTr}{iSubTr}{iDir + 1}.spd_sp_out(2) = corr(spd_cross_second(:, 2), spd_cross_second(:, 1), 'type', 'Pearson');
            
            % average parameters for second-order fit across ages
            out{iTr}{iSubTr}{iDir + 1}.lat_second_avparams = mean(lat_cross_second(:, 3:5));
            out{iTr}{iSubTr}{iDir + 1}.lat_second_CIinter = [quantile(lat_cross_second(:, 5), .025, 1) quantile(lat_cross_second(:, 5), .0975, 1)];
            
            out{iTr}{iSubTr}{iDir + 1}.spd_second_avparams = mean(spd_cross_second(:, 3:5));
            out{iTr}{iSubTr}{iDir + 1}.spd_second_CIinter = [quantile(spd_cross_second(:, 5), .025, 1) quantile(spd_cross_second(:, 5), .0975, 1)];
            
            % store the number of ages
            out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(3) = length(n1LatencyMeans);
            out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(3) = length(n1SpeedMeans);
            

            
            %
            % determine the best fit (1st or 2nd order polynomial)
            %

            x_age = 1:1:max([ccepData.age]);
            
            % check whether the first or the second polynomial is a better fit
            if out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1) > out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(2)
                % better fit with linear polynomial
                
                % prediction for every left out fit
                y_lat_n1 = lat_cross_linear(:, 3) * x_age + lat_cross_linear(:, 4);
                lat_cmap = [0.482352941176471, 0.341176470588235, 0.639215686274510];

                out{iTr}{iSubTr}{iDir + 1}.lat_fit      = 'linear';
                out{iTr}{iSubTr}{iDir + 1}.lat_delta    = out{iTr}{iSubTr}{iDir + 1}.lat_linear_avparams(1);
                out{iTr}{iSubTr}{iDir + 1}.lat_cod      = out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1);
                
            elseif out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1) < out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(2) 
                % better fit with second polynomial
                
                % prediction for every left out fit
                y_lat_n1 = lat_cross_second(:, 3) * x_age .^ 2 + lat_cross_second(:, 4) * x_age + lat_cross_second(:, 5);
                lat_cmap = [0.964705882352941, 0.674509803921569, 0.756862745098039];
                delta_ages = 0:10:50;

                out{iTr}{iSubTr}{iDir + 1}.lat_fit = 'second';
                out{iTr}{iSubTr}{iDir + 1}.lat_delta = diff(out{iTr}{iSubTr}{iDir + 1}.lat_second_avparams(1) * delta_ages .^ 2 + out{iTr}{iSubTr}{iDir + 1}.lat_second_avparams(2) * delta_ages + out{iTr}{iSubTr}{iDir + 1}.lat_second_avparams(3)) / 10;
                out{iTr}{iSubTr}{iDir + 1}.lat_cod = out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(2);
                
                clear delta_ages;
            end
            
            % check whether the first or the second polynomial is a better fit
            if out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1) > out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(2)
                % better fit with linear polynomial
                
                % prediction for every left out fit
                y_spd_n1 = spd_cross_linear(:, 3) * x_age + spd_cross_linear(:, 4);
                spd_cmap = [0.482352941176471 0.341176470588235 0.639215686274510];

                out{iTr}{iSubTr}{iDir + 1}.spd_fit      = 'linear';
                out{iTr}{iSubTr}{iDir + 1}.spd_delta    = out{iTr}{iSubTr}{iDir + 1}.spd_linear_avparams(1);
                out{iTr}{iSubTr}{iDir + 1}.spd_cod      = out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1);
                
            elseif out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1) < out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(2) 
                % better fit with second polynomial
                
                % prediction for every left out fit
                y_spd_n1 = spd_cross_second(:, 3) * x_age .^ 2 + spd_cross_second(:, 4) * x_age + spd_cross_second(:, 5);
                spd_cmap = [0.964705882352941, 0.674509803921569, 0.756862745098039];
                delta_ages = 0:10:50;

                out{iTr}{iSubTr}{iDir + 1}.spd_fit = 'second';
                out{iTr}{iSubTr}{iDir + 1}.spd_delta = diff(out{iTr}{iSubTr}{iDir + 1}.spd_second_avparams(1) * delta_ages .^ 2 + out{iTr}{iSubTr}{iDir + 1}.spd_second_avparams(2) * delta_ages + out{iTr}{iSubTr}{iDir + 1}.spd_second_avparams(3)) / 10;
                out{iTr}{iSubTr}{iDir + 1}.spd_cod = out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(2);
                
                clear delta_ages;
            end
            
            
            
            %
            % latency output figure
            %
            
            if strfind(rois(iTr).tract_name, '_U')
                figure('position',[0 0 600 400])
            else
                figure('position',[0 0 600 300])
            end
            hold on;
            
            % plot vertical histogram per subject in background
            % (shows every single subject and the effect of averaging within an age group)
            warning('off');
            for iSubj = 1:nsubs
                
                if ~isnan(subsValues(iSubj, 2))
                    distributionPlot(1000 * out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).latencies', ...
                                    'xValues', out{iTr}{iSubTr}{iDir + 1}.metrics(iSubj).age, ...
                                    'color', [.8 .8 .8], 'showMM', 0, 'histOpt', 2)
                end
                
                
            end
            warning('on');
            warning('backtrace', 'off')
                        
            % plot 95% CI
            low_ci = quantile(y_lat_n1, .025, 1);
            up_ci = quantile(y_lat_n1, .975, 1);
            fill([x_age x_age(end:-1:1)], [low_ci up_ci(end:-1:1)], lat_cmap, 'EdgeColor', lat_cmap)
            
            % plot fit trend (with confidence interval)
            plot(ages, n1LatencyMeans, 'k.', 'MarkerSize', 12, 'Color', [0 0 0]);

            % check if more than 20 subjects and 2nd order
            if out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(3) >= 20 && out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(2) > out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1)
                
                % plot 95% CI around age
                min_age = -lat_cross_second(:, 4) ./ (2 * lat_cross_second(:, 3));
                plot([quantile(min_age, 0.025, 1) quantile(min_age, 0.975, 1)], [5 5], 'Color', [.2 .7 .6], 'LineWidth', 10);
                disp([rois(iTr).tract_name, ' - ', strSubTitle, '  CI: ', num2str(quantile(min_age, 0.025, 1)), ' - ', num2str(quantile(min_age, 0.975, 1))]);
                
            end
            
            % plot 95% CI around latency
            if out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(2) > out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1)
                % 2nd order
                plot([0 0], [out{iTr}{iSubTr}{iDir + 1}.lat_second_CIinter], 'Color', [.2 .7 .6], 'LineWidth', 10);
                
            elseif out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(2) < out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1)
                % 1nd order
                plot([0 0], [out{iTr}{iSubTr}{iDir + 1}.lat_linear_CIinter], 'Color', [.2 .7 .6], 'LineWidth', 10);
                
            end
            
            % labels and axis
            title([rois(iTr).tract_name, ' - ', strSubTitle, ' - COD=' int2str(max(out{iTr}{iSubTr}{iDir + 1}.lat_cod_out(1:2)))]); % plot maximal COD (1st or 2nd order)
            xlim([0 60]), ylim([0 80]);
            set(gca, 'YTick', 20:20:100, 'FontName', 'Arial', 'FontSize', 12);
            set(gca, 'YGrid', 'on', 'XGrid', 'off');
            
            % save latency figure
            if ~exist(fullfile(myDataPath.output, 'derivatives', 'ageRobust'), 'dir')
                mkdir(fullfile(myDataPath.output, 'derivatives', 'ageRobust'));
            end
            figureName = fullfile(myDataPath.output, 'derivatives', 'ageRobust', ['ageVsLatency', '_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_')]);
            set(gcf,'PaperPositionMode', 'auto');
            print('-dpng', '-r300', figureName);
            print('-depsc', '-r300', figureName);
            close(gcf)

            
            
            %
            % speed output figure
            %
            
            if strfind(rois(iTr).tract_name, '_U')
                figure('position',[0 0 600 300])
            else
                figure('position',[0 0 600 200])
            end
            hold on;
            
            % plot fit trend (with confidence interval)
            plot(ages, n1SpeedMeans, 'k.', 'MarkerSize', 12, 'Color', [0 0 0]);

            % plot 95% CI
            low_ci = quantile(y_spd_n1, .025, 1);
            up_ci = quantile(y_spd_n1, .975, 1);
            fill([x_age x_age(end:-1:1)], [low_ci up_ci(end:-1:1)], spd_cmap, 'EdgeColor', spd_cmap)
            
            % check if more than 20 subjects and 2nd order
            if out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(3) >= 20 && out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(2) > out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1)
                
                % plot 95% CI around age
                min_age = -spd_cross_second(:, 4) ./ (2 * spd_cross_second(:, 3));
                plot([quantile(min_age, 0.025, 1) quantile(min_age, 0.975, 1)], [5 5], 'Color', [.2 .7 .6], 'LineWidth', 10);
                disp([rois(iTr).tract_name, ' - ', strSubTitle, '  CI: ', num2str(quantile(min_age, 0.025, 1)), ' - ', num2str(quantile(min_age, 0.975, 1))]);
                
            end

            % plot 95% CI around latency
            if out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(2) > out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1)
                % 2nd order
                plot([0 0], [out{iTr}{iSubTr}{iDir + 1}.spd_second_CIinter], 'Color', [.2 .7 .6], 'LineWidth', 10);
                
            elseif out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(2) < out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1)
                % 1nd order
                plot([0 0], [out{iTr}{iSubTr}{iDir + 1}.spd_linear_CIinter], 'Color', [.2 .7 .6], 'LineWidth', 10);
                
            end
            
            % labels and axis
            title([rois(iTr).tract_name, ' - ', strSubTitle, ' - COD=' int2str(max(out{iTr}{iSubTr}{iDir + 1}.spd_cod_out(1:2)))]); % plot maximal COD (1st or 2nd order)
            set(gca, 'yaxislocation', 'right');
            if strfind(rois(iTr).tract_name, '_U')
                xlim([0 60]), ylim([0 4]);
                set(gca, 'YTick', 0:1:4, 'FontName', 'Arial', 'FontSize', 12);
            else
                xlim([0 60]), ylim([0 12]);
                set(gca, 'YTick', 0:4:12, 'FontName', 'Arial', 'FontSize', 12);
            end
            set(gca, 'YGrid', 'on', 'XGrid', 'off');
            
            % save speed figure
            figureName = fullfile(myDataPath.output, 'derivatives', 'ageRobust', ['ageVsSpeed', '_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_')]);
            set(gcf,'PaperPositionMode', 'auto');
            print('-dpng', '-r300', figureName);
            print('-depsc', '-r300', figureName);
            close(gcf)

        end
    end
end

%% 
%  Display in command window the cod and delta for each subplot
%  this info is displayed in Figure 3 as well.

long_tracts = {'Y1065_TPAT - Parietal -> Temporal','Y1065_TPAT - Temporal -> Parietal',...
    'Y1065_AF - Frontal -> Temporal','Y1065_AF - Temporal -> Frontal',...
    'Y1065_SLF2 - Frontal -> Parietal','Y1065_SLF2 - Parietal -> Frontal',...
    'Y1065_SLF2 - Frontal -> Central','Y1065_SLF2 - Central -> Frontal'};
long_cod_lat = NaN(1,length(long_tracts));
long_cod_spd = NaN(1,length(long_tracts));

U_tracts = {'Y842_U - PreCentral -> PostCentral','Y842_U - PostCentral -> PreCentral',...
    'Y842_U - Frontal -> Frontal','Y842_U - Parietal -> Parietal'};
U_cod_lat = NaN(1,length(U_tracts));
U_cod_spd = NaN(1,length(U_tracts));

% latency
disp('Latency: ');
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]

            if strcmp(out{iTr}{iSubTr}{iDir + 1}.lat_fit, 'linear')
                fprintf(' - %s: best fit is linear, with COD = %2.0f, and delta %1.2f \n', ...
                        out{iTr}{iSubTr}{iDir + 1}.name, out{iTr}{iSubTr}{iDir + 1}.lat_cod, out{iTr}{iSubTr}{iDir + 1}.lat_delta)
                
            elseif strcmp(out{iTr}{iSubTr}{iDir + 1}.lat_fit, 'second')
                
                fprintf(' - %s: best fit is second, with COD = %2.0f, and delta: age0-10 = %1.2f, age10-20 = %1.2f, age20-30 = %1.2f, age30-40 = %1.2f, age40-50 = %1.2f \n', ...
                        out{iTr}{iSubTr}{iDir + 1}.name, out{iTr}{iSubTr}{iDir + 1}.lat_cod, out{iTr}{iSubTr}{iDir + 1}.lat_delta);
                
            end
            
            % get averages
            if ismember(out{iTr}{iSubTr}{iDir + 1}.name,U_tracts)
                U_cod_lat(find(ismember(U_tracts,out{iTr}{iSubTr}{iDir + 1}.name))) = out{iTr}{iSubTr}{iDir + 1}.lat_cod;
            elseif ismember(out{iTr}{iSubTr}{iDir + 1}.name,long_tracts)
                long_cod_lat(find(ismember(long_tracts,out{iTr}{iSubTr}{iDir + 1}.name))) = out{iTr}{iSubTr}{iDir + 1}.lat_cod;
            end

        end
    end
end

% speed
disp('speed: ');
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]

            if strcmp(out{iTr}{iSubTr}{iDir + 1}.spd_fit, 'linear')
                fprintf(' - %s: best fit is linear, with COD = %2.0f, and delta %1.2f \n', ...
                        out{iTr}{iSubTr}{iDir + 1}.name, out{iTr}{iSubTr}{iDir + 1}.spd_cod, out{iTr}{iSubTr}{iDir + 1}.spd_delta)
                
            elseif strcmp(out{iTr}{iSubTr}{iDir + 1}.spd_fit, 'second')
                
                fprintf(' - %s: best fit is second, with COD = %2.0f, and delta: age0-10 = %1.2f, age10-20 = %1.2f, age20-30 = %1.2f, age30-40 = %1.2f, age40-50 = %1.2f \n', ...
                        out{iTr}{iSubTr}{iDir + 1}.name, out{iTr}{iSubTr}{iDir + 1}.spd_cod, out{iTr}{iSubTr}{iDir + 1}.spd_delta);
                
            end
            % get averages
            if ismember(out{iTr}{iSubTr}{iDir + 1}.name,U_tracts)
                U_cod_spd(find(ismember(U_tracts,out{iTr}{iSubTr}{iDir + 1}.name))) = out{iTr}{iSubTr}{iDir + 1}.spd_cod;
            elseif ismember(out{iTr}{iSubTr}{iDir + 1}.name,long_tracts)
                long_cod_spd(find(ismember(long_tracts,out{iTr}{iSubTr}{iDir + 1}.name))) = out{iTr}{iSubTr}{iDir + 1}.spd_cod;
            end

        end
    end
end

mean(long_cod_lat)
mean(long_cod_spd)
mean(U_cod_lat)
mean(U_cod_spd)

