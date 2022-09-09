%
% This script produces Figure 3 of the article
%


%%
%  Load the n1Latencies from the derivatives

clear
close all

selectPat = input('Would you like to include all patients, or only the ones for whom it is certain that 8mA was applied (supplemental material)? [all/8] ','s');

if strcmp(selectPat, 'all')
    select_amplitude = 0; % make this 8 for only 8mA
elseif strcmp(selectPat, '8')
    select_amplitude = 8;
else
    error('Answer to previous question is not recognized.')
end
myDataPath = setLocalDataPath(1);

if select_amplitude == 0 
    if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'file')
        % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
        load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData')
    else
        disp('Run scripts ccep02_loadN1.m and ccep03_addtracts.m first')
    end
elseif select_amplitude == 8 % only 8 mA
    if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_8ma.mat'), 'file')
        % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
        load(fullfile(myDataPath.output,'derivatives', 'av_ccep', 'ccepData_8ma.mat'), 'n1Latencies8ma')

        n1Latencies = n1Latencies8ma;
        filename_averageCCEP = fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'average_ccep_age_8ma.mat');
    else
        disp('Run scripts ccep02_loadN1.m and ccep03_addtracts.m first')
    end
end



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
            
            % extract the latencies and number of N1s/CCEPs between the end-point ROIs for a specific (sub-)tract and direction
            temp = ccep_connectRegions( ccepData, ...
                                        rois(iTr).sub_tract(iSubTr).(['roi', num2str(iDir + 1)]), ...
                                        rois(iTr).sub_tract(iSubTr).(['roi', num2str(~iDir + 1)]));
            out{iTr}{iSubTr}{iDir + 1} = temp;
            
            % construct sub-tract string
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            if iDir == false
                strSubTitle = [subDir{1}, ' -> ', subDir{2}];
            else
                strSubTitle = [subDir{2}, ' -> ', subDir{1}];
            end
            
            % store the tract and ROIs designations
            out{iTr}{iSubTr}{iDir + 1}.name = [rois(iTr).tract_name, ' - ', strSubTitle];
            
            % add the native tract length 
            for iSubj = 1:length(out{iTr}{iSubTr}{iDir + 1}.sub)
                
                
                if any(contains(ccepData(iSubj).electrodes.jsonHemi, 'L')) && any(contains(ccepData(iSubj).electrodes.jsonHemi, 'R'))
                    % left and right
                    out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).nativeTractDist = mean(cell2mat(ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances), 'omitnan');
                    
                elseif any(contains(ccepData(iSubj).electrodes.jsonHemi, 'L'))
                    % left
                    out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).nativeTractDist = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{1};
                    
                else
                    % right
                    out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).nativeTractDist = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances{2};
                    
                end
                
            end
            
        end
    end
end



%%
%  ...

n1Type = 'Latency';
n1Type = 'Speed';

% loop over the (sub-)tracts and directions
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]

            % construct sub-tract string
            subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
            if iDir == false
                strSubTitle = [subDir{1}, ' -> ', subDir{2}];
            else
                strSubTitle = [subDir{2}, ' -> ', subDir{1}];
            end
            %%
            % initialize output: age, mean, variance in latency, and number of connections per subject
            nsubs = length(out{iTr}{iSubTr}{iDir + 1}.sub);
            subsValues = NaN(nsubs, 5);              % [subject, <age, mean latency, standard error, number of latency values, tract dist>]

            % retrieve the age, mean latency, standard error and number of latency values per subject
            for iSubj = 1:nsubs
                subsValues(iSubj, 1) = out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).age;
                subsValues(iSubj, 2) = mean(out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).latencies, 'omitnan');
                numLatencies = sum(~isnan(out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).latencies));
                subsValues(iSubj, 3) = std(out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).latencies, 'omitnan') ./ sqrt(numLatencies);
                subsValues(iSubj, 4) = numLatencies;
                
                % take the average distance of the tract in the left and right hemisphere
                subsValues(iSubj, 5) = out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).nativeTractDist;
                
                clear numLatencies nativeTractDist
            end

            % retrieve the unique ages and calculate mean across same age latencies
            ages = unique(sort(subsValues(~isnan(subsValues(:, 2)), 1)));
            n1Means = zeros(size(ages));
            for iAge = 1:length(ages)
                subjInclIndices = ismember(subsValues(:, 1), ages(iAge));
                if strcmpi(n1Type, 'speed')
                    % divide each subject by it's own length, then average
                    % across subjects of the same age in m/s
                    n1Means(iAge) = mean(.001 * subsValues(subjInclIndices, 5) ./ subsValues(subjInclIndices, 2), 'omitnan');
                    
                elseif strcmpi(n1Type, 'latency')
                    n1Means(iAge) = 1000 * mean(subsValues(subjInclIndices, 2), 'omitnan');
                    
                else
                    error('n1 type should either be set to ''latency'' or ''speed''');
                end
                
            end
            ages(isnan(n1Means)) = [];
            n1Means(isnan(n1Means)) = [];
            

            %
            % fit first order polynomial
            %
            
            % Test fitting a first order polynomial (leave 1 out cross validation)
            % y  =  w1*age + w2
            cross_val_linear = NaN(length(n1Means), 4);
            % size latency (ms) X prediction (ms) X p1 (slope) X p2 (intercept) of left out
            age_counter = 0;
            for iAge = 1:length(n1Means)
                age_counter = age_counter + 1;
                % leave out iAge
                subsTrain = ~ismember(1:length(n1Means), iAge)'; % leave out one age
                P = polyfit(ages(subsTrain), n1Means(subsTrain), 1);
                cross_val_linear(age_counter, 3:4) = P;                        % linear parameters
                cross_val_linear(age_counter, 1) = n1Means(iAge);              % measured N1 for iAge
                cross_val_linear(age_counter, 2) = P(1) * ages(iAge) + P(2);   % left out prediction for iAge
            end
            % coefficient of determination between prediction and left out
            % measurement
            out{iTr}{iSubTr}{iDir + 1}.cod_out(1) = calccod(cross_val_linear(:, 2), cross_val_linear(:, 1), 1);
            % correlation just for fun
            out{iTr}{iSubTr}{iDir + 1}.sp_out(1) = corr(cross_val_linear(:, 2), cross_val_linear(:, 1), 'type', 'Spearman');
            % average parameters for linear fit across ages
            out{iTr}{iSubTr}{iDir + 1}.linear_avparams = mean(cross_val_linear(:, 3:4));
            
            if isnan(out{iTr}{iSubTr}{iDir + 1}.cod_out(1))
                error('nan for cod');
            end
            
            %
            % fit second order polynomial
            %
            
            % Like Yeatman et al., for DTI fit a second order polynomial:
            cross_val_second = NaN(length(n1Means), 5);
            % size latency (ms) X prediction (ms) X p1 (age^2) X p2 (age) X p3 (intercept) of left out
            age_counter = 0;
            for iAge = 1:length(n1Means)
                age_counter = age_counter + 1;
                % leave out iAge
                subsTrain = ~ismember(1:length(n1Means), iAge)';
                P = polyfit(ages(subsTrain), n1Means(subsTrain), 2);
                cross_val_second(age_counter, 3:5) = P;
                cross_val_second(age_counter, 1) = n1Means(iAge);
                cross_val_second(age_counter, 2) = P(1) * ages(iAge) .^2 + P(2) * ages(iAge) + P(3);
            end
            out{iTr}{iSubTr}{iDir + 1}.cod_out(2) = calccod(cross_val_second(:, 2), cross_val_second(:, 1), 1);
            out{iTr}{iSubTr}{iDir + 1}.sp_out(2) = corr(cross_val_second(:, 2),cross_val_second(:, 1), 'type', 'Spearman');
            out{iTr}{iSubTr}{iDir + 1}.second_avparams = mean(cross_val_second(:, 3:5));

            out{iTr}{iSubTr}{iDir + 1}.cod_out(3) = length(n1Means); % number of ages
            
            
            %
            % make output figure
            %
            
            figure('position',[0 0 600 300])
            hold on;

            % plot vertical histogram per subject in background
            % add this if you want to see every single subject and the effect of averaging within an age group
            for iSubj = 1:nsubs
                
                if ~isnan(subsValues(iSubj, 2))
                    if strcmpi(n1Type, 'speed')
                        distributionPlot(.001 * out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).nativeTractDist ./ out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).latencies', ...
                                        'xValues', out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).age, ...
                                        'color', [.8 .8 .8], 'showMM', 0, 'histOpt', 2)
                    else

                        distributionPlot(1000 * out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).latencies', ...
                                        'xValues', out{iTr}{iSubTr}{iDir + 1}.sub(iSubj).age, ...
                                        'color', [.8 .8 .8], 'showMM', 0, 'histOpt', 2)
                    end
                end
                
                % plot mean + std err per subject
                %plot([my_output(kk, 1) my_output(kk, 1)], [1000 * (my_output(kk, 2) - my_output(kk, 3)) 1000 * (my_output(kk, 2)+my_output(kk, 3))], 'k', 'LineWidth', 1)
                
            end
            
            % plot mean per subject in a dot
            %plot(my_output(:, 1), 1000 * my_output(:, 2), 'ko', 'MarkerSize', 6)
            
            %
            % plot 1st or 2nd order polynomial
            %

            x_age = 1:1:max([ccepData.age]);
            
            % check whether the first or the second polynomial is a better fit
            if out{iTr}{iSubTr}{iDir + 1}.cod_out(1) > out{iTr}{iSubTr}{iDir + 1}.cod_out(2)
                % better fit with linear polynomial
                
                % prediction for every left out fit
                y_n1LatCross = cross_val_linear(:, 3) * x_age + cross_val_linear(:, 4);
                cmap = [0.6 0.2 1];

                out{iTr}{iSubTr}{iDir + 1}.fit = 'linear';
                out{iTr}{iSubTr}{iDir + 1}.delta = out{iTr}{iSubTr}{iDir + 1}.linear_avparams(1);
                out{iTr}{iSubTr}{iDir + 1}.cod = out{iTr}{iSubTr}{iDir + 1}.cod_out(1);
                
            elseif out{iTr}{iSubTr}{iDir + 1}.cod_out(1) < out{iTr}{iSubTr}{iDir + 1}.cod_out(2) 
                % better fit with second polynomial
                
                % prediction for every left out fit
                y_n1LatCross = cross_val_second(:, 3) * x_age .^ 2 + cross_val_second(:, 4) * x_age + cross_val_second(:, 5);
                cmap = [.2 0 1];
                delta_ages = 0:10:50;

                out{iTr}{iSubTr}{iDir + 1}.fit = 'second';
                out{iTr}{iSubTr}{iDir + 1}.delta = diff(out{iTr}{iSubTr}{iDir + 1}.second_avparams(1) * delta_ages .^ 2 + out{iTr}{iSubTr}{iDir + 1}.second_avparams(2) * delta_ages + out{iTr}{iSubTr}{iDir + 1}.second_avparams(3)) / 10;
                out{iTr}{iSubTr}{iDir + 1}.cod = out{iTr}{iSubTr}{iDir + 1}.cod_out(2);
            end

            
            %
            % plot confidence intervals and age
            %
            
            % get 95% confidence intervals
            low_ci = quantile(y_n1LatCross, .025, 1);
            up_ci = quantile(y_n1LatCross, .975, 1);
            if out{iTr}{iSubTr}{iDir + 1}.cod_out(3) < 20 % if less than 20 ages, plot light blue
                fill([x_age x_age(end:-1:1)], [low_ci up_ci(end:-1:1)], [.5 .7 1], 'EdgeColor', [.5 .7 1])
            else
                fill([x_age x_age(end:-1:1)], [low_ci up_ci(end:-1:1)], cmap, 'EdgeColor', cmap)
            end
            
            % if more than 20 subject and 2nd order, plot minimum 
            if out{iTr}{iSubTr}{iDir + 1}.cod_out(3) >= 20 && out{iTr}{iSubTr}{iDir + 1}.cod_out(2) > out{iTr}{iSubTr}{iDir + 1}.cod_out(1)
                
                % calculate minimum x
                min_age = -cross_val_second(:, 4) ./ (2 * cross_val_second(:, 3));
                if strcmpi(n1Type, 'speed')
                    plot([quantile(min_age, 0.025, 1) quantile(min_age, 0.975, 1)], [.05 .05], 'Color', [.2 .7 .6], 'LineWidth', 10);
                else
                    plot([quantile(min_age, 0.025, 1) quantile(min_age, 0.975, 1)], [5 5], 'Color', [.2 .7 .6], 'LineWidth', 10);
                end
                
            end
            
            
            % 
            title([rois(iTr).tract_name, ' - ', strSubTitle, ' - COD=' int2str(max(out{iTr}{iSubTr}{iDir + 1}.cod_out(1:2)))]); % plot maximal COD (1st or 2nd order)
            plot(ages, n1Means, 'k.', 'MarkerSize', 12);
            % xlabel('age (years)'), ylabel('mean dT (ms)')
            
            % 
            if strcmpi(n1Type, 'speed')
                xlim([0 60]), ylim([0 12]);
            else
                xlim([0 60]), ylim([0 80]);
                set(gca, 'XTick', 10:10:50, 'YTick', 20:20:100, 'FontName', 'Arial', 'FontSize', 12);
            end
            
            
            %
            % Save figure
            %
            
            if ~exist(fullfile(myDataPath.output, 'derivatives', 'age'), 'dir')
                mkdir(fullfile(myDataPath.output, 'derivatives', 'age'));
            end

            if select_amplitude == 0
                figureName = fullfile(myDataPath.output, 'derivatives', 'age', ['ageVs', n1Type, '_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_')]);
            elseif select_amplitude == 8
                figureName = fullfile(myDataPath.output, 'derivatives', 'age', ['ageVs', n1Type, '_', rois(iTr).tract_name, '_', strrep(strSubTitle, ' -> ', '_'), '_8mA']);
            end
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
clc 

% loop over the (sub-)tracts and directions
numOut = 0;
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]


            if strcmp(out{iTr}{iSubTr}{iDir + 1}.fit, 'linear')
                fprintf('%s: best fit is linear, with COD = %2.0f, and delta %1.2f \n', ...
                        out{iTr}{iSubTr}{iDir + 1}.name, out{iTr}{iSubTr}{iDir + 1}.cod, out{iTr}{iSubTr}{iDir + 1}.delta)
                
            elseif strcmp(out{iTr}{iSubTr}{iDir + 1}.fit, 'second')
                
                fprintf('%s: best fit is second, with COD = %2.0f, and delta: age0-10 = %1.2f, age10-20 = %1.2f, age20-30 = %1.2f, age30-40 = %1.2f, age40-50 = %1.2f \n', ...
                        out{iTr}{iSubTr}{iDir + 1}.name, out{iTr}{iSubTr}{iDir + 1}.cod, out{iTr}{iSubTr}{iDir + 1}.delta);
                
            end
            numOut = numOut + 1;
        end
    end
end


%%
%  Find average latencies
clc

delta_all = [];
y_lin = NaN(numOut, 3);
y_sec = NaN(numOut, 3);
min_age = NaN(numOut, 1);
connection = cell(numOut, 1);
fit = cell(numOut, 1);

% loop over the (sub-)tracts and directions
ii = 1;
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]

            if strcmp(out{iTr}{iSubTr}{iDir + 1}.fit, 'linear') && out{iTr}{iSubTr}{iDir + 1}.cod > 0
                delta_all = [delta_all, out{iTr}{iSubTr}{iDir + 1}.delta];

                age = [4, 25, 51];
                y_lin(ii, 1:3) = out{iTr}{iSubTr}{iDir + 1}.linear_avparams(1) * age + out{iTr}{iSubTr}{iDir + 1}.linear_avparams(2);
                connection{ii} = out{iTr}{iSubTr}{iDir + 1}.name;
                fit{ii} = out{iTr}{iSubTr}{iDir + 1}.fit;
                
            elseif strcmp(out{iTr}{iSubTr}{iDir + 1}.fit, 'second') && out{iTr}{iSubTr}{iDir + 1}.cod > 0
                min_age(ii) = -out{iTr}{iSubTr}{iDir + 1}.second_avparams(2) ./ (2 * out{iTr}{iSubTr}{iDir + 1}.second_avparams(1));
                connection{ii} = out{iTr}{iSubTr}{iDir + 1}.name;
                fit{ii} = out{iTr}{iSubTr}{iDir + 1}.fit;

                age = [4, min_age(ii), 51];
                y_sec(ii, 1:3) = out{iTr}{iSubTr}{iDir + 1}.second_avparams(1) * age .^ 2 + out{iTr}{iSubTr}{iDir + 1}.second_avparams(2) * age + out{iTr}{iSubTr}{iDir + 1}.second_avparams(3);

            else
                connection{ii} = out{iTr}{iSubTr}{iDir + 1}.name;        
            end
            
            ii = ii + 1;
        end
    end
end
fprintf('\n         LINEAR MODEL FIT \n')
fprintf('mean delta (min-max) = %1.2fms/year (%1.2f - %1.2f)\n', ...
        mean(delta_all), min(delta_all),max(delta_all))
fprintf('Mean latency at age 4 years: %1.2f ms \nMean latency at age 51 years: %1.2f ms\n \n',...
        mean(y_lin(:,1),'omitnan'), mean(y_lin(:,3),'omitnan'))

fprintf('         SECOND ORDER MODEL FIT \n')
delta_sec = diff(y_sec, [], 2) ./ diff([repmat(4, numOut, 1), min_age, repmat(51, numOut, 1)], [], 2);

fprintf('Mean minimal age (min-max) = %1.2f years (%1.2f - %1.2f)\n',...
    mean(min_age,'omitnan'), min(min_age),max(min_age))
fprintf('Mean delta until minimal latency (min-max) = %1.2fms/year (%1.2f - %1.2f)\n',...
    mean(delta_sec(:,1),'omitnan'), min(delta_sec(:,1)), max(delta_sec(:,1)))
fprintf('Mean delta after minimal latency (min-max) = %1.2fms/year (%1.2f - %1.2f)\n',...
    mean(delta_sec(:,2),'omitnan'), min(delta_sec(:,2)), max(delta_sec(:,2)))
fprintf('Minimal latency (min-max) = %1.2f ms (%1.2f - %1.2f)\n \n',...
    mean(y_sec(:,2),'omitnan'), min(y_sec(:,2)),max(y_sec(:,2)))

% show all latencies at age 4, minimal age/25years, 51 years. 
y = y_lin;
y(~isnan(y_sec(:, 1)), 1:3) = y_sec(~isnan(y_sec(:, 1)), 1:3);

disp([{'Connection'} ,{'Fit'}, {'Latency (4)'}, {'Latency(25/min_age)'}, {'Latency(51)'}, {'min_age'}; connection(:), fit(:), num2cell(y), num2cell(min_age)])

 
%% extra explained variance

%{
cod_out_check = cod_out;
cod_out_check(cod_out<0) = 0;
betterlinearfit = cod_out_check(:,1) - cod_out_check(:,2);
betterlinearfit(betterlinearfit<=0) = NaN;
mean(betterlinearfit,'omitnan')
%}