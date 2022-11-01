%
% This script produces supplementary figures 2-5
%
% Dora Hermes, Dorien van Blooijs, Max van den Boom, 2022

clear
close all
clc


%% 
%  Load the CCEP data from the derivatives

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'file')
    load(fullfile(myDataPath.output,'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'ccepData')
else
    disp('Run scripts ccep02_aggregateToStruct.m first')
end

stimStimElec_excludeDist = 18;     % the distance between the stimulated electrodes (in mm) above which N1s are excluded, 0 = not excluding
respStimElec_excludeDist = 13;     % the distance between a stimulated and response electrode (in mm) within which N1s are excluded, 0 = not excluding

% load tracts and their corresponding end-point ROIs
rois = ccep_categorizeAnatomicalRegions();



%% 
%  Define the connection to be displayed and prepare the data

conn_matrix = {[3 1 0], [4 1 0], [4 2 0], [1 1 0]; ...
               [3 1 1], [4 1 1], [4 2 1], [1 1 1]; ...
               };
           
all_varlat_p = [];
all_meanvarlat_p = [];
all_ratioN1s_p = [];
           
%
out = cell(size(conn_matrix));
for iRow = 1:size(conn_matrix, 1)
    for iCol = 1:size(conn_matrix, 2)
        iTr = conn_matrix{iRow, iCol}(1);
        iSubTr = conn_matrix{iRow, iCol}(2);
        iDir = conn_matrix{iRow, iCol}(3);
        
        % extract the connection name (includes tract-subtract and direction)
        strName = rois(iTr).tract_name;
        subDir = split(rois(iTr).sub_tract(iSubTr).name, '-');
        strName = [strName, ' - ',[subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}]];
        
        %
        disp(['Running - ', strName]);
        out{iRow, iCol}.name = strName;
        
        % extract the latencies and number of N1s between the end-point ROIs for a specific (sub-)tract and direction
        metrics = ccep_N1sBetweenRegions( ccepData, ...
                                          rois(iTr).sub_tract(iSubTr).(['roi', num2str(iDir + 1)]), ...
                                          rois(iTr).sub_tract(iSubTr).(['roi', num2str(~iDir + 1)]), ...
                                          stimStimElec_excludeDist, respStimElec_excludeDist);
        
        % retrieve metrics per subject, output format:
        % <subject> x <age, mean in latencies, variance in latencies, variance in (latencies * 1000), relative number of N1s>
        subsValues = NaN(length(metrics), 4);
        for iSubj = 1:length(metrics)
            subsValues(iSubj, 1) =      metrics(iSubj).age;
            subsValues(iSubj, 2) = mean(metrics(iSubj).latencies, 'omitnan');
            subsValues(iSubj, 3) =  var(metrics(iSubj).latencies, 'omitnan');
            subsValues(iSubj, 4) =  var(1000 * metrics(iSubj).latencies, 'omitnan');
            subsValues(iSubj, 5) = mean(metrics(iSubj).ratioN1s);
            
        end
        out{iRow, iCol}.subsValues = subsValues;
        
        % calculate and store the r and p for the age vs variance in latency
        [varlat_r, varlat_p] = corr(subsValues(~isnan(subsValues(:, 2)), 1), ...
                                    subsValues(~isnan(subsValues(:, 2)), 3), 'Type', 'Spearman');
        out{iRow, iCol}.varlat_r = varlat_r;
        out{iRow, iCol}.varlat_p = varlat_p;
        all_varlat_p(end + 1, :) = [iRow, iCol, varlat_p];
        
        % calculate and store the r and p for mean latency (x1000) vs variance in latency (base values already multiplied before by 1000)
        [meanvarlat_r, meanvarlat_p] = corr(subsValues(~isnan(subsValues(:, 2)), 2) * 1000, ...
                                            subsValues(~isnan(subsValues(:, 2)), 4), 'Type', 'Spearman');
        out{iRow, iCol}.meanvarlat_r = meanvarlat_r;
        out{iRow, iCol}.meanvarlat_p = meanvarlat_p;
        all_meanvarlat_p(end + 1, :) = [iRow, iCol, meanvarlat_p];

        % calculate and store the r and p for the age vs ratio of N1s
        [ratioN1s_r, ratioN1s_p] = corr(subsValues(~isnan(subsValues(:, 2)), 1), ...
                                        subsValues(~isnan(subsValues(:, 2)), 5), 'Type', 'Spearman');
        out{iRow, iCol}.ratioN1s_r = ratioN1s_r;
        out{iRow, iCol}.ratioN1s_p = ratioN1s_p;
        all_ratioN1s_p(end + 1, :) = [iRow, iCol, ratioN1s_p];
        
    end
end

% calculate the FDR corrected P-values
[~, ~, ~, all_varlat_p(:, 4)]  = fdr_bh(all_varlat_p(:, 3), 0.05, 'pdep');
[~, ~, ~, all_meanvarlat_p(:, 4)]  = fdr_bh(all_meanvarlat_p(:, 3), 0.05, 'pdep');
[~, ~, ~, all_ratioN1s_p(:, 4)]  = fdr_bh(all_ratioN1s_p(:, 3), 0.05, 'pdep');
for iP = 1:size(all_varlat_p, 1)
    out{all_varlat_p(iP, 1), all_varlat_p(iP, 2)}.varlat_p_fdr = all_varlat_p(iP, 4);
    out{all_meanvarlat_p(iP, 1), all_meanvarlat_p(iP, 2)}.meanvarlat_p_fdr = all_meanvarlat_p(iP, 4);
    out{all_ratioN1s_p(iP, 1), all_ratioN1s_p(iP, 2)}.ratioN1s_p_fdr = all_ratioN1s_p(iP, 4);
end



%% 
%  Generate supplementary figure 2 that displays the ratio of #N1s per #channels for each of the connections between the end-point areas

figure('position', [0 0 1200 600])
for iRow = 1:size(conn_matrix, 1)
    for iCol = 1:size(conn_matrix, 2)
        outInd = (iRow - 1) * size(conn_matrix, 2) + iCol;
        subsValues = out{iRow, iCol}.subsValues;
        
        % Plot age vs ratio of N1s
        subplot(size(conn_matrix, 1), size(conn_matrix, 2), outInd);    hold on;
        plot(subsValues(~isnan(subsValues(:, 2)), 1), subsValues(~isnan(subsValues(:, 2)), 5), 'k.', 'MarkerSize', 10);
        
        %
        title(strrep(out{iRow, iCol}.name, '_', '\_'));
        xlim([0 60]); ylim([0 1]);

        %
        text(40, 0.9, ['\rho=', num2str(round(out{iRow, iCol}.ratioN1s_r, 2))]);
        text(40, 0.8, ['P_f_d_r=', num2str(round(out{iRow, iCol}.ratioN1s_p_fdr, 2))]);

        hold off;
    end
end

%{
figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS2_AgeVsRatioN1s');
set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)
%}


%% 
%  Generate supplementary figure 3 that displays the variance in latencies for each of the connections between the end-point areas

figure('position', [0 0 1200 600])
for iRow = 1:size(conn_matrix, 1)
    for iCol = 1:size(conn_matrix, 2)
        outInd = (iRow - 1) * size(conn_matrix, 2) + iCol;
        subsValues = out{iRow, iCol}.subsValues;

        % age vs variance latencies
        subplot(size(conn_matrix, 1), size(conn_matrix, 2), outInd);    hold on;
        plot(subsValues(:, 1),1000 * subsValues(:, 3), 'k.', 'MarkerSize', 10);
        
        %
        title(strrep(out{iRow, iCol}.name, '_', '\_'));
        xlim([0 60]); ylim([0 1]);
        
        %
        text(40, 0.9, ['\rho=', num2str(round(out{iRow, iCol}.varlat_r, 2))]);
        text(40, 0.8, ['P_f_d_r=', num2str(round(out{iRow, iCol}.varlat_p_fdr, 2))]);

        hold off;
    end
end


figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS3_AgeVsLatVar_N1');
set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)



%% 
%  Generate supplementary figure 4 that displays the mean and variance in latencies for each of the connections between the end-point areas

figure('position', [0 0 1200 600])
for iRow = 1:size(conn_matrix, 1)
    for iCol = 1:size(conn_matrix, 2)
        outInd = (iRow - 1) * size(conn_matrix, 2) + iCol;
        subsValues = out{iRow, iCol}.subsValues;
        

        % mean latency (* 1000) vs latency variance (base values already multiplied before by 1000)
        subplot(size(conn_matrix, 1), size(conn_matrix, 2), outInd);    hold on;
        plot(subsValues(:, 2) * 1000, subsValues(:, 4), 'k.', 'MarkerSize', 10);
        
        %
        title(strrep(out{iRow, iCol}.name, '_', '\_'));
        xlim([0 100])
        ylim([0 (max(subsValues(:, 4)))])
        
        %
        text(75, max(subsValues(:, 4) * .95), ['\rho=', num2str(round(out{iRow, iCol}.meanvarlat_r, 2))]);
        if out{iRow, iCol}.meanvarlat_p_fdr < .001
            text(75, max(subsValues(:, 4) * .85), 'P_f_d_r < 0.001');
        elseif out{iRow, iCol}.meanvarlat_p_fdr < .01
            text(75, max(subsValues(:, 4) * .85), 'P_f_d_r < 0.01');
        else
            text(75, max(subsValues(:, 4) * .85), ['P_f_d_r=', num2str(round(out{iRow, iCol}.meanvarlat_p_fdr, 2))]);
        end
        
        %
        [P, S] = polyfit(subsValues(~isnan(subsValues(:, 2)), 2) * 1000, ...
                         subsValues(~isnan(subsValues(:, 2)), 4), 1);
        x_mean = 1:1:100;
        y_fit = P(1) * x_mean + P(2);
        if out{iRow, iCol}.meanvarlat_p_fdr < 0.05
            plot(x_mean, y_fit, 'r')
            %plot(100, 0, 'r*');
        end

        hold off;
    end
end
set(gca, 'XTick', 20:20:100)

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS4_MeanVsVariance_N1');
set(gcf,'PaperPositionMode', 'auto');
print('-dpng', '-r300', figureName);
print('-depsc', '-r300', figureName);
