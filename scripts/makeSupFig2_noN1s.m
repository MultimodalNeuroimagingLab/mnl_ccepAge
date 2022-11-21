%
% This script produces supplementary figures 2. It also checks how common
% N1s are detected on connections that are not governed by major fiber
% pathways according to Yeh et all., 2018
%
% Dora Hermes, Dorien van Blooijs, Max van den Boom, 2022

clear
close all
clc

%% 
%  Load the ccepData and ccepAverages from the derivatives

myDataPath = setLocalDataPath(1);
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

stimStimElec_excludeDist = 18;     % the distance between the stimulated electrodes (in mm) above which N1s are excluded, 0 = not excluding
respStimElec_excludeDist = 13;     % the distance between a stimulated and response electrode (in mm) within which N1s are excluded, 0 = not excluding

% load tracts and their corresponding end-point ROIs
rois = ccep_categorizeAnatomicalRegions();


%% 
%  Define the connection to be displayed and prepare the data
conn_matrix = {[2 1 0], [3 1 0], [3 2 0], [1 1 0]; ...
               [2 1 1], [3 1 1], [3 2 1], [1 1 1]; ...
               };

% get N1s for tested connections:
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
        
        
        %
        % latencies, number of N1s and ratio of N1s
        %
        
        % extract the latencies, number of N1s and ratio of N1s between the end-point ROIs for a specific (sub-)tract and direction
        metrics = ccep_N1sBetweenRegions( ccepData, ...
                                          rois(iTr).sub_tract(iSubTr).(['roi', num2str(iDir + 1)]), ...
                                          rois(iTr).sub_tract(iSubTr).(['roi', num2str(~iDir + 1)]), ...
                                          stimStimElec_excludeDist, respStimElec_excludeDist);
        
        % retrieve metrics per subject, output format:
        % <subject> x <age, mean in latencies, variance in latencies, variance in (latencies * 1000), relative number of N1s>
        subjectsN1Values = NaN(length(metrics), 4);
        for iSubj = 1:length(metrics)
            subjectsN1Values(iSubj, 1) = metrics(iSubj).age;
            subjectsN1Values(iSubj, 2) = mean(metrics(iSubj).latencies, 'omitnan');
            subjectsN1Values(iSubj, 3) = var(metrics(iSubj).latencies, 'omitnan');
            subjectsN1Values(iSubj, 4) = var(1000 * metrics(iSubj).latencies, 'omitnan');
            subjectsN1Values(iSubj, 5) = metrics(iSubj).numElecRespROI;
            subjectsN1Values(iSubj, 6) = metrics(iSubj).numElecStimROI;
            subjectsN1Values(iSubj, 7) = length(metrics(iSubj).latencies);
        end

        out{iRow, iCol}.subjectsN1Values = subjectsN1Values;
        
    end
end

%% 
%%  Generate supplementary figure 2 that displays the ratio of #N1s per #channels for each of the connections between the end-point areas
%%

p_all = [];
r_all = [];
n_all = [];

figure('position', [0 0 1200 600])
for iRow = 1:size(conn_matrix, 1)
    for iCol = 1:size(conn_matrix, 2)
        outInd = (iRow - 1) * size(conn_matrix, 2) + iCol;
        subjectsN1Values = out{iRow, iCol}.subjectsN1Values;
        
        % Plot age vs ratio of N1s
        subplot(size(conn_matrix, 1), size(conn_matrix, 2), outInd);    hold on;
        % plot(subjectsN1Values(~isnan(subjectsN1Values(:, 2)), 1), subjectsN1Values(~isnan(subjectsN1Values(:, 2)), 5), 'k.', 'MarkerSize', 10);
        
        % total possible N1s:
        total_possible_N1 = subjectsN1Values(:, 5) .* subjectsN1Values(:, 6);

        % number of N1s:
        number_N1 = subjectsN1Values(:, 6);
        
        plot(subjectsN1Values(total_possible_N1>0, 1),number_N1(total_possible_N1>0)./total_possible_N1(total_possible_N1>0), 'k.', 'MarkerSize', 10);
        
        [r,p] = corr(subjectsN1Values(total_possible_N1>0, 1),number_N1(total_possible_N1>0)./total_possible_N1(total_possible_N1>0),'Type','Spearman');
        
        r_all = [r_all r];
        p_all = [p_all p];
        n_all = [n_all outInd];
        
        title(strrep(out{iRow, iCol}.name, '_', '\_'));
        
        hold off;
    end
end

[~,~,~,p_all_fdr]  = fdr_bh(p_all, 0.05, 'pdep');

for iRow = 1:size(conn_matrix, 1)
    for iCol = 1:size(conn_matrix, 2)
        outInd = (iRow - 1) * size(conn_matrix, 2) + iCol;
        subplot(size(conn_matrix, 1), size(conn_matrix, 2), outInd);    hold on;
        xlim([0 60]); ylim([0 1]);
        if iRow == size(conn_matrix, 1),    xlabel('age'); end
        if iCol == 1,                       ylabel('ratio N1s'); end
        
        %
        text(40, 0.9, ['\rho=', num2str(r_all(outInd), 2) ]);
        text(40, 0.8, ['P_f_d_r=', num2str(p_all_fdr(outInd), 2)]);
    end
end


figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS2_AgeVsRatioN1s');
set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)














%% get N1s for non-present connection 
% Per Yeh et al., 2018) these connections should not be present
% Precentral <-> fusiform /// 29 G_precentral <-> 21 G_oc-temp_lat-fusifor
% Postcentral <-> inferior temporal gyrus /// 28 G_postcentral <-> 37 G_temporal_inf
% Supramarginal <-> fusiform /// 26 G_pariet_inf-Supramar <-> 21 G_oc-temp_lat-fusifor

nonrois(1).roi1 = [29];
nonrois(1).roi2 = [21];
nonrois(1).name = 'precentral-fusiform';
nonrois(2).roi1 = [28];
nonrois(2).roi2 = [37];
nonrois(2).name = 'postcentral-inferiortemporal';
nonrois(3).roi1 = [26];
nonrois(3).roi2 = [21];
nonrois(3).name = 'supramarginal-fusiform';

dir_rois = [0 1];

nonout = cell(length(nonrois),2);
for iRow = 1:length(nonrois)
    for iCol = 1:2
        
        iDir = dir_rois(iCol);

        % extract the connection name (includes tract-subtract and direction)
        subDir = split(nonrois(iRow).name, '-');
        strName = [subDir{iDir + 1}, ' -> ', subDir{~iDir + 1}];
        
        %
        disp(['Running - ', strName]);
        nonout{iRow, iCol}.name = strName;
        
        %
        % latencies, number of N1s and ratio of N1s
        %
        
        % extract the latencies, number of N1s and ratio of N1s between the end-point ROIs for a specific (sub-)tract and direction
        metrics = ccep_N1sBetweenRegions( ccepData, ...
                                          nonrois(iRow).(['roi', num2str(iDir + 1)]), ...
                                          nonrois(iRow).(['roi', num2str(~iDir + 1)]), ...
                                          stimStimElec_excludeDist, respStimElec_excludeDist);
        
        % retrieve metrics per subject, output format:
        % <subject> x <age, mean in latencies, variance in latencies, variance in (latencies * 1000), relative number of N1s>
        subjectsN1Values = NaN(length(metrics), 4);
        for iSubj = 1:length(metrics)
            subjectsN1Values(iSubj, 1) = metrics(iSubj).age;
            subjectsN1Values(iSubj, 2) = mean(metrics(iSubj).latencies, 'omitnan');
            subjectsN1Values(iSubj, 3) = var(metrics(iSubj).latencies, 'omitnan');
            subjectsN1Values(iSubj, 4) = var(1000 * metrics(iSubj).latencies, 'omitnan');
            subjectsN1Values(iSubj, 5) = metrics(iSubj).numElecRespROI;
            subjectsN1Values(iSubj, 6) = metrics(iSubj).numElecStimROI;
            subjectsN1Values(iSubj, 7) = length(metrics(iSubj).latencies);
        end

        nonout{iRow, iCol}.subjectsN1Values = subjectsN1Values;
        
    end
end
clear subjectsN1Values


%% distribution plot of the percentage of N1s in a connection

p_all = NaN(size(conn_matrix, 1),size(conn_matrix, 2));

figure('position', [0 0 200 300])
for iRow = 1:size(conn_matrix, 1)
    for iCol = 1:size(conn_matrix, 2)
        % outInd = (iRow - 1) * size(conn_matrix, 2) + iCol;
        subjectsN1Values = out{iRow, iCol}.subjectsN1Values;
        
        % Plot age vs ratio of N1s
        subplot(size(conn_matrix, 1), 1, iRow);    hold on;
        
        % total possible N1s:
        total_possible_N1 = subjectsN1Values(:, 5) .* subjectsN1Values(:, 6);

        % number of N1s:
        number_N1 = subjectsN1Values(:, 7);
        
        % boxplot(number_N1(total_possible_N1>0)./total_possible_N1(total_possible_N1>0));
        distributionPlot(number_N1(total_possible_N1>0)./total_possible_N1(total_possible_N1>0),...
            'xValues',iCol,'histOpt',2,'addSpread',1)

        [h,p,ci,stats] = ttest(number_N1(total_possible_N1>0)./total_possible_N1(total_possible_N1>0));
        
        p_all(iRow,iCol) = p;

        hold off;
    end
end

[~,~,~,p_all_fdr]  = fdr_bh(p_all, 0.05, 'pdep');
set(gca,'XTick',[1:4])

for iRow = 1:size(conn_matrix, 1)
    for iCol = 1:size(conn_matrix, 2)
        subplot(size(conn_matrix, 1), 1, iRow);    hold on;
        if p_all_fdr(iRow,iCol)<0.05
            plot(iCol,-0.1,'r*')
        end
        ylim([-0.15 1.02])
        xlim([0.5 size(conn_matrix, 2)+.5])
    end
end

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS2b_RatioN1sConn');
set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)



%% distribution plot of the percentage of N1s in a connection

p_all = NaN(length(nonrois),2);

figure('position', [0 0 160 300])
for iRow = 1:length(nonrois)
    for iCol = 1:2
        subjectsN1Values = nonout{iRow, iCol}.subjectsN1Values;
        
        % Plot age vs ratio of N1s
%         subplot(2,length(nonrois), iCol);    hold on;
        subplot(2,1, iCol);    hold on;

        % total possible N1s:
        total_possible_N1 = subjectsN1Values(:, 5) .* subjectsN1Values(:, 6);

        % number of N1s:
        number_N1 = subjectsN1Values(:, 7);
        
%         boxplot(number_N1(total_possible_N1>0)./total_possible_N1(total_possible_N1>0));
        distributionPlot(number_N1(total_possible_N1>0)./total_possible_N1(total_possible_N1>0),...
            'xValues',iRow,'histOpt',2,'addSpread',1)

        [h,p,ci,stats] = ttest(number_N1(total_possible_N1>0)./total_possible_N1(total_possible_N1>0));
        
        p_all(iRow,iCol) = p;
    end
end

[~,~,~,p_all_fdr]  = fdr_bh(p_all, 0.05, 'pdep');

for iRow = 1:length(nonrois)
    for iCol = 1:2        
        subplot(2,1, iCol);    hold on;
        if p_all_fdr(iRow,iCol)<0.05
            plot(iCol,-0.1,'r*')
        end
        xlim([0.5 3.5])
        ylim([-0.15 1.02])
    end
end

set(gca,'XTick',[1:3])

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS2c_RatioN1sNoconn');
set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)


