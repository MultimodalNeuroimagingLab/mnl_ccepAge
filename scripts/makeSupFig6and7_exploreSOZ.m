%
% Script for supplementary figure 6 and 7 - investigates the influence of seizure onset zones (SOZ) on latencies
%
% Max van den Boom. Dorien van Blooijs, Dora Hermes. 2022
%


clear
close all
clc

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'file')
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData')
else
    disp('Run scripts ccep02_aggregateToStruct.m and ccep03_addtracts.m first')
end

%%
% Gather all relevant ROI Destrieux codes

allROICodes = [];
rois = ccep_categorizeAnatomicalRegions();
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract) 
        allROICodes = [allROICodes, [rois(iTr).sub_tract(iSubTr).roi1, rois(iTr).sub_tract(iSubTr).roi2]];
    end
end
allROICodes = unique(allROICodes);


%%
%  Select patients in whom SOZ and/or RA are delineated
count = 1;
n1Latencies = struct;

for iSubj = 1:size(ccepData, 2)
    
    % include electrodes that are ECoG and on top of one of our ROIs
    % TODO: checking for 'strips' and 'grids' is not ideal, but electrodes.tsv does not have type (like channels does)
    %includedChannelIndices = ismember(lower(ccepData(iSubj).electrodes.group), {'grid', 'strip'});
    includedChannelIndices =    ismember(lower(ccepData(iSubj).electrodes.group), {'grid', 'strip'}) & ...
                                ismember(ccepData(iSubj).electrodes.Destrieux_label, allROICodes);
    
    
    
    if any(contains(ccepData(iSubj).electrodes.soz, 'yes')) || any(contains(ccepData(iSubj).electrodes.resected, 'yes'))
    
        n1Latencies(count).id         = ccepData(iSubj).id;
        n1Latencies(count).ses        = ccepData(iSubj).ses;
        n1Latencies(count).age        = ccepData(iSubj).age;
        n1Latencies(count).run        = ccepData(iSubj).run;
        
        % extract seizure onset zone (SOZ) and resected area (RA)
        n1Latencies(count).SOZ        = find(strcmpi(ccepData(iSubj).electrodes.soz, 'yes') == 1 & includedChannelIndices);
        n1Latencies(count).SOZ_names  = ccepData(iSubj).electrodes.name(n1Latencies(count).SOZ);
        n1Latencies(count).RA         = find(strcmpi(ccepData(iSubj).electrodes.resected, 'yes') == 1 & includedChannelIndices);
        n1Latencies(count).RA_names   = ccepData(iSubj).electrodes.name(n1Latencies(count).RA);
        
        % included channels (in one of the destrieux ROIs & ECoG)
        n1Latencies(count).included_channels = ccepData(iSubj).electrodes.name(includedChannelIndices);
        
        count = count + 1;
    end
end



%% 
%  Distinguish latencies SOZ and non-SOZ
clc
close all

% loop over the subjects
for iSubj = 1:size(n1Latencies, 2)
    
    % check if there are electrodes on the SOZs for this participant
    if ~isnan(n1Latencies(iSubj).SOZ)
        
        % loop over the runs
        for iRun = 1:size(n1Latencies(iSubj).run, 2)

            % retrieve the names of the electrodes that are on SOZs given the electrode tsv indexing
            SOZ_names = n1Latencies(iSubj).SOZ_names;
            
            % generate a variable that holds the indexes of the SOZ channels given the run's channel tsv indexing, and
            % a variable that holds the indexes of the stim-pairs of which either electrodes is on a SOZ (given the run's channel tsv indexing)
            respSOZ_channel_idxs = [];
            %respSOZ_channel_names = {};
            stimSOZ_stimpair_idxs = [];
            %stimSOZ_stimpair_names = {};
            for iName = 1:length(SOZ_names)
                
                % response channel (should match exactly)
                respIdx = find(ismember(lower(n1Latencies(iSubj).run(iRun).channel_names), lower(SOZ_names{iName})));
                if isempty(respIdx) || length(respIdx) ~= 1
                   error('Could not find electrode name in this run''s channel tsv'); 
                end
                respSOZ_channel_idxs(end + 1) = respIdx;
                %respSOZ_channel_names(end + 1) = n1LatenciesSOZ(iSubj).run(iRun).channel_names(respIdx);
                
                
                % stim-pairs (could be multiple)
                % TODO: counting switched polaries as seperate
                for iStimPair = 1:length(n1Latencies(iSubj).run(iRun).stimpair_names)
                    stimElecs = split(n1Latencies(iSubj).run(iRun).stimpair_names{iStimPair}, '-');
                    if strcmpi(SOZ_names{iName}, stimElecs{1}) || strcmpi(SOZ_names{iName}, stimElecs{2})
                        stimSOZ_stimpair_idxs(end + 1) = iStimPair;
                        %stimSOZ_stimpair_names(end + 1) = n1LatenciesSOZ(iSubj).run(iRun).stimpair_names(iStimPair);
                    end
                end
                
                clear respIdx stimIdxs stimElecs;
            end

            
            %
            % response electrodes
            %
            
            % retrieve the indices (and names) of the channels that were not on the SOZ (and good)
            nonSOZ_respChannel_indices = setdiff(n1Latencies(iSubj).run(iRun).good_channels, respSOZ_channel_idxs);
            nonSOZ_respChannel_names = n1Latencies(iSubj).run(iRun).channel_names(nonSOZ_respChannel_indices);
            
            % remove the response non-SOZ channels that are not included (not ECoG or on one of our ROIs)
            excludedRespChannels = ~ismember(lower(nonSOZ_respChannel_names), lower(n1Latencies(iSubj).included_channels));
            nonSOZ_respChannel_indices(excludedRespChannels) = [];
            nonSOZ_respChannel_names(excludedRespChannels) = [];
            
            % retrieve the N1 latencies (in ms) of the response electrodes that were on SOZs (and those which were not on SOZs)
            respSOZlatencies        = n1Latencies(iSubj).run(iRun).n1_peak_sample(respSOZ_channel_idxs, :);
            respSOZlatencies        = n1Latencies(iSubj).run(iRun).tt(respSOZlatencies(~isnan(respSOZlatencies)));
            respNonSOZlatencies     = n1Latencies(iSubj).run(iRun).n1_peak_sample(nonSOZ_respChannel_indices, :);
            respNonSOZlatencies     = n1Latencies(iSubj).run(iRun).tt(respNonSOZlatencies(~isnan(respNonSOZlatencies)));
            
            % store for run
            runRespNonSOZlat{iRun} = respNonSOZlatencies;
            runRespSOZlat{iRun} = respSOZlatencies;

            
            %
            % stimulated pairs
            %
            
            % retrieve the indices of the stim-pair of which none of the electrodes were are on the SOZ
            nonSOZ_stimpair_indices = setdiff(1:size(n1Latencies(iSubj).run(iRun).n1_peak_sample, 2), stimSOZ_stimpair_idxs);
            nonSOZ_stimpair_names = n1Latencies(iSubj).run(iRun).stimpair_names(nonSOZ_stimpair_indices);
            
            % remove the non-SOZ stim-pairs that have channels that are not included (not ECoG or on one of our ROIs)
            excludedStimpairs = [];
            for iStimPair = 1:length(nonSOZ_stimpair_indices)
                stimElecs = split(nonSOZ_stimpair_names{iStimPair}, '-');
                
                % check if either of the stim-pair electrodes should not be taken into account (because not ECoG or on one of our ROIs)
                if ~ismember(lower(stimElecs{1}), lower(n1Latencies(iSubj).included_channels)) || ~ismember(lower(stimElecs{2}), lower(n1Latencies(iSubj).included_channels))
                    excludedStimpairs(end + 1) = iStimPair;
                end
            end
            nonSOZ_stimpair_indices(excludedStimpairs) = [];
            nonSOZ_stimpair_names(excludedStimpairs) = [];
            
            % retrieve the N1 latencies (in ms) of the stim-pairs that were on SOZs (and those which were not on SOZs)
            stimSOZlatencies    = n1Latencies(iSubj).run(iRun).n1_peak_sample(:, stimSOZ_stimpair_idxs);
            stimSOZlatencies    = n1Latencies(iSubj).run(iRun).tt(stimSOZlatencies(~isnan(stimSOZlatencies)));
            stimNonSOZlatencies = n1Latencies(iSubj).run(iRun).n1_peak_sample(:, nonSOZ_stimpair_indices);
            stimNonSOZlatencies = n1Latencies(iSubj).run(iRun).tt(stimNonSOZlatencies(~isnan(stimNonSOZlatencies)));

            % store for run
            runStimNonSOZLatencies{iRun} = stimNonSOZlatencies;
            runStimSOZLatencies{iRun}    = stimSOZlatencies;

            clear respSOZlatencies respNonSOZlatencies stimSOZlatencies stimNonSOZlatencies
        end
        
        % concatenate response latencies over runs
        n1Latencies(iSubj).respSOZlatencies    = horzcat(runRespSOZlat{:});
        n1Latencies(iSubj).respNonSOZlatencies = horzcat(runRespNonSOZlat{:});

        if (any(~isnan(n1Latencies(iSubj).respSOZlatencies )) && any(~isnan(n1Latencies(iSubj).respNonSOZlatencies ))) || (~isempty(n1Latencies(iSubj).respSOZlatencies) && ~isempty(n1Latencies(iSubj).respNonSOZlatencies))
            
            % Wilcoxon rank sum/mann witney u test for responding electrode in either SOZ or non-SOZ
            n1Latencies(iSubj).respSOZ_p = ranksum(n1Latencies(iSubj).respNonSOZlatencies, n1Latencies(iSubj).respSOZlatencies);

            fprintf('-- %s: When comparing latency in response SOZ and response non-SOZ: p = %1.3f with median resp_SOZ = %1.3f sec and median resp_non-SOZ = %1.3f sec\n',...
                n1Latencies(iSubj).id, n1Latencies(iSubj).respSOZ_p, median([n1Latencies(iSubj).respSOZlatencies]),median([n1Latencies(iSubj).respNonSOZlatencies]))
        end

        % concatenate stim-pair latencies over runs
        n1Latencies(iSubj).stimSOZlatencies = horzcat(runStimSOZLatencies{:});
        n1Latencies(iSubj).stimNonSOZlatencies = horzcat(runStimNonSOZLatencies{:});
        
        
        if (any(~isnan(n1Latencies(iSubj).stimSOZlatencies )) && any(~isnan(n1Latencies(iSubj).stimNonSOZlatencies ))) || (~isempty(n1Latencies(iSubj).stimSOZlatencies) && ~isempty(n1Latencies(iSubj).stimNonSOZlatencies))

            % Wilcoxon rank sum/mann witney u test for stimulated electrode in either SOZ or non-SOZ
            n1Latencies(iSubj).stimSOZ_p = ranksum(n1Latencies(iSubj).stimNonSOZlatencies, n1Latencies(iSubj).stimSOZlatencies);

            fprintf('-- %s: When comparing latency in stimulated SOZ and stimulated nSOZ: p = %1.3f with median stim_SOZ = %1.3f sec and median stim_non-SOZ = %1.3f sec\n',...
                n1Latencies(iSubj).id, n1Latencies(iSubj).stimSOZ_p, median([n1Latencies(iSubj).stimSOZlatencies]),median([n1Latencies(iSubj).stimNonSOZlatencies]))
        end
        
        clear respnSOZlat respSOZlat stimnSOZlat stimSOZlat
    end
end

% calculate the FDR p-values
p_indices_resp = find(~cellfun(@isempty, {n1Latencies.respSOZ_p}));
p_indices_stim = find(~cellfun(@isempty, {n1Latencies.stimSOZ_p}));
[~, ~, ~, p_vals_resp_FDR]  = fdr_bh([n1Latencies(p_indices_resp).respSOZ_p], 0.05, 'pdep');
[~, ~, ~, p_vals_stim_FDR]  = fdr_bh([n1Latencies(p_indices_stim).stimSOZ_p], 0.05, 'pdep');

% put FDR p-values back into the table
for respIdx = 1:length(p_indices_resp)
    n1Latencies(p_indices_resp(respIdx)).respSOZ_pFDR = p_vals_resp_FDR(respIdx);
end
for stimIdx = 1:length(p_indices_stim)
    n1Latencies(p_indices_stim(stimIdx)).stimSOZ_pFDR = p_vals_stim_FDR(stimIdx);
end

% replace all empty p-values with nans
a = find(cellfun(@isempty, {n1Latencies.respSOZ_p}));
for aCount = 1:length(a),   n1Latencies(a(aCount)).respSOZ_p = nan;     end
a = find(cellfun(@isempty, {n1Latencies.stimSOZ_p}));
for aCount = 1:length(a),   n1Latencies(a(aCount)).stimSOZ_p = nan;     end
a = find(cellfun(@isempty, {n1Latencies.respSOZ_pFDR}));
for aCount = 1:length(a),   n1Latencies(a(aCount)).respSOZ_pFDR = nan;     end
a = find(cellfun(@isempty, {n1Latencies.stimSOZ_pFDR}));
for aCount = 1:length(a),   n1Latencies(a(aCount)).stimSOZ_pFDR = nan;     end


%%
%
for iSubj = 1:length(n1Latencies)
    n1Latencies(iSubj).meanRespSOZLatency = mean(n1Latencies(iSubj).respSOZlatencies);
    n1Latencies(iSubj).meanRespNonSOZLatency = mean(n1Latencies(iSubj).respNonSOZlatencies);
    
    n1Latencies(iSubj).meanStimSOZLatency = mean(n1Latencies(iSubj).stimSOZlatencies);
    n1Latencies(iSubj).meanStimNonSOZLatency = mean(n1Latencies(iSubj).stimNonSOZlatencies);
end


% Determine the number of subject that have significant differences between the SOZ and non-SOZ
respSignIncr = find( ([n1Latencies.meanRespSOZLatency] > [n1Latencies.meanRespNonSOZLatency]) & ...
                     (~isnan([n1Latencies.respSOZ_pFDR]) & [n1Latencies.respSOZ_pFDR] < .05));
respSignDecr = find( ([n1Latencies.meanRespSOZLatency] < [n1Latencies.meanRespNonSOZLatency]) & ...
                     (~isnan([n1Latencies.respSOZ_pFDR]) & [n1Latencies.respSOZ_pFDR] < .05));
stimSignIncr = find( ([n1Latencies.meanStimSOZLatency] > [n1Latencies.meanStimNonSOZLatency]) & ...
                     (~isnan([n1Latencies.stimSOZ_pFDR]) & [n1Latencies.stimSOZ_pFDR] < .05));
stimSignDecr = find( ([n1Latencies.meanStimSOZLatency] < [n1Latencies.meanStimNonSOZLatency]) & ...
                     (~isnan([n1Latencies.stimSOZ_pFDR]) & [n1Latencies.stimSOZ_pFDR] < .05));
disp(['- ', num2str(length(respSignIncr)), ' subjects have a significantly higher latency from response electrodes on the SOZ than on non-SOZ']);
disp(['- ', num2str(length(respSignDecr)), ' subjects have a significantly lower latency from response electrodes on the SOZ than on non-SOZ']);
disp(['- ', num2str(length(stimSignIncr)), ' subjects have a significantly higher latency from stim-pairs on the SOZ than on non-SOZ']);
disp(['- ', num2str(length(stimSignDecr)), ' subjects have a significantly lower latency from stim-pairs on the SOZ than on non-SOZ']);



%%
%  Generate supplement 6 - Violin plot responses in or outside SOZ

respall = cell(1);
SEs = cell(1);
names = cell(1);

count = 1;
for iSubj = 1:size(n1Latencies, 2)
    
    % check whether there are N1 latencies in electrodes that were on SOZs and electrodes that were not on SOZs
    if ~isempty(n1Latencies(iSubj).respSOZlatencies) && ~isempty(n1Latencies(iSubj).respNonSOZlatencies)
        
        % add N1 latencies for electrodes that were on SOZs to be plotted
        respall{count} = n1Latencies(iSubj).respSOZlatencies * 1000;
        SEs{count} = std(respall{count}) / sqrt(length(respall{count}));
        names{count} = cellstr(repmat([n1Latencies(iSubj).id ' soz'],size(respall{count}, 2), 1));
        count = count + 1;
        
        % add N1 latencies for electrodes that were on SOZs to be plotted
        respall{count} = n1Latencies(iSubj).respNonSOZlatencies * 1000;
        SEs{count} = std(respall{count}) / sqrt(length(respall{count}));
        names{count} = cellstr(repmat([n1Latencies(iSubj).id ' nsoz'],size(respall{count}, 2), 1));
        count = count + 1;
        
        if iSubj ~= size(n1Latencies, 2)
            respall{count} = -5 * ones(1, 1);
            names{count} = cellstr(repmat([n1Latencies(iSubj).id ' zempty'], 1, 1));
            count = count + 1;
        end
    
    end
       
end

% open the figure and draw the violins
close all
cmap = colormap('parula');
close(figure(1))
ymax = ceil(max(horzcat(respall{:})));
temp = horzcat(respall{:});
ymin = floor(min(temp(temp > 0)));
h = figure('position',[0 0 2400 1200]);
vs = violinplot(horzcat(respall{:}), vertcat(names{:}), 'ViolinColor', cmap(1, :), 'BoxColor', [0.4 0.4 0.6], 'Width', 0.4, 'BoxWidth', .08, 'ShowData', false);

% for each subject/violin plot pair (non-SOZ violin & SOZ violin)
for iViolin = 1:3:size(vs, 2)
    vs(iViolin).ViolinPlot.FaceColor = cmap(128, :);
    vs(iViolin).ScatterPlot.MarkerFaceColor = cmap(128, :);
    vs(iViolin).ScatterPlot.SizeData = 70;
    vs(iViolin).MedianPlot.SizeData = 70;
    vs(iViolin).ViolinPlot.FaceColor = cmap(128, :);
end

hold on
count = 1;
for iSubj = 1:size(n1Latencies, 2)
    
    % check whether there are N1 latencies in electrodes that were on SOZs and electrodes that were not on SOZs
    if ~isempty(n1Latencies(iSubj).respSOZ_pFDR) && ~isnan(n1Latencies(iSubj).respSOZ_pFDR)
        if n1Latencies(iSubj).respSOZ_pFDR < .05
            text(count + 0.3, 1.03 * ymax, '*', 'FontSize', 24)
            plot([count, count + 1], 1.02 * [ymax ymax], 'k')
        end
        count = count + 3;
    end
    
end
hold off

ylim([ymin 1.15 * ymax])
xlim([0 size(names, 2) + 1])
ax = gca;
ax.FontSize = 20;
ax.XTickLabelRotation = 0;
ax.XTick = 1.5:3:size(vs, 2);
ax.XTickLabel = 1:size([n1Latencies(:).respSOZ_p], 2);
ax.YLabel.String = 'N1 Peak Latency (ms)';
ax.XLabel.String = 'Subjects';
ax.Title.String = 'Latency of measured electrodes in or outside SOZ';
lgd = legend([vs(1).ViolinPlot, vs(2).ViolinPlot], 'Responses measured on non-SOZ', 'Responses measured on SOZ');
lgd.FontSize = 20;
set(gcf,'color', 'w');

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS6_LatencyRespSOZ');
set(gcf, 'renderer', 'Painters')
set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)



%%
%  Generate supplement 7 - violin plot stim-pairs in or outside SOZ

stimall = cell(1);
names = cell(1);

count = 1;
for iSubj = 1:size(n1Latencies, 2)
    
    % check whether there are N1 latencies in electrodes that were on SOZs and electrodes that were not on SOZs
    if ~isempty(n1Latencies(iSubj).stimSOZlatencies) && ~isempty(n1Latencies(iSubj).stimNonSOZlatencies)
        
        stimall{count} = n1Latencies(iSubj).stimSOZlatencies * 1000;
        names{count} = cellstr(repmat([n1Latencies(iSubj).id ' soz'],size(stimall{count}, 2), 1));
        count = count + 1;
        stimall{count} = n1Latencies(iSubj).stimNonSOZlatencies * 1000;
        names{count} = cellstr(repmat([n1Latencies(iSubj).id ' non-soz'],size(stimall{count}, 2), 1));
        count = count + 1;
        
        if iSubj ~= size(n1Latencies, 2)
            stimall{count} = -5 * ones(1, 1);
            names{count} = cellstr(repmat([n1Latencies(iSubj).id ' zempty'], 1, 1));
            count = count + 1;
        end
    
    end
       
end

% open the figure and draw the violins
close all
cmap = colormap('parula');
close(figure(1))
ymax = ceil(max(horzcat(stimall{:})));
temp = horzcat(stimall{:});
ymin = floor(min(temp(temp > 0)));
h = figure('position',[0 0 2000 1200]);
vs = violinplot(horzcat(stimall{:}), vertcat(names{:}), 'ViolinColor', cmap(1, :), 'BoxColor', [0.4 0.4 0.6], 'Width', 0.4, 'BoxWidth', .08, 'ShowData', false);

% for each subject/violin plot pair (non-SOZ violin & SOZ violin)
for iViolin = 1:3:size(vs, 2)
    vs(iViolin).ViolinPlot.FaceColor = cmap(128, :);
    vs(iViolin).ScatterPlot.MarkerFaceColor = cmap(128, :);
    vs(iViolin).ScatterPlot.SizeData = 70;
    
end

hold on
count = 1;
for iSubj = 1:size(n1Latencies, 2)
    
    % check whether there are N1 latencies in electrodes that were on SOZs and electrodes that were not on SOZs
    if ~isempty(n1Latencies(iSubj).stimSOZ_pFDR) && ~isnan(n1Latencies(iSubj).stimSOZ_pFDR)
        if n1Latencies(iSubj).stimSOZ_pFDR < .05
            text(count + 0.3, 1.03 * ymax, '*', 'FontSize', 24)
            plot([count, count + 1], 1.02 * [ymax ymax], 'k')
        end
        count = count + 3;
    end
    
end
hold off

ylim([ymin 1.15 * ymax])
xlim([0 size(names, 2) + 1])
ax = gca;
ax.FontSize = 20;
ax.XTickLabelRotation = 0;
ax.XTick = 1.5:3:size(vs, 2);
ax.XTickLabel = 1:size([n1Latencies(:).stimSOZ_p], 2);
ax.YLabel.String = 'N1 Peak Latency (ms)';
ax.XLabel.String = 'Subjects';
ax.Title.String = 'Latency of stim-pairs when stimulating in or outside SOZ';
lgd = legend([vs(1).ViolinPlot, vs(2).ViolinPlot], 'Responses when not stimulating SOZ', 'Responses when stimulating SOZ');
lgd.FontSize = 20;
set(gcf,'color', 'w');

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS7_LatencyStimSOZ');
set(gcf, 'renderer', 'Painters')
set(gcf, 'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)
