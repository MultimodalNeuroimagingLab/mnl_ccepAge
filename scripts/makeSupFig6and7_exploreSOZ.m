%
% This scripts investigates the influence of epilepsy on latencies
%
% Dorien van Blooijs, Dora Hermes, 2020
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
%  Select patients in whom SOZ and/or RA are delineated
count = 1;
n1LatenciesSOZ = struct;

for iSubj = 1:size(ccepData, 2)
    
    if any(contains(ccepData(iSubj).electrodes.soz, 'yes')) || any(contains(ccepData(iSubj).electrodes.resected, 'yes'))
    
        n1LatenciesSOZ(count).id         = ccepData(iSubj).id;
        n1LatenciesSOZ(count).ses        = ccepData(iSubj).ses;
        n1LatenciesSOZ(count).age        = ccepData(iSubj).age;
        n1LatenciesSOZ(count).run        = ccepData(iSubj).run;
        
        % extract seizure onset zone (SOZ) and resected area (RA)
        n1LatenciesSOZ(count).SOZ       = find(strcmp(ccepData(iSubj).electrodes.soz, 'yes') == 1);
        n1LatenciesSOZ(count).RA        = find(strcmp(ccepData(iSubj).electrodes.resected, 'yes') == 1);
        
        count = count + 1;
        
    end
end


%%
%  determine percentage SOZ electrodes of total electrodes

numSOZ = NaN(size(n1LatenciesSOZ));
numgrid = NaN(size(n1LatenciesSOZ));

for iSubj = 1:size(n1LatenciesSOZ, 2)
    
    if ~isempty(n1LatenciesSOZ(iSubj).SOZ)
        numSOZ(iSubj) = numel(n1LatenciesSOZ(iSubj).SOZ);
    end
   numgrid(iSubj) = numel(n1LatenciesSOZ(iSubj).run(1).good_channels);
    
end

median(numSOZ ./ numgrid, 'omitnan')
min(numSOZ ./ numgrid)
max(numSOZ ./ numgrid)


%% 
%  Distinguish latencies SOZ and nSOZ
clc
close all

% when SOZ is response electrode
for iSubj = 1:size(n1LatenciesSOZ, 2)
    if ~isnan(n1LatenciesSOZ(iSubj).SOZ)
        clear respnSOZlat respSOZlat stimnSOZlat stimSOZlat
        
        if ~isempty(n1LatenciesSOZ(iSubj).run)
            
            for ll = 1:size(n1LatenciesSOZ(iSubj).run, 2)
                clear respSOZlatencies respnSOZlatencies stimSOZlatencies stimnSOZlatencies
                                
                % when SOZ is response electrode
                SOZ = n1LatenciesSOZ(iSubj).SOZ;
                nSOZ = setdiff(n1LatenciesSOZ(iSubj).run(ll).good_channels,SOZ);
                respSOZlatencies  = n1LatenciesSOZ(iSubj).run(ll).n1_peak_sample(SOZ, :);
                respSOZlatencies  = n1LatenciesSOZ(iSubj).run(ll).tt(respSOZlatencies(~isnan(respSOZlatencies)));
                respnSOZlatencies = n1LatenciesSOZ(iSubj).run(ll).n1_peak_sample(nSOZ, :);
                respnSOZlatencies = n1LatenciesSOZ(iSubj).run(ll).tt(respnSOZlatencies(~isnan(respnSOZlatencies)));
                
                respnSOZlat{ll} = respnSOZlatencies;
                respSOZlat{ll} = respSOZlatencies;
                
                % when SOZ is stimulated
                stimSOZ = find(contains(n1LatenciesSOZ(iSubj).run(ll).stimpair_names, n1LatenciesSOZ(iSubj).run(ll).channel_names(SOZ)));
                stimnSOZ = setdiff(1:size(n1LatenciesSOZ(iSubj).run(ll).n1_peak_sample,2),stimSOZ);
                
                stimSOZlatencies  = n1LatenciesSOZ(iSubj).run(ll).n1_peak_sample(:, stimSOZ);
                stimSOZlatencies  = n1LatenciesSOZ(iSubj).run(ll).tt(stimSOZlatencies(~isnan(stimSOZlatencies)));
                stimnSOZlatencies = n1LatenciesSOZ(iSubj).run(ll).n1_peak_sample(:, stimnSOZ);
                stimnSOZlatencies = n1LatenciesSOZ(iSubj).run(ll).tt(stimnSOZlatencies(~isnan(stimnSOZlatencies)));
                
                stimnSOZlat{ll} = stimnSOZlatencies;
                stimSOZlat{ll} = stimSOZlatencies;
                
            end
            
            % mann witney u test for responding electrode in either SOZ or nSOZ
            n1LatenciesSOZ(iSubj).respSOZlatencies = horzcat(respSOZlat{:});
            n1LatenciesSOZ(iSubj).respnSOZlatencies = horzcat(respnSOZlat{:});
            if (any(~isnan(n1LatenciesSOZ(iSubj).respSOZlatencies )) && any(~isnan(n1LatenciesSOZ(iSubj).respnSOZlatencies ))) || (~isempty(n1LatenciesSOZ(iSubj).respSOZlatencies) && ~isempty(n1LatenciesSOZ(iSubj).respnSOZlatencies))
                n1LatenciesSOZ(iSubj).respSOZ_p = ranksum(n1LatenciesSOZ(iSubj).respnSOZlatencies,n1LatenciesSOZ(iSubj).respSOZlatencies);
                
                fprintf('-- %s: When comparing latency in response SOZ and response nSOZ: p = %1.3f with median resp_SOZ = %1.3f sec and median resp_nSOZ = %1.3f sec\n',...
                    n1LatenciesSOZ(iSubj).id, n1LatenciesSOZ(iSubj).respSOZ_p, median([n1LatenciesSOZ(iSubj).respSOZlatencies]),median([n1LatenciesSOZ(iSubj).respnSOZlatencies]))
            end
            
            % mann witney u test for stimulated electrode in either SOZ or nSOZ
            n1LatenciesSOZ(iSubj).stimSOZlatencies = horzcat(stimSOZlat{:});
            n1LatenciesSOZ(iSubj).stimnSOZlatencies = horzcat(stimnSOZlat{:});
            if (any(~isnan(n1LatenciesSOZ(iSubj).stimSOZlatencies )) && ...
                    any(~isnan(n1LatenciesSOZ(iSubj).stimnSOZlatencies ))) || ...
                    (~isempty(n1LatenciesSOZ(iSubj).stimSOZlatencies) && ...
                    ~isempty(n1LatenciesSOZ(iSubj).stimnSOZlatencies))
                
                n1LatenciesSOZ(iSubj).stimSOZ_p = ranksum(n1LatenciesSOZ(iSubj).stimnSOZlatencies,n1LatenciesSOZ(iSubj).stimSOZlatencies);
                
                fprintf('-- %s: When comparing latency in stimulated SOZ and stimulated nSOZ: p = %1.3f with median stim_SOZ = %1.3f sec and median stim_nSOZ = %1.3f sec\n',...
                    n1LatenciesSOZ(iSubj).id, n1LatenciesSOZ(iSubj).stimSOZ_p, median([n1LatenciesSOZ(iSubj).stimSOZlatencies]),median([n1LatenciesSOZ(iSubj).stimnSOZlatencies]))
            end
            
        end
    end
end

% calculate the FDR p-values
p_indices_resp = find(~cellfun(@isempty, {n1LatenciesSOZ.respSOZ_p}));
p_indices_stim = find(~cellfun(@isempty, {n1LatenciesSOZ.stimSOZ_p}));
[~, ~, ~, p_vals_resp_FDR]  = fdr_bh([n1LatenciesSOZ(p_indices_resp).respSOZ_p], 0.05, 'pdep');
[~, ~, ~, p_vals_stim_FDR]  = fdr_bh([n1LatenciesSOZ(p_indices_stim).stimSOZ_p], 0.05, 'pdep');


%{
figure
subplot(1,2,1), hold on, plot(thisVal), plot(p_sort_resp, 'r.'), title('Sig p-val in response SOZ/nSOZ after FDR correction')
subplot(1,2,2), hold on, plot(thisVal), plot(p_sort_stim, 'r.'), title('Sig p-val in stimulated SOZ/nSOZ after FDR correction')
%}

% drop the non-significant values
p_indices_resp(p_vals_resp_FDR >= .05) = [];
p_vals_resp_FDR(p_vals_resp_FDR >= .05) = [];

p_indices_stim(p_vals_stim_FDR >= .05) = [];
p_vals_stim_FDR(p_vals_stim_FDR >= .05) = [];

% assign to table
for idx = 1:length(p_indices_resp)
    n1LatenciesSOZ(p_indices_resp(idx)).respSOZ_pFDR = p_vals_resp_FDR(idx);
end
for idx = 1:length(p_indices_stim)
    n1LatenciesSOZ(p_indices_stim(idx)).stimSOZ_pFDR = p_vals_stim_FDR(idx);
end


%%
%  Violin plot responses in or outside SOZ
respall = cell(1);
names = cell(1);

count = 1;
for iSubj = 1:size(n1LatenciesSOZ, 2)
    
    if ~isempty(n1LatenciesSOZ(iSubj).respSOZlatencies) && ~isempty(n1LatenciesSOZ(iSubj).respnSOZlatencies)
        respall{count} = n1LatenciesSOZ(iSubj).respSOZlatencies * 1000;
        names{count} = cellstr(repmat([n1LatenciesSOZ(iSubj).id ' soz'],size(respall{count}, 2), 1));
        count = count + 1;
        respall{count} = n1LatenciesSOZ(iSubj).respnSOZlatencies * 1000;
        names{count} = cellstr(repmat([n1LatenciesSOZ(iSubj).id ' nsoz'],size(respall{count}, 2), 1));
        count = count + 1;
        
        if iSubj ~= size(n1LatenciesSOZ, 2)
            respall{count} = -5 * ones(1, 1);
            names{count} = cellstr(repmat([n1LatenciesSOZ(iSubj).id ' zempty'], 1, 1));
            count = count+1;
        end
    
    end
       
end

close all
cmap = colormap('parula');
ymax = ceil(max(horzcat(respall{:})));
temp = horzcat(respall{:});
ymin = floor(min(temp(temp > 0)));
h = figure(1);
vs = violinplot(horzcat(respall{:}), vertcat(names{:}), 'ViolinColor', cmap(1, :), 'Width', 0.3);
for iSubj = 1:3:size(vs, 2)
    vs(iSubj).ViolinPlot.FaceColor = cmap(128, :);
    vs(iSubj).ScatterPlot.MarkerFaceColor = cmap(128, :);
end

hold on
count = 1;
for iSubj = 1:size(n1LatenciesSOZ, 2)
    if ~isempty(n1LatenciesSOZ(iSubj).respSOZ_pFDR)
        if n1LatenciesSOZ(iSubj).respSOZ_pFDR == 1 
            text(count + 0.3, 1.03 * ymax, '*')
            plot([count, count + 1], 1.02 * [ymax ymax],'k')
        end
        count = count + 3;
    end
end
hold off

ylim([ymin 1.15 * ymax])
xlim([0 size(names, 2) + 1])
h.Units = 'normalized';
h.Position = [0.1 0.1 0.8 0.8];
ax = gca;
ax.FontSize = 13;
ax.XTickLabelRotation = 0;
ax.XTick = 1.5:3:size(vs, 2);
ax.XTickLabel = 1:size([n1LatenciesSOZ(:).respSOZ_p], 2);
ax.YLabel.String = 'Latency (ms)';
ax.XLabel.String = 'Subjects';
ax.Title.String = 'Latency of responses in or outside SOZ';
legend([vs(1).ViolinPlot,vs(2).ViolinPlot],'outside SOZ','inside SOZ')

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS6_LatencyRespSOZ');
set(gcf,'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)


%%
%  violin plot responses in or outside SOZ
stimall = cell(1);
names = cell(1);

count = 1;
for iSubj = 1:size(n1LatenciesSOZ,2)
    
    if ~isempty(n1LatenciesSOZ(iSubj).stimSOZlatencies) && ~isempty(n1LatenciesSOZ(iSubj).stimnSOZlatencies)
        stimall{count} = n1LatenciesSOZ(iSubj).stimSOZlatencies * 1000;
        names{count} = cellstr(repmat([n1LatenciesSOZ(iSubj).id ' soz'],size(stimall{count}, 2), 1));
        count = count + 1;
        stimall{count} = n1LatenciesSOZ(iSubj).stimnSOZlatencies * 1000;
        names{count} = cellstr(repmat([n1LatenciesSOZ(iSubj).id ' nsoz'],size(stimall{count}, 2), 1));
        count = count + 1;
        
        if iSubj ~= size(n1LatenciesSOZ, 2)
            stimall{count} = -5 * ones(1, 1);
            names{count} = cellstr(repmat([n1LatenciesSOZ(iSubj).id ' zempty'], 1, 1));
            count = count + 1;
        end
    
    end
       
end

close all
cmap = colormap('parula');
ymax = ceil(max(horzcat(stimall{:})));
temp = horzcat(stimall{:});
ymin = floor(min(temp(temp > 0)));
h = figure(1);
vs = violinplot(horzcat(stimall{:}), vertcat(names{:}), 'ViolinColor', cmap(1, :), 'Width', 0.3);
for iSubj = 1:3:size(vs, 2)
    vs(iSubj).ViolinPlot.FaceColor = cmap(128, :);
    vs(iSubj).ScatterPlot.MarkerFaceColor = cmap(128, :);
end

hold on
count = 1;
for iSubj = 1:size(n1LatenciesSOZ, 2)
    if ~isempty(n1LatenciesSOZ(iSubj).stimSOZ_pFDR)
        if n1LatenciesSOZ(iSubj).stimSOZ_pFDR ==1 
            text(count + 0.3, 1.03 * ymax,'*')
            plot([count, count + 1], 1.02 * [ymax ymax],'k')
        end
        count = count + 3;
    end
end
hold off

ylim([ymin 1.15 * ymax])
xlim([0 size(names, 2) + 1])
h.Units = 'normalized';
h.Position = [0.1 0.1 0.8 0.8];
ax = gca;
ax.FontSize = 13;
ax.XTickLabelRotation = 0;
ax.XTick = 1.5:3:size(vs, 2);
ax.XTickLabel = 1:size([n1LatenciesSOZ(:).stimSOZ_p], 2);
ax.YLabel.String = 'Latency (ms)';
ax.XLabel.String = 'Subjects';
ax.Title.String = 'Latency of responses when stimulating in or outside SOZ';
legend([vs(1).ViolinPlot, vs(2).ViolinPlot], 'outside SOZ', 'inside SOZ')

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS7_LatencyStimSOZ');
set(gcf, 'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)
