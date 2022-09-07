%% this scripts investigates the influence of epilepsy on latencies
% Dorien van Blooijs, Dora Hermes, 2020

clc
clear

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
    % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
else
    disp('Run first ccep02_loadN1.mat')
end

%% select patients in whom SOZ and/or RA are delineated
count = 1;
n1LatenciesSOZ = struct;

for subj = 1:size(n1Latencies,2)
    
    if any(contains(n1Latencies(subj).elecs_tsv.soz,'yes')) || any(contains(n1Latencies(subj).elecs_tsv.resected,'yes'))
    
        n1LatenciesSOZ(count).id = n1Latencies(subj).id;
        n1LatenciesSOZ(count).ses = n1Latencies(subj).ses;
        n1LatenciesSOZ(count).age = n1Latencies(subj).age;
        n1LatenciesSOZ(count).elecs_tsv = n1Latencies(subj).elecs_tsv;
        n1LatenciesSOZ(count).run = n1Latencies(subj).run;
        
        % extract seizure onset zone (SOZ) and resected area (RA)
        n1LatenciesSOZ(count).SOZ =  find(strcmp(n1Latencies(subj).elecs_tsv.soz,'yes')==1);
        n1LatenciesSOZ(count).RA =  find(strcmp(n1Latencies(subj).elecs_tsv.resected,'yes')==1);
        
        count = count+1;
        
    end
end

%% percentage SOZ electrodes of total electrodes
numSOZ = NaN(size(n1LatenciesSOZ));
numgrid = NaN(size(n1LatenciesSOZ));

for subj = 1:size(n1LatenciesSOZ,2)
    
    if ~isempty(n1LatenciesSOZ(subj).SOZ)
        numSOZ(subj) =  numel(n1LatenciesSOZ(subj).SOZ);
    end
   numgrid(subj) = numel(n1LatenciesSOZ(subj).run(1).good_channels);
    
end

median(numSOZ./numgrid,'omitnan')
min(numSOZ./numgrid)
max(numSOZ./numgrid)

%% distinguish latencies SOZ and nSOZ
clc
close all

% when SOZ is response electrode
for kk = 1:size(n1LatenciesSOZ,2)
    if ~isnan(n1LatenciesSOZ(kk).SOZ)
        clear respnSOZlat respSOZlat stimnSOZlat stimSOZlat
        
        if ~isempty(n1LatenciesSOZ(kk).run)
            
            for ll = 1:size(n1LatenciesSOZ(kk).run,2)
                clear respSOZlatencies respnSOZlatencies stimSOZlatencies stimnSOZlatencies
                                
                % when SOZ is response electrode
                SOZ = n1LatenciesSOZ(kk).SOZ;
                nSOZ = setdiff(n1LatenciesSOZ(kk).run(ll).good_channels,SOZ);
                respSOZlatencies = n1LatenciesSOZ(kk).run(ll).n1_peak_sample(SOZ,:);
                respSOZlatencies = n1LatenciesSOZ(kk).run(ll).tt(respSOZlatencies(~isnan(respSOZlatencies)));
                respnSOZlatencies = n1LatenciesSOZ(kk).run(ll).n1_peak_sample(nSOZ,:);
                respnSOZlatencies = n1LatenciesSOZ(kk).run(ll).tt(respnSOZlatencies(~isnan(respnSOZlatencies)));
                
                respnSOZlat{ll} = respnSOZlatencies; %#ok<SAGROW>
                respSOZlat{ll} = respSOZlatencies; %#ok<SAGROW>
                
                % when SOZ is stimulated
                stimSOZ = find(contains(n1LatenciesSOZ(kk).run(ll).average_ccep_names,...
                    n1LatenciesSOZ(kk).run(ll).channel_names(SOZ)));
                stimnSOZ = setdiff(1:size(n1LatenciesSOZ(kk).run(ll).n1_peak_sample,2),stimSOZ);
                
                stimSOZlatencies = n1LatenciesSOZ(kk).run(ll).n1_peak_sample(:,stimSOZ);
                stimSOZlatencies = n1LatenciesSOZ(kk).run(ll).tt(stimSOZlatencies(~isnan(stimSOZlatencies)));
                stimnSOZlatencies = n1LatenciesSOZ(kk).run(ll).n1_peak_sample(:,stimnSOZ);
                stimnSOZlatencies = n1LatenciesSOZ(kk).run(ll).tt(stimnSOZlatencies(~isnan(stimnSOZlatencies)));
                
                stimnSOZlat{ll} = stimnSOZlatencies; %#ok<SAGROW>
                stimSOZlat{ll} = stimSOZlatencies; %#ok<SAGROW>
                
            end
            
            % mann witney u test for responding electrode in either SOZ or
            % nSOZ
            n1LatenciesSOZ(kk).respSOZlatencies = horzcat(respSOZlat{:});
            n1LatenciesSOZ(kk).respnSOZlatencies = horzcat(respnSOZlat{:});
            if (any(~isnan(n1LatenciesSOZ(kk).respSOZlatencies )) && any(~isnan(n1LatenciesSOZ(kk).respnSOZlatencies ))) || (~isempty(n1LatenciesSOZ(kk).respSOZlatencies) && ~isempty(n1LatenciesSOZ(kk).respnSOZlatencies))
                n1LatenciesSOZ(kk).respSOZ_p = ranksum(n1LatenciesSOZ(kk).respnSOZlatencies,n1LatenciesSOZ(kk).respSOZlatencies);
                
                fprintf('-- %s: When comparing latency in response SOZ and response nSOZ: p = %1.3f with median resp_SOZ = %1.3f sec and median resp_nSOZ = %1.3f sec\n',...
                    n1LatenciesSOZ(kk).id, n1LatenciesSOZ(kk).respSOZ_p, median([n1LatenciesSOZ(kk).respSOZlatencies]),median([n1LatenciesSOZ(kk).respnSOZlatencies]))
            end  
            
            % mann witney u test for stimulated electrode in either SOZ or
            % nSOZ
            n1LatenciesSOZ(kk).stimSOZlatencies = horzcat(stimSOZlat{:});
            n1LatenciesSOZ(kk).stimnSOZlatencies = horzcat(stimnSOZlat{:});
            if (any(~isnan(n1LatenciesSOZ(kk).stimSOZlatencies )) && ...
                    any(~isnan(n1LatenciesSOZ(kk).stimnSOZlatencies ))) || ...
                    (~isempty(n1LatenciesSOZ(kk).stimSOZlatencies) && ...
                    ~isempty(n1LatenciesSOZ(kk).stimnSOZlatencies))
                
                n1LatenciesSOZ(kk).stimSOZ_p = ranksum(n1LatenciesSOZ(kk).stimnSOZlatencies,n1LatenciesSOZ(kk).stimSOZlatencies);
                
                fprintf('-- %s: When comparing latency in stimulated SOZ and stimulated nSOZ: p = %1.3f with median stim_SOZ = %1.3f sec and median stim_nSOZ = %1.3f sec\n',...
                    n1LatenciesSOZ(kk).id, n1LatenciesSOZ(kk).stimSOZ_p, median([n1LatenciesSOZ(kk).stimSOZlatencies]),median([n1LatenciesSOZ(kk).stimnSOZlatencies]))
            end  
        end
    end
end

% apply FDR correction
p_vals_resp = [n1LatenciesSOZ(:).respSOZ_p];
p_vals_stim = [n1LatenciesSOZ(:).stimSOZ_p];

[p_sort_resp,p_ind_resp] = sort(p_vals_resp(:));
[p_sort_stim,p_ind_stim] = sort(p_vals_stim(:));

m = length(p_vals_resp);
thisVal = NaN(size(p_sort_resp));
for kk = 1:length(p_sort_resp)
    thisVal(kk) = (kk/m)*0.05;
end

figure,subplot(1,2,1),hold on,plot(thisVal),plot(p_sort_resp,'r.'),title('Sig p-val in response SOZ/nSOZ after FDR correction')
subplot(1,2,2),hold on,plot(thisVal),plot(p_sort_stim,'r.'),title('Sig p-val in stimulated SOZ/nSOZ after FDR correction')

p_sig_resp = p_vals_resp;
p_sig_resp(p_ind_resp) = p_sort_resp<thisVal;
p_sig_stim = p_vals_stim;
p_sig_stim(p_ind_stim) = p_sort_stim<thisVal;
count = 1;
for subj = 1:size(n1LatenciesSOZ,2)
    
    if ~isempty(n1LatenciesSOZ(subj).SOZ)
        n1LatenciesSOZ(subj).respSOZ_pFDR = p_sig_resp(count);
        n1LatenciesSOZ(subj).stimSOZ_pFDR = p_sig_stim(count);
        count = count+1;
    end
end

%% violin plot responses in or outside SOZ
respall = cell(1);
names = cell(1);

count = 1;
for kk = 1:size(n1LatenciesSOZ,2)
    
    if ~isempty(n1LatenciesSOZ(kk).respSOZlatencies) && ~isempty(n1LatenciesSOZ(kk).respnSOZlatencies)
        respall{count} = n1LatenciesSOZ(kk).respSOZlatencies*1000;
        names{count} = cellstr(repmat([n1LatenciesSOZ(kk).id ' soz'],size(respall{count},2),1));
        count = count+1;
        respall{count} = n1LatenciesSOZ(kk).respnSOZlatencies*1000;
        names{count} = cellstr(repmat([n1LatenciesSOZ(kk).id ' nsoz'],size(respall{count},2),1));
        count = count+1;
        
        if kk ~= size(n1LatenciesSOZ,2)
            respall{count} = -5*ones(1,1);
            names{count} = cellstr(repmat([n1LatenciesSOZ(kk).id ' zempty'],1,1));
            count = count+1;
        end
    
    end
       
end

close all
cmap = colormap('parula');
ymax = ceil(max(horzcat(respall{:})));
temp = horzcat(respall{:});
ymin = floor(min(temp(temp>0)));
h = figure(1);
vs = violinplot(horzcat(respall{:}),vertcat(names{:}),'ViolinColor',cmap(1,:),'Width',0.3);
for subj = 1:3:size(vs,2)
    vs(subj).ViolinPlot.FaceColor = cmap(128,:);
    vs(subj).ScatterPlot.MarkerFaceColor = cmap(128,:);
end

hold on
count = 1;
for subj=1:size(n1LatenciesSOZ,2)
    if ~isempty(n1LatenciesSOZ(subj).respSOZ_pFDR)
        if n1LatenciesSOZ(subj).respSOZ_pFDR ==1 
            text(count+0.3,1.03*ymax,'*')
            plot([count, count+1],1.02*[ymax ymax],'k')
        end
        count = count+3;
    end
end
hold off

ylim([ymin 1.15*ymax])
xlim([0 size(names,2)+1])
h.Units = 'normalized';
h.Position = [0.1 0.1 0.8 0.8];
ax = gca;
ax.FontSize = 13;
ax.XTickLabelRotation = 0;
ax.XTick = 1.5:3:size(vs,2);
ax.XTickLabel = 1:size([n1LatenciesSOZ(:).respSOZ_p],2);
ax.YLabel.String = 'Latency (ms)';
ax.XLabel.String = 'Subjects';
ax.Title.String = 'Latency of responses in or outside SOZ';
legend([vs(1).ViolinPlot,vs(2).ViolinPlot],'outside SOZ','inside SOZ')

figureName = fullfile(myDataPath.output,'derivatives','age',...
    'latency_respSOZ');
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-depsc','-r300',figureName)

%% violin plot responses in or outside SOZ
stimall = cell(1);
names = cell(1);

count = 1;
for kk = 1:size(n1LatenciesSOZ,2)
    
    if ~isempty(n1LatenciesSOZ(kk).stimSOZlatencies) && ~isempty(n1LatenciesSOZ(kk).stimnSOZlatencies)
        stimall{count} = n1LatenciesSOZ(kk).stimSOZlatencies*1000;
        names{count} = cellstr(repmat([n1LatenciesSOZ(kk).id ' soz'],size(stimall{count},2),1));
        count = count+1;
        stimall{count} = n1LatenciesSOZ(kk).stimnSOZlatencies*1000;
        names{count} = cellstr(repmat([n1LatenciesSOZ(kk).id ' nsoz'],size(stimall{count},2),1));
        count = count+1;
        
        if kk ~= size(n1LatenciesSOZ,2)
            stimall{count} = -5*ones(1,1);
            names{count} = cellstr(repmat([n1LatenciesSOZ(kk).id ' zempty'],1,1));
            count = count+1;
        end
    
    end
       
end

close all
cmap = colormap('parula');
ymax = ceil(max(horzcat(stimall{:})));
temp = horzcat(stimall{:});
ymin = floor(min(temp(temp>0)));
h = figure(1);
vs = violinplot(horzcat(stimall{:}),vertcat(names{:}),'ViolinColor',cmap(1,:),'Width',0.3);
for subj = 1:3:size(vs,2)
    vs(subj).ViolinPlot.FaceColor = cmap(128,:);
    vs(subj).ScatterPlot.MarkerFaceColor = cmap(128,:);
end

hold on
count = 1;
for subj=1:size(n1LatenciesSOZ,2)
    if ~isempty(n1LatenciesSOZ(subj).stimSOZ_pFDR)
        if n1LatenciesSOZ(subj).stimSOZ_pFDR ==1 
            text(count+0.3,1.03*ymax,'*')
            plot([count, count+1],1.02*[ymax ymax],'k')
        end
        count = count+3;
    end
end
hold off

ylim([ymin 1.15*ymax])
xlim([0 size(names,2)+1])
h.Units = 'normalized';
h.Position = [0.1 0.1 0.8 0.8];
ax = gca;
ax.FontSize = 13;
ax.XTickLabelRotation = 0;
ax.XTick = 1.5:3:size(vs,2);
ax.XTickLabel = 1:size([n1LatenciesSOZ(:).stimSOZ_p],2);
ax.YLabel.String = 'Latency (ms)';
ax.XLabel.String = 'Subjects';
ax.Title.String = 'Latency of responses when stimulating in or outside SOZ';
legend([vs(1).ViolinPlot,vs(2).ViolinPlot],'outside SOZ','inside SOZ')

figureName = fullfile(myDataPath.output,'derivatives','age',...
    'latency_stimSOZ');
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-depsc','-r300',figureName)
