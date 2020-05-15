%% load the n1Latencies from the derivatives

% if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')

% these do not have the average CCEPs for areas or age groups yet, but that
% would get heavy on the memory

%%
% temporal areas:
% G_temporal_inf, G_temporal_middle
% roi = {'37','38'};
% G_temporal_inf, G_temporal_middle, G_temp_sup-Lateral,
% G_oc-temp_med-Parahip, G_oc-temp_lat-fusifor
roi{1} = {'37','38','34','23','21'};
roi_name{1} = 'temporal';
% frontal areas:
% G_front_inf-Triangul, G_front_middle, G_front_inf-Opercular
roi{2} = {'14','15','12'}; % maybe add 16: G_front_sup
roi_name{2} = 'frontal';% parietal areas:
% G_pariet_inf-Angular, G_pariet_inf-Supramar, G_parietal_sup
roi{3} = {'25','26','27'};
roi_name{3} = 'parietal';% occipital areas:
% G_occipital_middle G_oc-temp_med-Lingual Pole_occipital
roi{4} = {'19','22','42'};
roi_name{4} = 'occipital';% sensorimotor:
% G_postcentral G_precentral S_central
roi{5} = {'28','29','46'};
roi_name{5} = 'central';

average_ccep_age = cell(max([n1Latencies.age]),1);
for kk = 1:size(n1Latencies,2) % number of subjects
    
    age = n1Latencies(kk).age;
    average_ccep_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2),5*2048); % [roi_start, roi_end,run,tt]
    
    % get the session directory name
    sesDir = dir(fullfile(myDataPath.output,'derivatives','av_ccep',['sub-' n1Latencies(kk).id],'ses-*'));
    sesDir = sesDir.name;
    
    for ll = 1:size(n1Latencies(kk).run,2) % number of runs
        % get this run file name
        thisRun = fullfile(myDataPath.output,'derivatives','av_ccep',['sub-' n1Latencies(kk).id],sesDir,...
            theseSubs(kk).run{ll});
        thisData = load(thisRun);
        % we will resample all cceps to 2048Hz, no need to preselect
        for rr1 = 1%1:size(roi,2) % stimulation ROI
            % find stimulation pair within specific region
            chanPair = find(sum(contains(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr,roi{rr1}),2)>0);
            for rr2 = 3%1:size(roi,2) % response ROI
                % find response electrode within specific region
                chanSig = find(ismember(string(n1Latencies(kk).run(ll).channel_DestrieuxNr),roi{rr2})>0);
                
                % collect all signals with stimulation pair and response
                % electrode within specific region
                
                average_ccep_select = NaN(size(chanPair,1),size(chanSig,2),size(thisData.average_ccep,3));
                for cp = 1:size(chanPair,1)
                    for cs = 1:size(chanSig,1)
                        if ~isnan(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)))
                            average_ccep_select(cp,cs,:) = thisData.average_ccep(chanSig(cs),chanPair(cp),:);
                        end
                    end
                end
                thisResp = squeeze(nanmean(nanmean(average_ccep_select,1),2));
                % upsample if necessary  
                if length(thisResp)~=5*2048
                    thisResp = resample(thisResp,5*2048,length(thisResp));
                end
                average_ccep_run(rr1,rr2,ll,:) = thisResp;
            end
        end
        clear thisData thisRun
    end
    
    % average for this patient across these areas
    average_ccep_pat = squeeze(nanmean(average_ccep_run,3));
    
    if ~isempty(average_ccep_age{age})
        n = size(average_ccep_age{age},3);
        average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
    else
        n = 0;
        average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
    end
end

%% average per year

average_ccep_age_mean = cell(size(average_ccep_age));
for age = 1:max([n1Latencies.age])
    average_ccep_age_mean{age} = mean(average_ccep_age{age},3,'omitnan');
end

%%
%% make figure with all ccep signals
% define rois
rr1 = 1;
rr2 = 3;

tt = n1Latencies(74).run(1).tt;
ttmin = -0.02;
ttmax = 0.1;
amp = 500;
ymin = (min([n1Latencies.age])-1)*amp;
ymax = (max([n1Latencies.age])+1)*amp;

figure(2),
hold on
for age = 1:max([n1Latencies.age])
    plot(tt(tt>ttmin & tt< ttmax),zeros(size(tt(tt>ttmin & tt<ttmax)))+amp*age,'Color',[.8 .8 .8])
    if ~isempty(average_ccep_age_mean{age})
        plot(tt(tt>ttmin & tt<ttmax),squeeze(average_ccep_age_mean{age}(rr,rr2,:,tt>ttmin & tt<ttmax))+amp*age)
    end
end
hold off

xlabel('Time (s)')
ylabel('Age (years)')
ylim([ymin,ymax])
ax = gca;
ax.YTick = (1:max([n1Latencies.age]))*amp;
ax.YTickLabel = num2cell(1:max([n1Latencies.age]));
title(sprintf('N1 from %s to %s in increasing age',roi_name{rr},roi_name{rr2}))

%%
%% Add a normalization per subject

tt = n1Latencies(74).run(1).tt;

rr1 = 1;
rr2 = 3;

average_ccep_age_mean = cell(size(average_ccep_age));
for age = 1:max([n1Latencies.age])
    if size(average_ccep_age{age},1)>0 % there are some subjects at this age
        nr_subs = size(average_ccep_age{age},3);
        respMat = NaN(nr_subs,length(tt));
        for kk = 1:nr_subs
            thisResp = squeeze(average_ccep_age{age}(rr1,rr2,kk,:));
            % just devide by ccep amplitude:
            thisResp = -thisResp./min(thisResp(tt>.015 & tt<.1));
            respMat(kk,:) = thisResp;
        end
        average_ccep_age_mean{age}(rr1,rr2,1,:) = mean(respMat,1,'omitnan');
        clear thisResp respMat
    end
end

%%
%% make figure with all ccep signals
% define rois
rr1 = 1;
rr2 = 3;

tt = n1Latencies(74).run(1).tt;
ttmin = -0.02;
ttmax = 1;
amp = .25;
ymin = (min([n1Latencies.age])-5)*amp;
ymax = (max([n1Latencies.age])+5)*amp;

figure,
hold on
for age = 1:max([n1Latencies.age])
    plot(tt(tt>ttmin & tt< ttmax),zeros(size(tt(tt>ttmin & tt<ttmax)))+amp*age,'Color',[.8 .8 .8])
    if ~isempty(average_ccep_age_mean{age})
        plot(tt(tt>ttmin & tt<ttmax),squeeze(average_ccep_age_mean{age}(rr1,rr2,:,tt>ttmin & tt<ttmax))+amp*age)
    end
end
hold off

xlabel('Time (s)')
ylabel('Age (years)')
xlim([ttmin ttmax]),ylim([ymin,ymax])
ax = gca;
ax.YTick = (1:max([n1Latencies.age]))*amp;
ax.YTickLabel = num2cell(1:max([n1Latencies.age]));
title(sprintf('N1 from %s to %s in increasing age',roi_name{rr},roi_name{rr2}))

%% average over age groups
%% make figure with all ccep signals
% define rois
rr1 = 1;
rr2 = 3;

tt = n1Latencies(74).run(1).tt;
% ttmin = -0.02;
% ttmax = .5;
ttmin = 0.005;
ttmax = .4;
amp = .25;

age_groups = {[1:10],[11:20],[21:30],[31:40],[41:51]};
% age_groups = {[1:5],[6:10],[11:15],[16:20],[21:25],[26:30],[31:35],[36:40],[41:45],[46:51]};
% cm = hot(20);
% cm = [cm(1:4:16,:); cm(9:-4:1,:)];
cm = parula(length(age_groups)+1);

plotting_matrix = NaN(length(average_ccep_age_mean),length(tt));

for kk = 1:length(average_ccep_age_mean)
    if ~isempty(average_ccep_age_mean{kk})
        plotting_matrix(kk,:) = average_ccep_age_mean{kk}(rr1,rr2,1,:);
    end
end

figure,
hold on

for age = 1:length(age_groups)
    plot(tt(tt>ttmin & tt<ttmax),zeros(size(tt(tt>ttmin & tt<ttmax))),'Color',[.8 .8 .8]) 
    plot(tt(tt>ttmin & tt<ttmax),nanmean(plotting_matrix(age_groups{age},tt>ttmin & tt<ttmax),1),...
        'Color',cm(age,:),'LineWidth',2)
end


xlabel('Time (s)')
ylabel('Age (years)')
xlim([ttmin ttmax]),ylim([-1 .5])
ax = gca;
% title(sprintf('N1 from %s to %s in increasing age',roi_name{rr},roi_name{rr2}))

