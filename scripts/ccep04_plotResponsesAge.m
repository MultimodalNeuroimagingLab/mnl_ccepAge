%% load the n1Latencies from the derivatives

myDataPath = setLocalDataPath(1);

% if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')

% these do not have the average CCEPs for areas or age groups yet, but that
% would get heavy on the memory

% get the filename of all derivative datasets:
theseSubs = ccep_getSubFilenameInfo(myDataPath);

%% load ccep responses for each roi

% temporal areas:
% G_temporal_inf, G_temporal_middle
% roi = {'37','38'};
% G_temporal_inf, G_temporal_middle, G_temp_sup-Lateral,
% G_oc-temp_med-Parahip, G_oc-temp_lat-fusifor
roi{1} = {'37','38','34','23','21'};
roi_name{1} = 'temporal';

% central areas
% G_postcentral G_precentral S_central
roi{2} = {'28','29','46'};
roi_name{2} = 'central';

% parietal areas:
% G_pariet_inf-Angular, G_pariet_inf-Supramar, G_parietal_sup
roi{3} = {'25','26','27'};
roi_name{3} = 'parietal';

% frontal areas:
% G_front_inf-Triangul, G_front_middle, G_front_inf-Opercular
roi{4} = {'14','15','12'}; % maybe add 16: G_front_sup
roi_name{4} = 'frontal';

% occipital areas:
% % G_occipital_middle G_oc-temp_med-Lingual Pole_occipital
% roi{5} = {'19','22','42'};
% roi_name{5} = 'occipital';

tt = n1Latencies(74).run(1).tt;

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
        for rr1 = 1:size(roi,2) % stimulation ROI
            % find stimulation pair within specific region
            chanPair = find(sum(contains(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr,roi{rr1}),2)>0);
            for rr2 = 1:size(roi,2) % response ROI
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
                
                % 'normalize' response
                thisResp = -thisResp./min(thisResp(tt>.015 & tt<.1));
                
                average_ccep_run(rr1,rr2,ll,:) = thisResp;
            end
        end
        clear thisData thisRun
    end
    
    % average for this patient across these areas
    average_ccep_pat = squeeze(nanmean(average_ccep_run,3));
    clear average_ccep_run
    
    if ~isempty(average_ccep_age{age})
        n = size(average_ccep_age{age},3);
        average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
    else
        n = 0;
        average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
    end
    clear average_ccep_pat
end


save(fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age'),'average_ccep_age','tt','roi','roi_name');


%%
%% Average per year

average_ccep_age_mean = cell(size(average_ccep_age));
for rr1 = 1:4
    for rr2 = 1:4
        for age = 1:max([n1Latencies.age])
            if size(average_ccep_age{age},1)>0 % there are some subjects at this age
                nr_subs = size(average_ccep_age{age},3);
                respMat = NaN(nr_subs,length(tt));
                for kk = 1:nr_subs
                    thisResp = squeeze(average_ccep_age{age}(rr1,rr2,kk,:));
%                     % just devide by ccep amplitude:
%                     thisResp = -thisResp./min(thisResp(tt>.015 & tt<.1));
                    respMat(kk,:) = thisResp;
                end
                average_ccep_age_mean{age}(rr1,rr2,1,:) = mean(respMat,1,'omitnan');
                clear thisResp respMat nr_subs
            end
        end
    end
end

%%
%% make figure with all ccep signals underneath each other
% define rois
for rr1 = 1:4
    for rr2 = 1:4

        ttmin = 0.010;
        ttmax = .5;
        amp = .1;
        ymin = (min([n1Latencies.age])-12)*amp;
        ymax = (max([n1Latencies.age])+5)*amp;

        figure('Position',[0 0 350 600]),
        hold on
        cm = parula(max([n1Latencies.age]));

        % plot vertical grid in the back
        for kk = 20:20:100
            plot([kk kk],[ymin ymax],'Color',[.8 .8 .8])
        end
        for age = 1:max([n1Latencies.age])

            if ~isempty(average_ccep_age_mean{age})
                plot(tt(tt>0 & tt< ttmax)*1000,zeros(size(tt(tt>0 & tt<ttmax)))+amp*age,'Color',[.8 .8 .8])
                plot(tt(tt>ttmin & tt<ttmax)*1000,squeeze(average_ccep_age_mean{age}(rr1,rr2,:,tt>ttmin & tt<ttmax))+amp*age,...
                    'Color',cm(age,:),'LineWidth',2)
            end
        end
        hold off

        xlabel('Time (ms)')
        ylabel('Age (years)')
        xlim([0 ttmax]*1000),ylim([ymin,ymax])
        ax = gca;
        % ax.YTick = (1:max([n1Latencies.age]))*amp;
        % ax.YTickLabel = num2cell(1:max([n1Latencies.age]));
        ax.YTick = (0:5:50)*amp;
        ax.YTickLabel = num2cell(0:5:50);
        title(sprintf('N1 from %s to %s in increasing age',roi_name{rr1},roi_name{rr2}))

        figureName = fullfile(myDataPath.output,'derivatives','age',...
            ['CCEPvsAge_vertical_' roi_name{rr1} '_' roi_name{rr2} '_tmax' int2str(ttmax*1000)]);

        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',figureName)
        print('-depsc','-r300',figureName)
    end
end

%% average over age groups
%% make figure with all ccep signals

% ttmin = -0.02;
% ttmax = .5;
ttmin = 0.010;
ttmax = 0.300;
amp = 0.25;

% age_groups = {[1:10],[11:20],[21:30],[31:40],[41:51]};
age_groups = {[1:5],[6:10],[11:15],[16:20],[21:25],[26:30],[31:35],[36:40],[41:45],[46:51]};
% cm = hot(20);
% cm = [cm(1:4:16,:); cm(9:-4:1,:)];
cm = parula(length(age_groups)+1);

plotting_matrix = NaN(length(average_ccep_age_mean),length(tt));

figure('Position',[0 0 600 600]),

% define rois
for rr1 = 1:4
    for rr2 = 1:4
        subplot(4,4,(rr2-1)*4+rr1),hold on
        
        % plot vertical grid in the back
        for kk = 25:25:100
            plot([kk kk],[-3 12],'Color',[.8 .8 .8])
        end
        
        for kk = 1:length(average_ccep_age_mean)
            if ~isempty(average_ccep_age_mean{kk})
                plotting_matrix(kk,:) = average_ccep_age_mean{kk}(rr1,rr2,1,:);
            end
        end
        
        
        for age = 1:length(age_groups)
            plot(tt(tt>ttmin & tt<ttmax)*1000,zeros(size(tt(tt>ttmin & tt<ttmax))),'Color',[.8 .8 .8])
            plot(tt(tt>ttmin & tt<ttmax)*1000,nanmean(plotting_matrix(age_groups{age},tt>ttmin & tt<ttmax),1),...
                'Color',cm(age,:),'LineWidth',2)
        end
        
        xlim([0 ttmax]*1000),ylim([-1.2 .5])
        ax = gca;
    end
end

% title(sprintf('N1 from %s to %s in increasing age',roi_name{rr},roi_name{rr2}))

