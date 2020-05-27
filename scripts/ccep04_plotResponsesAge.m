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

tt = n1Latencies(74).run(1).tt;

average_ccep_age = cell(max([n1Latencies.age]),1);
average_n1_age = cell(max([n1Latencies.age]),1);
for kk = 1:size(n1Latencies,2) % number of subjects
    
    age = n1Latencies(kk).age;
    average_ccep_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2),5*2048); % [roi_start, roi_end,run,tt]
    average_N1_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2)); % [roi_start, roi_end,run,tt]
    
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
                average_ccep_select = NaN(size(chanPair,1),size(chanSig,1),size(thisData.average_ccep,3));
                average_N1_select = NaN(size(chanPair,1),size(chanSig,1)); % N1latency
                for cp = 1:size(chanPair,1)
                    for cs = 1:size(chanSig,1)
                        if ~isnan(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)))
                            if thisData.tt(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)))==0, disp('what???'),end
                            tempResp = squeeze(thisData.average_ccep(chanSig(cs),chanPair(cp),:));
                            average_ccep_select(cp,cs,:) = tempResp;
                            average_N1_select(cp,cs) = thisData.tt(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)));
                            clear tempResp
                        end
                    end
                end
                thisResp = squeeze(nanmean(nanmean(average_ccep_select,1),2));
                thisN1 = squeeze(nanmean(nanmean(average_N1_select,1),2));
                % upsample if necessary  
                if length(thisResp)~=5*2048
                    thisResp = resample(thisResp,5*2048,length(thisResp));
                end
                
                % 'normalize' 0.15-0.100 by unit length
                thisResp = thisResp./sqrt(sum(thisResp(tt>.015 & tt<0.100).^2));
                
                average_ccep_run(rr1,rr2,ll,:) = thisResp;
                average_N1_run(rr1,rr2,ll) = thisN1;
            end
        end
        clear thisData thisRun thisN1
    end
    
    % average for this patient across these areas
    average_ccep_pat = squeeze(nanmean(average_ccep_run,3));
    average_n1_pat = squeeze(nanmean(average_N1_run,3));
    clear average_ccep_run average_N1_run
    
    if ~isempty(average_ccep_age{age})
        n = size(average_ccep_age{age},3);
        average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
        average_n1_age{age}(:,:,n+1) = average_n1_pat;
    else
        n = 0;
        average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
        average_n1_age{age}(:,:,n+1) = average_n1_pat;
    end
    clear average_ccep_pat average_n1_pat
end


save(fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age'),'average_ccep_age','average_n1_age','tt','roi','roi_name');

% load(fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age'),'average_ccep_age','average_n1_age','tt','roi','roi_name');

%%
%% one approach is to sort according to age
%%

clear sortage

for rr1 = 1:4
    for rr2 = 1:4
        sortage(rr1,rr2).average_ccep = [];
        sortage(rr1,rr2).age_ind = [];
        sortage(rr1,rr2).average_n1 = [];

        for age = 1:max([n1Latencies.age])
            if size(average_ccep_age{age},1)>0 % there are some subjects at this age
                addThis = squeeze(average_ccep_age{age}(rr1,rr2,:,:));
                addThis(isnan(addThis(:,1)),:) = [];
                addN1 = squeeze(average_n1_age{age}(rr1,rr2,:));
                addN1(isnan(addN1)) = [];
                if ~isempty(addThis) % there are subjects with electrodes on ROI
                    if size(addThis,1)>size(addThis,2), addThis = addThis'; end
                    nr_subs = size(addThis,1);
                    sortage(rr1,rr2).average_ccep = [sortage(rr1,rr2).average_ccep; addThis];
                    sortage(rr1,rr2).average_n1 = [sortage(rr1,rr2).average_n1; addN1];
                    sortage(rr1,rr2).age_ind = [sortage(rr1,rr2).age_ind; zeros(nr_subs,1)+age];
    %                 clear addThis
                end
            end
        end
    end
end



%% figure of CCEPs + N1 sorted by age

ttmin = 0.010;
ttmax = .100;

figure('Position',[0 0 600 300])
for rr1 = 1:4
    for rr2 = 1:4
        subplot(4,4,(rr1-1)*4+rr2),hold on
        imagesc(1000*tt(tt>ttmin & tt< ttmax),1:length(sortage(rr1,rr2).age_ind),-sortage(rr1,rr2).average_ccep(:,tt>ttmin & tt< ttmax),...
            [-0.1 0.1])
        plot(1000*sortage(rr1,rr2).average_n1,1:length(sortage(rr1,rr2).age_ind),'k.')
        colormap(parula)
        hold on
%         plot([20 20],[1 length(sortage(rr1,rr2).age_ind)],'k')
%         plot([30 30],[1 length(sortage(rr1,rr2).age_ind)],'k')
%         plot([40 40],[1 length(sortage(rr1,rr2).age_ind)],'k')      

%         set(gca,'YTick',1:length(sortage(rr1,rr2).age_ind),'YTickLabel',sortage(rr1,rr2).age_ind)
        set(gca,'YTick',[])
        axis tight
    end
end

figureName = fullfile(myDataPath.output,'derivatives','age',...
            ['AllSortAge_tmax' int2str(ttmax*1000)]);

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)

figure('Position',[0 0 150 40])
imagesc(1:100)
colormap(parula)
axis off
figureName = fullfile(myDataPath.output,'derivatives','age',...
            ['AllSortAge_tmax' int2str(ttmax*1000) '_cm']);
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)

%%
%% approach to average per year
%%

average_ccep_age_mean = cell(size(average_ccep_age));
average_n1_age_mean = cell(size(average_n1_age));
for rr1 = 1:4
    for rr2 = 1:4
        for age = 1:max([n1Latencies.age])
            if size(average_ccep_age{age},1)>0 % there are some subjects at this age
                nr_subs = size(average_ccep_age{age},3);
                respMat = NaN(nr_subs,length(tt));
                for kk = 1:nr_subs
                    thisResp = squeeze(average_ccep_age{age}(rr1,rr2,kk,:));
                    respMat(kk,:) = thisResp;
                end
                average_ccep_age_mean{age}(rr1,rr2,1,:) = mean(respMat,1,'omitnan');
                average_n1_age_mean{age}(rr1,rr2) = mean(average_n1_age{age}(rr1,rr2,:),3,'omitnan');
                clear thisResp respMat nr_subs
            end
        end
    end
end


%%
%% make figure with all ccep signals underneath each other averaged per age
% define rois
for rr1 = 1%1:4
    for rr2 = 3%1:4

        ttmin = 0.010;
        ttmax = .1;
        amp = .01;
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

%         set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',figureName)
%         print('-depsc','-r300',figureName)
    end
end

%%
%% We can use this for an intro figure, using non-normalized CCEPs or so
% todo: copy code from above with this part
%%
%% average over age groups
%% make figure with all ccep signals

% ttmin = -0.02;
% ttmax = .5;
ttmin = 0.015;
ttmax = 0.200;
amp = 0.00001;

age_groups = {[1:20],[21:51]};
% age_groups = {[1:10],[11:20],[21:30],[31:40],[41:51]};
% age_groups = {[1:5],[6:10],[11:15],[16:20],[21:25],[26:30],[31:35],[36:40],[41:45],[46:51]};
% cm = hot(20);
% cm = [cm(1:4:16,:); cm(9:-4:1,:)];
cm = parula(length(age_groups)+1);
cm = cm(end:-1:1,:);
cm = cm.^1.5;
cm(1,:) = cm(1,:)*.9;
cm(1,1) = 1;
plotting_matrix = NaN(length(average_ccep_age_mean),length(tt));

figure('Position',[0 0 600 600]),

% define rois
for rr1 = 1:4
    for rr2 = 1:4
        subplot(4,4,(rr1-1)*4+rr2),hold on
        
        % plot vertical grid in the back
%         for kk = 25:25:100
%             plot([kk kk],[-3 12],'Color',[.8 .8 .8])
%         end
        
        for kk = 1:length(average_ccep_age_mean)
            if ~isempty(average_ccep_age_mean{kk})
                plotting_matrix(kk,:) = average_ccep_age_mean{kk}(rr1,rr2,1,:);
            end
        end
        
        
        for age = length(age_groups):-1:1
            plot(tt(tt>ttmin & tt<ttmax)*1000,zeros(size(tt(tt>ttmin & tt<ttmax))),'Color',[.8 .8 .8])
            plot(tt(tt>ttmin & tt<ttmax)*1000,nanmean(plotting_matrix(age_groups{age},tt>ttmin & tt<ttmax),1),...
                'Color',cm(age,:),'LineWidth',2)
        end
        
        xlim([ttmin ttmax]*1000)%,ylim([-0.2 .05])
        set(gca,'YTick',[])
        ax = gca;
    end
end

% title(sprintf('N1 from %s to %s in increasing age',roi_name{rr},roi_name{rr2}))

figureName = fullfile(myDataPath.output,'derivatives','age',...
            ['CCEPvsAgeGroups_tmax' int2str(ttmax*1000)]);

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)