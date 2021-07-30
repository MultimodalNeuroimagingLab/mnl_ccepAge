%% load the n1Latencies from the derivatives

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
    % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
else
    disp('Run first ccep02_loadN1.mat')
end


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

tt = n1Latencies(2).run(1).tt;

average_ccep_age = cell(max([n1Latencies.age]),1);
average_ccep_age_nonnorm = cell(max([n1Latencies.age]),1);
average_n1_age = cell(max([n1Latencies.age]),1);
for kk = 1:size(n1Latencies,2) % number of subjects
    
    age = n1Latencies(kk).age;
    average_ccep_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2),5*2048); % [roi_start, roi_end,run,tt]
    average_ccep_run_nonnorm = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2),5*2048); % [roi_start, roi_end,run,tt]
    average_N1_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2)); % [roi_start, roi_end,run,tt]
    
    % get the session directory name
    sesDir = dir(fullfile(myDataPath.output,'derivatives','av_ccep',n1Latencies(kk).id,'ses-*'));
    sesDir = sesDir.name;
    
    for ll = 1:size(n1Latencies(kk).run,2) % number of runs
        % get this run file name
        thisRun = fullfile(myDataPath.output,'derivatives','av_ccep',n1Latencies(kk).id,sesDir,...
            n1Latencies(kk).run(ll).runName);
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
                
                % if responses are added to average_ccep_select
                if any(~isnan(average_ccep_select(:)))
                    thisResp = squeeze(nanmean(nanmean(average_ccep_select,1),2));
                    thisN1 = squeeze(nanmean(nanmean(average_N1_select,1),2));
                    % upsample if necessary
                    if length(thisResp)~=5*2048
                        thisResp = resample(thisResp,5*2048,length(thisResp));
                    end
                    
                    thisResp2 = thisResp; % non-normalized
                    
                    % 'normalize' 0.15-0.100 by unit length
                    thisResp = thisResp./sqrt(sum(thisResp(tt>.015 & tt<0.100).^2));
                    
                    average_ccep_run_nonnorm(rr1,rr2,ll,:) = thisResp2; %[roi1,roi2,run,samples]
                    average_ccep_run(rr1,rr2,ll,:) = thisResp; %[roi1,roi2,run,samples]
                    average_N1_run(rr1,rr2,ll) = thisN1;
                else
                    average_ccep_run_nonnorm(rr1,rr2,ll,:) = NaN(5*2048,1);
                    average_ccep_run(rr1,rr2,ll,:) = NaN(5*2048,1);
                    average_N1_run(rr1,rr2,ll) = NaN;
                end
            end
        end
        clear thisData thisRun thisN1
    end
    
    % average for this patient across these areas
    average_ccep_pat_nonnorm = squeeze(nanmean(average_ccep_run_nonnorm,3));%[roi1, roi2,samples]
    average_ccep_pat = squeeze(nanmean(average_ccep_run,3));%[roi1, roi2,samples]
    average_n1_pat = squeeze(nanmean(average_N1_run,3));
    clear average_ccep_run average_N1_run
    
    if ~isempty(average_ccep_age{age})
        n = size(average_ccep_age{age},3);
        average_ccep_age_nonnorm{age}(:,:,n+1,:) = average_ccep_pat_nonnorm;
        average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
        average_n1_age{age}(:,:,n+1) = average_n1_pat;
    else
        n = 0;
        average_ccep_age_nonnorm{age}(:,:,n+1,:) = average_ccep_pat_nonnorm;
        average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
        average_n1_age{age}(:,:,n+1) = average_n1_pat;
    end
    clear average_ccep_pat average_ccep_pat_nonnorm average_n1_pat
end


% save(fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age'),'average_ccep_age','average_n1_age','tt','roi','roi_name');

% load(fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age'),'average_ccep_age','average_n1_age','tt','roi','roi_name');

%%
%% one approach is to sort according to age
%%

clear sortage

for rr1 = 1:4
    for rr2 = 1:4
        sortage(rr1,rr2).average_ccep = [];
        sortage(rr1,rr2).average_ccep_nonnorm = [];
        sortage(rr1,rr2).age_ind = [];
        sortage(rr1,rr2).average_n1 = [];

        for age = 1:max([n1Latencies.age])
            if size(average_ccep_age{age},1)>0 % there are some subjects at this age
                addThis = squeeze(average_ccep_age{age}(rr1,rr2,:,:));
                addThis(isnan(addThis(:,1)),:) = [];
                addThisNonnorm = squeeze(average_ccep_age_nonnorm{age}(rr1,rr2,:,:));
                addThisNonnorm(isnan(addThisNonnorm(:,1)),:) = [];
                addN1 = squeeze(average_n1_age{age}(rr1,rr2,:));
                addN1(isnan(addN1)) = [];
                if ~isempty(addThis) % there are subjects with electrodes on ROI
                    if size(addThis,1)>size(addThis,2), addThis = addThis'; addThisNonnorm = addThisNonnorm'; end
                    nr_subs = size(addThis,1);
                    sortage(rr1,rr2).average_ccep = [sortage(rr1,rr2).average_ccep; addThis];
                    sortage(rr1,rr2).average_ccep_nonnorm = [sortage(rr1,rr2).average_ccep_nonnorm; addThisNonnorm];
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
        set(gca,'XTick',[20:20:80],'YTick',[])
        axis tight
    end
end

% calculate correlation and p
p_all = zeros(4,4);
r_all = zeros(4,4);
for rr1 = 1:4
    for rr2 = 1:4
        n1_latency = 1000*sortage(rr1,rr2).average_n1;
        age = sortage(rr1,rr2).age_ind;
        [r,p] = corr(n1_latency,age);
        p_all(rr1,rr2) = p;
        r_all(rr1,rr2) = r;
    end
end

% FDR correction
p_vals = p_all(:);

m = length(p_vals);
[p_sort,p_ind] = sort(p_vals(:));
clear thisVal
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end
% figure,hold on,plot(thisVal),plot(p_sort,'r.')
p_sig = p_all;
p_sig(p_ind) = p_sort<thisVal';
for rr1 = 1:4
    for rr2 = 1:4
        subplot(4,4,(rr1-1)*4+rr2),hold on
        if p_sig(rr1,rr2)==1 % significant!
            plot(100,0,'r*')
        end
    end
end

figureName = fullfile(myDataPath.output,'derivatives','age',...
            ['AllSortAge_tmax' int2str(ttmax*1000)]);

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)

%%

figure('Position',[0 0 150 40])
imagesc(1:100)
colormap(parula)
axis off
figureName = fullfile(myDataPath.output,'derivatives','age',...
            ['AllSortAge_tmax' int2str(ttmax*1000) '_cm']);
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)

