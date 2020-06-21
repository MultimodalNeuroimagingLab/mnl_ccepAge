%% load the n1Latencies from the derivatives

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
    % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
else
    disp('Run first ccep02_loadN1.mat')
end

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

% s_nr = 4; % subj 4 age 4, subj 73 age 50 coverage parietal --> frontal
all_age  = zeros(size(n1Latencies,2),1); % number of subjects)
% check if there are N1s for this pair for this subject
rr1 = 3;%1:size(roi,2) % stimulation ROI
rr2 = 4;%1:size(roi,2) % response ROI 
all_age_pairN1 = zeros(size(n1Latencies,2),1); % number of subjects)
for kk = 1:size(n1Latencies,2) % subject number
    all_age(kk) = n1Latencies(kk).age;
    
    % find stimulation pair within specific region
    for ll = 1:size(n1Latencies(kk).run,2) % number of runs
        chanPair = find(sum(contains(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr,roi{rr1}),2)>0);
        chanSig = find(ismember(string(n1Latencies(kk).run(ll).channel_DestrieuxNr),roi{rr2})>0);
        % was there any N1 for this pair?
        for cp = 1:size(chanPair,1)
            for cs = 1:size(chanSig,1)
                if ~isnan(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)))
                    all_age_pairN1(kk) = 1;
                end
            end
        end 
    end
end

subjects_plot = [4 69]; % age 4 38

ccep_age = zeros(length(subjects_plot),1);
average_ccep_age = NaN(length(subjects_plot),5*2048);
average_ccep_age_nonnorm = NaN(length(subjects_plot),5*2048);
sterr_ccep_age_nonnorm = NaN(length(subjects_plot),5*2048);
average_n1_age = NaN(length(subjects_plot),1);

for s_nr = 1:length(subjects_plot) % subject number
    kk = subjects_plot(s_nr);
    ccep_age(s_nr) = n1Latencies(kk).age;

    % get the session directory name
    sesDir = dir(fullfile(myDataPath.output,'derivatives','av_ccep',['sub-' n1Latencies(kk).id],'ses-*'));
    sesDir = sesDir.name;

    av_ccep_allruns = [];
    av_n1_allruns = [];

    for ll = 1:size(n1Latencies(kk).run,2) % number of runs
        % get this run file name
        thisRun = fullfile(myDataPath.output,'derivatives','av_ccep',['sub-' n1Latencies(kk).id],sesDir,...
            theseSubs(kk).run{ll});
        thisData = load(thisRun);

        % find stimulation pair within specific region
        rr1 = 3;%1:size(roi,2) % stimulation ROI
        chanPair = find(sum(contains(n1Latencies(kk).run(ll).average_ccep_DestrieuxNr,roi{rr1}),2)>0);

        % find response electrode within specific region
        rr2 = 4;%1:size(roi,2) % response ROI
        chanSig = find(ismember(string(n1Latencies(kk).run(ll).channel_DestrieuxNr),roi{rr2})>0);

        % collect all signals with stimulation pair and response electrode
        % within specific region
        for cp = 1:size(chanPair,1)
            for cs = 1:size(chanSig,1)
                if ~isnan(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)))
                    % disp(['sub nr ' int2str(kk) ' age ' int2str(age) ' has ROIs'])
                    tempResp = squeeze(thisData.average_ccep(chanSig(cs),chanPair(cp),:));
                    tempN1 = thisData.tt(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)));
                    
                    if length(tempResp)~=5*2048
                        tempResp = resample(tempResp,5*2048,length(tempResp));
                    end                    
                    av_ccep_allruns = [av_ccep_allruns, tempResp];
                    av_n1_allruns = [av_n1_allruns; tempN1];
                    clear tempResp tempN1
                end
            end
        end


        clear thisData thisRun thisN1
    end % end run loop

    % average for this patient across these areas
    average_ccep_age_nonnorm(s_nr,:) = mean(av_ccep_allruns,2);%[roi1/roi2 pairs,samples]
    sterr_ccep_age_nonnorm(s_nr,:) = std(av_ccep_allruns,[],2)./sqrt(size(av_ccep_allruns,2));%[roi1/roi2 pairs,samples]
    % average_ccep_age(s_nr,:) = mean(av_ccep_allruns,2);%[roi1/roi2 pairs,samples]
    average_n1_age(s_nr) = mean(av_n1_allruns);
    clear av_ccep_allruns av_n1_allruns

end

%% plot an example for parietal --> frontal

ttmin = 0.010;
ttmax = .500;
% ttmin = -2;
% ttmax = 3;

cm = winter(size(average_ccep_age_nonnorm,1));

figure('Position',[0 0 600 300]),hold on        
plot(1000*tt(tt>ttmin & tt<ttmax),zeros(size(tt(tt>ttmin & tt<ttmax))),'k')
for kk = 1:size(average_ccep_age_nonnorm,1)
    plot(1000*tt(tt>ttmin & tt<ttmax),average_ccep_age_nonnorm(kk,tt>ttmin & tt< ttmax),'Color',cm(kk,:))
end

figureName = fullfile(myDataPath.output,'derivatives','age',...
            ['AgeExamples_tmax' int2str(ttmax*1000)]);

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)

%% plot rendering 
