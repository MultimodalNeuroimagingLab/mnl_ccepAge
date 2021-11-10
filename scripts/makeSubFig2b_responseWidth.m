%% load the n1Latencies from the derivatives
% This code is used to plot the normalized ccep of all patients in order of
% age. This figure is displayed as Figure 2 in the article.


%% 1. load the n1Latencies from the derivatives
% we use this code both for analysis in the main script, and for checks
% with only subjects in whom it is certain that 8mA is used for
% stimulation.

% in the script makeSubFig3_only8maSubs, mode is defined as follows:
% mode = {'8mA'}; --> if this is defined, n1Latencies from derivatives are
% not loaded.
clear
close all

selectPat = input('Would you like to include all patients, or only the ones for whom it is certain that 8mA was applied (supplemental material)? [all/8] ','s');

if strcmp(selectPat,'all')
    select_amplitude = 0; % make this 8 for only 8mA
elseif strcmp(selectPat,'8')
    select_amplitude = 8;
else
    error('Answer to previous question is not recognized.')
end

if select_amplitude==0
    myDataPath = setLocalDataPath(1);
    if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
        % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
        load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')

        filename_averageCCEP_width = fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age_width.mat');

    else
        disp('Run first ccep02_loadN1.mat')
    end
elseif select_amplitude==8 % only 8 mA
    myDataPath = setLocalDataPath(1);
    if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_8ma.mat'),'file')
        % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
        load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_8ma.mat'),'n1Latencies8ma')

        n1Latencies = n1Latencies8ma;
        filename_averageCCEP_width = fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age_width_8ma_.mat');

    else
        disp('Run first script ccep02_loadN1.m')
    end
end

%% 2. load ccep responses and categorize into connections from stimulated region to responding regions
%% skips automatically to 3. if you ran and saved output before (takes ~5 mins to run)
%
% rois: temporal, central, parietal, frontal

% the CCEPs are averaged for each run, and then averaged CCEPs per patient
% are collected for all subjects. Multiple subjects with the same age are
% collected for each age (average_ccep_age_nonnorm) and normalized
% (average_ccep_age)

if ~exist(filename_averageCCEP_width,'file')

    % categorize anatomical regions
    ccep_categorizeAnatomicalRegions

    tt = n1Latencies(2).run(1).tt; % this patient had fs=2048, so take this tt

    average_ccep_age = cell(max([n1Latencies.age]),1);
    average_ccep_age_nonnorm = cell(max([n1Latencies.age]),1);
    average_n1_age = cell(max([n1Latencies.age]),1);
    average_width_age = cell(max([n1Latencies.age]),1);

    for kk = 1:size(n1Latencies,2) % number of subjects

        fprintf('Load subj %d of %d \n',kk,size(n1Latencies,2))

        age = n1Latencies(kk).age;
        average_ccep_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2),5*2048); % [roi_start, roi_end,run,tt]
        average_ccep_run_nonnorm = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2),5*2048); % [roi_start, roi_end,run,tt]
        average_N1_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2)); % [roi_start, roi_end,run,tt]
        average_width_run = NaN(size(roi,2),size(roi,2),size(n1Latencies(kk).run,2)); % [roi_start, roi_end,run,tt]

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
                    average_width_select = NaN(size(chanPair,1),size(chanSig,1)); % N1width
                    for cp = 1:size(chanPair,1)
                        for cs = 1:size(chanSig,1)
                            if ~isnan(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)))
                                if thisData.tt(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)))==0, disp('what???'),end
                                tempResp = squeeze(thisData.average_ccep(chanSig(cs),chanPair(cp),:));
                                average_ccep_select(cp,cs,:) = tempResp;
                                average_N1_select(cp,cs) = thisData.tt(n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp)));

                                thisSample = n1Latencies(kk).run(ll).n1_peak_sample(chanSig(cs),chanPair(cp));
                                thisPeak = tempResp(thisSample);
                                this_th = .5*thisPeak;
                                steps_back = find(tempResp(thisSample:-1:1)>this_th,1)-1; % walk back from min
                                steps_forw = find(tempResp(thisSample:1:end)>this_th,1)-1; % walk forward from min
                                if ~isempty(steps_back) && ~isempty(steps_forw) % offset in response
                                    if thisData.tt(thisSample-steps_back)>0 && ... % can't onset <0, may occur if not at baseline at 9 ms
                                            (thisData.tt(thisSample+steps_forw)-thisData.tt(thisSample-steps_back))<.1 % width must be smaller than 100ms, also offset issue
                                        average_width_select(cp,cs) = thisData.tt(thisSample+steps_forw)-thisData.tt(thisSample-steps_back);
                                    end
                                end
                                clear tempResp
                            end
                        end
                    end

                    % if responses are added to average_ccep_select
                    if any(~isnan(average_ccep_select(:)))
                        thisResp = squeeze(mean(mean(average_ccep_select,1,'omitnan'),2,'omitnan'));
                        thisN1 = squeeze(mean(mean(average_N1_select,1,'omitnan'),2,'omitnan'));
                        thisWidth = squeeze(mean(mean(average_width_select,1,'omitnan'),2,'omitnan'));
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
                        average_width_run(rr1,rr2,ll) = thisWidth;
                    else
                        average_ccep_run_nonnorm(rr1,rr2,ll,:) = NaN(5*2048,1);
                        average_ccep_run(rr1,rr2,ll,:) = NaN(5*2048,1);
                        average_N1_run(rr1,rr2,ll) = NaN;
                        average_width_run(rr1,rr2,ll) = NaN;
                    end
                end
            end
            clear thisData thisRun thisN1
        end

        % average for this patient across these areas
        average_ccep_pat_nonnorm = squeeze(mean(average_ccep_run_nonnorm,3,'omitnan'));%[roi1, roi2,samples]
        average_ccep_pat = squeeze(mean(average_ccep_run,3,'omitnan'));%[roi1, roi2,samples]
        average_n1_pat = squeeze(mean(average_N1_run,3,'omitnan'));
        average_width_pat = squeeze(mean(average_width_run,3,'omitnan'));
        clear average_ccep_run average_N1_run

        if ~isempty(average_ccep_age{age})
            n = size(average_ccep_age{age},3);
            average_ccep_age_nonnorm{age}(:,:,n+1,:) = average_ccep_pat_nonnorm; % [roi1, roi2, subjects, samples]
            average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat; % [roi1, roi2, subjects, samples]
            average_n1_age{age}(:,:,n+1) = average_n1_pat;
            average_width_age{age}(:,:,n+1) = average_width_pat;
        else
            n = 0;
            average_ccep_age_nonnorm{age}(:,:,n+1,:) = average_ccep_pat_nonnorm;
            average_ccep_age{age}(:,:,n+1,:) = average_ccep_pat;
            average_n1_age{age}(:,:,n+1) = average_n1_pat;
            average_width_age{age}(:,:,n+1) = average_width_pat;
        end
        clear average_ccep_pat average_ccep_pat_nonnorm average_n1_pat
    end

    save(filename_averageCCEP_width,'average_ccep_age','average_ccep_age_nonnorm','average_n1_age','average_width_age','tt','roi','roi_name');
else
    load(filename_averageCCEP_width,'average_ccep_age','average_ccep_age_nonnorm','average_n1_age','average_width_age','tt','roi','roi_name');
end

%% 3. sort all cceps for each region to another region according to age

% this does not average any cceps when there are multiple subjects at the
% same age.

sortage = struct();

for rr1 = 1:4 % for each stimulation region
    for rr2 = 1:4 % for each response region
        sortage(rr1,rr2).average_ccep = [];
        sortage(rr1,rr2).average_ccep_nonnorm = [];
        sortage(rr1,rr2).age_ind = [];
        sortage(rr1,rr2).average_n1 = [];
        sortage(rr1,rr2).average_width = [];

        for age = 1:max([n1Latencies.age])
            if size(average_ccep_age{age},1)>0 % there are some subjects at this age
                addThis = squeeze(average_ccep_age{age}(rr1,rr2,:,:)); % normalized CCEP
                addThis(isnan(addThis(:,1)),:) = [];
                addThisNonnorm = squeeze(average_ccep_age_nonnorm{age}(rr1,rr2,:,:)); % nonnormalized CCEP
                addThisNonnorm(isnan(addThisNonnorm(:,1)),:) = [];
                addN1 = squeeze(average_n1_age{age}(rr1,rr2,:));
                addWidth = squeeze(average_width_age{age}(rr1,rr2,:));
                addWidth(isnan(addN1)) = [];
                addN1(isnan(addN1)) = [];
                if ~isempty(addThis) % there are subjects with electrodes on ROI
                    if size(addThis,1)>size(addThis,2), addThis = addThis'; addThisNonnorm = addThisNonnorm'; end
                    nr_subs = size(addThis,1);
                    sortage(rr1,rr2).average_ccep = [sortage(rr1,rr2).average_ccep; addThis];
                    sortage(rr1,rr2).average_ccep_nonnorm = [sortage(rr1,rr2).average_ccep_nonnorm; addThisNonnorm];
                    sortage(rr1,rr2).average_n1 = [sortage(rr1,rr2).average_n1; addN1];
                    sortage(rr1,rr2).average_width = [sortage(rr1,rr2).average_width; addWidth];
                    sortage(rr1,rr2).age_ind = [sortage(rr1,rr2).age_ind; zeros(nr_subs,1)+age];
                    clear addThis
                end
            end
        end
    end
end


%% figure of subplots for each stimulated and responding region:
% temporal, central, parietal, frontal
% with normalized CCEPs + N1 sorted by age

ttmin = 0.010;
ttmax = .100;

figure('Position',[0 0 600 300])
for rr1 = 1:4
    for rr2 = 1:4
        subplot(4,4,(rr1-1)*4+rr2),hold on
        plot(sortage(rr1,rr2).average_n1,sortage(rr1,rr2).average_width,'k.')
        hold on
        axis tight
    end
end

% calculate correlation and p
p_all = zeros(4,4);
r_all = zeros(4,4);
b_all = zeros(4,4);
for rr1 = 1:4
    for rr2 = 1:4
        n1_width = sortage(rr1,rr2).average_width;
        n1_latency = sortage(rr1,rr2).average_n1;
        [r,p] = corr(n1_width(~isnan(n1_width)),n1_latency(~isnan(n1_width)),'type','Spearman');
        p_all(rr1,rr2) = p;
        r_all(rr1,rr2) = r;
        b = regress(n1_width(~isnan(n1_width)),n1_latency(~isnan(n1_width)));
        b_all(rr1,rr2) = b;
    end
end

% FDR correction
p_vals = p_all(:);

m = length(p_vals);
[p_sort,p_ind] = sort(p_vals(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end
% figure,hold on,plot(thisVal),plot(p_sort,'r.'),title('Significant p-values after FDR correction')

% add significant 1:1 line indicating which subplots showed significant
% results after FDR corection
p_sig = p_all;
p_sig(p_ind) = p_sort<thisVal;
for rr1 = 1:4
    for rr2 = 1:4
        subplot(4,4,(rr1-1)*4+rr2),hold on
        if p_sig(rr1,rr2)==1 % significant!
            plot([0 .05],[0 .05],'r')
            % plot regression
            x = 0:.01:.05;
            plot(x,b_all(rr1,rr2)*x,'b')
        end
    end
end
figureName = fullfile(myDataPath.output,'derivatives','age',...
    ['All_widthVSlatency']);

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
