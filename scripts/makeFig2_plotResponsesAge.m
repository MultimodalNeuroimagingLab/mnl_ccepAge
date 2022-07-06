%% load the n1Latencies from the derivatives
% This code is used to plot the normalized ccep of all patients in order of
% age. This figure is displayed as Figure 2 in the article.


%% 1. load the n1Latencies from the derivatives
% we use this code both for analysis in the main script, and for checks
% with only subjects in whom it is certain that 8mA is used for
% stimulation. 

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

myDataPath = setLocalDataPath(1);

if select_amplitude == 0 
    if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
        % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
        load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')

        filename_averageCCEP = fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age.mat');
    else
        disp('Run first script ccep02_loadN1.m')
    end
elseif select_amplitude == 8 % only 8 mA
    if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_8ma.mat'),'file')
        % if the n1Latencies_8ma.mat was saved after ccep02_loadN1?, load the n1Latencies structure here
        load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_8ma.mat'),'n1Latencies8ma')

        n1Latencies = n1Latencies8ma;
        filename_averageCCEP = fullfile(myDataPath.output,'derivatives','av_ccep','average_ccep_age_8ma.mat');
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


if ~exist(filename_averageCCEP, 'file')

    % categorize anatomical regions
    ccep_categorizeAnatomicalRegions

    % retrieve a time vector in ms (this patient has fs-2048)
    tt = n1Latencies(2).run(1).tt;

    % initialize variables to store ...
    % (each cell represents a year of age, upto the highest age of all participants)
    average_ccep_age            = cell(max([n1Latencies.age]), 1);
    average_ccep_age_nonnorm    = cell(max([n1Latencies.age]), 1);
    average_n1_age              = cell(max([n1Latencies.age]), 1);
    
    % loop over the subjects
    for kk = 1:size(n1Latencies, 2)
        fprintf('Load subj %d of %d \n', kk, size(n1Latencies, 2))
        age = n1Latencies(kk).age;
        
        % variables to store the distance and included MNI tracks (file & line indices)
        % between a stimulation end-point ROI and a response end-point ROI
        ROIsDist                    = cell(size(roi, 2), size(roi, 2));          % [roi_start, roi_end]
        ROIsTrcs                    = cell(size(roi, 2), size(roi, 2));          % [roi_start, roi_end]

        % variables to collect the seperate CCEP responses and eventually contain the average
        average_ccep_run         = NaN(size(roi, 2), size(roi, 2), size(n1Latencies(kk).run, 2), 5 * 2048); % [roi_start, roi_end, run, tt]
        average_ccep_run_nonnorm = NaN(size(roi, 2), size(roi, 2), size(n1Latencies(kk).run, 2), 5 * 2048); % [roi_start, roi_end, run, tt]
        average_N1_run           = NaN(size(roi, 2), size(roi, 2), size(n1Latencies(kk).run, 2));           % [roi_start, roi_end, run]
        
        % get the session directory name
        sesDir = dir(fullfile(myDataPath.output, 'derivatives', 'av_ccep', n1Latencies(kk).id, 'ses-*'));
        sesDir = sesDir.name;

        % loop over the runs
        for iRun = 1:size(n1Latencies(kk).run, 2)
            
            % get this run file name
            thisRun = fullfile(myDataPath.output, 'derivatives', 'av_ccep', n1Latencies(kk).id, sesDir, n1Latencies(kk).run(iRun).runName);
            thisData = load(thisRun);
            
            % we will resample all cceps to 2048Hz, no need to preselect
            
            % loop over the "end-point" areas, as stimulation ROIs
            for rr1 = 1:size(roi, 2)
                
                % find stimulation pair within specific region
                chanPair = find(sum(contains(n1Latencies(kk).run(iRun).average_ccep_DestrieuxNr, roi{rr1}), 2) > 0);
                
                % loop over the "end-point" areas, as response ROIs
                for rr2 = 1:size(roi, 2)
                    
                    % find response electrode within specific region
                    chanSig = find(ismember(string(n1Latencies(kk).run(iRun).channel_DestrieuxNr), roi{rr2}) > 0);

                    
                    %
                    % collect all signals with stimulation pair and response electrode within specific region
                    %
                    
                    average_ccep_select  = NaN(size(chanPair, 1), size(chanSig, 1), size(thisData.average_ccep, 3));
                    average_N1_select    = NaN(size(chanPair, 1), size(chanSig, 1)); % N1latency
                    for cp = 1:size(chanPair, 1)
                        for cs = 1:size(chanSig, 1)
                            if ~isnan(n1Latencies(kk).run(iRun).n1_peak_sample(chanSig(cs),chanPair(cp)))
                                if thisData.tt(n1Latencies(kk).run(iRun).n1_peak_sample(chanSig(cs), chanPair(cp))) == 0, disp('what???'), end
                                tempResp = squeeze(thisData.average_ccep(chanSig(cs),chanPair(cp), :));
                                average_ccep_select(cp, cs, :) = tempResp;
                                average_N1_select(cp, cs) = thisData.tt(n1Latencies(kk).run(iRun).n1_peak_sample(chanSig(cs), chanPair(cp)));
                                clear tempResp
                            end
                        end
                    end

                    % if responses are added to average_ccep_select
                    if any(~isnan(average_ccep_select(:)))
                        thisResp = squeeze(mean(mean(average_ccep_select, 1, 'omitnan'), 2, 'omitnan'));
                        thisN1   = squeeze(mean(mean(average_N1_select,   1, 'omitnan'), 2, 'omitnan'));
                        
                        % upsample if necessary
                        if length(thisResp) ~= 5 * 2048
                            thisResp = resample(thisResp, 5 * 2048, length(thisResp));
                        end

                        thisResp2 = thisResp; % non-normalized

                        % 'normalize' 0.15-0.100 by unit length
                        thisResp = thisResp ./ sqrt(sum(thisResp(tt > .015 & tt < 0.100) .^ 2));
                        
                        % store the ...
                        average_ccep_run_nonnorm(rr1, rr2, iRun, :) = thisResp2;  %[roi1, roi2, run, samples]
                        average_ccep_run(rr1, rr2, iRun, :)         = thisResp;   %[roi1, roi2, run, samples]
                        average_N1_run(rr1, rr2, iRun)              = thisN1;     %[roi1, roi2, run]
                        
                    else
                        
                        average_ccep_run_nonnorm(rr1, rr2, iRun, :) = NaN(5 * 2048, 1);
                        average_ccep_run(rr1, rr2, iRun, :)         = NaN(5 * 2048, 1);
                        average_N1_run(rr1, rr2, iRun)              = NaN;
                        
                    end
                    
                    
                    %
                    % collect tract distance
                    %
                    
                    % equal for all runs, so only perform on the first run
                    if iRun == 1

                        % check if the stimulus and response end-point area concern the same tract
                        % and the stimulation and response end-point ROIs are not the same
                        if strcmpi(roi_track{rr1}, roi_track{rr2}) && rr1 ~= rr2
                            
                            trkFolder = '/Users/m218483/Documents/leadDBS/tracks/';
                            trkSets = {fullfile(trkFolder, [roi_track{rr1}, '_L.trk.gz']), fullfile(trkFolder, [roi_track{rr1}, '_R.trk.gz'])};
                            disp('Retrieving tract distance');
                            
                            % retrieve the distance between the stimulation and response end-point areas
                            % for this particular patient
                            [thisDist, thisTrcs] = ccep_retrieveInterROIDistance( ...
                                            trkSets, ...
                                            fullfile(myDataPath.input, 'derivatives', 'coreg_ANTs', n1Latencies(kk).id), ...
                                            roi{rr1}, roi{rr2});
                            
                            % TODO: store only if found. So if only LH
                            % electrodes, than RH tract file should not
                            % find distance/tract-lines
                                        
                            % store the distances and included MNI tract-lines
                            ROIsDist{rr1, rr2}            = thisDist;
                            ROIsTrcs{rr1, rr2}            = {trkSets, thisTrcs};
                            
                        end 
                    end
                      
                    
                end
            end
            clear thisData thisRun thisN1 thisLen;
        end

        % average for this patient across these areas
        average_ccep_pat_nonnorm  = squeeze(mean(average_ccep_run_nonnorm, 3, 'omitnan'));  %[roi1, roi2,samples]
        average_ccep_pat          = squeeze(mean(average_ccep_run, 3, 'omitnan'));          %[roi1, roi2,samples]
        average_n1_pat            = squeeze(mean(average_N1_run, 3, 'omitnan'));
        clear average_ccep_run_nonnorm average_ccep_run average_N1_run

        if ~isempty(average_ccep_age{age})
            n = size(average_ccep_age{age}, 3);
            average_ccep_age_nonnorm{age}(:, :, n + 1, :) = average_ccep_pat_nonnorm;       % [roi1, roi2, subjects, samples]
            average_ccep_age{age}(:, :, n + 1, :)         = average_ccep_pat;               % [roi1, roi2, subjects, samples]
            average_n1_age{age}(:, :, n + 1)              = average_n1_pat;
        else
            n = 0;
            average_ccep_age_nonnorm{age}(:, :, n + 1, :) = average_ccep_pat_nonnorm;
            average_ccep_age{age}(:, : ,n + 1, :)         = average_ccep_pat;
            average_n1_age{age}(:, : ,n + 1)              = average_n1_pat;
        end
        clear average_ccep_pat average_ccep_pat_nonnorm average_n1_pat
    end

    save(filename_averageCCEP, 'average_ccep_age', 'average_ccep_age_nonnorm', 'average_n1_age', 'tt', 'roi', 'roi_name');
else
    load(filename_averageCCEP, 'average_ccep_age', 'average_ccep_age_nonnorm', 'average_n1_age', 'tt', 'roi', 'roi_name');
end

%% 3. sort all cceps for each region to another region according to age

% this does not average any cceps when there are multiple subjects at the
% same age.

sortage = struct();

for rr1 = 1:4 % for each stimulation region
    for rr2 = 1:4 % for each response region
        sortage(rr1, rr2).average_ccep = [];
        sortage(rr1, rr2).average_ccep_nonnorm = [];
        sortage(rr1, rr2).age_ind = [];
        sortage(rr1, rr2).average_n1 = [];
        
        for age = 1:max([n1Latencies.age])
            if size(average_ccep_age{age},1)>0 % there are some subjects at this age
                addThis = squeeze(average_ccep_age{age}(rr1,rr2,:,:)); % normalized CCEP
                addThis(isnan(addThis(:,1)),:) = [];
                addThisNonnorm = squeeze(average_ccep_age_nonnorm{age}(rr1,rr2,:,:)); % nonnormalized CCEP
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
        imagesc(1000*tt(tt>ttmin & tt< ttmax),1:length(sortage(rr1,rr2).age_ind),-sortage(rr1,rr2).average_ccep(:,tt>ttmin & tt< ttmax),...
            [-0.1 0.1])
        plot(1000*sortage(rr1,rr2).average_n1,1:length(sortage(rr1,rr2).age_ind),'k.')
        colormap(parula)
        hold on
        set(gca,'XTick',20:20:80,'YTick',[])
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
        [r,p] = corr(n1_latency,age,'type','Spearman');
        p_all(rr1,rr2) = p;
        r_all(rr1,rr2) = r;
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

% add significant stars indicating which subplots showed significant
% results after FDR corection
p_sig = p_all;
p_sig(p_ind) = p_sort<thisVal;
for rr1 = 1:4
    for rr2 = 1:4
        subplot(4,4,(rr1-1)*4+rr2),hold on
        if p_sig(rr1,rr2)==1 % significant!
            plot(100,0,'r*')
        end
    end
end

if select_amplitude==0
    figureName = fullfile(myDataPath.output,'derivatives','age',...
        ['AllSortAge_tmax' int2str(ttmax*1000)]);
elseif select_amplitude==8
    figureName = fullfile(myDataPath.output,'derivatives','age',...
        ['AllSortAge_tmax' int2str(ttmax*1000), '_8mA']);
end

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-depsc','-r300',figureName)


%% figure with colormap

figure('Position',[0 0 150 40])
imagesc(1:100)
colormap(parula)
axis off
figureName = fullfile(myDataPath.output,'derivatives','age',...
    ['AllSortAge_tmax' int2str(ttmax*1000) '_cm']);

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)

