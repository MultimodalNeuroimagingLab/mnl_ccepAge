% 
%  Load the n1Latencies from the derivatives
%  This code is used to plot the normalized ccep of all patients in order of age. 
%  This figure is displayed as Figure 2 in the article.
%


%% 
%  1. Load the n1Latencies from the derivatives
%  We use this code both for analysis in the main script, and for checks with only
%  subjects in whom it is certain that 8mA is used for stimulation. 

clear
close all

selectPat = input('Would you like to include all patients, or only the ones for whom it is certain that 8mA was applied (supplemental material)? [all/8] ','s');

if strcmp(selectPat, 'all')
    select_amplitude = 0; % make this 8 for only 8mA
elseif strcmp(selectPat, '8')
    select_amplitude = 8;
else
    error('Answer to previous question is not recognized.')
end

myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

if select_amplitude == 0 
    if exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'n1Latencies_V2.mat'), 'file')
        % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
        load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'n1Latencies_V2.mat'), 'n1Latencies')

        filename_averageCCEP = fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'average_ccep_age.mat');
    else
        disp('Run scripts ccep02_loadN1.m and ccep03_addtracts.m first')
    end
elseif select_amplitude == 8 % only 8 mA
    if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_8ma.mat'),'file')
        % if the n1Latencies_8ma.mat was saved after ccep02_loadN1?, load the n1Latencies structure here
        load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_8ma.mat'),'n1Latencies8ma')

        n1Latencies = n1Latencies8ma;
        filename_averageCCEP = fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'average_ccep_age_8ma.mat');
    else
        disp('Run scripts ccep02_loadN1.m and ccep03_addtracts.m first')
    end
end


%% 
%  2. Load ccep responses and categorize into connections from stimulated region to responding regions
%  skips automatically to 3. if you ran and saved output before (takes ~5 mins to run)
%
%  the CCEPs are averaged for each run, and then averaged CCEPs per patient
%  are collected for all subjects. Multiple subjects with the same age are
%  collected for each age (average_ccep_age_nonnorm) and normalized
%  (average_ccep_age)

if ~exist(filename_averageCCEP, 'file')

    % categorize anatomical regions 
    %[roi_track,roi_name,roi] = ccep_categorizeAnatomicalRegions();
    rois = ccep_categorizeAnatomicalRegions();

    % retrieve a time vector in ms (this patient has fs-2048)
    tt = n1Latencies(2).run(1).tt;

    % initialize variables to store ...
    % (each cell represents a year of age, upto the highest age of all participants)
    for iTr = 1:length(rois)
        for iSubTr = 1:length(rois(iTr).sub_tract)
            average_ccep_age{iTr}{iSubTr}           = cell(max([n1Latencies.age]), 1);
            average_ccep_age_nonnorm{iTr}{iSubTr}   = cell(max([n1Latencies.age]), 1);
            average_n1_age{iTr}{iSubTr}             = cell(max([n1Latencies.age]), 1);
        end
    end
    
    % loop over the subjects
    for iSubj = 1:size(n1Latencies, 2)
        fprintf('Load subj %d of %d \n', iSubj, size(n1Latencies, 2))
        age = n1Latencies(iSubj).age;
        
        % create the cells and empty matrices to store CCEP traces and N1 quantifications in, and eventually contain the average
        for iTr = 1:length(rois)
            for iSubTr = 1:length(rois(iTr).sub_tract)
                average_ccep_run{iTr}{iSubTr}           = NaN(2, size(n1Latencies(iSubj).run, 2), 5 * 2048);    % [roiDir, run, tt]
                average_ccep_run_nonnorm{iTr}{iSubTr}   = NaN(2, size(n1Latencies(iSubj).run, 2), 5 * 2048);    % [roiDir, run, tt]
                average_N1_run{iTr}{iSubTr}             = NaN(2, size(n1Latencies(iSubj).run, 2));              % [roiDir, run]
            end
        end
        
        % get the session directory name
        sesDir = dir(fullfile(myDataPath.output, 'derivatives', 'av_ccep', n1Latencies(iSubj).id, 'ses-*'));
        sesDir = sesDir.name;

        % loop over the runs
        for iRun = 1:size(n1Latencies(iSubj).run, 2)
            
            % load the run data
            runData = load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', n1Latencies(iSubj).id, sesDir, n1Latencies(iSubj).run(iRun).runName));
            
            % we will resample all cceps to 2048Hz, no need to preselect
            
            % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
            for iTr = 1:length(rois)
                for iSubTr = 1:length(rois(iTr).sub_tract)

                    % use both directions (consider stim on roi1 and response on roi2, and vice versa)
                    for iDir = [false true]
                        stimRoi = rois(iTr).sub_tract(iSubTr).(['roi', num2str(iDir + 1)]);
                        respRoi = rois(iTr).sub_tract(iSubTr).(['roi', num2str(~iDir + 1)]);
                        
                        % find stimulation pairs within specific region
                        stimPairs = find(sum(ismember(str2double(n1Latencies(iSubj).run(iRun).average_ccep_DestrieuxNr), stimRoi), 2) > 0);
                        
                        % find response electrodes within specific region
                        respChan = find(ismember(str2double(n1Latencies(iSubj).run(iRun).channel_DestrieuxNr), respRoi) > 0);
                        
                        
                        %
                        % collect all signals with stimulation pair and response electrode within specific region
                        %

                        average_ccep_select  = NaN(size(stimPairs, 1), size(respChan, 1), size(runData.average_ccep, 3));
                        average_N1_select    = NaN(size(stimPairs, 1), size(respChan, 1)); % N1latency
                        for cp = 1:size(stimPairs, 1)
                            for cs = 1:size(respChan, 1)
                                if ~isnan(n1Latencies(iSubj).run(iRun).n1_peak_sample(respChan(cs),stimPairs(cp)))
                                    if runData.tt(n1Latencies(iSubj).run(iRun).n1_peak_sample(respChan(cs), stimPairs(cp))) == 0, disp('what???'), end
                                    tempResp = squeeze(runData.average_ccep(respChan(cs),stimPairs(cp), :));
                                    average_ccep_select(cp, cs, :) = tempResp;
                                    average_N1_select(cp, cs) = runData.tt(n1Latencies(iSubj).run(iRun).n1_peak_sample(respChan(cs), stimPairs(cp)));
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
                            average_ccep_run{iTr}{iSubTr}(iDir + 1, iRun, :)           = thisResp;   %[roiDir, run, samples]
                            average_ccep_run_nonnorm{iTr}{iSubTr}(iDir + 1, iRun, :)   = thisResp2;  %[roiDir, run, samples]
                            average_N1_run{iTr}{iSubTr}(iDir + 1, iRun)                = thisN1;     %[roiDir, run]

                        else
                            
                            average_ccep_run{iTr}{iSubTr}(iDir + 1, iRun, :)           = NaN(5 * 2048, 1);
                            average_ccep_run_nonnorm{iTr}{iSubTr}(iDir + 1, iRun, :)   = NaN(5 * 2048, 1);
                            average_N1_run{iTr}{iSubTr}(iDir + 1, iRun)                = NaN;

                        end
                        clear average_ccep_select average_N1_select;
                    
                    end
                    
                end
            end
            
            clear runData thisN1;
        end

        
        %
        % Average over the runs and store in an 'age' specific cell
        %
        
        % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
        for iTr = 1:length(rois)
            %if ~iscell(average_ccep_age{age}),  average_ccep_age{age} = cell(1, length(rois));  end
            
            for iSubTr = 1:length(rois(iTr).sub_tract)
                %if ~iscell(average_ccep_age{age}{iTr}),  average_ccep_age{age}{iTr} = cell(1, length(rois(iTr).sub_tract));  end
                
                % average over runs
                average_ccep_pat          = squeeze(mean(average_ccep_run{iTr}{iSubTr}, 2, 'omitnan'));         %[roiDir, samples]
                average_ccep_pat_nonnorm  = squeeze(mean(average_ccep_run_nonnorm{iTr}{iSubTr}, 2, 'omitnan')); %[roiDir, samples]
                average_n1_pat            = squeeze(mean(average_N1_run{iTr}{iSubTr}, 2, 'omitnan'));           %[roiDir]
                
                % 
                n = 0;
                if ~isempty(average_ccep_age{iTr}{iSubTr}{age})
                    n = size(average_ccep_age{iTr}{iSubTr}{age}, 2);
                end
                average_ccep_age{iTr}{iSubTr}{age}(:, n + 1, :)         = average_ccep_pat;          % [roiDir, subjects, samples]
                average_ccep_age_nonnorm{iTr}{iSubTr}{age}(:, n + 1, :) = average_ccep_pat_nonnorm;  % [roiDir, subjects, samples]
                average_n1_age{iTr}{iSubTr}{age}(:, n + 1)              = average_n1_pat;            % [roiDir, subjects]
                
            end
        end
        clear average_ccep_run average_ccep_run_nonnorm average_N1_run
        clear average_ccep_pat average_ccep_pat_nonnorm average_n1_pat ROIsDist ROIsTrcs;
    end

    save(filename_averageCCEP, 'average_ccep_age', 'average_ccep_age_nonnorm', 'average_n1_age', 'tt', 'rois');
else
    load(filename_averageCCEP, 'average_ccep_age', 'average_ccep_age_nonnorm', 'average_n1_age', 'tt', 'rois');
end



%%
%  3. sort all cceps for each region to another region according to age
%
%   this does not average any cceps when there are multiple subjects at the same age.
%

sortage = {};

for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        sortage{iTr}{iSubTr} = struct();
        
        for age = 1:max([n1Latencies.age])
            if ~isempty(average_ccep_age{iTr}{iSubTr}{age})
                
                for iDir = [false true]
                    
                    sortage{iTr}{iSubTr}(iDir + 1).average_ccep = [];
                    sortage{iTr}{iSubTr}(iDir + 1).average_ccep_nonnorm = [];
                    sortage{iTr}{iSubTr}(iDir + 1).age_ind = [];
                    sortage{iTr}{iSubTr}(iDir + 1).average_n1 = [];

                    addThis = squeeze(average_ccep_age{iTr}{iSubTr}{age}(iDir + 1, :, :))'; % normalized CCEP
                    addThis(isnan(addThis(:,1)),:) = [];
                    addThisNonnorm = squeeze(average_ccep_age_nonnorm{iTr}{iSubTr}{age}(iDir + 1, :, :))'; % nonnormalized CCEP
                    addThisNonnorm(isnan(addThisNonnorm(:,1)),:) = [];
                    addN1 = squeeze(average_n1_age{iTr}{iSubTr}{age}(iDir + 1, :))';
                    addN1(isnan(addN1)) = [];

                    if ~isempty(addThis) % there are subjects with electrodes on ROI
                        if size(addThis, 1) > size(addThis, 2), addThis = addThis'; addThisNonnorm = addThisNonnorm'; end
                        nr_subs = size(addThis, 1);

                        sortage{iTr}{iSubTr}(iDir + 1).average_ccep = [sortage{iTr}{iSubTr}(iDir + 1).average_ccep; addThis];
                        sortage{iTr}{iSubTr}(iDir + 1).average_ccep_nonnorm = [sortage{iTr}{iSubTr}(iDir + 1).average_ccep_nonnorm; addThisNonnorm];
                        sortage{iTr}{iSubTr}(iDir + 1).average_n1 = [sortage{iTr}{iSubTr}(iDir + 1).average_n1; addN1];
                        sortage{iTr}{iSubTr}(iDir + 1).age_ind = [sortage{iTr}{iSubTr}(iDir + 1).age_ind; zeros(nr_subs,1)+age];
                        clear addThis addThisNonnorm addN1
                    end

                end
            end
        end
    end
end


%% 
%  figure of subplots for each stimulated and responding region:
%  temporal, central, parietal, frontal
%  with normalized CCEPs + N1 sorted by age

for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        p_all{iTr}{iSubTr} = zeros(1, 2);
        r_all{iTr}{iSubTr} = zeros(1, 2);
        
        for iDir = [false true]
            
            n1_latency = 1000 * sortage{iTr}{iSubTr}(iDir + 1).average_n1;
            age = sortage{iTr}{iSubTr}(iDir + 1).age_ind;
            [r, p] = corr(n1_latency, age, 'type', 'Spearman');
            p_all{iTr}{iSubTr}(iDir + 1) = p;
            r_all{iTr}{iSubTr}(iDir + 1) = r;
            
        end
    end
end


ttmin = 0.010;
ttmax = .100;

for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        for iDir = [false true]
            
            %
            figure('Position',[0 0 600 300])
            hold on;
            imagesc(1000 * tt(tt > ttmin & tt < ttmax), 1:length(sortage{iTr}{iSubTr}(iDir + 1).age_ind), ...
                    -sortage{iTr}{iSubTr}(iDir + 1).average_ccep(:, tt > ttmin & tt < ttmax), ...
                    [-0.1 0.1]);
            
            plot(1000 * sortage{iTr}{iSubTr}(iDir + 1).average_n1, 1:length(sortage{iTr}{iSubTr}(iDir + 1).age_ind), 'k.')
            colormap(parula)
            hold on
            set(gca,'XTick',20:20:80,'YTick',[])
            axis tight
  
        end
    end
end

%%

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
for iSubj = 1:length(p_sort)
    thisVal(iSubj) = (iSubj/m)*0.05;
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

