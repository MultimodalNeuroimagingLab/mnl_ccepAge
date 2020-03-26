% This file contains the scripts to create the figures used in the final
% report. 

% Jaap van der Aar, aug. 2019 
% Age-related latency in signal propagation in children and adults 


% Matrices that can be used for plotting:
% for all latencies: 
% - matrix_reshape_all
% for all within ROI latencies:
% - ROI_plot_matrix_15
% - ROI_plot_matrix_26
% - ROI_plot_matrix_38
% - ROI_within_plot_matrix_all
% or for between ROI latencies:
% - ROI_between_plot_matrix

% - age_vector


addpath(genpath('/Fridge/users/jaap/github/ccep/'))

%% Figure 1: visualisation of how the detection algorithm works

subj = 1;
runs = 1;
ii = 3;
jj = 1;
minSD = 50;
n1_samples_start = find(tt>0.009,1);
n1_samples_end = find(tt>0.1,1);
thresh = 3.6;

baseline_tt = tt>-1 & tt<-.01;
signal_median = median(database(subj).metadata(runs).epoched_data_avg(ii,jj,baseline_tt),3);

% subtract median baseline from signal
new_signal = squeeze(database(subj).metadata(runs).epoched_data_avg(ii,jj,:)) - signal_median;
% testplot new signal: plot(tt,squeeze(new_signal))

% take area before the stimulation of the new signal and calculate its SD
pre_stim_sd_orig = std(new_signal(baseline_tt));

% if the pre_stim_sd is smaller that the minimally needed SD, use this
% the minSD as pre_stim_sd
if pre_stim_sd_orig < minSD
    pre_stim_sd = minSD * thresh;
else 
    pre_stim_sd =  pre_stim_sd_orig * thresh;
end

% display figure in specific size
fig10 = figure(10); set(fig10,'Position',[0 0 1000 400]), hold on


% plot relevant part of the signal (1s before till 500ms after stimulation)

plot(tt(tt>-1 & tt<.5),squeeze(new_signal(tt>-1 & tt<.5)),'Color',[0 0 0],'LineWidth',2)

xlabel('time (in sec)')
ylabel('amplitude (in uV)')
title('Visualisation of detection algorithm')
set(gca, 'FontName','arial','FontSize',20)

ylim([-1200 1200])
xlim([-1 .5])
box off

hold on
% plot stimulation artifact (first 20 samples, +/- 9ms)
plot(tt(5120:5139),squeeze(new_signal(5120:5139)),'r','LineWidth',2)
% plot found peaks and onsets
% plot(tt(output_ER_all(ii,jj,1)),output_ER_all(ii,jj,2),'b*')
% plot(tt(output_ER_all(ii,jj,3)),output_ER_all(ii,jj,4),'r*')
% plot(tt(output_ER_all(ii,jj,5)),output_ER_all(ii,jj,6),'g*')
% plot(tt(output_ER_all(ii,jj,7)),output_ER_all(ii,jj,8),'y*')

% plot calculated baseline standard deviation
plot(tt(baseline_tt), pre_stim_sd_orig+zeros(size(tt(baseline_tt))), 'm-','LineWidth',2)
plot(tt(baseline_tt), -pre_stim_sd_orig+zeros(size(tt(baseline_tt))), 'm-','LineWidth',2)

% plot adjusted baseline (when calculated < minSD)
plot(tt(baseline_tt), (pre_stim_sd)+zeros(size(tt(baseline_tt))), 'g-','LineWidth',2)
plot(tt(baseline_tt), (-pre_stim_sd)+zeros(size(tt(baseline_tt))), 'g-','LineWidth',2)

% plot lines with timeframe of n1 peak
plot([tt(n1_samples_start) tt(n1_samples_start)],[-2000 -180], 'b-','LineWidth',2)
plot([tt(n1_samples_end) tt(n1_samples_end)],[-2000 -180], 'b-','LineWidth',2)
plot(tt(n1_samples_start:n1_samples_end), (-pre_stim_sd)+zeros(size(tt(n1_samples_start:n1_samples_end))), 'b-','LineWidth',2)


plot(tt(database(subj).metadata(runs).n1_peak_sample(ii,jj)), database(subj).metadata(runs).n1_peak_amplitude(ii,jj), ...
    'o', 'LineWidth',4, 'MarkerSize',8,'MarkerEdgeColor', 'c');


hold off

%% Figure 2: Matrix with number of connections in data - finished

% First, create a connectivity matrix by running over all recordings to see
% which stimulations are applied, and on which regions there are recordings

% how many connections to each electrode would be possible 
connectivity_mat = zeros(75,75);

% iterate over all subjects in database
for subj = 1:length(database)
    % iterate over all their runs
    for runs = 1:length(database(subj).metadata)
        % for the SPESclin runs
        if strcmp(database(subj).metadata(runs).task, 'SPESclin')
            % run through all stimulations
            for stims = 1:length(database(subj).metadata(runs).stimulated_pairs)

                % add one to every recorded channel except for the ones that
                % are stimulated (e.g. add +1 to channel 3 till 64 one the
                % places (stim1,3:64) and (stim2,3:64). But then converted
                % to the brain re_jgions

                % find locations of the stimulated pairs
                stimnum_1 = database(subj).metadata(runs).stimulated_pairs(stims,1);
                stimnum_2 = database(subj).metadata(runs).stimulated_pairs(stims,1);
                
                if stimnum_1 <= length(database(subj).total_grid) && stimnum_2 <= length(database(subj).total_grid)
               
                    stimloc_1 = str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label(stimnum_1));            
                    stimloc_2 = str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label(stimnum_2));

                    % for all recorded channels
                    for zz = 1:length(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)
                        % if the row does not correspond with the stimulated
                        % electrodes, add +1 to the matrix
                        if zz ~= stimnum_1 || zz ~= stimnum_2

                            rec_loc = str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label(zz));
                            
                            if ~isnan(rec_loc) && ~isnan(stimloc_1) && ~isnan(stimloc_2) 
                                connectivity_mat(stimloc_1,rec_loc) = connectivity_mat(stimloc_1,rec_loc) + 1; 
                                connectivity_mat(stimloc_2,rec_loc) = connectivity_mat(stimloc_2,rec_loc) + 1; 
                            end

                        end
                    end              
                end               
            end
        end
    end
end

% currently, there is still a mistake in the code. The atlas starts
% labeling 1 later then how it is assigned in the electrodes tsv.
% in other code, there is correction for this problem, but in the figure
% making not, therefore, delete first row (destrieux should have 74 labels)
connectivity_mat = connectivity_mat(2:75,2:75);

% visualize matrix
fig61 = figure(61); set(fig61,'Position',[0 0 800 600])
imagesc(connectivity_mat,[0 max(connectivity_mat(:))])
axis square
title('Number of averaged epochs')
xlabel('Recorded region')
ylabel('Stimulated region')
set(gca,'FontName','arial','FontSize',20, 'XTick',[10:10:70],...
    'XDir','normal'); xtickangle(90)

% use color matrix 'hot', but make all combinations with > 1000 avg. epochs
% all white
cm = hot(1000);
cm = [0 0 0; cm];
cm(1000:max(connectivity_mat(:)),1) = [1];
cm(1000:max(connectivity_mat(:)),2) = [1];
cm(1000:max(connectivity_mat(:)),3) = [1];
colormap(cm)
hcb = colorbar;
set(gcf,'PaperPositionMode','auto')

hold off


%% Figure 3: Rendering of ROIs - finished

dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';
% add vistasoft for read_annotation
addpath('/home/jaap/vistasoft/external/freesurfer');
addpath('/Fridge/users/jaap/github/ecogBasicCode/render');

subjects = {'RESP0768'};
sessions = {'1'};
hemi_cap = {'R'}; 
hemi_small = {'r'};

v_dirs = [90 0];
s = 1;

% subject code
subj = subjects{s};
ses_label = sessions{s};

% gifti file name:
dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
% surface labels
surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'label',...
    [hemi_small{s} 'h.aparc.a2009s.annot']);
% surface_labels = MRIread(surface_labels_name);
[vertices, label, colortable] = read_annotation(surface_labels_name);
vert_label = label; 

for kk = 1:size(colortable.table,1) 
    vert_label(label==colortable.table(kk,5)) = kk;
end

% make manually a colormap which only highlights the selected ROIs
cmap = [repmat([.3 .3 .3],76,1)];
cmap(16,1:3) = [.6 .1 .1];
cmap(27,1:3) = [.1 .6 .1];
cmap(39,1:3) = [.1 .1 .6];

% electrode locations name:
dataLocName = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_electrodes.tsv'];
% load electrode locations
loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
elecmatrix = [loc_info.x loc_info.y loc_info.z];

% load gifti:
g = gifti(dataGiiName);

% figure with rendering for different viewing angles
for k = 1:size(v_dirs,1) % loop across viewing angles
    v_d = v_dirs(k,:);
    
    figure(71)
    ecog_RenderGiftiLabels(g,vert_label,cmap,colortable.struct_names)
    ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle
    

    set(gcf,'PaperPositionMode','auto')

end

title('Locations of the three ROIs') 
set(gca,'FontName','arial','FontSize',16)


%% Figure 4: The ROC-curves - finished

% First, load the validation matrices of all validated subjects,
% and calculate the averaged scores (variable: averaged_parameter_scores) 

% The validation_matrix and parameters_optimalize_mat of all validated data
% are loaded (the variables in the mat-files are renamed so they do not overwrite)
load(fullfile(top_path,'validation', 'sub-RESP0458_ses-1_run-011714_parameters_optimalization.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0458_ses-1_run-011714_validation_matrix.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0468_ses-1_run-031729_parameters_optimalization.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0468_ses-1_run-031729_validation_matrix.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0621_ses-1_run-021147_parameters_optimalization.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0621_ses-1_run-021147_validation_matrix.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0706_ses-1_run-041501_parameters_optimalization.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0706_ses-1_run-041501_validation_matrix.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0733_ses-1b_run-050941_parameters_optimalization.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0733_ses-1b_run-050941_validation_matrix.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0768_ses-1_run-021704_parameters_optimalization.mat'))
load(fullfile(top_path,'validation', 'sub-RESP0768_ses-1_run-021704_validation_matrix.mat'))


% The scores of the different datasets are based on a different amount of averaged epochs
% To calculate the average scores, it first needs recalculation.
total_validated_epochs = (size(validation_matrix_0458,1)*size(validation_matrix_0458,2)) ...
    + (size(validation_matrix_0468,1)*size(validation_matrix_0468,2)) ...
    + (size(validation_matrix_0621,1)*size(validation_matrix_0621,2)) ...
    + (size(validation_matrix_0706,1)*size(validation_matrix_0706,2)) ...
    + (size(validation_matrix_0733,1)*size(validation_matrix_0733,2)) ... 
    + (size(validation_matrix_0768,1)*size(validation_matrix_0768,2));

% multiply each parameters_optimalize_mat by the fraction on total epochs
parameters_optimalize_mat_0458 = parameters_optimalize_mat_0458 * ...
    (size(validation_matrix_0458,1)*size(validation_matrix_0458,2)/total_validated_epochs);
parameters_optimalize_mat_0468 = parameters_optimalize_mat_0468 * ...
    (size(validation_matrix_0468,1)*size(validation_matrix_0468,2)/total_validated_epochs);
parameters_optimalize_mat_0621 = parameters_optimalize_mat_0621 * ...
    (size(validation_matrix_0621,1)*size(validation_matrix_0621,2)/total_validated_epochs);
parameters_optimalize_mat_0706 = parameters_optimalize_mat_0706 * ...
    (size(validation_matrix_0706,1)*size(validation_matrix_0706,2)/total_validated_epochs);
parameters_optimalize_mat_0733 = parameters_optimalize_mat_0733 * ...
    (size(validation_matrix_0733,1)*size(validation_matrix_0733,2)/total_validated_epochs); 
parameters_optimalize_mat_0768 = parameters_optimalize_mat_0768 * ...
    (size(validation_matrix_0768,1)*size(validation_matrix_0768,2)/total_validated_epochs); 

% add scores to get 1 averaged matrix with parameter scores
averaged_parameter_scores = parameters_optimalize_mat_0458 + parameters_optimalize_mat_0468 + parameters_optimalize_mat_0621 ...
    + parameters_optimalize_mat_0706 + parameters_optimalize_mat_0733 + parameters_optimalize_mat_0768;

% load detection range and other variables 
thresh = [1:.2:5];   
subj = 1;
runs = 1;
tt = [1:database(subj).metadata(runs).epoch_length*database(subj).metadata(runs).data_hdr.Fs] / ...
    database(subj).metadata(runs).data_hdr.Fs - database(subj).metadata(runs).epoch_prestim_length;
n1_samples_end = [find(tt>0.04,1), find(tt>0.05,1), find(tt>0.06,1), find(tt>0.07,1), find(tt>0.08,1), find(tt>0.09,1), find(tt>0.1,1), find(tt>0.11,1)]; 


fig21 = figure(21); set(fig21,'Position',[0 0 1200 620])

% plot averaged ROCs
subplot(1,4,1:3)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(averaged_parameter_scores,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(averaged_parameter_scores,2)

    sens_plot = averaged_parameter_scores(:,time_th,1);
    spes_plot = averaged_parameter_scores(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1 - specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
axis square
title('ROC-curves of algorithm performance')
set(gca,'XTick',[0:.5:1],'YTick',[0:.5:1],'FontName','arial','FontSize',18)

% plot star at optimum
plot(((100-95.51)/100),(75.40/100),'p', 'MarkerFaceColor','black','MarkerSize',15,'MarkerEdgeColor', 'black');

% plot legenda 
time_end = round(tt(n1_samples_end),2) * 1000;

subplot(1,4,4),hold on
for time_th = 1:size(averaged_parameter_scores,2)
    for clrbar=1:75
        
        plot(clrbar,time_end(time_th),'.','MarkerSize',20,'Color',my_colors(time_th,:))
    end
end

xlim([-4 80]), set(gca,'XTick',[])
ylim([35,115]), ylabel('time end (ms)')
ttl = title('Detection range'); 
box on
set(gca,'FontName','arial','FontSize',18)
set(gcf,'PaperPositionMode','auto')


% plot figure with zoomed on part of the ROC-curves 
fig23 = figure(23);
plot([0 1],[0 1],'k'),hold on
% use jet as colors
my_colors = jet(size(averaged_parameter_scores,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(averaged_parameter_scores,2)

    sens_plot = averaged_parameter_scores(:,time_th,1);
    spes_plot = averaged_parameter_scores(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:),'LineWidth', 1.5)
end

%plot(1-averaged_parameter_scores_corrected(optimum_mat_x,optimum_mat_y,2),...
 %   averaged_parameter_scores_corrected(optimum_mat_x,optimum_mat_y,1),'k*','MarkerSize',10)
% title('Averaged ROC-curves')
% 
% xlabel('1 - specificity')
% ylabel('sensitivity')
xlim([0 .3]),ylim([0.7 1])
axis square

set(gcf,'PaperPositionMode','auto')
set(gca,'XTick',[0:.05:.3],'YTick',[.7:.05:1],'FontName','arial','FontSize',14)


%% Figure 5: Boxplot/scatter of all latencies 

% create age_vector
for subj = 1:length(database)
    
    age_vector(1,subj) = database(subj).age_ses1;
end


% manipulate age_vector for plotting (if two 17yr olds, change one to 16.9
% and 1 to 17.1)
age_vector_scale = age_vector;
age_vector_scale(6)= 16.8;
age_vector_scale(18) = 17.2;
age_vector_scale(7) = 9.8;
age_vector_scale(11) = 10.2;
age_vector_scale(15) = 10.8;
age_vector_scale(20) = 11.2;
age_vector_scale(9) = 12.8;
age_vector_scale(21) = 13.2;
age_vector_scale(8) = 14.4;
age_vector_scale(12) = 14.7;
age_vector_scale(24) = 15.3;
age_vector_scale(25) = 15.7;



% create boxplot of all latencies for each subject 
fig31 = figure(31); set(fig31,'Position',[0 0 2000 700]), hold on

boxplot(matrix_reshape_all, 'positions', age_vector_scale, 'labels', age_vector_scale,'Plotstyle', 'compact', ...
    'OutlierSize', 1, 'Symbol', 'x')


xlabel('Age subject (years)')
ylabel('Latency (ms)')
%title('Scatterplot of all detected latencies for each subject')

set(gca,'XTick',[0:10:50],'FontName','arial','FontSize',28)
xlim([0 55]),ylim([1 109])

% plot trendline
x = [6:1:50];
y = (-0.6547 * x) + 39.4192; 
plot(x,y, 'r')

hold off



%% Figure 6: Variation in latency - finished

% calculate SD of latencies
subj_std = nanstd(matrix_reshape_all,[],1);

% plot latency variation SD per subject
figure(41), hold on

plot(age_vector,subj_std,'.','MarkerSize',20)

xlim([0 55]),ylim([0 30])
set(gca,'XTick',[0:10:50], 'YTick', [0:10:30],'FontName','arial','FontSize',18)
xlabel('Age subject (years)')
ylabel('Variance SD (ms)')
%title('Variance in latency across subjects')

hold off

%% Figure 7: Violin + boxplot group differences all latencies - finished

% clear variables that will be used
clear group_1_cceps
clear group_2_cceps

% create two groups, an 18- and an 18+ group with the corresponding data
% also create age_group which corresponds with the group number
for gg = 1:length(age_vector)
    if age_vector(gg) < 18        
        
        group_1_cceps(:,gg) = matrix_reshape_all(:,gg);
        group_2_cceps(:,gg) = nan(length(matrix_reshape_all),1);
        
        age_group(gg) = 1;
                
    elseif age_vector(gg) >= 18
       
        group_1_cceps(:,gg) = nan(length(matrix_reshape_all),1);
        group_2_cceps(:,gg) = matrix_reshape_all(:,gg);        
        
        age_group(gg) = 2;
        
    end
end

% reshape for right format for violin function
group_1_cceps = reshape(group_1_cceps(:), (size(group_1_cceps,1) * size(group_1_cceps,2)),1);
group_2_cceps = reshape(group_2_cceps(:), (size(group_1_cceps,1) * size(group_1_cceps,2)),1);

% delete nan's 
group_1_cceps = group_1_cceps(~isnan(group_1_cceps));
group_2_cceps = group_2_cceps(~isnan(group_2_cceps));

% because they have different sizes, calculate difference and add the
% difference as number of nan's to the shortest group. This is necessary
% for the violin function
find_smallest = abs(size(group_2_cceps,1) - size(group_1_cceps,1));

if size(group_2_cceps,1) < size(group_1_cceps,1)
    
    group_2_cceps = [group_2_cceps;  nan(find_smallest,1)];
    
elseif size(group_2_cceps,1) >= size(group_1_cceps,1)
    
    group_1_cceps = [group_1_cceps;  nan(find_smallest,1)];
end

figure(51), hold on

% plot violin
violin([group_1_cceps group_2_cceps],'facecolor',[[1 0 0];[0 0 1]], 'edgecolor', 'none', 'medc','','mc','')
%set(gca, 'Box', 'off')

% plot boxplot within violin
boxplot(matrix_reshape_all,age_group,'Labels',{'Children','Adults'},'LabelOrientation', 'horizontal',...
    'Plotstyle', 'compact','Symbol','','Colors',[0 0 0])

ylim([9 109])
ylabel('Latency (ms)')
set(gca,'FontName','arial','FontSize',18)
set(gca, 'Box', 'off')
%title('Latency distribution in children and adults')

hold off


%% Figure 8: Within ROI scatterplots

% create cells with all within measures to put in loop
within_plots = {ROI_plot_matrix_15, ROI_plot_matrix_26,  ROI_plot_matrix_38, ROI_within_plot_matrix_all}; 
titles = {'Frontal Middle Gyrus', 'Supramarginal Gyrus', 'Temporal Middle Gyrus', 'Three ROIs combined'};

fig81 = figure(81); set(fig81,'Position',[0 0 800 600])

for qq = 1:length(within_plots)
    % for all three ROIs combined
    subplot(2,2,qq), hold on
    
    within_mat = cell2mat(within_plots(qq));
    
    % iterate over all subjects in database
    for subj = 1:length(database)

        % scatterplot
        scatter(repelem(age_vector(1,subj),length(within_mat)),within_mat(:,subj))

    end

    xlim([0 55]),ylim([1 109])
    
    if qq == 1 || qq == 3
        ylabel('Latency (ms)')       
    end
    if qq == 3 || qq == 4
        xlabel('Age subject (years)')
    end
    title(titles(qq))
    set(gca,'FontName','arial','FontSize',16, 'XTick', [10:10:50])
    
    % calculate means
    trendline_ROI = nanmean(within_mat);

    % do regress analysis to find the linear regression
    [B,BINT,R,RINT,STATS] = regress(age_vector',[trendline_ROI' ones(length(trendline_ROI),1)]);

    % plot regression trendine
    x = [6:1:50];
    y = B(1) * x + B(2); 
    plot(x,y, 'r')
    hold off

end

%% Figure 9: Within ROI violin plots

% create cells as input for loop
within_plots = {ROI_plot_matrix_15, ROI_plot_matrix_26,  ROI_plot_matrix_38, ROI_within_plot_matrix_all}; 
titles = {'Fro. Mid.', 'Sup. Mar.', 'Tem. Mid. ', 'Combined'};


fig91 = figure(91); set(fig91,'Position',[0 0 1100 250])

for yy = 1:length(within_plots)
    
    % clear variables that will be used
    clear group_1_cceps
    clear group_2_cceps
    
    
    % create two groups, an 18- and an 18+ group with the corresponding data
    % also create age_group which corresponds with the group number
    for gg = 1:length(age_vector)

        within_mat = cell2mat(within_plots(yy));

        if age_vector(gg) < 18        

            group_1_cceps(:,gg) = within_mat(:,gg);
            group_2_cceps(:,gg) = nan(length(within_mat),1);

            age_group(gg) = 1;

        elseif age_vector(gg) >= 18

            group_1_cceps(:,gg) = nan(length(within_mat),1);
            group_2_cceps(:,gg) = within_mat(:,gg);        

            age_group(gg) = 2;

        end
    end

    % reshape for right format for violin function
    group_1_cceps = reshape(group_1_cceps(:), (size(group_1_cceps,1) * size(group_1_cceps,2)),1);
    group_2_cceps = reshape(group_2_cceps(:), (size(group_1_cceps,1) * size(group_1_cceps,2)),1);

    % delete nan's 
    group_1_cceps = group_1_cceps(~isnan(group_1_cceps));
    group_2_cceps = group_2_cceps(~isnan(group_2_cceps));

    % because they have different sizes, calculate difference and add the
    % difference as number of nan's to the shortest group. This is necessary
    % for the violin function
    find_smallest = abs(size(group_2_cceps,1) - size(group_1_cceps,1));

    if size(group_2_cceps,1) < size(group_1_cceps,1)

        group_2_cceps = [group_2_cceps;  nan(find_smallest,1)];

    elseif size(group_2_cceps,1) >= size(group_1_cceps,1)

        group_1_cceps = [group_1_cceps;  nan(find_smallest,1)];
    end


    % for all three ROIs combined
    subplot(1,4,yy), hold on


    % plot violin
    violin([group_1_cceps group_2_cceps],'facecolor',[[1 0 0];[0 0 1]], 'edgecolor', 'none', 'medc','','mc','')

    % plot boxplot within violin
    boxplot(within_mat,age_group,'Labels',{'18-','18+'},'LabelOrientation', 'horizontal',...
        'Plotstyle', 'compact','Symbol','','Colors',[0 0 0])
    set(gca,'Box', 'off')

    ylim([9 109])
    
    if yy ==1 
        ylabel('Latency (ms)')
    end
    set(gca,'FontName','arial','FontSize',16)
    title(titles(yy))

    hold off
end


%% Figure 10: Between ROIs scatter

fig101 = figure(101); set(fig101,'Position',[0 0 450 350]), hold on

    
% iterate over all subjects in database
for subj = 1:length(database)
    
    % scatterplot
    scatter(repelem(age_vector(1,subj),length(ROI_between_plot_matrix)),ROI_between_plot_matrix(:,subj))
    
end

xlim([0 55]),ylim([1 109])
    

ylabel('Latency (ms)')
xlabel('Age subject (years)')

set(gca,'FontName','arial','FontSize',18, 'XTick', [10:10:50])

% calculate means
trendline_ROI = nanmean(ROI_between_plot_matrix);

% do regress analysis to find the linear regression
[B,BINT,R,RINT,STATS] = regress(age_vector',[trendline_ROI' ones(length(trendline_ROI),1)]);

% plot regression trendine
x = [6:1:50];
y = B(1) * x + B(2);
plot(x,y, 'r')
hold off


%% Figure 11: Between ROIs violin

% clear variables that will be used
clear group_1_cceps
clear group_2_cceps

% create two groups, an 18- and an 18+ group with the corresponding data
% also create age_group which corresponds with the group number
for gg = 1:length(age_vector)
    if age_vector(gg) < 18        
        
        group_1_cceps(:,gg) = ROI_between_plot_matrix(:,gg);
        group_2_cceps(:,gg) = nan(length(ROI_between_plot_matrix),1);
        
        age_group(gg) = 1;
                
    elseif age_vector(gg) >= 18
       
        group_1_cceps(:,gg) = nan(length(ROI_between_plot_matrix),1);
        group_2_cceps(:,gg) = ROI_between_plot_matrix(:,gg);        
        
        age_group(gg) = 2;
        
    end
end

% reshape for right format for violin function
group_1_cceps = reshape(group_1_cceps(:), (size(group_1_cceps,1) * size(group_1_cceps,2)),1);
group_2_cceps = reshape(group_2_cceps(:), (size(group_1_cceps,1) * size(group_1_cceps,2)),1);

% delete nan's 
group_1_cceps = group_1_cceps(~isnan(group_1_cceps));
group_2_cceps = group_2_cceps(~isnan(group_2_cceps));

% because they have different sizes, calculate difference and add the
% difference as number of nan's to the shortest group. This is necessary
% for the violin function
find_smallest = abs(size(group_2_cceps,1) - size(group_1_cceps,1));

if size(group_2_cceps,1) < size(group_1_cceps,1)
    
    group_2_cceps = [group_2_cceps;  nan(find_smallest,1)];
    
elseif size(group_2_cceps,1) >= size(group_1_cceps,1)
    
    group_1_cceps = [group_1_cceps;  nan(find_smallest,1)];
end

figure(111), hold on

% plot violin
violin([group_1_cceps group_2_cceps],'facecolor',[[1 0 0];[0 0 1]], 'edgecolor', 'none', 'medc','','mc','')

% plot boxplot within violin
boxplot(ROI_between_plot_matrix,age_group,'Labels',{'Children','Adults'},'LabelOrientation', 'horizontal',...
    'Plotstyle', 'compact','Symbol','','Colors',[0 0 0])
set(gca,'Box', 'off')

ylim([9 109])
ylabel('Latency (ms)')
set(gca,'FontName','arial','FontSize',18)
%title('Latency distribution in children and adults')

hold off

%% Figure 12: relative number of CCEPs scatter

% be sure to have database(subj).amount_cceps and database(subj).total_stims
% if not, first run analysis_percentage_CCEPs.m

% do some recalculation to add within data and between data 
for subj = 1:length(database)
    
     if ~isempty(database(subj).ROI_between_all)
         
         total_relative_cceps(subj) = ((database(subj).amount_cceps) + (sum(~isnan(database(subj).ROI_between_all(:,3))))) / ...
             ((database(subj).total_stims) + (length(database(subj).ROI_between_all(:,3))));

     elseif isempty(database(subj).ROI_between_all) 
         
         total_relative_cceps(subj) = (database(subj).amount_cceps) / (database(subj).total_stims);
     end
end


figure(121),hold on

for subj = 1:length(database)
    
    scatter(age_vector(1,subj),(rel_perc_cceps_within(subj)*100),50,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1])
      
end


for subj = 1:length(database)
    
    scatter(age_vector(1,subj),(rel_perc_cceps_between(subj)*100),50,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0])
    
end

xlim([1 55]), ylim([0 100])
xlabel('Age subject (years)')
ylabel('Percentage CCEPs')
set(gca,'FontName','arial','FontSize',18)
%title('Relative number of CCEPs')

hold off


%% Figure 13: Rendering of CCEP latencies

addpath('/Fridge/users/jaap/github/ecogBasicCode/render/')
addpath('/Fridge/users/jaap/github/ccep/oldCode/2019/')
top_path = '/Fridge/users/jaap/ccep/dataBIDS/';
%%
% Rendering brain of RESP0621 (ss = 1), stimpair 21
subjects = {'RESP0621'};
hemi_cap = {'R'};

% subject number that correspons with the RESP
ss = 1; 

% pick a viewing angle:
v_dirs = [90 0]; %;90 0;90 -60;270 -60;0 0];

% set stimulated pair you want to render
stim_pair = 21;

% select significant peaks in the other channels
n1_plot_sample = database(ss).metadata(1).n1_peak_sample(:,stim_pair);
n1_plot_amplitude =  database(ss).metadata(1).n1_peak_amplitude(:,stim_pair);

%%% Loop for rendering
for s = 1%:3 %1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(top_path,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);
    
    % electrode locations name:
    dataLocName = dir(fullfile(top_path,['sub-' subj],'ses-1','ieeg',...
        ['sub-' subj '_ses-1_electrodes.tsv']));
    dataLocName = fullfile(dataLocName(1).folder,dataLocName(1).name);
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure(131)
        ecog_RenderGifti(g) % render
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
        
        % ecog_Label(els,30,12) % add electrode positions
        % add electrodes
        ccep_el_add(els,[0.1 0.1 0.1],20)
        % give stimulated electrodes other color
        ccep_el_add(els(database(ss).metadata(1).stimulated_pairs(stim_pair,1:2),:),[0 0 0],40)
        % set size and color of channels with significant peak 
        % based on sample (from stimulation on, so -5120) and the amplitude
        % color = latency, size = amplitude
        ccep_el_add_size_and_color(els,n1_plot_amplitude,(n1_plot_sample-5120),500,200)

        set(gcf,'PaperPositionMode','auto')
        % print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
        % ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end




% Rendering brain of RESP0768 (ss = 4), stimpair 29
subjects = {'RESP0768'};
hemi_cap = {'R'};

% subject number that correspons with the RESP
ss = 4; 

% pick a viewing angle:
v_dirs = [90 0]; %;90 0;90 -60;270 -60;0 0];

% set stimulated pair you want to render
stim_pair = 29;

% select significant peaks in the other channels
n1_plot_sample = database(ss).metadata(1).n1_peak_sample(:,stim_pair);
n1_plot_amplitude =  database(ss).metadata(1).n1_peak_amplitude(:,stim_pair);

%%% Loop for rendering
for s = 1%:3 %1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(top_path,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);
    
    % electrode locations name:
    dataLocName = dir(fullfile(top_path,['sub-' subj],'ses-1','ieeg',...
        ['sub-' subj '_ses-1_electrodes.tsv']));
    dataLocName = fullfile(dataLocName(1).folder,dataLocName(1).name);
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure(132)
        ecog_RenderGifti(g) % render
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
        
        % ecog_Label(els,30,12) % add electrode positions
        % add electrodes
        ccep_el_add(els,[0.1 0.1 0.1],20)
        % give stimulated electrodes other color
        ccep_el_add(els(database(ss).metadata(1).stimulated_pairs(stim_pair,1:2),:),[0 0 0],40)
        % set size and color of channels with significant peak 
        % based on sample (from stimulation on, so -5120) and the amplitude
        % color = latency, size = amplitude
        ccep_el_add_size_and_color(els,n1_plot_amplitude,(n1_plot_sample-5120),500,200)

        set(gcf,'PaperPositionMode','auto')
        title('Latencies of detected CCEPs')
        set(gca,'FontName','arial','FontSize',24)
        % print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
        % ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end


% to create legenda
cm1=[repmat([0 0 0],120,1)];
cm1(1:60,2)=[1];
cm1(61:120,1)=[1];
cm1(61:120,2)=flip([0:1/59:1]);
cm1(1:60,1) = [0:1/59:1];


% elsize=[15:(45-15)/(100-1):45];
elsize = ones(1,120)*50;

% to create legenda
fig133 = figure(133); set(fig133, 'Color',[1 1 1],'Position',[30 50 500 150]),hold on
for k=1:120
    plot(k,1,'.','MarkerSize',elsize(k),'Color',cm1(k,:))
end
title('Electrode Color')
xlim([10 100])
set(gca,'YColor','w')
set(gca,'XTick', [10:10:100],'YTick',[], 'FontName','arial','FontSize',18)
xlabel('N1 peak latency (in ms)')
set(gcf, 'PaperPositionMode', 'auto');

%% Apendix II - individual ROC curves

%% figure with multiple ROC-curves for multiple datasets

fig_individual_rocs = figure('Position',[0 0 1200 2000]);

% The validation_matrix and parameters_optimalize_mat of all validated data
% are loaded (the variables in the mat-files are renamed so they do not overwrite)
load(fullfile(top_path,'validation', 'sub-RESP0458_ses-1_run-011714_parameters_optimalization.mat'))

load(fullfile(top_path,'validation', 'sub-RESP0468_ses-1_run-031729_parameters_optimalization.mat'))

load(fullfile(top_path,'validation', 'sub-RESP0621_ses-1_run-021147_parameters_optimalization.mat'))

load(fullfile(top_path,'validation', 'sub-RESP0706_ses-1_run-041501_parameters_optimalization.mat'))

load(fullfile(top_path,'validation', 'sub-RESP0733_ses-1b_run-050941_parameters_optimalization.mat'))

load(fullfile(top_path,'validation', 'sub-RESP0768_ses-1_run-021704_parameters_optimalization.mat'))


roc_curves_list = {parameters_optimalize_mat_0458, parameters_optimalize_mat_0468, parameters_optimalize_mat_0621, ...
    parameters_optimalize_mat_0706, parameters_optimalize_mat_0733, parameters_optimalize_mat_0768}; 

titles = {'RESP0458 - 9yr', 'RESP0468 - 16yr','RESP0621 - 18yr', 'RESP0706 - 50yr', 'RESP0733 - 8yr', 'RESP0768 - 6yr'};


for aa = 1:length(roc_curves_list)
    
    fig_individual_rocs = figure(aa);  set(fig_individual_rocs,'Position',[0 0 800 400]);

    % plot change level line
    plot([0 1],[0 1],'k'),hold on
    
    
    
    individual_roc = cell2mat(roc_curves_list(aa)); 
    
    % use jet as colors
    my_colors = jet(size(individual_roc,2));
    
    % for all different ranges plot a ROC- curve
    for time_th = 1:size(individual_roc,2)
        
        sens_plot = individual_roc(:,time_th,1);
        spes_plot = individual_roc(:,time_th,2);
        
        plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
    end
    xlabel('1 - specificity')
    ylabel('sensitivity')
    xlim([0 1]),ylim([0 1])
    axis square
    title(titles(aa))
    set(gca,'XTick',[0:.5:1],'YTick',[0:.5:1],'FontName','arial','FontSize',16)
    
end
