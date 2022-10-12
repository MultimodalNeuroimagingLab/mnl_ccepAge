%
% This script outputs supplementary figure 1, ...
%
% Dorien van Blooijs, Max van den Boom, 2022

clear
close all
clc

minDataPoints = 10;     % minimum data-points to correlate and linear-fit

%% 
%  Load the CCEP data from the derivatives

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'file')
    load(fullfile(myDataPath.output,'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'ccepData')
else
    disp('Run scripts ccep02_loadN1.m first')
end



%% 
%  ...

all_unique_currents =[];
all_currents = struct();

for iSubj = 1:length(ccepData)
    subject_unique_currents = [];
    subject_currents = struct();

    for iRun = 1:length(ccepData(iSubj).run)
        
        % variables to store the currents and latencies for each run
        run_unique_currents = [];
        run_currents = struct();
        
        % loop over the stim-pairs
        for iStimpair = 1:length(ccepData(iSubj).run(iRun).stimpair_currents)
            current = ccepData(iSubj).run(iRun).stimpair_currents{iStimpair} * 1000;
            
            if length(current) > 1
                warning('backtrace', 'off')
                warning(['Multiple currents on single stim-pair (subj ', ccepData(iSubj).id, ', run ', num2str(iRun), ', stimpair ', num2str(iStimpair), ' - ' , ccepData(iSubj).run(iRun).stimpair_names{iStimpair}, '), skipping']);
                continue;
            else
                
                % ensurure the existence of a field for the this current in this run
                run_unique_currents = [run_unique_currents; current];
                if ~isfield(run_currents, ['mA', num2str(current)])
                    run_currents.(['mA', num2str(current)]) = [];
                end
                
                % retrieve the (non-NaN) latencies for this pair (in samples) and convert them to seconds (from stim onset)
                lat = ccepData(iSubj).run(iRun).n1_peak_sample(~isnan(ccepData(iSubj).run(iRun).n1_peak_sample(:, iStimpair)), iStimpair);
                lat = ccepData(iSubj).run(iRun).tt(lat);
                
                % store the latencies for this stim-pair
                run_currents.(['mA', num2str(current)]) = [run_currents.(['mA', num2str(current)]), lat];
                %run_currents.(['mA', num2str(current)]) = [run_currents.(['mA', num2str(current)]), mean(lat)];
                
                clear lat
            end
        end
        run_unique_currents = unique(run_unique_currents);
        
        % add the run to the subject
        for iCurrent = 1:length(run_unique_currents)
            
            % ensure the existence of a field for this current for this subject
            subject_unique_currents = [subject_unique_currents; run_unique_currents(iCurrent)];
            if ~isfield(subject_currents, ['mA', num2str(run_unique_currents(iCurrent))])
                subject_currents.(['mA', num2str(run_unique_currents(iCurrent))]) = [];
            end
            
            % store the latencies for this subject
            subject_currents.(['mA', num2str(run_unique_currents(iCurrent))]) = [subject_currents.(['mA', num2str(run_unique_currents(iCurrent))]), ...
                                                            mean(run_currents.(['mA', num2str(run_unique_currents(iCurrent))]))];
            
        end
        
        clear run_unique_currents run_currents;   
        
    end
    subject_unique_currents = unique(subject_unique_currents);
    

    % take the mean latency over each current
    for iCurrent = 1:length(subject_unique_currents)
        subject_currents.(['mA', num2str(subject_unique_currents(iCurrent))]) = mean(subject_currents.(['mA', num2str(subject_unique_currents(iCurrent))]));
    end
    
    % add the subject to the global collection
    for iCurrent = 1:length(subject_unique_currents)
        %disp(['subj ', num2str(iSubj), ' - current ', num2str(subject_unique_currents(iCurrent)), ' - lat ', num2str(subject_currents.(['mA', num2str(subject_unique_currents(iCurrent))]), 2)]);
        
        % ensurure the existence of a field for the this current in this run
        all_unique_currents = [all_unique_currents; subject_unique_currents(iCurrent)];
        if ~isfield(all_currents, ['mA', num2str(subject_unique_currents(iCurrent))])
            all_currents.(['mA', num2str(subject_unique_currents(iCurrent))]) = [];
        end
        
        % for this sbuject and current, add the age and latency
        ind = length(all_currents.(['mA', num2str(subject_unique_currents(iCurrent))]));
        all_currents.(['mA', num2str(subject_unique_currents(iCurrent))])(ind + 1).age = ccepData(iSubj).age;
        all_currents.(['mA', num2str(subject_unique_currents(iCurrent))])(ind + 1).latency = subject_currents.(['mA', num2str(subject_unique_currents(iCurrent))]);
        
    end
    
    clear subject_unique_currents subject_currents;   
end
all_unique_currents = unique(all_unique_currents);


%%
% Calculate the correlation values, p-values (if over 10 data-points) and FDR corrected p-values

all_unique_currents_r = nan(1, length(all_unique_currents));
all_unique_currents_p = nan(1, length(all_unique_currents));

% for each unique current
for iCurrent = 1:length(all_unique_currents)
    ages = [all_currents.(['mA', num2str(all_unique_currents(iCurrent))]).age]';
    latencies = [all_currents.(['mA', num2str(all_unique_currents(iCurrent))]).latency]';
    
    if length(latencies) >= minDataPoints
        [r, p] = corr(ages, latencies * 1000, 'Type', 'Spearman');
        all_unique_currents_r(iCurrent) = r;
        all_unique_currents_p(iCurrent) = p;
        
    end
    
end

% FDR correct
pvals = all_unique_currents_p(~isnan(all_unique_currents_p));
[~, ~, ~, pvals_fdr]  = fdr_bh(pvals, 0.05, 'pdep');
all_unique_currents_p_fdr = all_unique_currents_p;
all_unique_currents_p_fdr(~isnan(all_unique_currents_p_fdr)) = pvals_fdr;
clear pvals pvals_fdr;


%%
%  Plot the figure

defaultColors = get(0, 'DefaultAxesColorOrder');


figure('Position', [0 0 800 800]);
hold on;

% for each unique current
for iCurrent = 1:length(all_unique_currents)
    ages = [all_currents.(['mA', num2str(all_unique_currents(iCurrent))]).age]';
    latencies = [all_currents.(['mA', num2str(all_unique_currents(iCurrent))]).latency]';
    
    % determine the color
    curColor = [0 0 0];
    if all_unique_currents(iCurrent) == 2
        curColor = defaultColors(1, :);
    elseif all_unique_currents(iCurrent) == 4
        curColor = defaultColors(2, :);
    elseif all_unique_currents(iCurrent) == 6
        curColor = defaultColors(3, :);
    elseif all_unique_currents(iCurrent) == 7
        curColor = defaultColors(4, :);
    end
    
    % plot the points
    plot(ages, latencies * 1000, '.', 'Color', curColor, 'MarkerSize', 12, 'DisplayName', [num2str(all_unique_currents(iCurrent)) ' mA (n=', num2str(length(latencies)), ')']);

    if length(latencies) >= minDataPoints
        
        [P, S] = polyfit(ages, latencies * 1000, 1);
        [y_fit, ~] = polyval(P, ages, S);

        % Plot polyfit throught data points
        plot(ages, y_fit, 'Color', curColor, 'LineWidth', 2, 'DisplayName', ...
             ['r=' num2str(all_unique_currents_r(iCurrent), 3) ' p=' num2str(all_unique_currents_p_fdr(iCurrent), 3)]);
    end
    
end
hold off;

xlabel('age (years)');
ylabel('mean latency (ms)');
xlim([0 50]), ylim([10 55]);
legend()

% save
if ~exist(fullfile(myDataPath.output, 'derivatives', 'age'), 'dir')
    mkdir(fullfile(myDataPath.output, 'derivatives', 'age'));
end

figureName = fullfile(myDataPath.output, 'derivatives', 'age', 'SupFigS1_corrAgeVslatency_8mA');
set(gcf, 'PaperPositionMode', 'auto')
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)

