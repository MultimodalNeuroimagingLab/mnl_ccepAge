%
%   Calculates average across condition numbers for ccep data 
%   [average_ccep, average_ccep_names, stimpair_currents, ccep, tt] = ccep_averageConditions(data, srate, events_table, channel_names, stim_pair_nr, stim_pair_name, params, verbose)
%
%       data                        = seeg or ecog data in electrodes X time
%       srate                       = sampling frequency
%       events_table                = loaded table with bids events
%       channel_names               = 
%       stim_pair_nr                = ccep condition number, starts at 1
%       stim_pair_name              = ccep stim pair name (e.g. F01-F02)
%       verbose                     = 
%
%       params                      = the configuration struct (pass [] for defaults)
%       params.epoch_length         = total epoch length in sec, default = 5
%       params.epoch_prestim_length = prestimulus epoch length in sec, default = 2
%       params.baseline_subtract    = subtract median baseline from each trial
%
% 
%   Returns: 
%       average_ccep                = the averaged and baselined CCEPs (format: electrodes X conditions/stim-pair X time), the stimulated
%                                     electrodes are set to NaN for each seperate stim-pair
%       stimpair_names              = conditions/stim-pair names (corresponds to the second dimension of the average_ccep variable)
%       stimpair_currents           = currents that were found in the trials for each specific condition/stim-pair
%       ccep                        = the CCEP for each trial (format: electrodes X conditions/stim-pair X trials X time)
%       tt                          = time for each epoch
%
% 
% 
%   Dora Hermes, Dorien van Blooijs, Max van den Boom, 2022
%
function [average_ccep, stimpair_names, stimpair_currents, ccep, tt] = ccep_averageConditions(data, srate, events_table, channel_names, stim_pair_nr, stim_pair_name, params, verbose)

    % defaults/optional parameters
    if ~exist('params', 'var') || isempty(params)
        % epochs of -2:3 seconds
        epoch_length = 5; 
        epoch_prestim_length = 2;
        baseline_subtract = 1;
    else
        epoch_length = params.epoch_length; 
        epoch_prestim_length = params.epoch_prestim_length; 
        baseline_subtract = params.baseline_subtract;
    end
    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 1;
    end

    nr_channels = size(data, 1);
    max_nr_epochs = max(histcounts(stim_pair_nr, 'BinMethod', 'integers'));
    % set epoch parameters
    tt = (1:epoch_length * srate) / srate - epoch_prestim_length;

    % initialize output variables
    average_ccep        = NaN(nr_channels, max(stim_pair_nr), round(epoch_length * srate));
    stimpair_names      = cell(max(stim_pair_nr), 1);
    stimpair_currents   = cell(max(stim_pair_nr), 1);
    ccep                = NaN(nr_channels, max(stim_pair_nr), max_nr_epochs, round(epoch_length * srate));
    

    %%
    %  
    
    for iCond = 1:max(stim_pair_nr)
        if verbose == 1
            disp(['loading data for condition/stimpair: ' int2str(iCond) ' out of ' int2str(max(stim_pair_nr))])
        end
        
        % trials/epochs belonging to this condition/stim-pair
        trialIndices = find(stim_pair_nr == iCond);

        % save name of the current stimulation pair
        stimpair_names{iCond} = stim_pair_name{trialIndices(1)};

        % store the currents for this condition/stim-pair
        stimpair_currents{iCond} = events_table.electrical_stimulation_current(trialIndices);
        stimpair_currents{iCond}(isnan(stimpair_currents{iCond})) = [];
        stimpair_currents{iCond} = unique(stimpair_currents{iCond});
        if length(stimpair_currents{iCond}) > 1
            warning('Same stim-pair was stimulated with different currents over the span of the run');
        end
        
        % initialize matrix with this epoch type (channels X individual trials X time)
        stimpair_epochedData = NaN(nr_channels, length(trialIndices), round(epoch_length * srate));
        
        % for this stimulation pair number (iStimPair), load each epoch (iEpoch)
        for iTrial = 1:length(trialIndices)

            % the sample_start is correct
            epoch_start = round(events_table.sample_start(trialIndices(iTrial)) - epoch_prestim_length * srate);
            epoch_end = round(events_table.sample_start(trialIndices(iTrial)) + (epoch_length - epoch_prestim_length) * srate);

            % load data
            if ~isnan(epoch_start) && epoch_end < size(data, 2) && epoch_start > 0
                stimpair_epochedData(:, iTrial, :) = data(:, epoch_start + 1:epoch_end);
            end

        end    

        % baseline subtract for each individual epoch (iEpoch)
        if baseline_subtract == 1
            samples_base = find(tt > -1 & tt < -0.1);
            stimpair_epochedData = ccep_baselinesubtract(stimpair_epochedData, samples_base, 'median');
        end

        if size(stimpair_epochedData, 2) == max_nr_epochs
            % save all cceps
            ccep(:, iCond, :, :) = stimpair_epochedData;
        else
            ccep(:, iCond, 1:size(stimpair_epochedData, 2), :) = stimpair_epochedData;
        end

        % average CCEPs for this condition/stim-pair
        stimpair_average_ccep = squeeze(mean(stimpair_epochedData, 2, 'omitnan'));   % [channels x samples]

        % since each trial has had baseline subtraction between -1:-0.1, we
        % average all trials, there could be a baseline shift from 0. So here,
        % again, we apply a baseline correction (this was first part of the
        % function ccep_detect_n1peak_ECoG.m, but we relocated this here.
        baseline_tt = tt > -2 & tt < -.1;
        signal_median = median(stimpair_average_ccep(:, baseline_tt), 2);

        % save average CCEPs (for this condition/stim-pair) into a global variable after subtracting the median (for each channel)
        average_ccep(:, iCond, :) = stimpair_average_ccep - signal_median;

        clear stimpair_epochedData epoch_start epoch_end stimpair_average_ccep
    end

    
    %% 
    %  Set stimulated electrodes to NaN 
    
    for iCond = 1:size(average_ccep, 2) % epochs
        
        % stimulated electrode names
        el1 = extractBefore(stimpair_names{iCond}, '-');
        if contains(extractAfter(stimpair_names{iCond}, '-'), '-')     % in case el1 - el2 - XmA is present in average_ccep_names
            el2 = extractBefore(extractAfter(stimpair_names{iCond}, '-'), '-');
        else
            el2 = extractAfter(stimpair_names{iCond}, '-');
        end
        
        % stimulated electrode index in the data
        el1_nr = ismember(channel_names, el1);
        el2_nr = ismember(channel_names, el2);

        % set to NaN
        average_ccep(el1_nr == 1, iCond, :) = NaN;
        average_ccep(el2_nr == 1, iCond, :) = NaN;

        clear el1 el2 el1_nr el2_nr
    end

end
