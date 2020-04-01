function [stim_pair_nr,stim_pair_name] = ccep_bidsEvents2conditions(events_table,events_include,params)
% function [stim_pair_nr,stim_pair_name] = bidsEvents2conditions(events_table,events_include)
% makes a stim_pair vector to average epochs
%
% input 
%   events_table: loaded table with bids events
%   events_include: optional vector with included trials, set to [] if
%       using all events
%   params: 
%   	[]: defaults
%       params.mergePlusMin 
%           1: merge F01-F02 and F02-F01 (default)
%           0: see F01-F02 and F02-F01 as separate conditions
%       params.mergeAmp 
%           0: merge stimulations of all amplitudes
%           1: test whether one pair has been stimulated with different
%               amplitudes and then separate conditions
%
% output
%   stim_pair_nr
%   stim_pair_name
% 
%   - stim_pair_nr & stim_pair_name have length of the events table
%   - stim_pair_nr indicates a condition number, starts at 1
%   - stim_pair_name contains the stim pair name (if mergePlusMin F02-F01 is switched to F01-F02 in the name)
%
% dhermes, 2020, Multimodal Neuroimaging Lab

if isempty(params)
    params.mergePlusMin = 1;
end

% include all events if input of events_include is empty
if isempty(events_include)
    events_include = ones(height(events_table),1);
end

% initialize stim_pair name and number (conditions vectors)
stim_pair_nr = NaN(height(events_table),1); % number of pair (e.g. 1)
stim_pair_name = cell(height(events_table),1); % name of pair (e.g. LTG1-LTG2) 

% get all stim + electrodes
stimEl1 = extractBefore(events_table.electrical_stimulation_site,'-');
% get all stim - electrodes
stimEl2 = extractAfter(events_table.electrical_stimulation_site,'-');

condition_type_counter = 0;
for kk = 1:height(events_table)

    % which electrodes are stimulated
    el1 = stimEl1{kk};
    el2 = stimEl2{kk};
        
    % is this trial a stimulation trial: do el1 & el2 have content & can
    % the event be included
    if ~isempty(el1) && ~isempty(el2) && events_include(kk)==1
        
        if params.mergePlusMin==1 % merge el1-el with el2-el1
            % if this trial type does not exist yet & is a stimulation trial
            if sum(strcmp(stim_pair_name,[el1 '-' el2]))==0 && ... % does el1-el2 already exist?
                    sum(strcmp(stim_pair_name,[el2 '-' el1]))==0 % group el2-el1 with el1-el2
                condition_type_counter = condition_type_counter+1;

                % find all trials with el1 & el2 | el2 & el1
                theseTrials = strcmp(stimEl1,el1) & strcmp(stimEl2,el2) | ...
                    strcmp(stimEl2,el1) & strcmp(stimEl1,el2);
                trial_nrs = find(theseTrials==1); % number of trials of this type 
                for ll = 1:sum(theseTrials)
                    stim_pair_name{trial_nrs(ll),1} = [el1 '-' el2];   
                end
                stim_pair_nr(theseTrials) = condition_type_counter;   
            end
            
        elseif params.mergePlusMin==0 % separate el1-el with el2-el1
            
            % if this trial type does not exist yet & is a stimulation trial
            if sum(strcmp(stim_pair_name,[el1 '-' el2]))==0 % el1-el2
                condition_type_counter = condition_type_counter+1;

                % find all trials with el1 & el2 
                theseTrials = strcmp(stimEl1,el1) & strcmp(stimEl2,el2);
                trial_nrs = find(theseTrials==1); % number of trials of this type 
                for ll = 1:sum(theseTrials)
                    stim_pair_name{trial_nrs(ll),1} = [el1 '-' el2];   
                end
                stim_pair_nr(theseTrials) = condition_type_counter;   
            end
        end
    end
end

clear el1 el2 trial_nrs theseTrials epoch_type_counter

% should conditions of different amplitudes be separated?
if params.mergeAmp ==0
    % TODO
end