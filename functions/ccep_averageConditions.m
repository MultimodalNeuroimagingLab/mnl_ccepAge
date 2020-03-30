function [average_ccep,average_ccep_names] = ccep_averageConditions(data,srate,events_table,stim_pair_nr,stim_pair_name,params)
%
% function [average_ccep,average_ccep_names] = ccep_averageConditions(data,srate,events_table,stim_pair_nr,stim_pair_name,params)
% calculates average across condition numbers for ccep data 
%
% input 
%   data: seeg or ecog data in electrodes X time
%   srate: sampling frequency
%   events_table: loaded table with bids events
%   stim_pair_nr: ccep condition number, starts at 1
%   stim_pair_name: ccep stim pair name (e.g. F01-F02)
%   params: parameters for epoch length
%       params.epoch_length: total epoch length in sec, default = 5
%       params.epoch_prestim_length: prestimulus epoch length in sec, default = 2
%       params.baseline_subtract: subtract median baseline from each trial
%
% output
%   average_ccep: electrodes X condition (stim pair) X time
%   average_ccep_names: ccep condition (stim pair) names
%
% Dora Hermes, 2020, Multimodal Neuroimaging Lab
% Dorien van Blooijs, 2020, UMC Utrecht

if isempty(params)
    % epochs of -2:3 seconds
    epoch_length = 5; 
    epoch_prestim_length = 2;
    baseline_subtract = 1;
else
    epoch_length = params.epoch_length; 
    epoch_prestim_length = params.epoch_prestim_length; 
    baseline_subtract = params.baseline_subtract;
end

nr_channels = size(data,1);
    
% set epoch parameters
tt = (1:epoch_length*srate)/srate - epoch_prestim_length;

% initialize output
average_ccep = NaN(nr_channels,max(stim_pair_nr),epoch_length*srate);
average_ccep_names = cell(max(stim_pair_nr),1);

for kk = 1:max(stim_pair_nr) % condition number
    disp(['loading data for condition ' int2str(kk) ' out of ' int2str(max(stim_pair_nr))])
    
    % epochs of this condition
    these_epochs = find(stim_pair_nr==kk);
    
    % save name of the current epoch
    average_ccep_names{kk} = stim_pair_name{these_epochs(1)};
    
    % initialize matrix with this epoch type (channels X conditions X time)
    these_epochs_data = NaN(nr_channels,length(these_epochs),epoch_length*srate);
    
    % for this condition number (kk), load each epoch (ll)
    
    for ll = 1:length(these_epochs)
        
        ll_start = round((events_table.onset(these_epochs(ll))-epoch_prestim_length)*srate);
        ll_end = round((events_table.onset(these_epochs(ll))+epoch_length-epoch_prestim_length)*srate);
        
        % load data
        if ~isnan(ll_start)
            these_epochs_data(:,ll,:) = data(:,ll_start+1:ll_end);
        end
        
    end    
        
    % baseline subtract
    if baseline_subtract==1
        samples_base = find(tt>-1 & tt<-0.1);
        these_epochs_data = ccep_baselinesubtract(these_epochs_data,samples_base,'median');
    end    
    
    % save average
    average_ccep(:,kk,:) = squeeze(nanmean(these_epochs_data,2));
    
    clear these_epochs_data ll_start ll_end
end