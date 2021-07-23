function data_epoch = ccep_baselinesubtract(data_epoch,samples_base,varargin)
%
% function data_epoch = ieeg_baselinesubtract(data_epoch,samples_base,baseline_function)
% subtracts the mean signal during samples_base from each epoch
%
% input
%   data_epoch: data with electrodes X epoch X t
%   samples_base: samples for the baseline calculation
%   baseline_function: OPTIONAL mean or median, calculated mean or median from the
%   baseline to subtract from each trial, default is mean
% 
% output
%   data_epoch
% 
%
% dhermes, multimodal neurimaging lab, 2020
% dvanblooijs, umcutrecht, 2021 - removed for-loop over channels to make it
% more time-efficient

if isempty(varargin)
    baseline_function = 'mean';
elseif ~isempty(varargin)
    baseline_function = varargin{1};
end

% baseline correct
for mm = 1:size(data_epoch,2)%epochs
    x = squeeze(data_epoch(:,mm,:));
    
    if strcmp(baseline_function,'mean')
        x = x-mean(x(:,samples_base),2);
    elseif strcmp(baseline_function,'median')
        x = x-median(x(:,samples_base),2);
    end
    data_epoch(:,mm,:) = x;
end

