function data_epoch = ieeg_baselinesubtract(data_epoch,samples_base,varargin)
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

if isempty(varargin)
    baseline_function = 'mean';
elseif ~isempty(varargin)
    baseline_function = varargin{1};
end

% baseline correct
for kk = 1:size(data_epoch,1)%channels
    for mm = 1:size(data_epoch,2)%epochs
        x = squeeze(data_epoch(kk,mm,:));
        
        if strcmp(baseline_function,'mean')
            x = x-mean(x(samples_base));
        elseif strcmp(baseline_function,'median')
            x = x-median(x(samples_base));
        end
        data_epoch(kk,mm,:) = x;
    end
end
