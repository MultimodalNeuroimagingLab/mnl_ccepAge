function [signalOut,spatfiltmatrix] = ieeg_car(signal, chans2incl)

% This function performs Common Avregae Reference (CAR) filtering on a
% signal
%
% Inputs:
% signal = time x channels
% chans2inc = channel indiced to include in common average reference
%
% Outputs:
% signalOut=signal*spatfiltmatrix;
%
% Adapted from BCI2000 code from Joshua Fialkoff, clarified inputs/outputs
% and added option to exclude channels.
%
% DH 2010

if size(signal,1) < size(signal,2) % signal samples X electrodes
    disp('transpose signal to be samples X electrodes')
    return
end

num_chans = length(chans2incl);

% spatfiltmatrix=[];
% create a CAR spatial filter
spatfiltmatrix=-1/num_chans*ones(size(signal,2)); % ischanged
for k = 1:size(signal,2)
    spatfiltmatrix(k,k)=1-1/num_chans;% was num_chans-1;
end

% start exclude channels
% signalOut1=double(signal)*spatfiltmatrix;
for k = 1:size(signal,2)
    if ismember(k,chans2incl)==0
        spatfiltmatrix(k,:)=0;
        spatfiltmatrix(:,k)=0;
        spatfiltmatrix(k,k)=1;
    end
end
% end exclude channels

signalOut=signal*spatfiltmatrix;

% clear spatfiltmatrix;