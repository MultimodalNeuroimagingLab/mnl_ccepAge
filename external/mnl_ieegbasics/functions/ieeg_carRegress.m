function [signal] = ieeg_carRegress(signal, chans2incl)
% common average referencing flunction with regressing out the mean, rather
% than subtracting out the mean
%
% [signal] = ieeg_CarRegress(signal, chans2incl)
%
% signal = electrodes X samples
% chans2incl: channels to include in CAR
%
% dh - Oct 2010


if size(signal,2) < size(signal,1) % signal samples X electrodes
    disp('transpose signal to be electrodes X samples')
    return
end

ca_signal = mean(signal(chans2incl,:),1);

% regress off the mean signal
for kk = 1:size(signal,1) % elecs
%     disp(['elec ' int2str(kk)])
    if ismember(kk,chans2incl)
        [~,~,R] = regress(signal(kk,:)',ca_signal');
        signal(kk,:) = R';
    end
end

