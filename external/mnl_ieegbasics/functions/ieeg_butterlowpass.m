function [band_sig]=butterlowpass_eeglabdata(signal,band,srate)
% function to do lowpass filtering for specified band
% input:
% signal = time X channels;
% band = [stop];

%automatic stuff:
num_chans=size(signal,2); % number of channels
Rp = 3; Rs = 60; % third order Butterworth
delta = 0.001*2/srate;
high_p = band(1)*2/srate;
high_s = high_p + 0.01;

[n_band, wn_band] = buttord(high_p, high_s, Rp, Rs);
[bf_b bf_a] = butter(n_band, wn_band,'low');

% fband=filtfilt(bf_b, bf_a, signal(:,84)); %band pass
band_sig=zeros(size(signal));
for k=1:num_chans
    % just for nice disp:
    if mod(k,5)==0,disp(strcat(num2str(k),'/',num2str(num_chans))),end %this is to tell us our progress as the program runs
    band_sig(:,k)=filtfilt(bf_b,bf_a, signal(:,k)); %band pass
end




