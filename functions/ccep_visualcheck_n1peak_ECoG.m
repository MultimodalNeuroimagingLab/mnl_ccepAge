function [n1_peak_amplitude_check, n1_peak_sample_check ] = ccep_visualcheck_n1peak_ECoG(average_ccep, ccep,average_ccep_names,channel_names,tt,n1_peak_amplitude,n1_peak_sample)

%
% function ccep_visualcheck_n1peak_ECoG(average_ccep, average_ccep_names,...
%   ccep,tt,n1_peak_amplitude,n1_peak_sample)
%
% Function plots average CCEPs and all CCEPs to visually check detected N1s
%
% input
%   average_ccep: electrodes X condition (stim pair) X time
%   tt: time
%   n1_peak_sample: electrodes X condition, [] to not show
%   n1_peak_amplitude: [] or electrodes X condition, [] to not show
%   average_ccep_names: ccep condition (stim pair) names
%   ccep: electrodes X condition (stim pair) X trial X time
%   channel_names: names of the channels (size electrodes)
%
% Dora Hermes, 2020, Multimodal Neuroimaging Lab, Mayo Clinic
% Dorien van Blooijs, 2020, UMC Utrecht

% INSTRUCTIONS
% select new point: select new point > enter or y
% correct: y
% incorrect: n or enter

close all
clc

n1_peak_amplitude_check = NaN(size(n1_peak_amplitude));
n1_peak_sample_check = NaN(size(n1_peak_sample));

for stimp = 1:size(average_ccep,2)
    for chan = 1 :size(average_ccep,1)
        
        if ~isnan(n1_peak_sample(chan,stimp))
            % figure with left the epoch, and right zoomed in
            H = figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.11 0.77 0.8];
            this_plot = squeeze(ccep(chan,stimp,:,:));
            %             this_plot(:,tt>0 &tt<0.01) = NaN;
            
            this_plot_avg = squeeze(average_ccep(chan,stimp,:));
            this_plot_avg(tt>0 & tt<0.01) = NaN;
            
            subplot(1,2,1)
            plot(tt,this_plot,':r','linewidth',1);
            hold on
            plot(tt,this_plot_avg,'k','linewidth',2);
            plot(tt(n1_peak_sample(chan,stimp)),n1_peak_amplitude(chan,stimp),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3)
            hold off
            xlim([-2 2])
            ylim([-2000 2000])
            xlabel('time(s)')
            ylabel('amplitude(uV)')
            title(sprintf('Electrode %s, stimulated by %s',channel_names{chan},average_ccep_names{stimp}))
            
            subplot(1,2,2)
            plot(tt,this_plot,':r','linewidth',1);
            hold on
            plot(tt,this_plot_avg,'k','linewidth',2);
            plot(tt(n1_peak_sample(chan,stimp)),n1_peak_amplitude(chan,stimp),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
            hold off
            xlim([-0.2 0.5])
            ylim([-750 750])
            title('Zoomed average signal')
            xlabel('Time (s)')
            ylabel('Voltage (uV)')
            
            currkey = 0;
            fprintf('N1 [y/n], if incorrect N1, select correct N1 \n')
            
            % select new N1 or categorize as good N1 or no N1
            while ~strcmp(currkey,{'y','n',char(13)})
                
                w = waitforbuttonpress;
                if w == 0
                    % draw correct N1
                    cp = get(gca,'CurrentPoint');
                    
                    % find sample number closest to the selected point
                    [~,sampnum] = min(abs(tt-cp(1,1)));
                    
                    % find nearby peak
                    [~,locs] = findpeaks(-1*this_plot_avg(sampnum-50:sampnum+50),...
                        'NPeaks',1,'SortStr','descend');
                    
                    % find x-position of nearby peak
                    locsamp = sampnum-50+locs-1;
                    
                    hold on
                    plot(tt(locsamp),this_plot_avg(locsamp),'bo','MarkerFaceColor','b','MarkerSize',4); drawnow;
                    hold off
                    
                    n1_peak_sample_check(chan,stimp) = locsamp ;
                    n1_peak_amplitude_check(chan,stimp) = this_plot_avg(locsamp) ;
                    
                elseif w == 1
                    currkey = get(gcf,'CurrentCharacter');
                    
                    if strcmp(currkey,'y') && isempty(cp)
                        n1_peak_amplitude_check(chan,stimp) = n1_peak_amplitude(chan,stimp) ;
                        n1_peak_sample_check(chan,stimp) = n1_peak_sample(chan,stimp) ;
                    elseif strcmp(currkey,'n')
                        n1_peak_amplitude_check(chan,stimp) = NaN ;
                        n1_peak_sample_check(chan,stimp) = NaN ;
                    end
                end
            end
            
        end
    end
end
disp('Finished visual check!')








