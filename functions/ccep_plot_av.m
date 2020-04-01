function ccep_plot_av(average_ccep,tt,n1_peak_sample, n1_peak_amplitude,average_ccep_names,...
    channel_names,good_channels,myDataPath,bids_sub,bids_ses,bids_task,bids_runs,save_fig)
%
% function ccep_plot_av(average_ccep,tt,n1_peak_sample, n1_peak_amplitude,average_ccep_names,...
%     channel_names,good_channels,myDataPath,bids_sub,bids_ses,bids_task,bids_runs,save_fig)
%
% Function plots average CCEPs across conditions per electrode.
%
% input 
%   average_ccep: electrodes X condition (stim pair) X time
%   tt: time
%   n1_peak_sample: electrodes X condition, [] to not show
%   n1_peak_amplitude: [] or electrodes X condition, [] to not show
%   average_ccep_names: ccep condition (stim pair) names
%   channel_names: names of the channels (size electrodes)
%   good_channels: electrode indices to plot
%   myDataPath:
%   bids_sub:
%   bids_ses:
%   bids_ses:
%   bids_task:
%   bids_runs:
%   save_fig: 1 to save, 0 not to save
%
% Dora Hermes, 2020, Multimodal Neuroimaging Lab, Mayo Clinic
% Dorien van Blooijs, 2020, UMC Utrecht

if isempty(n1_peak_sample)
    n1_peak_sample = NaN(size(average_ccep,1),size(average_ccep,2));
end
if isempty(n1_peak_amplitude)
    n1_peak_amplitude = NaN(size(average_ccep,1),size(average_ccep,2));
end
    
elnrs_plot = good_channels;

for ll = 20%1:length(elnrs_plot)
    el_plot = elnrs_plot(ll);
    figure('Position',[0 0 700 700]),hold on
    for kk = 1:length(average_ccep_names)
        this_ccep_plot = squeeze(average_ccep(el_plot,kk,:));
        %         this_ccep_plot(tt>-0.010 & tt<0.010) = NaN;
        
        plot(tt,kk*500+zeros(size(tt)),'Color',[.8 .8 .8])
        plot(tt,kk*500+this_ccep_plot)
        if ~isnan(n1_peak_sample(el_plot,kk))
            plot(tt(n1_peak_sample(el_plot,kk)),n1_peak_amplitude(el_plot,kk)+kk*500,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',2)
        end
    end
    xlim([-.2 1])
    ylim([-500,(kk+2)*500])
    set(gca,'YTick',500*(1:length(average_ccep_names)),'YTickLabel',average_ccep_names)
    title([channel_names{el_plot}])
    
    ylabel('stimulated electrodes')
    xlabel('time(s)')
    
    % add amplitude bar
    plot([0.9 0.9],[1000 1500],'k','LineWidth',2)
    text(0.91,1250,['500 ' native2unicode(181,'latin1') 'V'])
    
    if save_fig==1
        % create folder to save figures
        if ~ exist(fullfile(myDataPath.output,'derivatives','av_ccep_figures',bids_sub,bids_ses,bids_runs),'dir')
            mkdir(fullfile(myDataPath.output,'derivatives','av_ccep_figures',bids_sub,bids_ses,bids_runs));
        end

        % filename
        figureName = fullfile(myDataPath.output,'derivatives','av_ccep_figures',bids_sub,bids_ses,bids_runs,...
            [bids_sub '_' bids_ses '_' bids_task '_' bids_runs '_incomingCCEP_el' channel_names{el_plot}]);
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',figureName)
        print('-depsc','-r300',figureName)
    end
%     close all
end

end