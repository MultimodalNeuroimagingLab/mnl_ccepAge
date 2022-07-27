%
% [roi_track, roi_name, roi] = ccep_categorizeAnatomicalRegions()
% 
% Here we define the areas (using one or more Destrieux labels from the _electrodes.tsv) that will be used as 
% "end-point" areas given specific tract(-segments) and areas from the Yeh HCP YA842 tractography atlas.
% Note: The areas defined here will - in some scripts - be compared all-to-all. Which does not always makes sense (inter-tract), but we will be selecting the required tract connections after
%
% Also see: Yeh, F. C., Panesar, S., Fernandes, D., Meola, A., Yoshino, M., Fernandez-Miranda, J. C., ... & Verstynen, T. (2018). Population-averaged atlas of the macroscale human structural connectome and its network topology. Neuroimage, 178, 57-68.
%
% DH, MvdB, DvB, 2022
%
function [rois] = ccep_categorizeAnatomicalRegions()

    %  SLF tract-segments
    rois(1).tract_name = 'SLF';
    
    rois(1).sub_tract(1).name       = 'Frontal-Parietal';
    rois(1).sub_tract(1).roi1       = [15, 53];               % Yeh: Fmid
    rois(1).sub_tract(1).roi2       = [25, 26, 27];           % Yeh: Pinf, PSp, Ag
    rois(1).sub_tract(1).interHemi  = 0;
    
    rois(1).sub_tract(2).name       = 'Frontal-Central';
    rois(1).sub_tract(2).roi1       = [12, 13, 14, 15, 16, 52, 53, 54];     % Yeh: FSp, Finf, Fmid
    rois(1).sub_tract(2).roi2       = [3, 4, 29, 68, 69];     % Yeh: PreC, ParaC
    rois(1).sub_tract(2).interHemi  = 0;
    
    %  Arcuate Fasciculus tract-segments
    
    rois(2).tract_name = 'AF';
    
    rois(2).sub_tract(1).name       = 'Frontal-Temporal';
    rois(2).sub_tract(1).roi1       = [12, 13, 14, 15, 52, 53];    % Yeh: Finf, Fmid
    rois(2).sub_tract(1).roi2       = [34, 36, 37, 38, 72, 73];    % Yeh: Tinf, TMd, TSp
    rois(2).sub_tract(1).interHemi  = 0;

    




    %{

    %%
    %  Regions on initial article submission

    % temporal areas:
    % G_temporal_inf, G_temporal_middle, G_temp_sup-Lateral,
    % G_oc-temp_med-Parahip, G_oc-temp_lat-fusifor
    roi{1} = {'37','38','34','23','21'};
    roi_name{1} = 'temporal';
    roi_temporal = [37 38 34 23 21];

    % central areas
    % G_postcentral G_precentral S_central
    roi{2} = {'28','29','45'};
    roi_name{2} = 'central';
    roi_central = [28 29 45];

    % parietal areas:
    % G_pariet_inf-Angular, G_pariet_inf-Supramar, G_parietal_sup
    roi{3} = {'25','26','27'};
    roi_name{3} = 'parietal';
    roi_parietal = [25 26 27];

    % frontal areas:
    % G_front_inf-Triangul, G_front_middle, G_front_inf-Opercular
    roi{4} = {'14','15','12'}; 
    roi_name{4} = 'frontal';
    roi_frontal = [14 15 12]; 

    %}
    
end
