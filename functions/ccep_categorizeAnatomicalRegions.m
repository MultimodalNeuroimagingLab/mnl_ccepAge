%
%   Define the tracts, sub-tracts their end-point ROIs using (one or more Destrieux labels)
%   [roi_track, roi_name, roi] = ccep_categorizeAnatomicalRegions()
%
%   Returns: 
%       roi_track    = 
%       roi_name     = 
%       roi          = 
%
%
%   Note:     Tracts are based on the Yeh HCP YA842 tractography atlas, while their end-point areas were matched based on the Destrieux atlas
%   Note 2:   Subtracts are definedsubtracts defines SLF parts with endpoint parietal or central, since these have different lengths
%
%   Also see: Yeh, F. C., Panesar, S., Fernandes, D., Meola, A., Yoshino, M., Fernandez-Miranda, J. C., ... & Verstynen, T. (2018). Population-averaged atlas of the macroscale human structural connectome and its network topology. Neuroimage, 178, 57-68.
%
%   Dora Hermes, Max van den Boom, Dorien van Blooijs, 2022
%
function [rois] = ccep_categorizeAnatomicalRegions()

    %  Yeh 1065 - TPAT (Temporo-Parietal Aslant Tract) segments
    rois(1).tract_name                  = 'Y1065_TPAT';
    rois(1).interHemi                   = 0;
    
    rois(1).sub_tract(1).name           = 'Parietal-Temporal';
    rois(1).sub_tract(1).roi1           = [25, 26, 27, 56];
    rois(1).sub_tract(1).roi2           = [33, 37, 38, 60, 73, 74];
    rois(1).sub_tract(1).allowIntraROI  = 0;
    
    
    %  Yeh 1065 - VOF tract-segments
    rois(2).tract_name                  = 'Y1065_VOF';
    rois(2).interHemi                   = 0;
    
    rois(2).sub_tract(1).name           = 'Dorsal-Ventral';
    rois(2).sub_tract(1).roi1           = [19, 20, 58];
    rois(2).sub_tract(1).roi2           = [2, 21, 59];
    rois(2).sub_tract(1).allowIntraROI  = 0;

    
    %  Yeh 1065 - AF (Arcuate Fasciculus) tract-segments
    rois(3).tract_name                  = 'Y1065_AF';
    rois(3).interHemi                   = 0;
    
    rois(3).sub_tract(1).name           = 'Frontal-Temporal';
    rois(3).sub_tract(1).roi1           = [12, 14, 15, 52, 53];
    rois(3).sub_tract(1).roi2           = [34, 36, 37, 38, 72, 73];
    rois(3).sub_tract(1).allowIntraROI  = 0;

    
    %  Yeh 1065 - SLF2 tract-segments
    rois(4).tract_name                  = 'Y1065_SLF2';
    rois(4).interHemi                   = 0;
    
    rois(4).sub_tract(1).name           = 'Frontal-Parietal';
    rois(4).sub_tract(1).roi1           = [12, 14, 15, 52, 53];
    rois(4).sub_tract(1).roi2           = [25, 26, 27, 56];
    rois(4).sub_tract(1).allowIntraROI  = 0;
    
    rois(4).sub_tract(2).name           = 'Frontal-Central';
    rois(4).sub_tract(2).roi1           = [12, 14, 15, 52, 53, 54]; 
    rois(4).sub_tract(2).roi2           = [3, 4, 28, 29, 45];
    rois(4).sub_tract(2).allowIntraROI  = 0;


    %  Yeh 842 - U tract-segments
    rois(5).tract_name                  = 'Y842_U';
    rois(5).interHemi                   = 0;

    rois(5).sub_tract(1).name           = 'PreCentral-PostCentral';
    rois(5).sub_tract(1).roi1           = [29];                         % Yeh: PreC
    rois(5).sub_tract(1).roi2           = [28];                         % Yeh: PostC
    rois(5).sub_tract(1).allowIntraROI  = 0;
    
    rois(5).sub_tract(2).name           = 'Central-Central';
    rois(5).sub_tract(2).roi1           = [28, 29, 4, 45];
    rois(5).sub_tract(2).roi2           = [28, 29, 4, 45];
    rois(5).sub_tract(2).allowIntraROI  = 1;
    
    rois(5).sub_tract(3).name           = 'Frontal-Frontal';
    rois(5).sub_tract(3).roi1           = [12, 14, 15];
    rois(5).sub_tract(3).roi2           = [12, 14, 15];
    rois(5).sub_tract(3).allowIntraROI  = 1;
    
    rois(5).sub_tract(4).name           = 'Parietal-Parietal';
    rois(5).sub_tract(4).roi1           = [25, 26];
    rois(5).sub_tract(4).roi2           = [25, 26];
    rois(5).sub_tract(4).allowIntraROI  = 1;
    
    %{
    %  Yeh 842 - Superior Longitudinal Fasciculus tract-segments
    rois(6).tract_name                  = 'Y842_SLF';
    rois(6).interHemi                   = 0;
    
    rois(6).sub_tract(1).name           = 'Frontal-Parietal';
    rois(6).sub_tract(1).roi1           = [15, 53];                     % Yeh: Fmid
    rois(6).sub_tract(1).roi2           = [25, 26, 27, 56];             % Yeh: Pinf, PSp, Ag
    rois(6).sub_tract(1).allowIntraROI  = 0;
    
    rois(6).sub_tract(2).name           = 'Frontal-Central';
    rois(6).sub_tract(2).roi1           = [12, 13, 14, 15, 16, 52, 53, 54];     % Yeh: FSp, Finf, Fmid
    rois(6).sub_tract(2).roi2           = [3, 4, 29, 68, 69];           % Yeh: PreC, ParaC
    rois(6).sub_tract(2).allowIntraROI  = 0;
    
    % wang 2015
    rois(6).sub_tract(3).name           = 'WFrontal-Wparietal';
    rois(6).sub_tract(3).roi1           = [12, 13, 14, 15, 52, 53];
    rois(6).sub_tract(3).roi2           = [25, 26, 56];
    rois(6).sub_tract(3).allowIntraROI  = 0;
    %}
    
    
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
