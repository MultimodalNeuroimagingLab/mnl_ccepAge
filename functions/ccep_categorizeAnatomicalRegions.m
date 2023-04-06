%
%   Define the tracts, sub-tracts their end-point ROIs using (one or more Destrieux labels)
%   [tracts] = ccep_categorizeAnatomicalRegions()
%
%   Returns: 
%       tracts    = 
%
%
%   Note:     Tracts are based on the Yeh HCP YA842 tractography atlas, while their end-point areas were matched based on the Destrieux atlas
%   Note 2:   Subtracts are defined ...
%
%   Also see: Yeh, F. C., Panesar, S., Fernandes, D., Meola, A., Yoshino, M., Fernandez-Miranda, J. C., ... & Verstynen, T. (2018). Population-averaged atlas of the macroscale human structural connectome and its network topology. Neuroimage, 178, 57-68.
%
%   Dora Hermes, Max van den Boom, Dorien van Blooijs, 2022
%
function [tracts] = ccep_categorizeAnatomicalRegions()

    %  Yeh 1065 - TPAT (Temporo-Parietal Aslant Tract) segments
    tracts(1).tract_name                  = 'Y1065_TPAT';
    
    tracts(1).sub_tract(1).name           = 'Parietal-Temporal';
    tracts(1).sub_tract(1).roi1           = [25, 26, 27, 56];
    tracts(1).sub_tract(1).roi2           = [33, 37, 38, 60, 73, 74];
    tracts(1).sub_tract(1).allowIntraROI  = 0;
    
    
    %  Yeh 1065 - AF (Arcuate Fasciculus) tract-segments
    tracts(2).tract_name                  = 'Y1065_AF';
    
    tracts(2).sub_tract(1).name           = 'Frontal-Temporal';
    tracts(2).sub_tract(1).roi1           = [12, 14, 15, 52, 53];
    tracts(2).sub_tract(1).roi2           = [34, 36, 37, 38, 72, 73];
    tracts(2).sub_tract(1).allowIntraROI  = 0;

    
    %  Yeh 1065 - SLF2 tract-segments
    tracts(3).tract_name                  = 'Y1065_SLF2';
    
    tracts(3).sub_tract(1).name           = 'Frontal-Parietal';
    tracts(3).sub_tract(1).roi1           = [12, 14, 15, 52, 53];
    tracts(3).sub_tract(1).roi2           = [25, 26, 27, 56];
    tracts(3).sub_tract(1).allowIntraROI  = 0;
    
    tracts(3).sub_tract(2).name           = 'Frontal-Central';
    tracts(3).sub_tract(2).roi1           = [12, 14, 15, 52, 53, 54]; 
    tracts(3).sub_tract(2).roi2           = [3, 4, 28, 29, 45];
    tracts(3).sub_tract(2).allowIntraROI  = 0;


    %  Yeh 842 - U tract-segments
    tracts(4).tract_name                  = 'Y842_U';

    tracts(4).sub_tract(1).name           = 'PreCentral-PostCentral';
    tracts(4).sub_tract(1).roi1           = [29];                         % Yeh: PreC
    tracts(4).sub_tract(1).roi2           = [28];                         % Yeh: PostC
    tracts(4).sub_tract(1).allowIntraROI  = 0;
    
    tracts(4).sub_tract(2).name           = 'Central-Central';
    tracts(4).sub_tract(2).roi1           = [28, 29, 4, 45];
    tracts(4).sub_tract(2).roi2           = [28, 29, 4, 45];
    tracts(4).sub_tract(2).allowIntraROI  = 1;
    
    tracts(4).sub_tract(3).name           = 'Frontal-Frontal';
    tracts(4).sub_tract(3).roi1           = [12, 14, 15];
    tracts(4).sub_tract(3).roi2           = [12, 14, 15];
    tracts(4).sub_tract(3).allowIntraROI  = 1;
    
    tracts(4).sub_tract(4).name           = 'Parietal-Parietal';
    tracts(4).sub_tract(4).roi1           = [25, 26];
    tracts(4).sub_tract(4).roi2           = [25, 26];
    tracts(4).sub_tract(4).allowIntraROI  = 1;

    %{
    %  Yeh 1065 - VOF tract-segments
    tracts(5).tract_name                  = 'Y1065_VOF';
    
    tracts(5).sub_tract(1).name           = 'Dorsal-Ventral';
    tracts(5).sub_tract(1).roi1           = [19, 20, 58];
    tracts(5).sub_tract(1).roi2           = [2, 21, 59];
    tracts(5).sub_tract(1).allowIntraROI  = 0;
    %}
    
end
