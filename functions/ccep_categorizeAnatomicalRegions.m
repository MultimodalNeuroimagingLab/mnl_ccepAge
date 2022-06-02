% ccep_categorizeAnatomicalRegions
%
% Here we define the areas (using one or more Destrieux labels from the _electrodes.tsv) that will be used as 
% "end-point" areas given specific tract(-segments) and areas from the Yeh HCP YA842 tractography atlas.
% Note: The areas defined here will - in some scripts - be compared all-to-all. Which does not always makes sense (inter-tract), but we will be selecting the required tract connections after
%
% Also see: Yeh, F. C., Panesar, S., Fernandes, D., Meola, A., Yoshino, M., Fernandez-Miranda, J. C., ... & Verstynen, T. (2018). Population-averaged atlas of the macroscale human structural connectome and its network topology. Neuroimage, 178, 57-68.
%


%%
%  SLF tract-segments


% Frontal <----> Parietal

% Yeh: Fmid
roi_name{1} = 'SLF_FP_front';
roi{1} = {'15', '53'};

% Yeh: Pinf, PSp, Ag
roi_name{2} = 'SLF_FP_parietal';
roi{2} = {'25', '26', '27'};


% Frontal <----> Central

% Yeh: FSp, Finf, Fmid
roi_name{3} = 'SLF_FC_frontal';
roi{3} = {'12', '13', '14', '15', '16', '52', '53', '54'};

% Yeh: PreC, ParaC
roi_name{4} = 'SLF_FC_central';
roi{4} = {'3', '4', '29', '68', '69'};





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
