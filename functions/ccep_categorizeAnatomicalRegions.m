% ccep_categorizeAnatomicalRegions
%
% Regions are defined by eind-point labels of Yeh's HCP YA842 tractography atlas.
% Here the Destrieux labels from the _electrodes.tsv are mapped (e.g. 'G_front_middle') to tract end-point labels (e.g. 'Fmid')
%
%
% Also see: Yeh, F. C., Panesar, S., Fernandes, D., Meola, A., Yoshino, M., Fernandez-Miranda, J. C., ... & Verstynen, T. (2018). Population-averaged atlas of the macroscale human structural connectome and its network topology. Neuroimage, 178, 57-68.
%


%[19;37;38;34;25;36;4;26;15;14;68;52;12;0;23;43;10;21;72;9;22;8;16;3;30;27;42;2;67;28;29;54;71;59;50;57;20;56;44;73;11;46;24;65;58;5;1;13;7;31;53;45;64;61;6;60;66;55;62;47;35;49;63;69;32;33;40;18;41;74;17]
%{'G_occipital_middle';'G_temporal_inf';'G_temporal_middle';'G_temp_sup-Lateral';'G_pariet_inf-Angular';'G_temp_sup-Plan_tempo';'G&S_subcentral';'G_pariet_inf-Supramar';'G_front_middle';'G_front_inf-Triangul';'S_precentral-inf-part';'S_front_inf';'G_front_inf-Opercular';'Unknown';'G_oc-temp_med-Parahip';'Pole_temporal';'G_cingul-Post-ventral';'G_oc-temp_lat-fusifor';'S_temporal_inf';'G_cingul-Post-dorsal';'G_oc-temp_med-Lingual';'G&S_cingul-Mid-Post';'G_front_sup';'G&S_paracentral';'G_precuneus';'G_parietal_sup';'Pole_occipital';'G&S_occipital_inf';'S_postcentral';'G_postcentral';'G_precentral';'S_front_sup';'S_subparietal';'S_occipital_ant';'S_collat_transv_ant';'S_oc_middle&Lunatus';'G_occipital_sup';'S_intrapariet&P_trans';'S_calcarine';'S_temporal_sup';'G_cuneus';'S_cingul-Marginalis';'G_orbital';'S_parieto_occipital';'S_oc_sup&transversal';'G&S_transv_frontopol';'G&S_frontomargin';'G_front_inf-Orbital';'G&S_cingul-Mid-Ant';'G_rectus';'S_front_middle';'S_central';'S_orbital-H_Shaped';'S_oc-temp_med&Lingual';'G&S_cingul-Ant';'S_oc-temp_lat';'S_pericallosal';'S_interm_prim-Jensen';'S_orbital_lateral';'S_circular_insula_ant';'G_temp_sup-Plan_polar';'S_circular_insula_sup';'S_orbital_med-olfact';'S_precentral-sup-part';'G_subcallosal';'G_temp_sup-G_T_transv';'Lat_Fis-ant-Vertical';'G_insular_short';'Lat_Fis-post';'S_temporal_transverse';'G_Ins_lg&S_cent_ins'}


%{
a = [65;58;5;1;31;64;61;6;60;66;55;35;63;32;33;40;41;74];
b = {'S_parieto_occipital';'S_oc_sup&transversal';'G&S_transv_frontopol';'G&S_frontomargin';'G_rectus';'S_orbital-H_Shaped';'S_oc-temp_med&Lingual';'S_oc-temp_lat';'S_pericallosal';'S_interm_prim-Jensen';'G_temp_sup-Plan_polar';'S_orbital_med-olfact';'G_subcallosal';'G_temp_sup-G_T_transv';'Lat_Fis-ant-Vertical';'Lat_Fis-post';'S_temporal_transverse'};
find(contains(b, 'S_orbital_lateral'))
%}



%%
% Frontal end-point areas:

% Yeh:              Frontal Interior (Finf)
% Destrieux:        G_front_inf-Opercular (12), G_front_inf-Triangul (14), S_front_inf (52), S_orbital_lateral (62), G_orbital (24), G_front_inf-Orbital (13)
roi{1} = {'12', '13', '14', '24', '52', '62'};
roi_name{1} = 'F_inf';

% Yeh:              Frontal Middle (Fmid)
% Destrieux:        G_front_middle (15), S_front_middle (53)
roi{2} = {'15', '53'};
roi_name{2} = 'F_mid';

% Yeh:              Frontal Superior (FSp)
% Destrieux:        G_front_sup (16), S_front_sup (54)
roi{1} = {'16', '54'};
roi_name{1} = 'F_Sp';



%%
% Central end-point areas

% Yeh:              Precentral (PreC)
% Destrieux:        G_precentral (29), S_precentral-inf-part (68), S_central (45), S_precentral-sup-part (69)
roi{1} = {'29', '45', '68', '69'};
roi_name{1} = 'C_PreC';

% Yeh:              Postcentral (PostC)
% Destrieux:        G&S_subcentral (4), G_postcentral (28), S_postcentral (67)
roi{1} = {'4', '28', '67'};
roi_name{1} = 'C_PostC';

% Yeh:              Paracentral (ParaC)
% Destrieux:        G&S_paracentral (3)
roi{1} = {'3'};
roi_name{1} = 'C_ParaC';



%%
% Parietal end-point areas

% Yeh:              Parietal Superior (PSp)
% Destrieux:        G_parietal_sup (27)
roi{1} = {'27'};
roi_name{1} = 'P_Sp';

% Yeh:              Angular (Ag)
% Destrieux:        G_pariet_inf-Angular (25)
roi{1} = {'25'};
roi_name{1} = 'P_Ag';

% Yeh:              SupraMarginal (SpMar)
% Destrieux:        G_pariet_inf-Supramar (26)
roi{1} = {'26'};
roi_name{1} = 'P_SpMar';

% Yeh:              Precuneus (PreCun)
% Destrieux:        G_precuneus (30)
roi{1} = {'30'};
roi_name{1} = 'P_PreCun';



%%
% Occipital end-point areas

% Yeh:              Occipital Inferior (Oinf)
% Destrieux:        G&S_occipital_inf (2), S_occipital_ant (59)
roi{1} = {'2', '59'};
roi_name{1} = 'O_inf';

% Yeh:              Occipital Middle (Omd)
% Destrieux:        G_occipital_middle (19), S_oc_middle&Lunatus (57)
roi{1} = {'19', '57'};
roi_name{1} = 'O_md';

% Yeh:              Occipital Superior (Osp)
% Destrieux:        G_occipital_sup (20), S_intrapariet&P_trans (56)
roi{1} = {'20', '56'};
roi_name{1} = 'O_sp';

% Yeh:              Lingual (Lin)
% Destrieux:        G_oc-temp_med-Lingual (22)
roi{1} = {'22'};
roi_name{1} = 'O_Lin';

% Yeh:              Calcarine (Cal)
% Destrieux:        S_calcarine (44), Pole_occipital (42)
roi{1} = {'42', '44'};
roi_name{1} = 'O_Cal';

% Yeh:              Cuneus (Cun)
% Destrieux:        G_cuneus (11)
roi{1} = {'11'};
roi_name{1} = 'O_Cun';



%%
% Temporal end-point areas

% Yeh:              Temporal Inferior (Tinf)
% Destrieux:        G_temporal_inf (37), Pole_temporal (43), S_collat_transv_ant (50) 
roi{1} = {'37', '43', '50'};
roi_name{1} = 'T_inf';

% Yeh:              Temporal Middle (TMd)
% Destrieux:        G_temporal_middle (38), S_temporal_inf (72), S_temporal_sup (73)
roi{1} = {'38', '72', '73'};
roi_name{1} = 'T_Md';

% Yeh:              Temporal Superior (TSp)
% Destrieux:        G_temp_sup-Lateral (34), G_temp_sup-Plan_tempo (36)
roi{1} = {'34', '36'};
roi_name{1} = 'T_Sp';

% Yeh:              Fusiform (Fu)
% Destrieux:        G_oc-temp_lat-fusifor (21)
roi{1} = {'21'};
roi_name{1} = 'T_Fu';

% Yeh:              Hippocampus (Hipp)
% Destrieux:        G_oc-temp_med-Parahip (23)
roi{1} = {'23'};
roi_name{1} = 'T_Hipp';



%%
% Internal/Other end-point areas

% Yeh:              Cingulum (C)
% Destrieux:        G&S_cingul-Ant (6), G&S_cingul-Mid-Ant (7), G&S_cingul-Mid-Post (8), G_cingul-Post-dorsal (9), G_cingul-Post-ventral (10), S_cingul-Marginalis (46)
roi{1} = {'6', '7', '8', '9', '10', '46'};
roi_name{1} = 'I_Cing';

% Yeh:              Insula (Insula)
% Destrieux:        S_circular_insula_ant (47), S_circular_insula_sup (49), G_insular_short (18), G_Ins_lg&S_cent_ins (17)
roi{1} = {'17', '18', '47', '49'};
roi_name{1} = 'I_Insula';



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
