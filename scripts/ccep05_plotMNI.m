%
% This script can be used as workflow script to create average CCEPs
% for the CCEP data in the RESPect database.
%
% Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, 2019
%

%% Set paths
clc
clear
myDataPath = setLocalDataPath(1);

%% get a list of datasets

theseSubs = ccep_getSubFilenameInfo(myDataPath);

%% get standardized electrodes from surface based registration (through sphere)

% Freesurfer subjects directory
FSsubjectsdir = fullfile(localDataPath,'derivatives/freesurfer');

elec_coords = [];

for kk = 1:length(theseSubs) 
    disp(['subj ' int2str(kk) ' of ' int2str(length(theseSubs))])
    
    % subject freesurfer dir
    FSdir = fullfile(myDataPath.input,'derivatives','freesurfer',theseSubs(kk).name);
    
    % get electrodes info
    elec_coords(kk).elecs_tsv = read_tsv(fullfile(myDataPath.input,theseSubs(kk).name,theseSubs(kk).ses,'ieeg',...
        [theseSubs(kk).name,'_',theseSubs(kk).ses,'_electrodes.tsv']));
    
    elecmatrix = [elec_coords(kk).elecs_tsv.x elec_coords(kk).elecs_tsv.y elec_coords(kk).elecs_tsv.z];
    nElec = size(elecmatrix,1);
    
    % hemi (we need to get this, probably from json file or so)
    hemi = loc_info.hemisphere;

    % convert to MNI using surface
    mni_coords = ieeg_mni305ThroughFsSphere(elecmatrix,hemi,FSdir,FSsubjectsdir);
    
end


%%

% get the mni sphere index in the mni pial
% mni305 pial
[mnipial_vert,mnipial_face] = read_surf(fullfile(FSsubjects,'fsaverage','surf','lh.pial'));

% figure to check electrodes in mni space
figure
g.faces = mnipial_face+1;
g.vertices = mnipial_vert;
g = gifti(g);
tH = ieeg_RenderGifti(g);
ieeg_label(mni_coords)
set(tH,'FaceAlpha',.5) % make transparent
