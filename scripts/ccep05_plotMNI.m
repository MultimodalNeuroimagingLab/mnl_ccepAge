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
FSsubjectsdir = fullfile(myDataPath.input,'derivatives','freesurfer');

elec_coords = [];

for kk = 1:length(theseSubs) 
    disp(['subj ' int2str(kk) ' of ' int2str(length(theseSubs))])
    
    % subject freesurfer dir
    FSdir = fullfile(myDataPath.input,'derivatives','freesurfer',theseSubs(kk).name,theseSubs(kk).ses);
    
    % get electrodes info
    elec_coords(kk).elecs_tsv = readtable(fullfile(myDataPath.input,theseSubs(kk).name,theseSubs(kk).ses,'ieeg',...
        [theseSubs(kk).name,'_',theseSubs(kk).ses,'_electrodes.tsv']),'FileType','text','Delimiter','\t');
    if iscell(elec_coords(kk).elecs_tsv.x)
        elecmatrix = NaN(size(elec_coords(kk).elecs_tsv,1),3);
        for ll = 1:size(elec_coords(kk).elecs_tsv,1)
            if ~isequal(elec_coords(kk).elecs_tsv.x{ll},'n/a')
                elecmatrix(ll,:) = [str2double(elec_coords(kk).elecs_tsv.x{ll}) str2double(elec_coords(kk).elecs_tsv.y{ll}) str2double(elec_coords(kk).elecs_tsv.z{ll})];
            end
        end
    else
        elecmatrix = [elec_coords(kk).elecs_tsv.x elec_coords(kk).elecs_tsv.y elec_coords(kk).elecs_tsv.z];
    end
    nElec = size(elecmatrix,1);
    
    % get hemisphere for each electrode
    these_json = dir(fullfile(myDataPath.input,theseSubs(kk).name,theseSubs(kk).ses,'ieeg',[theseSubs(kk).name,'_',theseSubs(kk).ses,'_task-SPESclin*_ieeg.json']));
    ieeg_json = jsonread(fullfile(these_json(1).folder,these_json(1).name));
    if isequal(ieeg_json.iEEGPlacementScheme,'left')
        hemi = num2cell(repmat('L',nElec,1));
    elseif isequal(ieeg_json.iEEGPlacementScheme,'right')
        hemi = num2cell(repmat('R',nElec,1));
    elseif isequal(ieeg_json.iEEGPlacementScheme,'left,right')
        hemi = num2cell(repmat('n/a',nElec,1));
        for ll = 1:nElec
            if elecmatrix(kk,1)<0
                hemi{ll} = 'L';
            else
                hemi{ll} = 'R';
            end
        end
    end
    elec_coords(kk).hemi = hemi;
    % convert to MNI using surface
    elec_coords(kk).mni_coords = ccep_mni305ThroughFsSphere(elecmatrix,hemi,FSdir,FSsubjectsdir);
    
end

save(fullfile(myDataPath.output,'derivatives','elec_coordinatesMNI305.mat'),'elec_coords')


%%

% get the mni sphere index in the mni pial
% mni305 pial
[mnipial_vert,mnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));

% figure to check electrodes in mni space
figure
g.faces = mnipial_face+1;
g.vertices = mnipial_vert;
g = gifti(g);
tH = ieeg_RenderGifti(g);
ieeg_label(mni_coords)
set(tH,'FaceAlpha',.5) % make transparent
