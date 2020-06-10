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

%% Get standardized electrodes from surface based registration (through sphere)

% get a list of datasets
theseSubs = ccep_getSubFilenameInfo(myDataPath);

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


%% Start here to make figures

load(fullfile(myDataPath.output,'derivatives','elec_coordinatesMNI305.mat'),'elec_coords')

% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.input,'derivatives','freesurfer');

% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));

% load mni305 inflated
[Lmniinfl_vert,Lmniinfl_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.inflated'));
[Rmniinfl_vert,Rmniinfl_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.inflated'));


%% add all electrodes labels and left or right hemisphere
allmni_coords = [];
allmni_coords_infl = [];

allmni_labels = [];
all_hemi = [];
for kk = 1:length(elec_coords)
    Destrieux_label = elec_coords(kk).elecs_tsv.Destrieux_label;
    if iscell(Destrieux_label)
        for ll = 1:size(Destrieux_label,1)
            if ischar(Destrieux_label{ll})
                if isequal(Destrieux_label,'n/a') 
                    Destrieux_label{ll} = NaN;
                else
                    Destrieux_label{ll} = str2double(Destrieux_label{ll});
                end
            end
        end
        Destrieux_label = cell2mat(Destrieux_label);
    end
    
    allmni_coords = [allmni_coords; elec_coords(kk).mni_coords];
    allmni_labels = [allmni_labels; Destrieux_label];
    all_hemi = [all_hemi; elec_coords(kk).hemi];
    
    % run through all coordinates and find the inflated points
    temp_inflated = NaN(size(elec_coords(kk).mni_coords));
    for ll = 1:size(Destrieux_label,1)
        if isequal(elec_coords(kk).hemi{ll},'L')
            [~,min_ind] = min(sqrt(sum((Lmniinfl_vert-elec_coords(kk).mni_coords(ll,:)).^2,2)));
            temp_inflated(ll,:) = Lmniinfl_vert(min_ind,:);
        elseif isequal(elec_coords(kk).hemi{ll},'R')
            [~,min_ind] = min(sqrt(sum((Rmniinfl_vert-elec_coords(kk).mni_coords(ll,:)).^2,2)));
            temp_inflated(ll,:) = Rmniinfl_vert(min_ind,:);
        end
    end
    
    allmni_coords_infl = [allmni_coords_infl; temp_inflated];
end

%% labels for electrode areas we want to color

% G_temporal_inf, G_temporal_middle, G_temp_sup-Lateral,
% G_oc-temp_med-Parahip, G_oc-temp_lat-fusifor
roi_temporal = [37 38 34 23 21];
% G_front_inf-Triangul, G_front_middle, G_front_inf-Opercular
roi_frontal = [14 15 12]; 
% G_pariet_inf-Angular, G_pariet_inf-Supramar, G_parietal_sup
roi_parietal = [25 26 27];
% G_postcentral G_precentral S_central
roi_central = [28 29 46];

%% Left pial with electrodes in mni space

v_d = [270 0];

figure
gl.faces = Lmnipial_face+1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);
tH = ieeg_RenderGifti(gl);

% make sure electrodes pop out
a_offset = .5*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      
% els = allmni_coords;

ieeg_elAdd(els(ismember(all_hemi,'L'),:),'k',10)
% set(tH,'FaceAlpha',.5) % make transparent
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_temporal),:),[0 0 .8],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_frontal),:),[1 .8 0],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_central),:),[0 .5 .5],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_parietal),:),[0 .5 0],15)
ieeg_viewLight(v_d(1),v_d(2))

%%
v_d = [96 6];

figure
gr.faces = Rmnipial_face+1;
gr.vertices = Rmnipial_vert;
gr = gifti(gr);
tH = ieeg_RenderGifti(gr);

% make sure electrodes pop out
a_offset = .5*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

ieeg_elAdd(els(ismember(all_hemi,'R'),:),'k',10)
% set(tH,'FaceAlpha',.5) % make transparent
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_temporal),:),[0 0 .8],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_frontal),:),[1 .8 0],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_central),:),[0 .5 .5],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_parietal),:),[0 .5 0],15)
ieeg_viewLight(v_d(1),v_d(2))


%% Left inflated with electrodes in mni space

v_d = [270 0];

figure
gl.faces = Lmniinfl_face+1;
gl.vertices = Lmniinfl_vert;
gl = gifti(gl);
tH = ieeg_RenderGifti(gl);

% make sure electrodes pop out
% a_offset = .1*max(abs(allmni_coords_infl(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
% els = allmni_coords_infl+repmat(a_offset,size(allmni_coords_infl,1),1);      
els = allmni_coords_infl;

ieeg_elAdd(els(ismember(all_hemi,'L'),:),'k',10)
% set(tH,'FaceAlpha',.5) % make transparent
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_temporal),:),[0 0 .8],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_frontal),:),[1 .8 0],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_central),:),[0 .5 .5],15)
ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels,roi_parietal),:),[0 .5 0],15)
ieeg_viewLight(v_d(1),v_d(2))

%%
v_d = [96 6];

figure
gr.faces = Rmnipial_face+1;
gr.vertices = Rmnipial_vert;
gr = gifti(gr);
tH = ieeg_RenderGifti(gr);

% make sure electrodes pop out
a_offset = .5*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

ieeg_elAdd(els(ismember(all_hemi,'R'),:),'k',10)
% set(tH,'FaceAlpha',.5) % make transparent
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_temporal),:),[0 0 .8],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_frontal),:),[1 .8 0],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_central),:),[0 .5 .5],15)
ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels,roi_parietal),:),[0 .5 0],15)
ieeg_viewLight(v_d(1),v_d(2))