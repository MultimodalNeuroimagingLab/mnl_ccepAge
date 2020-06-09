function mni_coords = ccep_mni305ThroughFsSphere(elecmatrix,hemi,FSdir,FSsubjectsdir)
%
% Convert electrodes in native space to MNI through surface based
% registration (through Freesurfer sphere). Since it is surface based,
% it should ony be used for ECoG.
%
% Input:
%   elecmatrix: nElec x 3 in the original T1 space (before freesurfer)
%   hemi: hemisphere to load for each electrode (l/r)
%   FSdir: subjects freesurfer director
%   FSsubjectsdir: freesurfer subjects directory with fsaverage
%
% Example:
%   mni_coords = ccep_mni305ThroughFsSphere(elecmatrix,hemi,FSdir,FSsubjectsdir);
% 
%
% Figure to check: get the mni sphere index in the mni pial
% load the freesurfer MNI surface
%   [mnipial_vert,mnipial_face] = read_surf(fullfile(FSsubjects,'fsaverage','surf','lh.pial'));
% 
%   figure
%   g.faces = mnipial_face+1;
%   g.vertices = mnipial_vert;
%   g = gifti(g);
%   tH = ieeg_RenderGifti(g);
%   ieeg_label(mni_coords)
%
%
%
% Dora Hermes, 2020
% Mutimodal Neuroimaging Lab
% Mayo Clinic

if isempty(FSsubjectsdir)
    disp('select freesurfer subjects directory with fsaverage dir')
    [FSsubjectsdir] = uigetdir(pwd,'select freesurfer directory');
end

if isempty(FSdir)
    disp('select this subjects freesurfer directory')
    [FSdir] = uigetdir(pwd,'select this subjects freesurfer directory');
end

% number of electrodes
nElec = size(elecmatrix,1);
disp(['getting MNI305 coordinates for ' int2str(nElec) ' electrodes'])

% load mri orig header
origName = fullfile(FSdir,'mri','orig.mgz');
orig = MRIread(origName,'true');
Norig = orig.vox2ras; 
Torig = orig.tkrvox2ras;

% electrodes to freesurfer space
freeSurfer2T1 = inv(Norig*inv(Torig));
elCoords = freeSurfer2T1*([elecmatrix'; ones(1, nElec)]);
elCoords = elCoords(1:3,:)';

% subject pial
[Lsubpial_vert,Lsubpial_face] = read_surf(fullfile(FSdir,'surf','lh.pial'));
[Rsubpial_vert,Rsubpial_face] = read_surf(fullfile(FSdir,'surf','rh.pial'));

% figure to check electrodes in freesurfer space
% figure
% g.faces = pial_face+1;
% g.vertices = pial_vert;
% tH = ieeg_RenderGifti(g);
% ieeg_label(elCoords)
% set(tH,'FaceAlpha',.5) % make transparent

% subject sphere
[Lsubsphere_vert] = read_surf(fullfile(FSdir,'surf','lh.sphere'));
[Rsubsphere_vert] = read_surf(fullfile(FSdir,'surf','rh.sphere'));

% mni305 sphere
[Lmnisphere_vert] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.sphere'));
[Rmnisphere_vert] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.sphere'));

% mni305 pial
[Lmnipial_vert] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));

% index for closest point to subject's pial
s_pial_ind = zeros(nElec,1);
% subject sphere coords at these indices
sphere_coords = zeros(nElec,3);

mni_coords = NaN(nElec,3);

for kk = 1:nElec
    xyz = elCoords(kk,:);
    if ~isnan(xyz(1))
        if isequal(upper(hemi{kk}),'L')

            % index for closest point to subject's pial
            [~,min_ind] = min(sqrt(sum((Lsubpial_vert-xyz).^2,2)));
            s_pial_ind(kk) = min_ind;

            % get the same index on subjects sphere
            sphere_coords(kk,:) = Lsubsphere_vert(min_ind,:);

            % closest point subjects sphere to mni sphere
            xyz_sphere = sphere_coords(kk,:);
            [~,min_ind_sphere] = min(sqrt(sum((Lmnisphere_vert-xyz_sphere).^2,2)));

            % get the  mni pial at the mni sphere index
            mni_coords(kk,:) = Lmnipial_vert(min_ind_sphere,:);


        elseif isequal(upper(hemi{kk}),'R')
            % index for closest point to subject's pial
            [~,min_ind] = min(sqrt(sum((Rsubpial_vert-xyz).^2,2)));
            s_pial_ind(kk) = min_ind;

            % get the same index on subjects sphere
            sphere_coords(kk,:) = Rsubsphere_vert(min_ind,:);

            % closest point subjects sphere to mni sphere
            xyz_sphere = sphere_coords(kk,:);
            [~,min_ind_sphere] = min(sqrt(sum((Rmnisphere_vert-xyz_sphere).^2,2)));

            % get the  mni pial at the mni sphere index
            mni_coords(kk,:) = Rmnipial_vert(min_ind_sphere,:);    
        end
    end
end

