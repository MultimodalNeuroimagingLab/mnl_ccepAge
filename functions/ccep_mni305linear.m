function mni_coords = ccep_mni305linear(elecmatrix,FSdir)
%
% Convert electrodes in native space to MNI using linear transformation
%
% Input:
%   elecmatrix: nElec x 3 in the original T1 space (before freesurfer)
%   hemi: hemisphere to load for each electrode (l/r)
%   FSdir: subjects freesurfer director
%
% Example:
%   mni_coords = ccep_mni305linear(elecmatrix,FSdir);
% 
%
% Figure to check: get the mni sphere index in the mni pial
% load the freesurfer MNI surface from freesurfer subjects directory with fsaverage
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

% transformation to mni305 to be applied
targ2mni305 = xfm_read(fullfile(FSdir,'mri','transforms','talairach.xfm'));

% do the transformation
mni_coords = (targ2mni305*Norig*(Torig\[elecmatrix'; ones(1, nElec)]))';
mni_coords = mni_coords(:,1:3);



