function [outdir] = ieeg_FSsurf2T1space(FSdir,outdir)

% Function ieeg_FSsurf2T1space(FSdir,T1file,outdir) converts Freesurfer
% surface to a gifti file in freesurfer space and then converts this gifti
% to a gifti file in the space of the T1file
%
% Input:  
%   FSdir: freesurfer directory, if [] will provide a popup
%   output: output directory, if [] will save in freesurfer root
%
% Output: 
%   output directory
% 
% Need vistasoft in the path for Gifti and Freesurfers Matlab functions.
%
% Dora Hermes, Multimodal Neuroimaging Lab, 2020

if isempty(FSdir)
    disp('select freesurfer directory')
    [FSdir] = uigetdir(pwd,'select freesurfer directory');
end

if isempty(outdir) 
   outdir = FSdir; 
end

%% 
%% Safe pial surface as a gifti
%% Write gifti file in derivatives/surfaces with original MRI coordinates from freesurfer


for kk = 1:2
    clear g
    if kk==1
        hemi = 'l'; 
    elseif kk==2
        hemi = 'r';
    end

    % name to save the surface
    giftiSafeName = fullfile(outdir,['pial.' upper(hemi) '.surf.gii']);

    % load the Freesurfer surface
    [vertex_coords, faces] = read_surf(fullfile(FSdir,'surf',[hemi 'h.pial']));
    
    % load the Freesurfer mri_orig for the coordinates
    mri_orig = fullfile(FSdir,'mri','orig.mgz');
    
    % convert from freesurfer space to original space
    orig = MRIread(mri_orig);
    Torig = orig.tkrvox2ras;
    Norig = orig.vox2ras;
    freeSurfer2T1 = Norig*inv(Torig);

    % convert vertices to original space
    vert_mat = double(([vertex_coords ones(size(vertex_coords,1),1)])');
    vert_mat = freeSurfer2T1*vert_mat;
    vert_mat(4,:) = [];
    vert_mat = vert_mat';
    g.vertices = vert_mat; clear vert_mat
    g.faces = faces+1;
    g = gifti(g);
    
    disp(['saving ' giftiSafeName])
    save(g,giftiSafeName,'Base64Binary')
end

%% Safe inflated surface as a gifti
%% Write gifti file in derivatives/surfaces with original MRI coordinates from freesurfer

for kk = 1:2
    clear g
    if kk==1
        hemi = 'l'; 
    elseif kk==2
        hemi = 'r';
    end

    % name to save the surface
    giftiSafeName = fullfile(outdir,['inflated.' upper(hemi) '.surf.gii']);

    % load the Freesurfer surface
    [vertex_coords, faces] = read_surf(fullfile(FSdir,'surf',[hemi 'h.inflated']));
    
    % load the Freesurfer mri_orig for the coordinates
    mri_orig = fullfile(FSdir,'mri','orig.mgz');
    
    % convert from freesurfer space to original space
    orig = MRIread(mri_orig);
    Torig = orig.tkrvox2ras;
    Norig = orig.vox2ras;
    freeSurfer2T1 = Norig*inv(Torig);

    % convert vertices to original space
    vert_mat = double(([vertex_coords ones(size(vertex_coords,1),1)])');
    vert_mat = freeSurfer2T1*vert_mat;
    vert_mat(4,:) = [];
    vert_mat = vert_mat';
    g.vertices = vert_mat; clear vert_mat
    g.faces = faces+1;
    g = gifti(g);

    disp(['saving ' giftiSafeName])
    save(g,giftiSafeName,'Base64Binary')
end
