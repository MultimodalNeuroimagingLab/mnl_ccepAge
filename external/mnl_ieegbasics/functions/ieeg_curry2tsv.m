function [t] = ieeg_curry2tsv(t1File,pomFile)
%
% Function writes a tsv file with electrodes names and xyz coordinates that
% should match the T1. 
%
% Inputs
%   t1File: name of the t1File that was used to localize electrodes
%   pomFile: name of the .pom file created by Curry with electrodes in XYNo coordinates
%
% Outputs
%   t: table with name, x, y, z, size of electrodes
%       size is n/a 
%
%
% [t] = ieeg_curry2tsv(t1File,pomFile)
%
% dhermes, Multimodal Neuroimaging Lab, Mayo Clinic, 2020

if isempty(t1File)
    disp('select T1 nifti file')
    [FILENAME, PATHNAME] = ...
        uigetfile({'*.nii;*.nii.gz','nifti-files (*.nii, *.nii.gz)'},'select T1 nifti file');
    t1File = [PATHNAME,FILENAME];
end

if isempty(pomFile)
    disp('select .pom file in XYNo coordinates')
    [FILENAME, PATHNAME] = ...
        uigetfile('*.pom','select .pom file in XYNo coordinates');
    pomFile = [PATHNAME,FILENAME];
end


% load the nifti file
ni = niftiRead(t1File);

% read the pom file and put coordinates in xyz_cu matrix
aa = ieeg_inchannel_curry_pom(pomFile);
nr_channels = length(aa.Channel);
xyz_cu = [];
for kk = 1:nr_channels
    xyz_cu(kk,1) = aa.Channel(kk).Loc(1)+1; 
    xyz_cu(kk,2) = ni.dim(2)-aa.Channel(kk).Loc(2);%ni.dim(2)-aa.Channel(kk).Loc(2); 
    xyz_cu(kk,3) = aa.Channel(kk).Loc(3);%aa.Channel(kk).Loc(3)+1; % add 1 because Curry may start counting at zero? 
end

% convert to xyz coordinates in the T1 frame
xyz = ni.qto_xyz * [xyz_cu(:,1) xyz_cu(:,2) xyz_cu(:,3) ones(length(xyz_cu),1)]';
xyz = xyz(1:3,:)';

%%

% Fields for sub-XX_ses-XX_electrodes.tsv
% Required fields:
name = cell(nr_channels,1);
for kk = 1:nr_channels
    name{kk}        = aa.Channel(kk).Name;% label of the channel
end

hemisphere = cell(nr_channels,1);
for kk = 1:nr_channels
    if ismember('R',name{kk})
        hemisphere{kk}  = 'R';
    elseif ismember('L',name{kk})
        hemisphere{kk}  = 'L';
    else
        hemisphere{kk}  = 'n/a';
    end
end
x               = xyz(:,1);
y               = xyz(:,2);
z               = xyz(:,3);
size            = NaN(nr_channels,1); % Surface area of electrode in mm
seizure_zone    = NaN(nr_channels,1); % Surface area of electrode in mm

t = table(name,x,y,z,size,hemisphere,seizure_zone);
clear size x y z name seizure_zone% housekeeping

pathForSave = fileparts(pomFile);

for kk = 1:1000
    electrodes_tsv_name = fullfile(pathForSave,...
        ['electrodes_' int2str(kk) '.tsv']);
    if exist(electrodes_tsv_name,'file')
        continue
    else
        break
    end
end

t = bids_tsv_nan2na(t);

disp(['saving ' electrodes_tsv_name])
writetable(t,electrodes_tsv_name,'FileType','text','Delimiter','\t');

