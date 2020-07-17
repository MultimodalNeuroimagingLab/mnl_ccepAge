function [t_new] = ieeg_labelElectrodesDestrieux(FSdir,electrodes_tsv_name,saveNew,circleradius)
%
% this script labels electrodes based on the Destrieux Atlas from
% Freesurfer
%
% requires vistasoft in the path
%
%%% Preperation step: 
%
% Tn the terminal go to the subjects freesurfer folder and type: 
%
%  mri_convert mri/aparc.a2009s+aseg.mgz mri/aparc.a2009s+aseg.nii.gz -rt nearest
%
% This converts the mri to the original T1 space, using nearest neighbor
% reslicing to keep one label for each voxel (no interpolation)
%
%% 

% load electrode positions
loc_info = readtable(electrodes_tsv_name,...
    'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

elecmatrix = [loc_info.x loc_info.y loc_info.z];

% load nifti file with labels
niDestrieux = niftiRead(fullfile(FSdir,'mri',...
    'aparc.a2009s+aseg.nii.gz'));

% load color table 
surface_labels_Destrieux = fullfile(FSdir,'label',...
    'lh.aparc.a2009s.annot');
[~, ~, colortable_Destrieux] = read_annotation(surface_labels_Destrieux);

% load labels for subcortical areas
subcortical_labels = readtable(fullfile(FSdir,'mri',...
    'aseg.auto_noCCseg.label_intensities.txt'),...
    'FileType','text','Delimiter',' ','ReadVariableNames',false,'HeaderLines',0);


%% Define output

Destrieux_label = zeros(size(elecmatrix,1),1);
Destrieux_label_text = cell(size(elecmatrix,1),1);

%%% LOOP THROUGH ELECTRODES AND ASSIGN LABELS
% for every electrode, look up the most common label within 2 mm radius
if isempty(circleradius)
    % use default
    circleradius = 3; % millimeters
end
voxelsize = niDestrieux.pixdim;

% electrodes xyz to indices
els_ind = niDestrieux.qto_ijk * [elecmatrix ones(size(elecmatrix,1),1)]';
els_ind = round(els_ind(1:3,:)');

for elec = 1:size(elecmatrix,1) % loop across electrodes
    if isnan(elecmatrix(elec,1)) % electrode no position
        Destrieux_label(elec) = NaN;
        Destrieux_label_text{elec} = 'n/a';
    elseif ~isnan(elecmatrix(elec,1))
        temp.electrode = zeros(size(niDestrieux.data));

        % box with ones
        minidata.elecbox = NaN(round(abs(circleradius*4/voxelsize(1))),...
            round(abs(circleradius*4/voxelsize(2))),...
            round(abs(circleradius*4/voxelsize(3))));
        % xyz
        [minidata.x,minidata.y,minidata.z] = ...
            ind2sub(size(minidata.elecbox),find(isnan(minidata.elecbox)));
        minidata.x = minidata.x*voxelsize(1);
        minidata.y = minidata.y*voxelsize(2);
        minidata.z = minidata.z*voxelsize(3);
        minidata.mean = [mean(minidata.x) mean(minidata.y) mean(minidata.z)];
        minidata.straal = sqrt((minidata.x-minidata.mean(1)).^2+... 
            (minidata.y-minidata.mean(2)).^2+... 
            (minidata.z-minidata.mean(3)).^2);
        minidata.elecbox(minidata.straal<circleradius) = 1;
        minidata.elecbox(minidata.straal==min(minidata.straal)) = 1;

        % set size box 
        temp.xsize = length(minidata.elecbox(:,1,1));
        temp.ysize = length(minidata.elecbox(1,:,1));
        temp.zsize = length(minidata.elecbox(1,1,:));

        if mod(temp.xsize,2)==0 % even
            temp.xmindefine = floor((length(minidata.elecbox(:,1,1))-1)/2);
            temp.xmaxdefine = ceil((length(minidata.elecbox(:,1,1))-1)/2);
        else
            temp.xmindefine = floor(length(minidata.elecbox(:,1,1))/2);
            temp.xmaxdefine = floor(length(minidata.elecbox(:,1,1))/2);
        end 
        if mod(temp.ysize,2)==0 % even
            temp.ymindefine = floor((length(minidata.elecbox(1,:,1))-1)/2);
            temp.ymaxdefine = ceil((length(minidata.elecbox(1,:,1))-1)/2);
        else
            temp.ymindefine = floor(length(minidata.elecbox(1,:,1))/2);
            temp.ymaxdefine = floor(length(minidata.elecbox(1,:,1))/2);
        end 
        if mod(temp.zsize,2)==0 % even
            temp.zmindefine = floor((length(minidata.elecbox(1,1,:))-1)/2);
            temp.zmaxdefine = ceil((length(minidata.elecbox(1,1,:))-1)/2);
        else
            temp.zmindefine = floor(length(minidata.elecbox(1,1,:))/2);
            temp.zmaxdefine = floor(length(minidata.elecbox(1,1,:))/2);
        end

        diffxmin=1;
        diffxmax=0;
        diffymin=1;
        diffymax=0;
        diffzmin=1;
        diffzmax=0;

        if els_ind(elec,1)-temp.xmindefine<=0
            disp('ERROR: electrode x < minimal x of image')
            diffxmin=abs(els_ind(elec,1)-temp.xmindefine)+1;
            temp.xmindefine=temp.xmindefine-diffxmin;
        end
        if els_ind(elec,2)-temp.ymindefine<=0
            disp('ERROR: electrode y < minimal y of image')
            diffymin=abs(els_ind(elec,2)-temp.ymindefine)+1;
            temp.ymindefine=temp.ymindefine-diffymin;%return
        end
        if els_ind(elec,3)-temp.zmindefine<=0
            disp('ERROR: electrode z < minimal z of image')
            diffzmin=abs(els_ind(elec,3)-temp.zmindefine)+1;
            temp.zmindefine=temp.zmindefine-diffzmin;%return
        end
        if els_ind(elec,1)+temp.xmaxdefine>length(temp.electrode(:,1,1))
            disp('ERROR: electrode x > maximal x of image')
            diffxmax=els_ind(elec,1)+temp.xmaxdefine-length(temp.electrode(:,1,1));
            temp.xmaxdefine=temp.xmaxdefine-(els_ind(elec,1)+temp.xmaxdefine-length(temp.electrode(:,1,1)));
        end
        if els_ind(elec,2)+temp.ymaxdefine>length(temp.electrode(1,:,1))
            disp('ERROR: electrode y > maximal y of image')
            diffymax=els_ind(elec,2)+temp.ymaxdefine-length(temp.electrode(1,:,1));
            temp.ymaxdefine=temp.ymaxdefine-(els_ind(elec,2)+temp.ymaxdefine-length(temp.electrode(1,:,1)));
        end
        if els_ind(elec,3)+temp.zmaxdefine>length(temp.electrode(1,1,:))
            disp(['ERROR: for electrode ' int2str(elec) ' z > maximal z of image'])
            diffzmax=els_ind(elec,3)+temp.zmaxdefine-length(temp.electrode(1,1,:));
            temp.zmaxdefine=temp.zmaxdefine-(els_ind(elec,3)+temp.zmaxdefine-length(temp.electrode(1,1,:)));
            %return
            %break
        end  

        minidata.elecbox = minidata.elecbox(diffxmin:end-diffxmax,diffymin:end-diffymax,diffzmin:end-diffzmax);

        % for each electrode draw circle
        temp.electrode(els_ind(elec,1)-temp.xmindefine:els_ind(elec,1)+temp.xmaxdefine,...
            els_ind(elec,2)-temp.ymindefine:els_ind(elec,2)+temp.ymaxdefine,...
            els_ind(elec,3)-temp.zmindefine:els_ind(elec,3)+temp.zmaxdefine) = ...
            minidata.elecbox;

        % get the label using the mode:
        thisLabel = double(mode(niDestrieux.data(temp.electrode==1)));

        % if <60% white matter, use the most common label of other voxels
        % 2 is label for left WM, 41 is label for right WM
        if thisLabel == 2
            wm_fraction = length(find(niDestrieux.data(temp.electrode==1)==2))./...
                length(niDestrieux.data(temp.electrode==1));
            if wm_fraction<.7 % close to WM, but use other label
                temp_labels = niDestrieux.data(temp.electrode==1);
                temp_labels(temp_labels==2) = [];
                thisLabel = mode(temp_labels);
            end
        elseif thisLabel == 41
            wm_fraction = length(find(niDestrieux.data(temp.electrode==1)==41))./...
                length(niDestrieux.data(temp.electrode==1));
            if wm_fraction<.7 % close to WM, but use other label
                temp_labels = niDestrieux.data(temp.electrode==1);
                temp_labels(temp_labels==41) = [];
                thisLabel = mode(temp_labels);
            end
        elseif thisLabel == 0 % if CSF, use other label if 10% of voxels have another label
            other_fraction = length(find(niDestrieux.data(temp.electrode==1)==0))./...
                length(niDestrieux.data(temp.electrode==1));
            if other_fraction<.90 % 
                temp_labels = niDestrieux.data(temp.electrode==1);
                temp_labels(temp_labels==0) = []; % remove zero labels
                thisLabel = mode(temp_labels);
            end
        end

        Destrieux_label(elec) = double(thisLabel);

        % get the text
        if Destrieux_label(elec)>0 && Destrieux_label(elec)<200
            Destrieux_label_text{elec} = subcortical_labels.Var2{subcortical_labels.Var1==thisLabel};
        elseif Destrieux_label(elec)>11099 && Destrieux_label(elec)<12100
            Destrieux_label_text{elec} = ['lh_' colortable_Destrieux.struct_names{Destrieux_label(elec)-11099}];
        elseif Destrieux_label(elec)>12099 
            Destrieux_label_text{elec} = ['rh_' colortable_Destrieux.struct_names{Destrieux_label(elec)-12099}];
        else
            Destrieux_label_text{elec} = 'n/a';
        end
    end
end

% append Destrieux columns to loc_info
if ~ismember('Destrieux_label',loc_info.Properties.VariableNames)
    tDes = table(Destrieux_label,Destrieux_label_text);
    t_new = [loc_info tDes(:,{'Destrieux_label','Destrieux_label_text'})];
else
    t_new = loc_info;
    t_new.Destrieux_label = Destrieux_label;
    t_new.Destrieux_label_text = Destrieux_label_text;
end

if saveNew==1
    writetable(t_new,electrodes_tsv_name,'FileType','text','Delimiter','\t'); 
end