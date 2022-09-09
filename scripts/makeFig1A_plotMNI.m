%
% This script can be used to create an MNI cortex (inflated) with
% electrodes in different colors for different locations for all patients
% used in this study. 
%
% Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, Max van den Boom, 2022
%



%%
%  Set paths and load data

clc
clear
close all
myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.input, 'derivatives', 'freesurfer');

% load the data struct (for electrode information and (sub-)tracts that run between different the end-point ROIs)
load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData');

% load the tract/ROIs end-point structure to store the concatenated tract-lines and vertex colors in
rois = ccep_categorizeAnatomicalRegions();

% check leadDBS availability and setup
if exist('ea_getearoot') ~= 2 || exist('ea_prefs') ~= 2
   error('Could not find LeadDBS functions. LeadDBS is an external dependency, make sure it installed and added as a MATLAB path.');
end
addpath(genpath(ea_getearoot()));
prefs = ea_prefs;
if exist('ea_trk2ftr') ~= 2 || exist('ea_ants_apply_transforms_to_points') ~= 2
    error('Could not find LeadDBS functions. LeadDBS is an external dependency, make sure it installed and added as a MATLAB path.');
end



%%
%  Retrieve the included line-tracts from all subjects

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        % init fields
        rois(iTr).sub_tract(iSubTr).allMNIlineIndices = cell(1, 2);
        rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount = cell(1, 2);
        if rois(iTr).interHemi == 1
            error('interhemi not supported');
        end
        
        % loop over the subjects
        for iSubj = 1:size(ccepData, 2)

            %nativeDistances = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances;
            %MNIfiles = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).MNIfiles;
            MNIlineIndices = ccepData(iSubj).rois(iTr).sub_tract(iSubTr).MNIlineIndices;

            
            % for each hemisphere, concatenate the MNI line indices (for either the inter, or individual hemispheres)
            for iHemi = 1:2
                rois(iTr).sub_tract(iSubTr).allMNIlineIndices{iHemi} = unique([rois(iTr).sub_tract(iSubTr).allMNIlineIndices{iHemi}, MNIlineIndices{iHemi}]);

                if ~isempty(MNIlineIndices{iHemi})
                    if isempty(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi})
                        rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi} = zeros(1, max(MNIlineIndices{iHemi}));
                    else
                        if max(MNIlineIndices{iHemi}) > length(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi})
                            rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi} = [rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi}, ...
                                                                                         zeros(1, max(MNIlineIndices{iHemi}) - length(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi}))];
                        end
                    end
                    rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi}(MNIlineIndices{iHemi}) = rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi}(MNIlineIndices{iHemi}) + 1;
                end
            end
            clear MNIlineIndices;
            
        end
        
        %
        for iHemi = 1:2
            
            % determine the 80th percentile of the number of subjects per tract-line (given only the relevant tract-lines)
            percentile80 = prctile(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi}(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi} > 0), 80);
            
            % take the mean over the number of subjects for the tract-lines that are above the 80th percentile
            mm = mean(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi}(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi} > percentile80));
            
            % allow tract-lines that have at least 1/3 of the "average" of the subjects
            mmlim = round(mm/3);
            
            % hist
            %{
            figure; bar(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi});
            title([strrep(rois(iTr).tract_name, '_', '\_'), ' - ',  rois(iTr).sub_tract(iSubTr).name, ' - ', num2str(iHemi), ' - top average: ', num2str(mm), ' - 1/3 rnd: ', num2str(mmlim)]);
            %}
            
            %
            rois(iTr).sub_tract(iSubTr).allMNIlineIndices{iHemi} = find(rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount{iHemi} > mmlim);
            
        end
    end
end
    


%%
%  load the MNI pial, inflated, and surface labels


% load mni305 pial
[Lmnipial_vert, Lmnipial_face] = read_surf(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'lh.pial'));
[Rmnipial_vert, Rmnipial_face] = read_surf(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'rh.pial'));

% load mni305 inflated
[Lmniinfl_vert, Lmniinfl_face] = read_surf(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'lh.inflated'));
[Rmniinfl_vert, Rmniinfl_face] = read_surf(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'rh.inflated'));

% read the annotations for the pial brain (a2009s = Destrieux) and relabel the vertices to the Destrieux atlas codes
dstrxCodeToLabel = {'G_and_S_frontomargin';'G_and_S_occipital_inf';'G_and_S_paracentral';'G_and_S_subcentral';'G_and_S_transv_frontopol';'G_and_S_cingul-Ant';'G_and_S_cingul-Mid-Ant';'G_and_S_cingul-Mid-Post';'G_cingul-Post-dorsal';'G_cingul-Post-ventral';'G_cuneus';'G_front_inf-Opercular';'G_front_inf-Orbital';'G_front_inf-Triangul';'G_front_middle';'G_front_sup';'G_Ins_lg_and_S_cent_ins';'G_insular_short';'G_occipital_middle';'G_occipital_sup';'G_oc-temp_lat-fusifor';'G_oc-temp_med-Lingual';'G_oc-temp_med-Parahip';'G_orbital';'G_pariet_inf-Angular';'G_pariet_inf-Supramar';'G_parietal_sup';'G_postcentral';'G_precentral';'G_precuneus';'G_rectus';'G_subcallosal';'G_temp_sup-G_T_transv';'G_temp_sup-Lateral';'G_temp_sup-Plan_polar';'G_temp_sup-Plan_tempo';'G_temporal_inf';'G_temporal_middle';'Lat_Fis-ant-Horizont';'Lat_Fis-ant-Vertical';'Lat_Fis-post';'Pole_occipital';'Pole_temporal';'S_calcarine';'S_central';'S_cingul-Marginalis';'S_circular_insula_ant';'S_circular_insula_inf';'S_circular_insula_sup';'S_collat_transv_ant';'S_collat_transv_post';'S_front_inf';'S_front_middle';'S_front_sup';'S_interm_prim-Jensen';'S_intrapariet_and_P_trans';'S_oc_middle_and_Lunatus';'S_oc_sup_and_transversal';'S_occipital_ant';'S_oc-temp_lat';'S_oc-temp_med_and_Lingual';'S_orbital_lateral';'S_orbital_med-olfact';'S_orbital-H_Shaped';'S_parieto_occipital';'S_pericallosal';'S_postcentral';'S_precentral-inf-part';'S_precentral-sup-part';'S_suborbital';'S_subparietal';'S_temporal_inf';'S_temporal_sup';'S_temporal_transverse'};
[~, Llabel, Lcolortable] = read_annotation(fullfile(FSsubjectsdir, 'fsaverage', 'label', 'lh.aparc.a2009s.annot'));
Lvert_label = nan(1, size(Llabel, 1));
for iLbl = 2:size(Lcolortable.table, 1)
    a = find(ismember(dstrxCodeToLabel, Lcolortable.struct_names{iLbl}));
    if ~isempty(a)
        Lvert_label(Llabel == Lcolortable.table(iLbl, 5)) = a;
    end
end

[~, Rlabel, Rcolortable] = read_annotation(fullfile(FSsubjectsdir, 'fsaverage', 'label', 'rh.aparc.a2009s.annot'));
Rvert_label = nan(1, size(Rlabel, 1));
for iLbl = 2:size(Rcolortable.table, 1)
    a = find(ismember(dstrxCodeToLabel, Rcolortable.struct_names{iLbl}));
    if ~isempty(a)
        Rvert_label(Rlabel == Rcolortable.table(iLbl, 5)) = a;
    end
end
clear Llabel Rlabel a iLbl;


      
%% 
%  Add all electrodes labels and left or right hemisphere into variables: allmni305_coords and allmni305_coords_infl

allmni305_coords        = [];
allmni305_coords_infl   = [];
allmni305_labels        = [];
allmni305_Destrlabels   = [];
allmni305_hemi          = [];

for iSubj = 1:length(ccepData)

    elecs = ccepData(iSubj).electrodes;
    
    Destrieux_label = elecs.Destrieux_label;
    if iscell(Destrieux_label)                  % TODO: this is because bad BIDS store/loading (unnecessary)
        for iElec = 1:size(Destrieux_label,1)
            if ischar(Destrieux_label{iElec})
                if isequal(Destrieux_label, 'n/a') 
                    Destrieux_label{iElec} = NaN;
                else
                    Destrieux_label{iElec} = str2double(Destrieux_label{iElec});
                end
            end
        end
        Destrieux_label = cell2mat(Destrieux_label);
    end
    
    mni305_coords           = [elecs.x, elecs.y, elecs.z];
    
    % 
    allmni305_labels        = [allmni305_labels; strcat(['s', ccepData(iSubj).id(end - 1:end), '-'], elecs.name)];
    allmni305_coords        = [allmni305_coords; mni305_coords];
    allmni305_Destrlabels   = [allmni305_Destrlabels; Destrieux_label];
    allmni305_hemi          = [allmni305_hemi; elecs.jsonHemi];
    
    % run through all coordinates and find the inflated points
    temp_inflated = NaN(size(mni305_coords));
    for iElec = 1:size(Destrieux_label, 1)
        
        if isequal(elecs.jsonHemi{iElec}, 'L')
            [~, min_ind] = min(sqrt(sum((Lmnipial_vert - mni305_coords(iElec, :)) .^ 2, 2)));
            temp_inflated(iElec, :) = Lmniinfl_vert(min_ind, :);
            
        elseif isequal(elecs.jsonHemi{iElec}, 'R')
            [~, min_ind] = min(sqrt(sum((Rmnipial_vert - mni305_coords(iElec, :)) .^ 2, 2)));
            temp_inflated(iElec, :) = Rmniinfl_vert(min_ind, :);
            
        end
        
    end
    allmni305_coords_infl = [allmni305_coords_infl; temp_inflated];
    
end
clear elecs temp_inflated;

%%
%  Label, for each (sub)tract, the ROIs for display (in color)

% loop over the tracts and sub-tracts
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        % re-label the vertex labels (0 = no ROIs, 1 = ROI1, and 2 = ROI2)
        lroi_label = nan(1, numel(Lvert_label));
        lroi_label(ismember(Lvert_label, rois(iTr).sub_tract(iSubTr).roi1)) = 1;
        lroi_label(ismember(Lvert_label, rois(iTr).sub_tract(iSubTr).roi2)) = 2;
        rois(iTr).sub_tract(iSubTr).Lvert_labels = lroi_label;
        
    end
    
end



%%
%  Plot figure with left pial with electrodes in mni space

v_d = [270 0];

% pour the faces and vertices into a gifti
gl.faces = Lmnipial_face + 1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);

% move the electrodes a bit away from the origin (0, 0, 0), to make sure electrodes pop out
a_offset = .1 * max(abs(allmni305_coords(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = allmni305_coords + repmat(a_offset, size(allmni305_coords, 1), 1);

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
%for iTr = 1:1
for iTr = 1:length(rois)
    %for iSubTr = 1:1
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        
        
        % load the (sub)tract file (is in MNI152 space)
        trkFile = fullfile(track_path, [rois(iTr).tract_name, '_L.trk.gz']);
        [fibers, idx] = ea_trk2ftr(trkFile, 1);
        
        %
        roi1elecs = ismember(allmni305_hemi, 'L') & ismember(allmni305_Destrlabels, rois(iTr).sub_tract(iSubTr).roi1);
        roi2elecs = ismember(allmni305_hemi, 'L') & ismember(allmni305_Destrlabels, rois(iTr).sub_tract(iSubTr).roi2);
        
        % open the MNI pial
        tic
        hFig = figure;
        set(hFig, 'Visible', 'off');
        %tH = ieeg_RenderGifti(gl);
        vLabels = rois(iTr).sub_tract(iSubTr).Lvert_labels;
        vLabels(isnan(vLabels)) = 0;
        tH = ieeg_RenderGiftiLabels(gl, vLabels', 'jet');
        set(tH,'FaceAlpha', .2) % make transparent
        
        %{
        toolConfig = {};
        toolConfig.hideToolWindow           = 1;
        toolConfig.yokeCam                  = 1;
        
        toolConfig.('overlay1')             = rois(iTr).sub_tract(iSubTr).Lvert_labels;
        toolConfig.('overlay1PosEnabled')   = 1;
        toolConfig.('overlay1PosColormap')  = 'green';
        toolConfig.('overlay1PosMin')       = 1;
        toolConfig.('overlay1PosMax')       = max(rois(iTr).sub_tract(iSubTr).Lvert_labels);
        toolConfig.('overlay1NegEnabled')   = 0;
        %toolConfig.pointSet1                = els(roi2elecs, :);
        %toolConfig.pointSet1Text            = string(allmni305_Destrlabels(roi2elecs));
        toolConfig.defaultBackgroundAlpha   = .6;
        mx.three_dimensional.giftiTools(gl, toolConfig);
        %}
        
        % add the (sub)tracts
        roisTrkLines = rois(iTr).sub_tract(iSubTr).allMNIlineIndices{1, 1}; % cell 1 = left
        hold on;
        startV = 1;
        for i = roisTrkLines
            if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
            endV = startV + idx(i) - 1;
            plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
        end
        hold off;
        
        % plot electrodes not part of the end-point ROIs 
        %ieeg_elAdd(els(ismember(allmni305_hemi, 'L') & ~ismember(allmni305_Destrlabels, [rois(iTr).sub_tract(iSubTr).roi1 rois(iTr).sub_tract(iSubTr).roi2]), :), 'k', 10)

        % set the electrode colormap
        colorMap                = [];
        colorMap.temporal       = [0 0 .8];
        colorMap.frontal        = [1 .8 0];
        colorMap.central        = [.8 .3 0];
        colorMap.parietal       = [0 .5 0];
        
        colorMap.ventral        = [0 0 .8];
        colorMap.dorsal         = [0 .5 0];
        colorMap.precentral     = [.8 .7 0];
        colorMap.postcentral    = [.8 .4 0];
        
        % plot the electrodes for the end-point ROIs
        subDir = lower(split(rois(iTr).sub_tract(iSubTr).name, '-'));
        ieeg_elAdd(els(roi1elecs, :), colorMap.(subDir{1}), 7)
        ieeg_elAdd(els(roi2elecs, :), colorMap.(subDir{2}), 7)
        ieeg_viewLight(v_d(1), v_d(2))

        % save the image
        figureName = fullfile(myDataPath.output, 'derivatives', 'render', ['leftMNIpial_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '.png']);
        set(gcf, 'PaperPositionMode', 'auto')
        set(hFig, 'Visible', 'on');
        print('-dpng', '-r300', figureName)
        close(hFig)
        toc
        
    end
end

return;


%% 
%   Plot figure with right pial with electrodes in mni space

v_d = [96 6];

% pour the faces and vertices into a gifti
gr.faces = Rmnipial_face + 1;
gr.vertices = Rmnipial_vert;
gr = gifti(gr);

% make sure electrodes pop out
a_offset = .5 * max(abs(allmni305_coords(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = allmni305_coords + repmat(a_offset, size(allmni305_coords, 1), 1);      

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)

        % load the (sub)tract file
        trkFile = fullfile(track_path, [rois(iTr).tract_name, '_R.trk.gz']);
        [fibers, idx] = ea_trk2ftr(trkFile, 1);
        
        % open the MNI pial
        figure
        tH = ieeg_RenderGifti(gr);
        set(tH,'FaceAlpha', .2) % make transparent
        %viewGii(gr, 'Trans.2')

        % add the (sub)tracts
        if rois(iTr).sub_tract(iSubTr).interHemi == 0
            roisTrkLines = rois(iTr).sub_tract(iSubTr).allMNIlineIndices{1, 2}; % cell 1 = left
        else
            roisTrkLines = rois(iTr).sub_tract(iSubTr).allMNIlineIndices{1, 1};
        end
        hold on;
        startV = 1;
        for i = roisTrkLines
            if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
            endV = startV + idx(i) - 1;
            plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
        end
        hold off;
        
        % plot electrodes not part of the end-point ROIs 
        %ieeg_elAdd(els(ismember(allmni305_hemi, 'R') & ~ismember(allmni305_Destrlabels,[roi_temporal roi_frontal roi_central roi_parietal]),:),'k',10)
        
        % plot the electrodes for the end-point ROIs
        ieeg_elAdd(els(ismember(allmni305_hemi, 'R') & ismember(allmni305_Destrlabels, rois(iTr).sub_tract(iSubTr).roi1), :), [0 0 .8], 15);
        ieeg_elAdd(els(ismember(allmni305_hemi, 'R') & ismember(allmni305_Destrlabels, rois(iTr).sub_tract(iSubTr).roi2), :), [0 .5 0], 15);
        ieeg_viewLight(v_d(1), v_d(2))

        % save the image
        figureName = fullfile(myDataPath.output, 'derivatives', 'render', ['rightMNIpial_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '.png']);
        set(gcf,'PaperPositionMode', 'auto')
        print('-dpng', '-r300', figureName);
        
    end
end


%%
%  Plot left inflated brain surface with electrodes in mni space

Lsulcal_labels = read_curv(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'lh.sulc'));

% make a colormap for the labels
cmap = Lcolortable.table(:, 1:3) ./ 256;

v_d = [270 0];

% pour the faces and vertices into a gifti
gl.faces = Lmniinfl_face + 1;
gl.vertices = Lmniinfl_vert;
gl = gifti(gl);

% make sure electrodes pop out
a_offset = .1 * max(abs(allmni305_coords_infl(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = allmni305_coords_infl+repmat(a_offset, size(allmni305_coords_infl, 1), 1);      
% els = allmni305_coords_infl;

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)

        figure
        tH = ieeg_RenderGifti(gl);
        set(tH,'FaceAlpha', .5);
        
        % with Destrieux labels:
        % tH = ieeg_RenderGiftiLabels(gl, Lvert_label,cmap, Lcolortable.struct_names);
        tH = ieeg_RenderGiftiLabels(gl, Lsulcal_labels, [.5 .5 .5; .8 .8 .8]);
        % sulci_rois = Lsulcal_labels;
        % sulci_rois(lroi_label == 1) = 3;
        % sulci_rois(lroi_label == 2) = 4;
        % tH = ieeg_RenderGiftiLabels(gl, sulci_rois, [.5 .5 .5; .8 .8 .8; 1 0 0; 0 1 0; 0 0 1]);

        % plot electrodes not part of the end-point ROIs
        %ieeg_elAdd(els(ismember(allmni305_hemi,'L') & ~ismember(allmni305_Destrlabels,[roi_temporal roi_frontal roi_central roi_parietal]),:),'k',10)
        
        % plot the electrodes for the end-point ROIs
        ieeg_elAdd(els(ismember(allmni305_hemi,'L') & ismember(allmni305_Destrlabels, rois(iTr).sub_tract(iSubTr).roi1), :), [0 0 .8], 15);
        ieeg_elAdd(els(ismember(allmni305_hemi,'L') & ismember(allmni305_Destrlabels, rois(iTr).sub_tract(iSubTr).roi2), :), [.8 .3 0], 15);
        ieeg_viewLight(v_d(1), v_d(2));
        
        %figureName = fullfile(myDataPath.output, 'derivatives', 'render', ['leftMNIinflated_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '.png']);
        set(gcf,'PaperPositionMode', 'auto');
        % print('-dpng', '-r300', figureName)
        
    end
end



%%
%  Plot right inflated brain surface with electrodes in mni space

v_d = [96 6];
Rsulcal_labels = read_curv(fullfile(FSsubjectsdir,'fsaverage', 'surf', 'rh.sulc'));

% make a colormap for the labels
cmap = Rcolortable.table(:, 1:3) ./ 256;

% pour the faces and vertices into a gifti
gr.faces = Rmniinfl_face + 1;
gr.vertices = Rmniinfl_vert;
gr = gifti(gr);

% make sure electrodes pop out
a_offset = .1 * max(abs(allmni305_coords_infl(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = allmni305_coords_infl + repmat(a_offset,size(allmni305_coords_infl, 1), 1);
% els = allmni305_coords_infl;

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)

        figure
        % tH = ieeg_RenderGifti(gl);
        % set(tH,'FaceAlpha',.5) % make transparent
        
        % with Destrieux labels:
        tH = ieeg_RenderGiftiLabels(gr,Rsulcal_labels,[.5 .5 .5;.8 .8 .8]);

        % plot electrodes not part of the end-point ROIs
        %ieeg_elAdd(els(ismember(allmni305_hemi,'R') & ~ismember(allmni305_Destrlabels,[roi_temporal roi_frontal roi_central roi_parietal]),:),'k',10)
        
        % plot the electrodes for the end-point ROIs
        ieeg_elAdd(els(ismember(allmni305_hemi,'R') & ismember(allmni305_Destrlabels, rois(iTr).sub_tract(iSubTr).roi1), :), [0 0 .8], 15);
        ieeg_elAdd(els(ismember(allmni305_hemi,'R') & ismember(allmni305_Destrlabels, rois(iTr).sub_tract(iSubTr).roi2), :), [.8 .3 0], 15);
        ieeg_viewLight(v_d(1), v_d(2));

        %figureName = fullfile(myDataPath.output, 'derivatives', 'render', ['rightMNIinflated_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '.png']);
        set(gcf,'PaperPositionMode','auto')
        % print('-dpng','-r300',figureName)

    end
end



%% 
%  Plot individual subjects rendering           
%
%  Note: this section cannot be run with the shared data alone, electrodes
%        coordinates in native space cannot be made public due to privacy
%        regulations

% TODO: perhaps rename rois here for this section to prevent conflict?
% load tracts and their corresponding end-point ROIs
rois = ccep_categorizeAnatomicalRegions();

elec_coords = [];

iSubj = 70; % in Fig1A of the article number 4 (4 years of age) and 70 (38 years of age) are used
%iSubj = 2;

disp(['sub ' int2str(iSubj)])

% subject freesurfer dir
FSdir = fullfile(myDataPath.input, 'derivatives', 'freesurfer', subjects(iSubj).name);

% get electrodes info
elec_coords(iSubj).elecs_tsv = readtable(fullfile(myDataPath.input, 'derivatives', 'native_electrodes', subjects(iSubj).name, ...
                                                [subjects(iSubj).name, '_', subjects(iSubj).ses, '_electrodes.tsv']), ...
                                                'FileType', 'text', 'Delimiter', '\t');
if iscell(elec_coords(iSubj).elecs_tsv.x)
    elecmatrix = NaN(size(elec_coords(iSubj).elecs_tsv, 1), 3);
    for iElec = 1:size(elec_coords(iSubj).elecs_tsv, 1)
        if ~isequal(elec_coords(iSubj).elecs_tsv.x{iElec}, 'n/a')
            elecmatrix(iElec, :) = [str2double(elec_coords(iSubj).elecs_tsv.x{iElec}) str2double(elec_coords(iSubj).elecs_tsv.y{iElec}) str2double(elec_coords(iSubj).elecs_tsv.z{iElec})];
        end
    end
else
    elecmatrix = [elec_coords(iSubj).elecs_tsv.x elec_coords(iSubj).elecs_tsv.y elec_coords(iSubj).elecs_tsv.z];
end
nElec = size(elecmatrix, 1);

% get hemisphere for each electrode
hemi = ccep_retrieveElecsHemisphere(fullfile(myDataPath.input, ccepData(iSubj).id, ccepData(iSubj).ses, 'ieeg', ...
                                             [ccepData(iSubj).id, '_', ccepData(iSubj).ses, '_task-SPESclin*_ieeg.json']), ...
                                    elecs_tsv);
                                
% TODO: check transformation steps below, from which to which
                                
% load mri orig header
origName = fullfile(FSdir, 'mri', 'orig.mgz');
orig = MRIread(origName, 'true');
Norig = orig.vox2ras; 
Torig = orig.tkrvox2ras;

% electrodes to freesurfer space
freeSurfer2T1 = inv(Norig * inv(Torig));
elCoords = freeSurfer2T1 * ([elecmatrix'; ones(1, nElec)]);
elCoords = elCoords(1:3, :)';

% subject pial
[Lsubpial_vert, Lsubpial_face] = read_surf(fullfile(FSdir, 'surf', 'lh.pial'));
[Rsubpial_vert, Rsubpial_face] = read_surf(fullfile(FSdir, 'surf', 'rh.pial'));

% set the view for the correct hemisphere
if isequal(hemi{1}, 'L')
    g.faces = Lsubpial_face + 1; % correct for zero index
    g.vertices = Lsubpial_vert;
    v_d = ([270 0]);
elseif isequal(hemi{1}, 'R')
    g.faces = Rsubpial_face + 1; % correct for zero index
    g.vertices = Rsubpial_vert;
    v_d = ([90 0]);
end

Destrieux_label = elec_coords(iSubj).elecs_tsv.Destrieux_label;
if iscell(Destrieux_label)
    for iElec = 1:size(Destrieux_label, 1)
        if ischar(Destrieux_label{iElec})
            if isequal(Destrieux_label, 'n/a') 
                Destrieux_label{iElec} = NaN;
            else
                Destrieux_label{iElec} = str2double(Destrieux_label{iElec});
            end
        end
    end
    Destrieux_label = cell2mat(Destrieux_label);
end

% make the electrodes be out of the brain cortex
a_offset = .1 * max(abs(elCoords(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = elCoords + repmat(a_offset, size(elCoords, 1), 1);

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)

        % load the (sub)tract file
        trkFile = fullfile(track_path, rois(iTr).tract_name);
        disp('Retrieving tracts in native space');

        % retrieve the distance between the stimulation and response end-point ROIs
        % for this particular patient given specific tracts
        [~, ~, trkLineIndices, trkNativeFibers] = ccep_retrieveInterROIDistance( ...
                                            rois(iTr).sub_tract(iSubTr).interHemi, ...
                                            trkFile, ...
                                            fullfile(myDataPath.input, 'derivatives', 'coreg_ANTs', ccepData(iSubj).id), ...
                                            fullfile(myDataPath.input, 'derivatives', 'freesurfer', ccepData(iSubj).id), ...
                                            rois(iTr).sub_tract(iSubTr).roi1, ...
                                            rois(iTr).sub_tract(iSubTr).roi2);


        
        % open the native pial
        figure
        tH = ieeg_RenderGifti(g);
        set(tH,'FaceAlpha', .2)

        % add the (sub)tracts
        if isequal(hemi{1}, 'L')
            roisTrkLines = trkLineIndices{1, 1};
            fibers = trkNativeFibers{1, 1}{1, 1};
            idx = trkNativeFibers{1, 1}{1, 2};
        else
            roisTrkLines = trkLineIndices{1, 2};
            fibers = trkNativeFibers{1, 2}{1, 1};
            idx = trkNativeFibers{1, 2}{1, 2};
        end
        hold on;
        startV = 1;
        for i = roisTrkLines
            if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
            endV = startV + idx(i) - 1;
            plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
        end
        hold off;
        
        %
        % ieeg_label(elCoords)
        
        % plot electrodes not part of the end-point ROIs
        %ieeg_elAdd(els(~ismember(Destrieux_label,[roi_temporal roi_frontal roi_central roi_parietal]), :), 'k', 20)
        
        % plot the electrodes for the end-point ROIs
        ieeg_elAdd(els(ismember(Destrieux_label, rois(iTr).sub_tract(iSubTr).roi1), :), [0 0 .8], 20)
        ieeg_elAdd(els(ismember(Destrieux_label, rois(iTr).sub_tract(iSubTr).roi2), :), [.8 .3 0], 20)

        
        ieeg_viewLight(v_d(1), v_d(2))

        pause(3)
        figureName = fullfile(myDataPath.output, 'derivatives', 'render', [subjects(iSubj).name  '_elsColors_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '.png']);
        set(gcf,'PaperPositionMode','auto');
        % print('-dpng','-r300',figureName)

    end
end
