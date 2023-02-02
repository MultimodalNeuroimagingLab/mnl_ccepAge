%
% This script produces the images for:
%    - Figure 1A - that depict a MNI surface with tracts, subtracts, ROIs and electrodes
%    - Supplemental Figure 9 - distribution of electrodes over different age groups
%
% Max van den Boom, Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, 2022
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

% loop over the tracts (SLF, AF, etc...) and sub-tracts
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        % init fields
        rois(iTr).sub_tract(iSubTr).allMNIlineIndices = cell(1, 2);
        rois(iTr).sub_tract(iSubTr).allMNIlineIndicesCount = cell(1, 2);
        
        % loop over the subjects
        for iSubj = 1:size(ccepData, 2)

            % retrieve the lines for this subject
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
        
        % each hemisphere
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

allmni305_coords            = [];
allmni305_subjIds           = [];
allmni305_labels            = [];
allmni305_destrLabels       = [];
allmni305_destrLabelTexts   = [];
allmni305_hemi              = [];
allmni305_age               = [];
allmni305_age2               = [];
short_ages                  = nan(1, length(ccepData));

for iSubj = 1:length(ccepData)

    elecs = ccepData(iSubj).electrodes;
    elecInclude = find(ismember(lower(elecs.group), {'strip', 'grid'}));
    elecs = elecs(elecInclude, :);
    
    % 
    allmni305_coords            = [allmni305_coords; [elecs.x, elecs.y, elecs.z]];
    allmni305_subjIds           = [allmni305_subjIds; repmat(str2double(ccepData(iSubj).id(end - 1:end)), size(elecs, 1), 1)];
    allmni305_labels            = [allmni305_labels; strcat(['s', ccepData(iSubj).id(end - 1:end), '-'], elecs.name)];
    allmni305_destrLabels       = [allmni305_destrLabels; elecs.Destrieux_label];
    allmni305_destrLabelTexts   = [allmni305_destrLabelTexts; elecs.Destrieux_label_text];
    allmni305_hemi              = [allmni305_hemi; elecs.jsonHemi];
    allmni305_age               = [allmni305_age; repmat(ccepData(iSubj).age, size(elecs, 1), 1)];
    allmni305_age2              = [allmni305_age2; strcat(['sub', num2str(iSubj), '-'], string(repmat(ccepData(iSubj).age, size(elecs, 1), 1)))];
    
    short_ages(iSubj)           = ccepData(iSubj).age;
    
end
clear elecs;


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
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)

        % load the (sub)tract file (is in MNI152 space)
        trkFile = fullfile(track_path, [rois(iTr).tract_name, '_L.trk.gz']);
        [fibers, idx] = ea_trk2ftr(trkFile, 1);
        
        %
        roi1elecs = ismember(allmni305_hemi, 'L') & ismember(allmni305_destrLabels, rois(iTr).sub_tract(iSubTr).roi1);
        roi2elecs = ismember(allmni305_hemi, 'L') & ismember(allmni305_destrLabels, rois(iTr).sub_tract(iSubTr).roi2);
        
        % open the MNI pial
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
        %toolConfig.pointSet1Text            = string(allmni305_destrLabels(roi2elecs));
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
        %ieeg_elAdd(els(ismember(allmni305_hemi, 'L') & ~ismember(allmni305_destrLabels, [rois(iTr).sub_tract(iSubTr).roi1 rois(iTr).sub_tract(iSubTr).roi2]), :), 'k', 10)

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
        if ~exist(fullfile(myDataPath.output, 'derivatives', 'images'), 'dir')
            mkdir(fullfile(myDataPath.output, 'derivatives', 'images'));
        end
        figureName = fullfile(myDataPath.output, 'derivatives', 'images', ['leftMNIpial_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '.png']);
        set(gcf, 'PaperPositionMode', 'auto')
        set(hFig, 'Visible', 'on');
        print('-dpng', '-r300', figureName)
        close(hFig)
        
    end
end



%%
%  Print number of electrodes and generate supplement Figure 9

edges = [0 10 20 30 40 50 60];

% count the total number of connections
numPlots = 0;
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        numPlots = numPlots + 1;
    end
end

% open figure
hFigElec = figure('Position', [0 0 2400 800]);
elecPlotCounter = 1;

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)

        %
        subNames = split(rois(iTr).sub_tract(iSubTr).name, '-');
        disp(' ');
        disp(['----   ', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '   ----']);
        
        %
        roi1Elecs = ismember(allmni305_destrLabels, rois(iTr).sub_tract(iSubTr).roi1);
        roi2Elecs = ismember(allmni305_destrLabels, rois(iTr).sub_tract(iSubTr).roi2);
        
        % determine the unique ages
        uniqueAges = unique([allmni305_age2(roi1Elecs); allmni305_age2(roi2Elecs)]);
        for iAge = 1:length(uniqueAges)
            uAge = split(uniqueAges(iAge), '-');
            uniqueAges(iAge) = str2double(uAge(2));
        end
        uniqueAges = str2double(uniqueAges);
        
        % calculate relative count
        [hstCountAges, ~] = histcounts(uniqueAges, edges);
        [hstCountRoi1, ~] = histcounts(allmni305_age(roi1Elecs), edges);
        [hstCountRoi2, ~] = histcounts(allmni305_age(roi2Elecs), edges);
        hstRelRoi1 = hstCountRoi1 ./ hstCountAges;
        hstRelRoi2 = hstCountRoi2 ./ hstCountAges;
        
        % display descriptives
        disp(['--     ', subNames{1}, '     --']);
        disp(['labels: ', strjoin(unique(allmni305_destrLabelTexts(roi1Elecs)), ', ')]);
        disp(['total # elec: ', num2str(length(allmni305_age(roi1Elecs)))]);
        disp(['total # subjects with coverage: ', num2str(length(unique(allmni305_subjIds(roi1Elecs))))]);
        %disp(['median elecs per subject: ', num2str(length(unique(allmni305_subjIds(roi1Elecs))))]);
        
        
        disp(['--     ', subNames{2}, '     --']);
        disp(['labels: ', strjoin(unique(allmni305_destrLabelTexts(roi2Elecs)), ', ')]);
        disp(['total # elec: ', num2str(length(allmni305_age(roi2Elecs)))]);
        disp(['total # subjects with coverage: ', num2str(length(unique(allmni305_subjIds(roi2Elecs))))]);
        
        % ROI 1 - electrode per age bin
        subplot(3, numPlots, elecPlotCounter, 'Parent', hFigElec);
        histogram(allmni305_age(ismember(allmni305_destrLabels, rois(iTr).sub_tract(iSubTr).roi1)), edges)
        title({strrep(rois(iTr).tract_name, '_', '\_'), ' ', subNames{1}});
        xlim([edges(1), edges(end)]);
        set(gca, 'XTick', edges)
        if elecPlotCounter == 1
           ylabel('#electrodes'); 
        end

        %{
        % ROI 1 - electrodes per subject
        yyaxis right
        for iEdge = 1:length(edges) - 1
            line([edges(iEdge) edges(iEdge + 1)], [hstRelRoi1(iEdge), hstRelRoi1(iEdge)]);
        end
        if elecPlotCounter == numPlots
           ylabel('#electrodes per subject'); 
        end
        %}
        
        % ROI 2 - electrode per age bin
        subplot(3, numPlots, numPlots + elecPlotCounter, 'Parent', hFigElec);
        histogram(allmni305_age(ismember(allmni305_destrLabels, rois(iTr).sub_tract(iSubTr).roi2)), edges)
        title(subNames{2});
        xlim([edges(1), edges(end)]);
        set(gca, 'XTick', edges)
        if elecPlotCounter == 1
           ylabel('# electrodes'); 
        end

        %{
        % ROI 2 - electrodes per subject
        yyaxis right
        for iEdge = 1:length(edges) - 1
            line([edges(iEdge) edges(iEdge + 1)], [hstRelRoi2(iEdge), hstRelRoi2(iEdge)]);
        end
        if elecPlotCounter == numPlots
           ylabel('#electrodes per subject'); 
        end
        %}
        
        %
        if iTr == 1 && iSubTr == 1
            subplot(3, numPlots, numPlots * 2 + elecPlotCounter, 'Parent', hFigElec);
            histogram(short_ages, edges, 'FaceColor', [.8 .8 0])
            %title({strrep(rois(iTr).tract_name, '_', '\_'), ' ', subNames{1}});
            xlim([edges(1), edges(end)]);
            set(gca, 'XTick', edges)
            if elecPlotCounter == 1
               ylabel('# subjects'); 
            end
        end
        
        
        elecPlotCounter = elecPlotCounter + 1;
        
    end
end

% save supplement 9 figure
figure(hFigElec);
set(hFigElec, 'renderer', 'Painters')
set(hFigElec, 'PaperPositionMode', 'auto')

if ~exist(fullfile(myDataPath.output, 'derivatives', 'images'), 'dir')
    mkdir(fullfile(myDataPath.output, 'derivatives', 'images'));
end
figureName = fullfile(myDataPath.output, 'derivatives', 'images', 'SupFigS9_ElecCoverage');
print('-dpng', '-r300', figureName)
print('-depsc', '-r300', figureName)
%close(hFigElec)
