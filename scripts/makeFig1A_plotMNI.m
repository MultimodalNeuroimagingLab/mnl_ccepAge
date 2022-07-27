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
myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

% get a list of datasets
subjects = ccep_getSubFilenameInfo(myDataPath);

% load the data struct to retrieve the (sub-)tracts that run between different the end-point ROIs
load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'n1Latencies_V2.mat'), 'n1Latencies');

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
%  Get standardized electrodes through surface based registration or linear
%  convert electrodes from patient's individual MRI to MNI305 space

% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.input, 'derivatives', 'freesurfer');

elec_coords = struct();

for iSub = 1:length(subjects) 
    disp(['subj ' int2str(iSub) ' of ' int2str(length(subjects))])
    
    % get electrodes info
    elec_coords(iSub).elecs_tsv = readtable(fullfile(myDataPath.input, subjects(iSub).name, subjects(iSub).ses, 'ieeg', ...
                                            [subjects(iSub).name, '_', subjects(iSub).ses, '_electrodes.tsv']), 'FileType', 'text', 'Delimiter', '\t');
    if iscell(elec_coords(iSub).elecs_tsv.x)
        elecmatrix = NaN(size(elec_coords(iSub).elecs_tsv, 1), 3);
        for ll = 1:size(elec_coords(iSub).elecs_tsv, 1)
            if ~isequal(elec_coords(iSub).elecs_tsv.x{ll}, 'n/a')
                elecmatrix(ll, :) = [str2double(elec_coords(iSub).elecs_tsv.x{ll}) str2double(elec_coords(iSub).elecs_tsv.y{ll}) str2double(elec_coords(iSub).elecs_tsv.z{ll})];
            end
        end
    else
        elecmatrix = [elec_coords(iSub).elecs_tsv.x elec_coords(iSub).elecs_tsv.y elec_coords(iSub).elecs_tsv.z];
    end
    nElec = size(elecmatrix, 1);
    
    % get hemisphere for each electrode
    these_json = dir(fullfile(myDataPath.input,subjects(iSub).name,subjects(iSub).ses, 'ieeg', [subjects(iSub).name, '_', subjects(iSub).ses, '_task-SPESclin*_ieeg.json']));
    ieeg_json = jsonread(fullfile(these_json(1).folder, these_json(1).name));
    if isequal(ieeg_json.iEEGPlacementScheme, 'left') || isequal(ieeg_json.iEEGPlacementScheme, 'left;')
        hemi = num2cell(repmat('L', nElec,1));
    elseif isequal(ieeg_json.iEEGPlacementScheme, 'right')|| isequal(ieeg_json.iEEGPlacementScheme, 'right;')
        hemi = num2cell(repmat('R', nElec,1));
    elseif contains(ieeg_json.iEEGPlacementScheme, {'left', 'right'}) % check with kk=17
        hemi = cell(nElec, 1);
        [hemi{:}] = deal('n/a');
        
        schemesplit = strsplit(ieeg_json.iEEGPlacementScheme, ';');
        rightcell = find(contains(schemesplit, 'right'));
        leftcell = find(contains(schemesplit, 'left'));
        
        if rightcell < leftcell
            leftcells = extractAfter(ieeg_json.iEEGPlacementScheme, 'left');
            rightcells = extractBetween(ieeg_json.iEEGPlacementScheme, 'right', 'left');
            rightcells = rightcells{:};
        else
            rightcells = extractAfter(ieeg_json.iEEGPlacementScheme, 'right');
            leftcells = extractBetween(ieeg_json.iEEGPlacementScheme, 'left', 'right');
            leftcells = leftcells{:};
        end
        
        leftelec = strsplit(leftcells, ';');
        leftelec =  leftelec(~cellfun('isempty', leftelec));
        rightelec = strsplit(rightcells, ';');
        rightelec = rightelec(~cellfun('isempty', rightelec));
        
        for elec=1:size(leftelec,2)
           C = strsplit(leftelec{elec}, {'[', ']'});
           elecInd = find(contains(elec_coords(iSub).elecs_tsv.name, C{1}));
           [hemi{elecInd}] = deal('L');
        end
        
        for elec=1:size(rightelec,2)
           C = strsplit(rightelec{elec}, {'[', ']'});
           elecInd = find(contains(elec_coords(iSub).elecs_tsv.name, C{1}));
           [hemi{elecInd}] = deal('R');
        end
    end
    elec_coords(iSub).hemi = hemi;
    
    % conversion to MNI already done and stored in tsv, just use
    elec_coords(iSub).mni_coords = [elec_coords(iSub).elecs_tsv.x, elec_coords(iSub).elecs_tsv.y, elec_coords(iSub).elecs_tsv.z];
    
end

save(fullfile(myDataPath.output, 'derivatives', 'elec_coordinatesMNI305.mat'), 'elec_coords')
% save(fullfile(myDataPath.output,'derivatives','elec_coordinatesMNI305lin.mat'),'elec_coords')



%%
%  Retrieve the tracts from all subjects
        
% loop over the subjects
for iSubj = 1:size(n1Latencies, 2)
    
    % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
    for iTr = 1:length(n1Latencies(iSubj).rois)
        for iSubTr = 1:length(n1Latencies(iSubj).rois(iTr).sub_tract)

            %nativeDistances = n1Latencies(iSubj).rois(iTr).sub_tract(iSubTr).nativeDistances;
            %MNIfiles = n1Latencies(iSubj).rois(iTr).sub_tract(iSubTr).MNIfiles;
            MNIlineIndices = n1Latencies(iSubj).rois(iTr).sub_tract(iSubTr).MNIlineIndices;
            
            % concatenate the MNI line indices (for either the inter, or individual hemispheres)
            if ~isfield(rois(iTr).sub_tract(iSubTr), 'allMNIlineIndices') || isempty(rois(iTr).sub_tract(iSubTr).allMNIlineIndices)
                rois(iTr).sub_tract(iSubTr).allMNIlineIndices = cell(1, length(MNIlineIndices));
            end
            for iHemi = 1:length(MNIlineIndices)
                rois(iTr).sub_tract(iSubTr).allMNIlineIndices{iHemi} = unique([rois(iTr).sub_tract(iSubTr).allMNIlineIndices{iHemi}, MNIlineIndices{iHemi}]);
            end
            clear MNIlineIndices;
        
        end
    end
    
end


%% 
%  Start here to make figures
%  load MNI electrode positions (saved in previous section), MNI sphere, pial, and surface labels

% linear is not so nice...
% load(fullfile(myDataPath.output,'derivatives','elec_coordinatesMNI305lin.mat'),'elec_coords')
% surface based is nice:
load(fullfile(myDataPath.output, 'derivatives', 'elec_coordinatesMNI305.mat'), 'elec_coords')

% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.input, 'derivatives', 'freesurfer');

% load mni305 pial
[Lmnipial_vert, Lmnipial_face] = read_surf(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'lh.pial'));
[Rmnipial_vert, Rmnipial_face] = read_surf(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'rh.pial'));

% load mni305 inflated
[Lmniinfl_vert, Lmniinfl_face] = read_surf(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'lh.inflated'));
[Rmniinfl_vert, Rmniinfl_face] = read_surf(fullfile(FSsubjectsdir, 'fsaverage', 'surf', 'rh.inflated'));

% surface labels
[Lvertices, Llabel, Lcolortable] = read_annotation(fullfile(FSsubjectsdir, 'fsaverage', 'label', 'lh.aparc.a2009s.annot'));
Lvert_label = Llabel; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for iSub = 1:size(Lcolortable.table, 1) % 76 are labels
    Lvert_label(Llabel == Lcolortable.table(iSub, 5)) = iSub;
end
[Rvertices, Rlabel, Rcolortable] = read_annotation(fullfile(FSsubjectsdir, 'fsaverage', 'label', 'rh.aparc.a2009s.annot'));
Rvert_label = Rlabel; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for iSub = 1:size(Rcolortable.table, 1) % 76 are labels
    Rvert_label(Rlabel == Rcolortable.table(iSub, 5)) = iSub;
end



%% 
%  Add all electrodes labels and left or right hemisphere into one  variable: allmni_coords and allmni_coords_infl

allmni_coords = [];
allmni_coords_infl = [];

allmni_labels = [];
all_hemi = [];
for iSub = 1:length(elec_coords)
    Destrieux_label = elec_coords(iSub).elecs_tsv.Destrieux_label;
    if iscell(Destrieux_label)
        for ll = 1:size(Destrieux_label,1)
            if ischar(Destrieux_label{ll})
                if isequal(Destrieux_label, 'n/a') 
                    Destrieux_label{ll} = NaN;
                else
                    Destrieux_label{ll} = str2double(Destrieux_label{ll});
                end
            end
        end
        Destrieux_label = cell2mat(Destrieux_label);
    end
    
    allmni_coords   = [allmni_coords; elec_coords(iSub).mni_coords];
    allmni_labels   = [allmni_labels; Destrieux_label];
    all_hemi        = [all_hemi; elec_coords(iSub).hemi];
    
    % run through all coordinates and find the inflated points
    temp_inflated = NaN(size(elec_coords(iSub).mni_coords));
    for ll = 1:size(Destrieux_label, 1)
        if isequal(elec_coords(iSub).hemi{ll}, 'L')
            [~, min_ind] = min(sqrt(sum((Lmnipial_vert-elec_coords(iSub).mni_coords(ll, :)) .^ 2, 2)));
            temp_inflated(ll, :) = Lmniinfl_vert(min_ind, :);
        elseif isequal(elec_coords(iSub).hemi{ll}, 'R')
            [~, min_ind] = min(sqrt(sum((Rmnipial_vert-elec_coords(iSub).mni_coords(ll, :)) .^ 2, 2)));
            temp_inflated(ll, :) = Rmniinfl_vert(min_ind, :);
        end
    end
    
    allmni_coords_infl = [allmni_coords_infl; temp_inflated];
end



%%
%  Label, for each (sub)tract, the ROIs for display (in color)

% loop over the tracts and sub-tracts
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        lroi_label = Lvert_label;
        lroi_label(ismember(lroi_label, rois(iTr).sub_tract(iSubTr).roi1 + 1)) = 200;
        lroi_label(ismember(lroi_label, rois(iTr).sub_tract(iSubTr).roi2 + 1)) = 300;
        lroi_label(lroi_label < 100) = 0;
        lroi_label = lroi_label / 100;
        rois(iTr).sub_tract(iSubTr).vert_labels = lroi_label;
        
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
a_offset = .1 * max(abs(allmni_coords(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords + repmat(a_offset, size(allmni_coords, 1), 1);

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)
        
        % load the (sub)tract file
        trkFile = fullfile(track_path, [rois(iTr).tract_name, '_L.trk.gz']);
        [fibers, idx] = ea_trk2ftr(trkFile, 1);

        % open the MNI pial
        figure
        tH = ieeg_RenderGifti(gl);
        set(tH,'FaceAlpha', .2) % make transparent
        %viewGii(gl, 'Trans.2')

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
        %ieeg_elAdd(els(ismember(all_hemi, 'L') & ~ismember(allmni_labels, [roi_temporal roi_frontal roi_central roi_parietal]), :), 'k', 10)
        
        % plot the electrodes for the end-point ROIs
        ieeg_elAdd(els(ismember(all_hemi, 'L') & ismember(allmni_labels, rois(iTr).sub_tract(iSubTr).roi1), :), [0 0 .8], 15)
        ieeg_elAdd(els(ismember(all_hemi, 'L') & ismember(allmni_labels, rois(iTr).sub_tract(iSubTr).roi2), :), [.8 .3 0], 15)
        ieeg_viewLight(v_d(1), v_d(2))

        % save the image
        figureName = fullfile(myDataPath.output, 'derivatives', 'render', ['leftMNIpial_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '.png']);
        set(gcf, 'PaperPositionMode', 'auto')
        print('-dpng', '-r300', figureName)

    end
end



%% 
%   Plot figure with right pial with electrodes in mni space

v_d = [96 6];

% pour the faces and vertices into a gifti
gr.faces = Rmnipial_face + 1;
gr.vertices = Rmnipial_vert;
gr = gifti(gr);

% make sure electrodes pop out
a_offset = .5 * max(abs(allmni_coords(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords + repmat(a_offset, size(allmni_coords, 1), 1);      

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
        %ieeg_elAdd(els(ismember(all_hemi, 'R') & ~ismember(allmni_labels,[roi_temporal roi_frontal roi_central roi_parietal]),:),'k',10)
        
        % plot the electrodes for the end-point ROIs
        ieeg_elAdd(els(ismember(all_hemi, 'R') & ismember(allmni_labels, rois(iTr).sub_tract(iSubTr).roi1), :), [0 0 .8], 15);
        ieeg_elAdd(els(ismember(all_hemi, 'R') & ismember(allmni_labels, rois(iTr).sub_tract(iSubTr).roi2), :), [0 .5 0], 15);
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
a_offset = .1 * max(abs(allmni_coords_infl(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords_infl+repmat(a_offset, size(allmni_coords_infl, 1), 1);      
% els = allmni_coords_infl;

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
        % sulci_rois(lroi_label == 2) = 3;
        % sulci_rois(lroi_label == 3) = 4;
        % tH = ieeg_RenderGiftiLabels(gl, sulci_rois, [.5 .5 .5; .8 .8 .8; 1 0 0; 0 1 0; 0 0 1]);

        % plot electrodes not part of the end-point ROIs
        %ieeg_elAdd(els(ismember(all_hemi,'L') & ~ismember(allmni_labels,[roi_temporal roi_frontal roi_central roi_parietal]),:),'k',10)
        
        % plot the electrodes for the end-point ROIs
        ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels, rois(iTr).sub_tract(iSubTr).roi1), :), [0 0 .8], 15);
        ieeg_elAdd(els(ismember(all_hemi,'L') & ismember(allmni_labels, rois(iTr).sub_tract(iSubTr).roi2), :), [.8 .3 0], 15);
        ieeg_viewLight(v_d(1), v_d(2));
        
        %figureName = fullfile(myDataPath.output, 'derivatives', 'render', ['leftMNIinflated_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', ''), '.png']);
        set(gcf,'PaperPositionMode', 'auto');
        % print('-dpng', '-r300', figureName)
        
    end
end



%%
%  Plot right inflated brain surface with electrodes in mni space

v_d = [96 6];
Rsulcal_labels = read_curv(fullfile(FSsubjectsdir,'fsaverage','surf','rh.sulc'));

% make a colormap for the labels
cmap = Rcolortable.table(:, 1:3) ./ 256;

% pour the faces and vertices into a gifti
gr.faces = Rmniinfl_face + 1;
gr.vertices = Rmniinfl_vert;
gr = gifti(gr);

% make sure electrodes pop out
a_offset = .1 * max(abs(allmni_coords_infl(:, 1))) * [cosd(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(1) - 90) * cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords_infl + repmat(a_offset,size(allmni_coords_infl, 1), 1);
% els = allmni_coords_infl;

% loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
for iTr = 1:length(rois)
    for iSubTr = 1:length(rois(iTr).sub_tract)

        figure
        % tH = ieeg_RenderGifti(gl);
        % set(tH,'FaceAlpha',.5) % make transparent
        
        % with Destrieux labels:
        tH = ieeg_RenderGiftiLabels(gr,Rsulcal_labels,[.5 .5 .5;.8 .8 .8]); %#ok<NASGU>

        % plot electrodes not part of the end-point ROIs
        %ieeg_elAdd(els(ismember(all_hemi,'R') & ~ismember(allmni_labels,[roi_temporal roi_frontal roi_central roi_parietal]),:),'k',10)
        
        % plot the electrodes for the end-point ROIs
        ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels, rois(iTr).sub_tract(iSubTr).roi1), :), [0 0 .8], 15);
        ieeg_elAdd(els(ismember(all_hemi,'R') & ismember(allmni_labels, rois(iTr).sub_tract(iSubTr).roi2), :), [.8 .3 0], 15);
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

% categorize anatomical regions
% TODO: perhaps rename rois here for this section to prevent conflict?
rois = ccep_categorizeAnatomicalRegions();

% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.input, 'derivatives', 'freesurfer');

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
    for ll = 1:size(elec_coords(iSubj).elecs_tsv, 1)
        if ~isequal(elec_coords(iSubj).elecs_tsv.x{ll}, 'n/a')
            elecmatrix(ll, :) = [str2double(elec_coords(iSubj).elecs_tsv.x{ll}) str2double(elec_coords(iSubj).elecs_tsv.y{ll}) str2double(elec_coords(iSubj).elecs_tsv.z{ll})];
        end
    end
else
    elecmatrix = [elec_coords(iSubj).elecs_tsv.x elec_coords(iSubj).elecs_tsv.y elec_coords(iSubj).elecs_tsv.z];
end
nElec = size(elecmatrix, 1);

% get hemisphere for each electrode
these_json = dir(fullfile(myDataPath.input,subjects(iSubj).name,subjects(iSubj).ses, 'ieeg', [subjects(iSubj).name, '_', subjects(iSubj).ses, '_task-SPESclin*_ieeg.json']));
ieeg_json = jsonread(fullfile(these_json(1).folder, these_json(1).name));
if isequal(ieeg_json.iEEGPlacementScheme, 'left') || isequal(ieeg_json.iEEGPlacementScheme, 'left;')
    hemi = num2cell(repmat('L', nElec, 1));
elseif isequal(ieeg_json.iEEGPlacementScheme, 'right') || isequal(ieeg_json.iEEGPlacementScheme, 'right;')
    hemi = num2cell(repmat('R', nElec, 1));
elseif contains(ieeg_json.iEEGPlacementScheme, {'left', 'right'}) % check with kk=17
    hemi = cell(nElec, 1);
    [hemi{:}] = deal('n/a');
    
    schemesplit = strsplit(ieeg_json.iEEGPlacementScheme,';'); 
    rightcell = find(contains(schemesplit, 'right'));
    leftcell = find(contains(schemesplit, 'left'));
    
    if rightcell < leftcell
        leftcells = extractAfter(ieeg_json.iEEGPlacementScheme, 'left');
        rightcells = extractBetween(ieeg_json.iEEGPlacementScheme, 'right', 'left');
        rightcells = rightcells{:};
    else
        rightcells = extractAfter(ieeg_json.iEEGPlacementScheme, 'right');
        leftcells = extractBetween(ieeg_json.iEEGPlacementScheme, 'left', 'right');
        leftcells = leftcells{:};
    end
    
    leftelec = strsplit(leftcells,';');
    leftelec =  leftelec(~cellfun('isempty',leftelec));
    rightelec = strsplit(rightcells,';');
    rightelec = rightelec(~cellfun('isempty',rightelec));
    
    % set L in variable hemi for electrodes in the left hemisphere
    for elec=1:size(leftelec,2)
        C = strsplit(leftelec{elec},{'[',']'});
        elecInd = find(contains(elec_coords(iSubj).elecs_tsv.name,C{1}));
        [hemi{elecInd}] = deal('L');
    end
    
    % set R in variable hemi for electrodes in the right hemisphere
    for elec=1:size(rightelec,2)
        C = strsplit(rightelec{elec},{'[',']'});
        elecInd = find(contains(elec_coords(iSubj).elecs_tsv.name,C{1}));
        [hemi{elecInd}] = deal('R');
    end
end

% number of electrodes
nElec = size(elecmatrix, 1);

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
    for ll = 1:size(Destrieux_label, 1)
        if ischar(Destrieux_label{ll})
            if isequal(Destrieux_label, 'n/a') 
                Destrieux_label{ll} = NaN;
            else
                Destrieux_label{ll} = str2double(Destrieux_label{ll});
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
                                            fullfile(myDataPath.input, 'derivatives', 'coreg_ANTs', n1Latencies(iSubj).id), ...
                                            fullfile(myDataPath.input, 'derivatives', 'freesurfer', n1Latencies(iSubj).id), ...
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
