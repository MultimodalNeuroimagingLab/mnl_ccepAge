%
% Script to add the (sub-)tracts lines and distance between the end-point ROIs
% 
% Max van den Boom, MultimodalNeuroimaging Lab (MNL), 2022
%


%% 
%  Set paths
clc
clear
myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');


%%
%  

if ~exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'file')
    error('Could not find/load _V1 file');
else
    
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'ccepData');
    
    % loop over the subjects
    %for iSubj = 1:length(ccepData)
    for iSubj = 2:2
        fprintf('Load subj %d of %d (%s)\n', iSubj, length(ccepData), ccepData(iSubj).id);
        
        subjFSDir   = fullfile(myDataPath.input, 'derivatives', 'freesurfer', ccepData(iSubj).id);
        subjElecDir = fullfile(myDataPath.input, 'derivatives', 'native_electrodes', ccepData(iSubj).id);
        
        % (re)load the tract/ROIs end-point structure
        rois = ccep_categorizeAnatomicalRegions();
        
        % retrieve channels from the first run and determine the bad electrodes
        jsonFiles = dir(fullfile(myDataPath.input, ccepData(iSubj).id, ccepData(iSubj).ses, 'ieeg', ...
                                                     [ccepData(iSubj).id, '_', ccepData(iSubj).ses, '_task-SPESclin*_channels.tsv']));
        channels = readtable(fullfile(jsonFiles(1).folder, jsonFiles(1).name), ...
                             'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);
        bad_channels = find(~strcmpi(channels.status, 'good'));
        good_channels = setdiff(1:size(channels, 1), bad_channels);
        clear jsonFiles;
        
        % retrieve the electrodes and sort the electrodes so they match the channels
        electrodes = ccepData(iSubj).elecs;
        electrodes = ccep_sortElectrodes(electrodes, channels, 0);

        %{
        % debug, brain & electrodes
        if any(contains(electrodes.jsonHemi, 'L')) && any(contains(electrodes.jsonHemi, 'R'))
            [vertexMatrix, facesMatrix] = mx.three_dimensional.merge3DObjects(gifti(fullfile(subjFSDir, 'pial.L.surf.gii')), gifti(fullfile(subjFSDir, 'pial.R.surf.gii')));
            gPial = gifti(struct('faces', facesMatrix, 'vertices', vertexMatrix));
            clear vertexMatrix facesMatrix;
        elseif any(contains(electrodes.jsonHemi, 'L'))
            gPial = gifti(fullfile(subjFSDir, 'pial.L.surf.gii'));
        elseif any(contains(electrodes.jsonHemi, 'R'))
            gPial = gifti(fullfile(subjFSDir, 'pial.R.surf.gii'));
        end
        viewGii(gPial, [electrodes(good_channels, :).x electrodes(good_channels, :).y electrodes(good_channels, :).z]);
        %}
        
        % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
        %for iTr = 1:length(rois)
        for iTr = 2:2
            for iSubTr = 1:length(rois(iTr).sub_tract)
                
                % 
                trkFile = fullfile(track_path, rois(iTr).tract_name);
                disp('Retrieving tract distance');

                % retrieve the distance between the stimulation and response end-point ROIs
                % for this particular patient given specific tracts
                [trkDist, trkFiles, trkLineIndices, trkNativeFibers, trkExtElecs] = ccep_retrieveInterROIDistance( ...
                                                                                                    rois(iTr).sub_tract(iSubTr).interHemi, ...
                                                                                                    trkFile, ...
                                                                                                    fullfile(myDataPath.input, 'derivatives', 'coreg_ANTs', ccepData(iSubj).id), ...
                                                                                                    subjFSDir, ...
                                                                                                    rois(iTr).sub_tract(iSubTr).roi1, ...
                                                                                                    rois(iTr).sub_tract(iSubTr).roi2, ...
                                                                                                    [electrodes.x electrodes.y electrodes.z], ...
                                                                                                    electrodes.jsonHemi);

                % store the distances, included MNI files and tract-lines, and elec
                rois(iTr).sub_tract(iSubTr).MNIfiles = trkFiles;
                rois(iTr).sub_tract(iSubTr).MNIlineIndices = trkLineIndices;
                rois(iTr).sub_tract(iSubTr).nativeDistances = trkDist;
                rois(iTr).sub_tract(iSubTr).extElecNames = cell(1, length(trkExtElecs));
                for iCell = 1:length(trkExtElecs)
                    rois(iTr).sub_tract(iSubTr).extElecNames{iCell} = electrodes(trkExtElecs{iCell}, :).name;
                end
                
                % store the native tract/end-point ROIs track-lines in output struct
                ccepData(iSubj).rois = rois;
                
            end
        end
        
    end
    
end


%%
%  (optionally) store the updated structure

s = input('Do you want to save the ccepData structure? [y/n]: ', 's');
if strcmp(s, 'y')
    save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData')
end

