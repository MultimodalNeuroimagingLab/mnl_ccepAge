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

if ~exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'n1Latencies_V1.mat'), 'file')
    error('Could not find/load _V1 file');
else
    
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'n1Latencies_V1.mat'), 'n1Latencies');
    
    % loop over the subjects
    for iSubj = 1:length(n1Latencies)
        fprintf('Load subj %d of %d (%s)\n', iSubj, length(n1Latencies), n1Latencies(iSubj).id);
        
        subjFSDir   = fullfile(myDataPath.input, 'derivatives', 'freesurfer', n1Latencies(iSubj).id);
        subjElecDir = fullfile(myDataPath.input, 'derivatives', 'native_electrodes', n1Latencies(iSubj).id);
        
        % (re)load the tract/ROIs end-point structure
        rois = ccep_categorizeAnatomicalRegions();
        
        % retrieve the electrodes and determine the hemisphere for each electrode
        electrodes = readtable(fullfile(subjElecDir, [n1Latencies(iSubj).id, '_', n1Latencies(iSubj).ses, '_electrodes.tsv']), ...
                               'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);
        elecHemi = ccep_retrieveElecsHemisphere(fullfile(myDataPath.input, n1Latencies(iSubj).id, n1Latencies(iSubj).ses, 'ieeg', ...
                                                     [n1Latencies(iSubj).id, '_', n1Latencies(iSubj).ses, '_task-SPESclin*_ieeg.json']), ...
                                                electrodes);
        
        % retrieve channels from the first run and determine the bad electrodes
        jsonFiles = dir(fullfile(myDataPath.input, n1Latencies(iSubj).id, n1Latencies(iSubj).ses, 'ieeg', ...
                                                     [n1Latencies(iSubj).id, '_', n1Latencies(iSubj).ses, '_task-SPESclin*_channels.tsv']));
        channels = readtable(fullfile(jsonFiles(1).folder, jsonFiles(1).name), ...
                             'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);
        bad_channels = find(~strcmpi(channels.status, 'good'));
        good_channels = setdiff(1:size(channels, 1), bad_channels);
        clear jsonFiles;

        % sort the electrodes so they match the channels
        electrodes = sortElectrodes(electrodes, channels, 0);

        %{
        % debug, brain & electrodes
        if any(contains(hemi, 'L')) && any(contains(hemi, 'R'))
            [vertexMatrix, facesMatrix] = mx.three_dimensional.merge3DObjects(gifti(fullfile(subjFSDir, 'pial.L.surf.gii')), gifti(fullfile(subjFSDir, 'pial.R.surf.gii')));
            gPial = gifti(struct('faces', facesMatrix, 'vertices', vertexMatrix));
            clear vertexMatrix facesMatrix;
        elseif any(contains(hemi, 'L'))
            gPial = gifti(fullfile(subjFSDir, 'pial.L.surf.gii'));
        elseif any(contains(hemi, 'R'))
            gPial = gifti(fullfile(subjFSDir, 'pial.R.surf.gii'));
        end
        viewGii(gPial, [electrodes(good_channels, :).x electrodes(good_channels, :).y electrodes(good_channels, :).z]);
        %}
        
        % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
        for iTr = 1:length(rois)
            for iSubTr = 1:length(rois(iTr).sub_tract)
                
                % 
                trkFile = fullfile(track_path, rois(iTr).tract_name);
                disp('Retrieving tract distance');

                % retrieve the distance between the stimulation and response end-point ROIs
                % for this particular patient given specific tracts
                [trkDist, trkFiles, trkLineIndices, trkNativeFibers, trkExtElecs] = ccep_retrieveInterROIDistance( ...
                                                                                                    rois(iTr).sub_tract(iSubTr).interHemi, ...
                                                                                                    trkFile, ...
                                                                                                    fullfile(myDataPath.input, 'derivatives', 'coreg_ANTs', n1Latencies(iSubj).id), ...
                                                                                                    subjFSDir, ...
                                                                                                    rois(iTr).sub_tract(iSubTr).roi1, ...
                                                                                                    rois(iTr).sub_tract(iSubTr).roi2, ...
                                                                                                    [electrodes.x electrodes.y electrodes.z], ...
                                                                                                    elecHemi);

                % store the distances, included MNI files and tract-lines, and elec
                rois(iTr).sub_tract(iSubTr).MNIfiles = trkFiles;
                rois(iTr).sub_tract(iSubTr).MNIlineIndices = trkLineIndices;
                rois(iTr).sub_tract(iSubTr).nativeDistances = trkDist;
                rois(iTr).sub_tract(iSubTr).extElecNames = cell(1, length(trkExtElecs));
                for iCell = 1:length(trkExtElecs)
                    rois(iTr).sub_tract(iSubTr).extElecNames{iCell} = electrodes(trkExtElecs{iCell}, :).name;
                end
                
                % store the native tract/end-point ROIs track-lines in output struct
                n1Latencies(iSubj).rois = rois;
                
            end
        end
        
    end
    
end


%%
%  (optionally) store the updated structure

s = input('Do you want to save the n1Latencies structure? [y/n]: ', 's');
if strcmp(s, 'y')
    save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'n1Latencies_V2.mat'), 'n1Latencies')
end

