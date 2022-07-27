%
% Script to add the (sub-)tracts lines and distance between the end-point ROIs
% 
% Max van den Boom, MultimodalNeuroimaging Lab (MN), 2022
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
    for iSubj = 1:size(n1Latencies, 2)
        fprintf('Load subj %d of %d \n', iSubj, size(n1Latencies, 2))
    
        % (re)load the tract/ROIs end-point structure
        rois = ccep_categorizeAnatomicalRegions();
        
        % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
        for iTr = 1:length(rois)
            for iSubTr = 1:length(rois(iTr).sub_tract)
                
                
                % Note: currently not considering the electrode positions to
                % select tract-lines, just using end-point ROIs
                

                trkFile = fullfile(track_path, rois(iTr).tract_name);
                disp('Retrieving tract distance');

                % retrieve the distance between the stimulation and response end-point ROIs
                % for this particular patient given specific tracts
                [trkDist, trkFiles, trkLineIndices] = ccep_retrieveInterROIDistance( ...
                                                    rois(iTr).sub_tract(iSubTr).interHemi, ...
                                                    trkFile, ...
                                                    fullfile(myDataPath.input, 'derivatives', 'coreg_ANTs', n1Latencies(iSubj).id), ...
                                                    fullfile(myDataPath.input, 'derivatives', 'freesurfer', n1Latencies(iSubj).id), ...
                                                    rois(iTr).sub_tract(iSubTr).roi1, ...
                                                    rois(iTr).sub_tract(iSubTr).roi2);

                % store the distances and included MNI files and tract-lines
                rois(iTr).sub_tract(iSubTr).MNIfiles = trkFiles;
                rois(iTr).sub_tract(iSubTr).MNIlineIndices = trkLineIndices;
                rois(iTr).sub_tract(iSubTr).nativeDistances = trkDist;
                
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

