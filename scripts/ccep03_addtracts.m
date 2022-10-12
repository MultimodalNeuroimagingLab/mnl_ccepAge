%
%  Script to add the (sub-)tracts lines and distance between the end-point ROIs
%
%
%   Note:   This script requires the freesurfer output of the individual subjects to run. Unfortunately, for privacy
%           reasons, we are not allowed to share the freesurfer output files. Instead, the file that this script would
%           produce, ccepData_V2.mat, is shared and can be used directly in the subsequent scripts
%
%   Note 2: For Linux and Mac, execute permissions might need to be set for
%           the leadDBS binary files ('/leadDBS/ext_libs/ANTs').
%           Navigate to that folder in the terminal and run 'chmod u+x antsApplyTransformsToPoints.*'
%
%   Max van den Boom, MultimodalNeuroimaging Lab (MNL), 2022
%


%% 
%  Set paths
clc
clear
warning('on');
warning('backtrace', 'off')

myDataPath = setLocalDataPath(1);
track_path = fullfile(myDataPath.input, 'sourcedata', 'tracks');

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
%  

if ~exist(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'file')
    error('Could not find/load _V1 file');
else
    
    load(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V1.mat'), 'ccepData');
    
    % loop over the subjects
    for iSubj = 1:length(ccepData)
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
        

        electrodes = [];
        elecPositions = [];
        if isfolder(subjElecDir)
            
            % load the electrodes (in native space)
            electrodes = readtable(fullfile(subjElecDir, [ccepData(iSubj).id, '_', ccepData(iSubj).ses, '_electrodes.tsv']), ...
                                   'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);

            % check if the native electrode order matches the MNI electrodes by name
            if length(electrodes.name) ~= length(ccepData(iSubj).electrodes.name) || sum(~strcmp(electrodes.name, ccepData(iSubj).electrodes.name)) > 0
                error('The native electrodes files does not match the MNI electrodes file');
            end
            
            % calculate and store the distances (in native space) between the electrodes
            dist = (electrodes.x' - electrodes.x) .^ 2 + ...
                   (electrodes.y' - electrodes.y) .^ 2 + ...
                   (electrodes.z' - electrodes.z) .^ 2;
            dist = sqrt(dist);
            ccepData(iSubj).nativeElecDistances = dist;

            % retrieve the electrodes and sort the electrodes so they match the channels
            % Note: this might remove electrode rows (if they are not in channels), so only do this after calculating the distance
            electrodes = ccep_sortElectrodes(electrodes, channels, 0);
            elecPositions = [electrodes.x electrodes.y electrodes.z];
            
        end
%{
        % loop over the tracts (SLF, AF, etc...) and sub-tracts (frontal, central, parietal, etc...)
        for iTr = 1:length(rois)

            %
            if rois(iTr).interHemi == 1
                error('Inter-hemisphere not implemented, only processing per hemisphere');
            end

            % for each hemisphere
            for iHemi = 1:2
                hemi = 'l';
                if iHemi == 2, hemi = 'r';  end

                % load the subject freesurfer hemisphere pial
                gPial = gifti(fullfile(subjFSDir, ['/pial.', upper(hemi), '.surf.gii']));
                %viewGii(gPial, 'trans.7', elecPositions, 'WireSpheres1');
                
                % 
                trkFile = fullfile(track_path, rois(iTr).tract_name);
                disp('Retrieving tract distance');
                

                %%
                %  Load the tract file and transform the line vertices from MNI-152 space to native space

                rois(iTr).MNIfiles{iHemi} = [trkFile, '_', upper(hemi), '.trk.gz'];
                [fibers, idx] = ea_trk2ftr(rois(iTr).MNIfiles{iHemi}, 1);
                fibers = fibers(:, 1:3);

                % RAS to LPS, ANTs (ITK) uses LPS coords
                fibers(:, 1:2) = -fibers(:, 1:2);

                % Apply transform (input should be Nx3 = <points> x <X,Y,Z>)
                fibers = ea_ants_apply_transforms_to_points(fullfile(myDataPath.input, 'derivatives', 'coreg_ANTs', ccepData(iSubj).id), fibers, 0);

                % LPS to RAS, restore to RAS coords
                fibers(:, 1:2) = -fibers(:, 1:2);

                % get the fiber end-point coordinates
                trcLineEnds = nan(length(idx) * 2, 6);              % the last two points of each end of a line (order: inner x,y,z + outer x,y,z)
                startV = 1;
                for iTrk = 1:length(idx)

                    % determine the start- and end-vertex indices
                    if iTrk > 1,   startV = sum(idx(1:iTrk - 1)) + 1;   end
                    endV = startV + idx(iTrk) - 1;

                    %trcLineEnds((iTrk - 1) * 2 + 1, :)              = [fibers(startV + 3, 1:3), fibers(startV, 1:3)];
                    %trcLineEnds((iTrk - 1) * 2 + 2, :)              = [fibers(endV - 3, 1:3), fibers(endV, 1:3)];

                    % start and extend from a part of the line further from the actual end of the line
                    trcLineEnds((iTrk - 1) * 2 + 1, :)              = [fibers(startV + 10, 1:3), fibers(startV + 4, 1:3)];
                    trcLineEnds((iTrk - 1) * 2 + 2, :)              = [fibers(endV - 10, 1:3), fibers(endV - 4, 1:3)];

                end

                % extend the tract-lines only forward
                %extLines = ((trcLineEnds(:, 4:6) - trcLineEnds(:, 1:3)) * 25) + trcLineEnds(:, 1:3);
                %extLines = [trcLineEnds(:, 1:3), extLines];

                % extend the tract-lines forward and backward
                extLines = [((trcLineEnds(:, 1:3) - trcLineEnds(:, 4:6)) * 3) + trcLineEnds(:, 1:3), ...
                            ((trcLineEnds(:, 4:6) - trcLineEnds(:, 1:3)) * 10) + trcLineEnds(:, 1:3)];

                %{
                % debug, check extended lines
                viewGii([], extLines);
                hold on;
                startV = 1;
                for i = 1:length(idx)    
                    if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
                    endV = startV + idx(i) - 1;
                    plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
                end
                hold off;


                viewGii([]);
                hold on;
                startV = 1;
                for i = 1:length(idx)    
                    if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
                    endV = startV + idx(i) - 1;
                    plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
                end
                hold off;
                %}


                
                %
                % 
                %

                
                % loop over the sub-tracts (frontal, central, parietal, etc...)
                for iSubTr = 1:length(rois(iTr).sub_tract)

                    % retrieve the distance between the stimulation and response end-point ROIs
                    % for this particular patient given specific tracts
                    [trkDist, trkLineIndices, trkExtElecs, gROIPial1, gROIPial2] = ccep_retrieveInterROIDistance( ...
                                                                                            gPial, hemi, ...
                                                                                            fibers, idx, ...
                                                                                            subjFSDir, ...
                                                                                            rois(iTr).sub_tract(iSubTr).roi1, ...
                                                                                            rois(iTr).sub_tract(iSubTr).roi2, ...
                                                                                            rois(iTr).sub_tract(iSubTr).allowIntraROI, ...
                                                                                            elecPositions, ...
                                                                                            trcLineEnds, ...
                                                                                            extLines);

                    % warn if the tract length could not be determined
                    if isempty(trkDist) || isnan(trkDist)
                        warning(['Could not determine tract-length: subject ', ccepData(iSubj).id, ' - hemisphere ', upper(hemi), ' - tract-roi', rois(iTr).tract_name, ' ', rois(iTr).sub_tract(iSubTr).name])
                    end
                                                                                        
                    % store the distances, included MNI files and tract-lines, and elec
                    rois(iTr).sub_tract(iSubTr).MNIlineIndices{iHemi} = trkLineIndices;
                    rois(iTr).sub_tract(iSubTr).nativeDistances{iHemi} = trkDist;
                    if ~isempty(trkExtElecs)
                        rois(iTr).sub_tract(iSubTr).extElecNames{iHemi} = electrodes(trkExtElecs, :).name;
                    else
                        rois(iTr).sub_tract(iSubTr).extElecNames{iHemi} = {};
                    end


                    %{
                    % debug, show all tracts in native
                    viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge');
                    hold on;
                    startV = 1;
                    for i = 1:length(idx)    
                        if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
                        endV = startV + idx(i) - 1;
                        plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
                    end
                    hold off;
                    %}

                    %{
                    % debug, show ROI tracts in native with relevant electrodes
                    if (iHemi == 1 && any(contains(ccepData(iSubj).electrodes.jsonHemi, 'L'))) || (iHemi == 2 && any(contains(ccepData(iSubj).electrodes.jsonHemi, 'R')))   % only on hemisphere that matters

                        excludedTrkElecs = 1:size(elecPositions, 1);
                        excludedTrkElecs(trkExtElecs) = [];
                        b = (trkLineIndices - 1) * 2 + 1;
                        %viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(trkExtElecs, :), 'WireSpheres3');
                        %viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(trkExtElecs, :), elecPositions(excludedTrkElecs, :), 'WireSpheres1');
                        %viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(trkExtElecs, :), 'WireSpheres1');
                        viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge', elecPositions(trkExtElecs, :), 'WireSpheres1');
                        set(gcf, 'Visible', 'off');
                        hold on;
                        startV = 1;
                        for iLine = trkLineIndices
                            if iLine > 1,   startV = sum(idx(1:iLine - 1)) + 1;   end
                            endV = startV + idx(iLine) - 1;
                            plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
                        end
                        hold off;
                        if iHemi == 1
                            campos([-1452.248, 40.763, 153.4652]);
                            camtarget([-59.8763, -11.4336, 25.2902]);
                            camup([0.092386, 0.019813, 0.99553]);
                            camva(7.4861);
                        else
                            campos([1284.2761, 19.1907, 389.4063]);
                            camtarget([55.8074, 1.6975, 30.1294]);
                            camup([-0.28023, -0.027494, 0.95954]);
                            camva(6.9761);
                        end
                        delete(findall(gcf, 'Type', 'light'));
                        camlight(gca, 'headlight');

                        %
                        if ~exist(fullfile(myDataPath.output, 'derivatives', 'render', 'tractsROIs'), 'dir')
                            mkdir(fullfile(myDataPath.output, 'derivatives', 'render', 'tractsROIs'));
                        end
                        figureName = fullfile(myDataPath.output, 'derivatives', 'render', 'tractsROIs', [ccepData(iSubj).id, '_', upper(hemi), '_', rois(iTr).tract_name, '_',  strrep(rois(iTr).sub_tract(iSubTr).name, '-', '_'), '.png']);
                        set(gcf,'PaperPositionMode', 'auto');
                        set(gcf, 'Visible', 'on');
                        print('-dpng', '-r300', figureName);
                        close(gcf)

                    end
                    %}

                    %{
                    % debug, show excluded tracts in native
                    excludedTrkLines = 1:length(idx);
                    excludedTrkLines(trkLineIndices) = [];
                    b = (excludedTrkLines - 1) * 2 + 1;
                    viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge');
                    hold on;
                    startV = 1;
                    for iLine = excludedTrkLines
                        if iLine > 1,   startV = sum(idx(1:iLine - 1)) + 1;   end
                        endV = startV + idx(iLine) - 1;
                        plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
                    end
                    hold off;
                    %}

                    %{
                    % debug, show electrodes detection regions (included and excluded
                    elecRadius = 1;
                    excludedTrkElecs = 1:size(elecPositions, 1);
                    excludedTrkElecs(trkExtElecs) = [];
                    viewGii(gROIPial1, gROIPial2, 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(trkExtElecs, :), ['WireSpheres', num2str(elecRadius)]);
                    viewGii(gROIPial1, gROIPial2, 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(excludedTrkElecs, :), ['WireSpheres', num2str(elecRadius)]);
                    %}


                    % store the native tract/end-point ROIs track-lines in output struct
                    ccepData(iSubj).rois = rois;

                end     % end of sub-tract loop

            end     % end hemisphere loop
            
        end     % end of tract loop
        %}
    end     % end subjects loop
    
end


%%
%  (optionally) store the updated structure

s = input('Do you want to save the ccepData structure? [y/n]: ', 's');
if strcmp(s, 'y')
    save(fullfile(myDataPath.output, 'derivatives', 'av_ccep', 'ccepData_V2.mat'), 'ccepData')
end
