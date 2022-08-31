%
%   Retrieve the distance in subject native space (mm) between end-point ROIs given specific tracts
%   [distance] = ccep_retrieveInterROIDistance(interHemi, trkFile, subjectANTsFolder, subjectFsFolder, roi1, roi2)
%
%       interHemi         = indicate if both hemispheres are processed together with one tract
%                           file and combined FS hemisphere surfaces (1), or seperately with different
%                           tract files and FS surfaces for left and right (0).
%       trkFile           = the tract file(s) to load the MNI-152 tracts from. If 'interHemi' is set
%                           to 0, tract files will be loaded seperately for the left and the right
%                           hemisphere and either '_L.trk.gz' '_R.trk.gz' will be appended (e.g. 
%                           when '/tracks/FA' is passed, both '/tracks/FA_L.trk.gz' and 
%                           '/tracks/FA_R.trk.gz' will be processed seperately with their respective
%                           hemispheres). If 'interHemi' is set to 1, only '.trk.gz' will be appended.
%       subjectANTsFolder = the subject specific folder with the ANTs transformation files, this is 
%                           used to transform the tracts from MNI-152 space to subject native space
%       subjectFsFolder   = the subject-specific freesurfer folder, ...
%       roi1              = ...
%       roi2              = ...
%       allowIntraROI     = Allow tracts to begin and end in the same ROI
%       elecPositions     = ...
%       elecHemi          = Indicates for the electrodes on which hemisphere they are
%
%   Returns: 
%       trkDistance       = the distance in mm between the ROI's end-points
%       trkFiles          = the tract files that were used to calculate the distance
%       trkIndices        = the indices of the tract-lines that were used to calculate the distance
%       trkNativeFibers   = the positions of the fiber vertices in native space
%
%
%   Note: For Linux and Mac, execute permissions might need to be set for
%         the leadDBS binary files ('/leadDBS/ext_libs/ANTs').
%         Navigate to that folder in the terminal and run 'chmod u+x antsApplyTransformsToPoints.*'
%
%
%
%   Max van den Boom, Multimodal Neuroimaging Lab (MNL), Mayo Clinic, 2022
%
function [trkDistance, trkFiles, trkIndices, trkNativeFibers, trkExtElecs] = ccep_retrieveInterROIDistance(interHemi, trkFile, subjectANTsFolder, subjectFsFolder, roi1, roi2, allowIntraROI, elecPositions, elecHemi)
    
    % check leadDBS availability and setup
    if exist('ea_getearoot') ~= 2 || exist('ea_prefs') ~= 2
       error('Could not find LeadDBS functions. LeadDBS is an external dependency, make sure it installed and added as a MATLAB path.');
    end
    addpath(genpath(ea_getearoot()));
    prefs = ea_prefs;
    if exist('ea_trk2ftr') ~= 2 || exist('ea_ants_apply_transforms_to_points') ~= 2
        error('Could not find LeadDBS functions. LeadDBS is an external dependency, make sure it installed and added as a MATLAB path.');
    end
    
    
    
    if interHemi == 1
        % process both hemispheres together
        
    else
        % process each hemisphere seperately

        % loop over the tract files (likely to be two, one for each hemisphere)
        for iHemi = 1:2
            hemi = 'l';
            if iHemi == 2, hemi = 'r';  end

            
            %%
            %  Load the tract file and transform the line vertices from MNI-152 space to native space

            trkFiles{iHemi} = [trkFile, '_', upper(hemi), '.trk.gz'];
            [fibers, idx] = ea_trk2ftr(trkFiles{iHemi}, 1);
            fibers = fibers(:, 1:3);

            % RAS to LPS, ANTs (ITK) uses LPS coords
            fibers(:, 1:2) = -fibers(:, 1:2);

            % Apply transform (input should be Nx3 = <points> x <X,Y,Z>)
            tic
            fibers = ea_ants_apply_transforms_to_points(subjectANTsFolder, fibers, 0);
            toc
            
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
                        ((trcLineEnds(:, 4:6) - trcLineEnds(:, 1:3)) * 14) + trcLineEnds(:, 1:3)];
            
            %{
            % debug, check extended lines
            gPial = gifti(fullfile(subjectFsFolder, ['/pial.', upper(hemi), '.surf.gii']));
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
                        
            
            %%
            %  Determine which tract-lines end in the ROIs

            % debug
            gPial = gifti(fullfile(subjectFsFolder, ['/pial.', upper(hemi), '.surf.gii']));

            % read the annotations for the pial brain (a2009s = Destrieux) and relabel 
            % the vertex annotation (color) labels to the ROI area indices
            [~, annotVertexLabels, annotColortable] = read_annotation(fullfile(subjectFsFolder, 'label', [hemi, 'h.aparc.a2009s.annot']));
            
            % for each ROI, retrieve the tract-lines with a forward search
            [proxTrkLines1, gROIPial1] = retrieveROILines(roi1, annotColortable, annotVertexLabels, trcLineEnds, extLines, gPial);
            [proxTrkLines2, gROIPial2] = retrieveROILines(roi2, annotColortable, annotVertexLabels, trcLineEnds, extLines, gPial);
            
            
            if allowIntraROI == 1

                % determine which tract-lines that touch upon both ROIs (with each end in one ROI) or that touch upon the same ROI (with both ends)
                roisTrkLines = find((proxTrkLines1(:, 1) & proxTrkLines2(:, 2)) | (proxTrkLines1(:, 2) & proxTrkLines2(:, 1)) | ...
                                     all(proxTrkLines1, 2) | all(proxTrkLines2, 2))';
                
            else
                
                % determine which tract-lines that touch upon both ROIs (with each end in one ROI)
                roisTrkLines = find((proxTrkLines1(:, 1) & proxTrkLines2(:, 2)) | (proxTrkLines1(:, 2) & proxTrkLines2(:, 1)))';
            
            end
            
            
            
            % calculate the length of each tract-line
            for iTrk = 1:length(roisTrkLines)

                % determine the start- and end-vertex indices
                startV = 1;
                if roisTrkLines(iTrk) > 1,   startV = sum(idx(1:roisTrkLines(iTrk) - 1)) + 1;   end
                endV = startV + idx(roisTrkLines(iTrk)) - 1;

                % calculate and store the length of the current tract line
                trkLengths(iTrk) = sqrt(sum(sum(diff(fibers(startV:endV, 1:3), 1, 1) .^ 2, 2)));

            end

            % return the tract-lines that run between the ROIs
            trkIndices{iHemi} = roisTrkLines;
            
            
            
            %%
            %  Determine which electrodes are at the end of the tract-lines

            % specify the point of detection along the extended lines of the line-tracts that pierce the ROIs
            steps = 20;
            extPoints = nan(length(roisTrkLines), 2, steps, 3);
            for iLine = 1:length(roisTrkLines)
                lineIdx = (roisTrkLines(iLine) - 1) * 2 + 1;

                extPoints(iLine, 1, :, :)  = [linspace(extLines(lineIdx, 1),     extLines(lineIdx, 4), steps); ...
                                              linspace(extLines(lineIdx, 2),     extLines(lineIdx, 5), steps); ...
                                              linspace(extLines(lineIdx, 3),     extLines(lineIdx, 6), steps)]';
                extPoints(iLine, 2, :, :)  = [linspace(extLines(lineIdx + 1, 4), extLines(lineIdx + 1, 1), steps); ...
                                              linspace(extLines(lineIdx + 1, 5), extLines(lineIdx + 1, 2), steps); ...
                                              linspace(extLines(lineIdx + 1, 6), extLines(lineIdx + 1, 3), steps)]';
            end
            extPoints = reshape(extPoints, [], 3);

            %{
            % debug, check detection points and lines
            b = (roisTrkLines - 1) * 2 + 1;
            viewGii(gROIPial1, gROIPial2, 'merge', extLines(b, :), extLines(b + 1, :), extPoints);
            %}

            % check which electrodes are at the end of the extended tract-lines
            
            %pInter = collPoints(elecPositions, extPoints, 2);
             
            
            elecRadius = 1;
            dist = (elecPositions(:, 1)' - extPoints(:, 1)) .^ 2 + ...
                   (elecPositions(:, 2)' - extPoints(:, 2)) .^ 2 + ...
                   (elecPositions(:, 3)' - extPoints(:, 3)) .^ 2;
            dist = dist < (elecRadius  .^ 2);
            roisTrkElecs = find(any(dist, 1));

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

            
            % debug, show ROI tracts in native with relevant electrodes
            if (iHemi == 1 && any(contains(elecHemi, 'L'))) || (iHemi == 2 && any(contains(elecHemi, 'R')))   % only on hemisphere that matters

                excludedTrkElecs = 1:size(elecPositions, 1);
                excludedTrkElecs(roisTrkElecs) = [];
                b = (roisTrkLines - 1) * 2 + 1;
                %viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(roisTrkElecs, :), 'WireSpheres3');
                viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(roisTrkElecs, :), elecPositions(excludedTrkElecs, :), 'WireSpheres3');
                %viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge');
                hold on;
                startV = 1;
                for iLine = roisTrkLines
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
                [~, trc] = fileparts(trkFile);
                [~, sub] = fileparts(subjectFsFolder);
                myDataPath = setLocalDataPath(1);
                if ~exist(fullfile(myDataPath.output, 'derivatives', 'render', 'tractsROIs'), 'dir')
                    mkdir(fullfile(myDataPath.output, 'derivatives', 'render', 'tractsROIs'));
                end
                figureName = fullfile(myDataPath.output, 'derivatives', 'render', 'tractsROIs', [sub, '_', upper(hemi), '_', trc, '_',  strrep(num2str(roi1), ' ', ''), '_', strrep(num2str(roi2), ' ', ''), '.png']);
                set(gcf,'PaperPositionMode', 'auto')
                print('-dpng', '-r300', figureName);
                close(gcf)
                
            end
            
            
            %{
            % debug, show excluded tracts in native
            excludedTrkLines = 1:length(idx);
            excludedTrkLines(roisTrkLines) = [];
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
            excludedTrkElecs = 1:size(elecPositions, 1);
            excludedTrkElecs(roisTrkElecs) = [];
            viewGii(gROIPial1, gROIPial2, 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(roisTrkElecs, :), ['WireSpheres', num2str(elecRadius)]);
            viewGii(gROIPial1, gROIPial2, 'merge', extLines(b, :), extLines(b + 1, :), elecPositions(excludedTrkElecs, :), ['WireSpheres', num2str(elecRadius)]);
            %}
            
            % return the electrodes that the extended tract-lines pierces
            trkExtElecs{iHemi} = roisTrkElecs;
            


            %%
            %  Calculate the average distance over the tract-lines

            % loop over the tract-lines
            trkLengths = nan(1, length(trkIndices{iHemi}));
            for iTrk = 1:length(trkIndices{iHemi})

                % determine the start- and end-vertex indices
                startV = 1;
                if trkIndices{iHemi}(iTrk) > 1,   startV = sum(idx(1:trkIndices{iHemi}(iTrk) - 1)) + 1;   end
                endV = startV + idx(trkIndices{iHemi}(iTrk)) - 1;

                % calculate and store the length of the current tract line
                trkLengths(iTrk) = sum(sqrt(sum(diff(fibers(startV:endV, 1:3), 1, 1) .^ 2, 2)));

            end
            
            % store(/return) the average over all tract-lines
            trkDistance{iHemi} = mean(trkLengths);
            
            if nargout > 3
                trkNativeFibers{iHemi} = {fibers, idx};
            end

        end
    end
    
end

function [proxTrkLines, gROIPial] = retrieveROILines(roiCodes, annotColortable, annotVertexLabels, trcLineEnds, trcExtLines, gPial)

    % convert the Destrieux codes to Destrieux labels
    dstrxCodeToLabel = {'G_and_S_frontomargin';'G_and_S_occipital_inf';'G_and_S_paracentral';'G_and_S_subcentral';'G_and_S_transv_frontopol';'G_and_S_cingul-Ant';'G_and_S_cingul-Mid-Ant';'G_and_S_cingul-Mid-Post';'G_cingul-Post-dorsal';'G_cingul-Post-ventral';'G_cuneus';'G_front_inf-Opercular';'G_front_inf-Orbital';'G_front_inf-Triangul';'G_front_middle';'G_front_sup';'G_Ins_lg_and_S_cent_ins';'G_insular_short';'G_occipital_middle';'G_occipital_sup';'G_oc-temp_lat-fusifor';'G_oc-temp_med-Lingual';'G_oc-temp_med-Parahip';'G_orbital';'G_pariet_inf-Angular';'G_pariet_inf-Supramar';'G_parietal_sup';'G_postcentral';'G_precentral';'G_precuneus';'G_rectus';'G_subcallosal';'G_temp_sup-G_T_transv';'G_temp_sup-Lateral';'G_temp_sup-Plan_polar';'G_temp_sup-Plan_tempo';'G_temporal_inf';'G_temporal_middle';'Lat_Fis-ant-Horizont';'Lat_Fis-ant-Vertical';'Lat_Fis-post';'Pole_occipital';'Pole_temporal';'S_calcarine';'S_central';'S_cingul-Marginalis';'S_circular_insula_ant';'S_circular_insula_inf';'S_circular_insula_sup';'S_collat_transv_ant';'S_collat_transv_post';'S_front_inf';'S_front_middle';'S_front_sup';'S_interm_prim-Jensen';'S_intrapariet_and_P_trans';'S_oc_middle_and_Lunatus';'S_oc_sup_and_transversal';'S_occipital_ant';'S_oc-temp_lat';'S_oc-temp_med_and_Lingual';'S_orbital_lateral';'S_orbital_med-olfact';'S_orbital-H_Shaped';'S_parieto_occipital';'S_pericallosal';'S_postcentral';'S_precentral-inf-part';'S_precentral-sup-part';'S_suborbital';'S_subparietal';'S_temporal_inf';'S_temporal_sup';'S_temporal_transverse'};
    roi_lbls = {};
    for roiCode = roiCodes,     roi_lbls{end + 1} = dstrxCodeToLabel{roiCode};      end

    % label the vertices according to the freesurfer labels
    roiVertexLabels = mx.freesurfer.fsRelabelToAreas(roi_lbls, annotColortable, annotVertexLabels);

    % extract the vertices/faces that belong to the ROI
    roiVertexIDs = find(~isnan(roiVertexLabels));
    roiFacesLabels = ismember(gPial.faces, roiVertexIDs);
    roiFacesLabels = all(roiFacesLabels, 2);
    [roiVertices, roiFaces, ~] = mx.three_dimensional.extract3DFaces(gPial, roiFacesLabels == 1);
    clear roiVertexIDs roiFacesLabels;
    
    % build output gifti
    gROIPial = gPial;
    gROIPial.vertices = roiVertices;
    gROIPial.faces = roiFaces;
    %viewGii(gROIPial, trcLineEnds(:, 1:3), trcLineEnds(:, 4:6), trcExtLines);
    
    % build points along the extended lines
    % TODO: probably could optimize this
    steps = 30;
    extPoints = nan(size(trcLineEnds, 1), steps, 3);
    for iLine = 1:2:size(trcLineEnds, 1)
        
        extPoints(iLine, :, :)      = [linspace(trcExtLines(iLine, 1),     trcExtLines(iLine, 4), steps); ...
                                       linspace(trcExtLines(iLine, 2),     trcExtLines(iLine, 5), steps); ...
                                       linspace(trcExtLines(iLine, 3),     trcExtLines(iLine, 6), steps)]';
        extPoints(iLine + 1, :, :)  = [linspace(trcExtLines(iLine + 1, 4), trcExtLines(iLine + 1, 1), steps); ...
                                       linspace(trcExtLines(iLine + 1, 5), trcExtLines(iLine + 1, 2), steps); ...
                                       linspace(trcExtLines(iLine + 1, 6), trcExtLines(iLine + 1, 3), steps)]';
        
    end
    extPoints = reshape(extPoints, [], 3);
    %viewGii(gROIPial, trcExtLines, reshape(extPoints, [], 3));

    % Check extPoints to roiVertices, this will allow for a "dilated" line piercing
    % Note: some ROIs the sulci might not be included so an extended tract-line can pierce just between the gyri
    pInter = collPoints(extPoints, roiVertices, 2);
    
    % return which tract-lines touch ROI <1st dim = tract-line, 2nd dim = which end of the tract line>
    proxTrkLines = any(reshape(pInter, [], steps, 1), 2);
    proxTrkLines = [proxTrkLines(1:2:end), proxTrkLines(2:2:end)];
    
end

function pInter = collPoints(extPoints, roiVertices, radius)

    % 
    pInter = false(1, size(extPoints, 1));

    % split into part of 70000 points (requires about 3GB of memory)
    maxPointPerLoop = 70000;
    numLoops = ceil(size(extPoints, 1) / maxPointPerLoop);

    % set the point ranges for each set
    iCounter = 1;
    for i = 1:numLoops
        range = [];
        if iCounter + maxPointPerLoop <= size(extPoints, 1)
            range = [iCounter, iCounter + maxPointPerLoop - 1];
        else
            range = [iCounter, size(extPoints, 1)];
        end
        
        % determine which points intersect with the ROI vertices
        dist = (extPoints(range(1):range(2), 1)' - roiVertices(:, 1)) .^ 2 + ...
               (extPoints(range(1):range(2), 2)' - roiVertices(:, 2)) .^ 2 + ...
               (extPoints(range(1):range(2), 3)' - roiVertices(:, 3)) .^ 2;
        dist = dist < (radius  .^ 2);
        pInter(range(1):range(2)) = any(dist, 1);
        clear dist;
        
        %
        iCounter = iCounter + maxPointPerLoop;
    end
    

end

