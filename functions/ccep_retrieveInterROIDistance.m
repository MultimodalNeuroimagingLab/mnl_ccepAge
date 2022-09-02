%
%   Retrieve the distance in subject native space (mm) between end-point ROIs given specific tracts
%   [distance] = ccep_retrieveInterROIDistance(gPial, hemi, fibers, idx, subjectFsFolder, roi1, roi2, allowIntraROI, elecPositions, trcLineEnds, extLines)
%
%       gPial             = ...
%       hemi              = which hemisphere the tract apply to. Used to determine whether to load the FS files for the left or right hemi
%       fibers            = ...
%       idx               = ...
%       subjectFsFolder   = the subject-specific freesurfer folder, ...
%       roi1              = ...
%       roi2              = ...
%       allowIntraROI     = Allow tracts to begin and end in the same ROI
%       elecPositions     = ...
%       trcLineEnds       = 
%       extLines          = 
%
%   Returns: 
%       trkDistance       = the distance in mm between the ROI's end-points
%       trkIndices        = the indices of the tract-lines that were used to calculate the distance
%       trkExtElecs       = the indices of the electrodes (as provided by the 'elecPositions' input argument) that the extended
%                           tract-lines intersect with (given a specific electrode radius)
%       gROIPial1         = [optional]
%       gROIPial2         = [optional]
%
%
%   Max van den Boom, Multimodal Neuroimaging Lab (MNL), Mayo Clinic, 2022
%
function [trkDistance, trkIndices, trkExtElecs, gROIPial1, gROIPial2] = ccep_retrieveInterROIDistance(gPial, hemi, fibers, idx, subjectFsFolder, roi1, roi2, allowIntraROI, elecPositions, trcLineEnds, extLines)
    

    %%
    %  Determine which tract-lines end in the ROIs

    % read the annotations for the pial brain (a2009s = Destrieux) and relabel 
    % the vertex annotation (color) labels to the ROI area indices
    [~, annotVertexLabels, annotColortable] = read_annotation(fullfile(subjectFsFolder, 'label', [hemi, 'h.aparc.a2009s.annot']));

    % for each ROI, retrieve the tract-lines with a forward search
    [proxTrkLines1, gROIPial1] = retrieveROILines(roi1, annotColortable, annotVertexLabels, trcLineEnds, extLines, gPial);
    [proxTrkLines2, gROIPial2] = retrieveROILines(roi2, annotColortable, annotVertexLabels, trcLineEnds, extLines, gPial);
    
    % determine and return the tract-lines that run between the ROIs
    if allowIntraROI == 1

        % determine which tract-lines that touch upon both ROIs (with each end in one ROI) or that touch upon the same ROI (with both ends)
        trkIndices = find((proxTrkLines1(:, 1) & proxTrkLines2(:, 2)) | (proxTrkLines1(:, 2) & proxTrkLines2(:, 1)) | ...
                             all(proxTrkLines1, 2) | all(proxTrkLines2, 2))';

    else

        % determine which tract-lines that touch upon both ROIs (with each end in one ROI)
        trkIndices = find((proxTrkLines1(:, 1) & proxTrkLines2(:, 2)) | (proxTrkLines1(:, 2) & proxTrkLines2(:, 1)))';

    end
    
    % make sure there are at least 5 tract-lines
    if length(trkIndices) < 5
        trkDistance = nan;
        trkIndices = [];
        trkExtElecs = [];
        return;
    end
                    
    % calculate the length of each tract-line
    for iTrk = 1:length(trkIndices)

        % determine the start- and end-vertex indices
        startV = 1;
        if trkIndices(iTrk) > 1,   startV = sum(idx(1:trkIndices(iTrk) - 1)) + 1;   end
        endV = startV + idx(trkIndices(iTrk)) - 1;

        % calculate and store the length of the current tract line
        trkLengths(iTrk) = sqrt(sum(sum(diff(fibers(startV:endV, 1:3), 1, 1) .^ 2, 2)));

    end


    %%
    %  Determine which electrodes are at the end of the tract-lines

    % specify the point of detection along the extended lines of the line-tracts that pierce the ROIs
    steps = 20;
    extPoints = nan(length(trkIndices), 2, steps, 3);
    for iLine = 1:length(trkIndices)
        lineIdx = (trkIndices(iLine) - 1) * 2 + 1;

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
    b = (trkIndices - 1) * 2 + 1;
    viewGii(gROIPial1, gROIPial2, 'merge', extLines(b, :), extLines(b + 1, :), extPoints);
    %}

    % check and return which electrodes are at the end of the extended tract-lines
    pInter = collPoints(elecPositions, extPoints, 1);
    trkExtElecs = find(pInter);


    %%
    %  Calculate the average distance over the tract-lines

    % loop over the tract-lines
    trkLengths = nan(1, length(trkIndices));
    for iTrk = 1:length(trkIndices)

        % determine the start- and end-vertex indices
        startV = 1;
        if trkIndices(iTrk) > 1,   startV = sum(idx(1:trkIndices(iTrk) - 1)) + 1;   end
        endV = startV + idx(trkIndices(iTrk)) - 1;

        % calculate and store the length of the current tract line
        trkLengths(iTrk) = sum(sqrt(sum(diff(fibers(startV:endV, 1:3), 1, 1) .^ 2, 2)));

    end

    % store(/return) the average over all tract-lines
    trkDistance = mean(trkLengths);
    
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
    steps = 40;
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
    % Note: for some some ROIs the sulci might not be included and can create small gaps between the gyri, so an extended tract-line
    %       can pierce just trought and not be include, therefore allow a bit of a radius around the line
    pInter = collPoints(extPoints, roiVertices, .5);
    
    % return which tract-lines touch ROI <1st dim = tract-line, 2nd dim = which end of the tract line>
    proxTrkLines = any(reshape(pInter, [], steps, 1), 2);
    proxTrkLines = [proxTrkLines(1:2:end), proxTrkLines(2:2:end)];
    
end

%
%  Find the 
%
%   checkPoints = The points around which (given the radius) to search for 'passivePoints'
%
function pInter = collPoints(checkPoints, passivePoints, radius)

    % 
    pInter = false(1, size(checkPoints, 1));

    % split into part of 70000 points (requires about 3GB of memory)
    maxPointPerLoop = 70000;
    numLoops = ceil(size(checkPoints, 1) / maxPointPerLoop);

    % set the point ranges for each set
    iCounter = 1;
    for i = 1:numLoops
        range = [];
        if iCounter + maxPointPerLoop <= size(checkPoints, 1)
            range = [iCounter, iCounter + maxPointPerLoop - 1];
        else
            range = [iCounter, size(checkPoints, 1)];
        end
        
        % determine which points intersect with the ROI vertices
        dist = (checkPoints(range(1):range(2), 1)' - passivePoints(:, 1)) .^ 2 + ...
               (checkPoints(range(1):range(2), 2)' - passivePoints(:, 2)) .^ 2 + ...
               (checkPoints(range(1):range(2), 3)' - passivePoints(:, 3)) .^ 2;
        dist = dist < (radius  .^ 2);
        pInter(range(1):range(2)) = any(dist, 1);
        clear dist;
        
        %
        iCounter = iCounter + maxPointPerLoop;
    end 

end

