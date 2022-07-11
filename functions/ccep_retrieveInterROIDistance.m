%
%   Retrieve the distance in subject native space (mm) between end-point ROIs given specific tracts
%   [distance] = ccep_retrieveInterROIDistance(interHemi, trkFile, subjectANTsFolder, subjectFsFolder, roi1, roi2)
%
%       interHemi         = indicate if both hemispheres are processed together with one tract
%                           file and combined FS hemisphere surfaces (1), or seperately with different
%                           tract files and FS surfaces for left and right (0).
%       trkFile           = the tract file(s) to load the MNI tracts from. If 'interHemi' is set
%                           to 0, tract files will be loaded seperately for the left and the right
%                           hemisphere and either '_L.trk.gz' '_R.trk.gz' will be appended (e.g. 
%                           when '/tracks/FA' is passed, both '/tracks/FA_L.trk.gz' and 
%                           '/tracks/FA_R.trk.gz' will be processed seperately with their respective
%                           hemispheres). If 'interHemi' is set to 1, only '.trk.gz' will be appended.
%       subjectANTsFolder = the subject specific folder with the ANTs transformation files, this is 
%                           used to transform the tracts from MNI space to subject native space
%       subjectFsFolder   = the subject-specific freesurfer folder, ...
%
%   Returns: 
%       trkDistance       = the distance in mm between the ROI's end-points
%       trkFiles          = the tract files that were used to calculate the distance
%       trkIndices        = the indices of the tract-lines that were used to calculate the distance
%
%
%   Note: For Linux and Mac, execute permissions might need to be set for
%         the leadDBS binary files ('/leadDBS/ext_libs/ANTs').
%         Navigate to that folder in the terminal and run 'chmod u+x antsApplyTransformsToPoints.*'
%
%
%
%   Max van den Boom, Multimodal Neuroimaging Lab, Mayo Clinic, 2022
%
function [trkDistance, trkFiles, trkIndices] = ccep_retrieveInterROIDistance(interHemi, trkFile, subjectANTsFolder, subjectFsFolder, roi1, roi2)
    
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
            %  Load the tract file and transform the line vertices from MNI space to native space

            % load
            trkFiles{iHemi} = [trkFile, '_', upper(hemi), '.trk.gz'];
            [fibers, idx] = ea_trk2ftr(trkFiles{iHemi}, 1);
            fibers = fibers(:, 1:3);

            % RAS to LPS, ANTs (ITK) uses LPS coords
            fibers(:, 1:2) = -fibers(:, 1:2);

            % Apply transform (input should be Nx3 = <points> x <X,Y,Z>)
            fibers = ea_ants_apply_transforms_to_points(subjectANTsFolder, fibers, 0);

            % LPS to RAS, restore to RAS coords
            fibers(:, 1:2) = -fibers(:, 1:2);

            % get the fiber end-point coordinates
            trcLineEndPoints = nan(length(idx) * 2, 3);
            trcLineEndPoints2 = nan(length(idx) * 2, 1);
            startV = 1;
            for iTrk = 1:length(idx)

                % determine the start- and end-vertex indices
                if iTrk > 1,   startV = sum(idx(1:iTrk - 1)) + 1;   end
                endV = startV + idx(iTrk) - 1;

                trcLineEndPoints2((iTrk - 1) * 2 + 1, :)   = startV;
                trcLineEndPoints2((iTrk - 1) * 2 + 2, :)   = endV;
                trcLineEndPoints((iTrk - 1) * 2 + 1, :)   = fibers(startV, 1:3);
                trcLineEndPoints((iTrk - 1) * 2 + 2, :)   = fibers(endV, 1:3);

            end
            
            b = reshape(trcLineEndPoints2, 2, []);
            

            
            %%
            %  Determine which tract-lines that are close enough to the ROIs

            % debug
            %gPial = gifti(fullfile(subjectFsFolder, ['/surf/', hemi, 'h.pial.gii']));
            gPial = gifti(fullfile(subjectFsFolder, ['/pial.', upper(hemi), '.surf.gii']));

            % read the annotations for the pial brain (a2009s = Destrieux) and relabel 
            % the vertex annotation (color) labels to the ROI area indices
            [~, annotVertexLabels, annotColortable] = read_annotation(fullfile(subjectFsFolder, 'label', [hemi, 'h.aparc.a2009s.annot']));
            
            
            radius = 5; % in mm
            
            % for each ROI, retrieve the tract-lines of which an end-point is within x radius of any of the ROI's vertices
            [proxTrkLines1, gROIPial1] = retrieveROI(roi1, annotColortable, annotVertexLabels, radius, trcLineEndPoints, gPial);
            [proxTrkLines2, gROIPial2] = retrieveROI(roi2, annotColortable, annotVertexLabels, radius, trcLineEndPoints, gPial);
            
            % determine whether there are tract-lines that touch upon both ROIs
            roisTrkLines = find(all([proxTrkLines1; proxTrkLines2], 1));

            % debug, show ROI pials
            %viewGii(gPial, gROIPial1, gROIPial2);

            % debug, show all tracts in native
            %viewGii(gPial, 'trans.1');
            viewGii(gROIPial1, gROIPial2, 'trans.7', 'merge');
            hold on;
            startV = 1;
            for i = 1:length(idx)    
                if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
                endV = startV + idx(i) - 1;

                plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));
            end
            plot3(trcLineEndPoints(:, 1), trcLineEndPoints(:, 2), trcLineEndPoints(:, 3), 'ob');
            hold off;

            % debug, show ROI tracts in native
            %viewGii(gPial, 'trans.1');
            viewGii(gROIPial1, gROIPial2, 'trans.8', 'merge');
            %g = b(:, roisTrkLines);
            %viewGii(gROIPial1, gROIPial2, 'trans.1', 'merge', fibers(g(:), :), 'WireSpheres5');
            hold on;
            for i = roisTrkLines
            
                startV = 1;
                if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
                endV = startV + idx(i) - 1;

                plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));

            end
            %plot3(trcLineEndPoints(:, 1), trcLineEndPoints(:, 2), trcLineEndPoints(:, 3), 'ob');
            hold off;
            

            % debug, show excluded tracts in native
            
            excludedTrkLines = 1:length(idx);
            excludedTrkLines(roisTrkLines) = [];
            %g = b(:, excludedTrkLines);
            %viewGii(gROIPial1, gROIPial2, 'trans.8', 'merge', fibers(g(:), :), ['WireSpheres', radius]);
            viewGii(gROIPial1, gROIPial2, 'trans.8', 'merge');
            hold on;
            for i = excludedTrkLines
            
                startV = 1;
                if i > 1,   startV = sum(idx(1:i - 1)) + 1;   end
                endV = startV + idx(i) - 1;

                plot3(fibers(startV:endV, 1), fibers(startV:endV, 2), fibers(startV:endV, 3));

            end
            %plot3(trcLineEndPoints(:, 1), trcLineEndPoints(:, 2), trcLineEndPoints(:, 3), 'ob');
            hold off;

            % return the tract-lines between the ROIs
            trkIndices{iHemi} = roisTrkLines;
            
            % debug, pick all
            %trkIndices{iSet} = 1:length(idx);


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
                trkLengths(iTrk) = sqrt(sum(sum(diff(fibers(startV:endV, 1:3), 1, 1) .^ 2, 2)));

            end

            % store(/return) the average over all tract-lines
            trkDistance{iHemi} = mean(trkLengths);

        end
    end
    
end

function [proxTrkLines, gROIPial] = retrieveROI(roiCodes, annotColortable, annotVertexLabels, radius, trcLineEndPoints, gPial)

    % convert the Destrieux codes to Destrieux labels
    dstrxCodeToLabel = {'G_and_S_frontomargin';'G_and_S_occipital_inf';'G_and_S_paracentral';'G_and_S_subcentral';'G_and_S_transv_frontopol';'G_and_S_cingul-Ant';'G_and_S_cingul-Mid-Ant';'G_and_S_cingul-Mid-Post';'G_cingul-Post-dorsal';'G_cingul-Post-ventral';'G_cuneus';'G_front_inf-Opercular';'G_front_inf-Orbital';'G_front_inf-Triangul';'G_front_middle';'G_front_sup';'G_Ins_lg_and_S_cent_ins';'G_insular_short';'G_occipital_middle';'G_occipital_sup';'G_oc-temp_lat-fusifor';'G_oc-temp_med-Lingual';'G_oc-temp_med-Parahip';'G_orbital';'G_pariet_inf-Angular';'G_pariet_inf-Supramar';'G_parietal_sup';'G_postcentral';'G_precentral';'G_precuneus';'G_rectus';'G_subcallosal';'G_temp_sup-G_T_transv';'G_temp_sup-Lateral';'G_temp_sup-Plan_polar';'G_temp_sup-Plan_tempo';'G_temporal_inf';'G_temporal_middle';'Lat_Fis-ant-Horizont';'Lat_Fis-ant-Vertical';'Lat_Fis-post';'Pole_occipital';'Pole_temporal';'S_calcarine';'S_central';'S_cingul-Marginalis';'S_circular_insula_ant';'S_circular_insula_inf';'S_circular_insula_sup';'S_collat_transv_ant';'S_collat_transv_post';'S_front_inf';'S_front_middle';'S_front_sup';'S_interm_prim-Jensen';'S_intrapariet_and_P_trans';'S_oc_middle_and_Lunatus';'S_oc_sup_and_transversal';'S_occipital_ant';'S_oc-temp_lat';'S_oc-temp_med_and_Lingual';'S_orbital_lateral';'S_orbital_med-olfact';'S_orbital-H_Shaped';'S_parieto_occipital';'S_pericallosal';'S_postcentral';'S_precentral-inf-part';'S_precentral-sup-part';'S_suborbital';'S_subparietal';'S_temporal_inf';'S_temporal_sup';'S_temporal_transverse'};
    roi_lbls = {};
    for roiCode = roiCodes,     roi_lbls{end + 1} = dstrxCodeToLabel{str2num(roiCode{1})};  end

    % label the vertices according to the freesurfer labels
    roiVertexLabels = mx.freesurfer.fsRelabelToAreas(roi_lbls, annotColortable, annotVertexLabels);

    % extract the vertices/faces that belong to the ROI
    roiVertexIDs = find(~isnan(roiVertexLabels));
    roiFacesLabels = ismember(gPial.faces, roiVertexIDs);
    roiFacesLabels = all(roiFacesLabels, 2);
    [vertexMatrix, facesMatrix, ~] = mx.three_dimensional.extract3DFaces(gPial, roiFacesLabels == 1);
    gROIPial = gPial;
    gROIPial.vertices = vertexMatrix;
    gROIPial.faces = facesMatrix;
    
    % calculate the distance from each search point to each retrieval point
    % and determine which end points are within radius distance
    dist = (trcLineEndPoints(:, 1)' - vertexMatrix(:, 1)) .^ 2 + ...
           (trcLineEndPoints(:, 2)' - vertexMatrix(:, 2)) .^ 2 + ...
           (trcLineEndPoints(:, 3)' - vertexMatrix(:, 3)) .^ 2;
    dist = dist < (radius  .^ 2);

    % resolve which tract lines have both end-points within the ROIs radius
    b2 = reshape(any(dist, 1), 2, []);
    proxTrkLines = any(b2, 1);

end