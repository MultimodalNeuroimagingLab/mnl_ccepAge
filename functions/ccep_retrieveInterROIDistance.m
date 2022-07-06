%
%   Retrieve the distance in subject native space (mm) between one ROI and another ROI
%   [distance] = ccep_retrieveInterROIDistance(track, roi1, roi2)
%
%       trkFiles          = the tract files to load the MNI tracts from. If a cell-array is
%                           passed (e.g. for one file for each hemisphere) then each
%                           then each cell is processed seperately and output is given per tract file
%       subjectANTsFolder = the subject specific folder with the ANTs transformation files, this is 
%                           used to transform the tracts from MNI space to subject native space
%
%   Returns: 
%       trkDistance       = the distance in mm between the ROI's end-points
%       trkIndices        = the indices of the tract-lines that were used to calculate the distance
%
%
%   Note: For Linux and Mac, execute permissions might need to be set for
%         the leadDBS binary files in '<mnl_ccepBids>/external/leadDBS/ext_libs/ANTs'.
%         Navigate to that folder in the terminal and run 'chmod u+x antsApplyTransformsToPoints.*'
%
%
%
%   Max van den Boom, Multimodal Neuroimaging Lab, Mayo Clinic
%
function [trkDistance, trkIndices] = ccep_retrieveInterROIDistance(trkFiles, subjectANTsFolder, roi1, roi2)
    
    % setup leadDBS
    prefs = ea_prefs;
    
    % loop over the tract files (likely to be two, one for each hemisphere)
    for iSet = 1:length(trkFiles)
        
        %%
        %  Load the tract file and transform the line vertices from MNI space to native space
        
        % load
        [fibers, idx] = ea_trk2ftr(trkFiles{iSet}, 2);
        fibers = fibers(:, 1:3);

        % RAS to LPS, ANTs (ITK) uses LPS coords
        fibers(:, 1:2) = -fibers(:, 1:2);

        % Apply transform (input should be Nx3 = <points> x <X,Y,Z>)
        fibers = ea_ants_apply_transforms_to_points(subjectANTsFolder, fibers, 0);

        % LPS to RAS, restore to RAS coords
        fibers(:, 1:2) = -fibers(:, 1:2);


        
        %%
        %  Determine which tract-lines that are close enough to the ROIs
        
        % TODO, for now pick all
        trkIndices{iSet} = 1:length(idx);
        
        
        %%
        %  Calculate the average distance over the tract-lines
        
        % loop over the tract-lines
        trkLengths = nan(1, length(trkIndices{iSet}));
        for iTrk = 1:length(trkIndices{iSet})
            
            % determine the start- and end-vertex indices
            startV = 1;
            if iTrk > 1,   startV = sum(idx(1:trkIndices{iSet}(iTrk) - 1)) + 1;   end
            endV = startV + idx(trkIndices{iSet}(iTrk)) - 1;
            
            % calculate and store the length of the current tract line
            trkLengths(iTrk) = sqrt(sum(sum(diff(fibers(startV:endV, 1:3), 1, 1) .^ 2, 2)));
            
        end
        
        % store(/return) the average over all tract-lines
        trkDistance{iSet} = mean(trkLengths);
        
    end
    
end
