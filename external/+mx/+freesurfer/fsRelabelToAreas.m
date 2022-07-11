%   
%   Relabel the annotation labels to match the indices of a given (string) list of annotation areas
%   reAnnotLabels = fsRelabelToAreas(fsAreas, annotColortable, annotVertexLabels)
%
%       fsAreas             = a string cell array of the areas to look up
%       annotColortable     = the color table, read in using freesurfer's read_annotation
%       annotVertexLabels   = the vector of vector annotation labels which should be relabeled
%
%
%   Returns: 
%       A vector of vertex annotation labels, where the values match the
%       indices of the given list of annotation areas and nan values for
%       the labels not in the list
%
%   Example:
%
%       fsAreas = { 'precentral', 'postcentral'}
%       [~, annotVertexLabels, annotColortable] = read_annotation(fsAnnotationFilepath);
%
%       reAnnotVertexLabels = fsRelabelToAreas(fsAreas, annotColortable, annotLabel);
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function reAnnotVertexLabels = fsRelabelToAreas(fsAreas, annotColortable, annotVertexLabels)
    
    % find the color codes for the required annotation areas
    areaColorCodes = mx.freesurfer.fsRetrieveColorCodes(fsAreas, annotColortable);

    % create a copy of the input list
    reAnnotVertexLabels = annotVertexLabels;
    
    % loop and replace the labels
    uniqueLabels = unique(reAnnotVertexLabels);
    for iLabel = 1:length(uniqueLabels)
        labelIndex = find(areaColorCodes == uniqueLabels(iLabel));
        if isempty(labelIndex)
            reAnnotVertexLabels(reAnnotVertexLabels == uniqueLabels(iLabel)) = nan;
        else
            reAnnotVertexLabels(reAnnotVertexLabels == uniqueLabels(iLabel)) = labelIndex;
        end
    end
    
end