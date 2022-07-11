% 
%   Retrieve the color codes for a given (string) list of annotation areas
%   fsAreaColorCodes = fsRetrieveColorCodes(fsAreas, annotColortable)
%
%       fsAreas           = a string cell array of the areas to look up
%       annotColortable   = the color table, read in using freesurfer's read_annotation
%
% 
%   Returns: 
%       An array with the color codes of each annotated area in the
%       input list. A cell will hold nan if the area or matching color
%       code was not found
%
%   Example:
%
%       fsAreas = { 'precentral', 'postcentral'}
%       [annotVertices, annotLabel, annotColortable] = read_annotation(fsAnnotationFilepath);
%
%       fsAreaColorCodes = fsRetrieveColorCodes(fsAreas, annotColortable);
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function fsAreaColorCodes = fsRetrieveColorCodes(fsAreas, annotColortable)
    fsAreaColorCodes = [];
    
    % loop through the areas
    for iFsArea = 1:length(fsAreas)

        % find the indices for the area names
        fsAreaIndex = find(strcmp(annotColortable.struct_names, fsAreas{iFsArea}));

        % retrieve the color code
        if isempty(fsAreaIndex)
            fsAreaColorCodes(end + 1) = nan;
        else
            fsAreaColorCodes(end + 1) = annotColortable.table(fsAreaIndex, 5);
        end

    end
    
end