%
%   Extract and reconstruct part of a 3D object to make up the new object
%
%   [vertexMatrix, facesMatrix, vertexConversion] = extract3DByFaces(inputGifti, extractFaceIndices)
%   [vertexMatrix, facesMatrix, vertexConversion] = extract3DByFaces(vertices, faces, extractFaceIndices)
%   
%       inputGifti           = the input gifti to extract the 3D object from
%       vertices, faces      = the input object as a vertex matrix (N x 3) and faces matrix (N x 3)
%       extractFaceIndices   = a vector with the indices of faces which should be extracted
%
%
%   Returns: 
%       vertexMatrix     = the new vertex matrix
%       facesMatrix      = the new faces matrix (with vertex indices that
%                          correspond to the vertexMatrix)
%       vertexConversion = the conversion table, the first column holds the
%                          old vertex indices, the second the new vertex indices
% 
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [vertexMatrix, facesMatrix, vertexConversion] = extract3DFaces(varargin)
    
    % check the input (gifti or vertices-and-faces)
    vertices = [];
    faces = [];
    nextArgIndex = 1;
    if isstruct(varargin{1}) || isobject(varargin{1})
            % if struct or object
            
            if isa(varargin{1}, 'gifti')
                % gifti object
                
                vertices = varargin{1}.vertices;
                faces = varargin{1}.faces;
                
            elseif isstruct(varargin{1})
                % struct
                
                if isfield(varargin{1}, 'vertices')
                    vertices = varargin{1}.vertices;
                elseif isfield(varargin{1}, 'vert')
                    vertices = varargin{1}.vert;
                end
                
                if isfield(varargin{1}, 'faces')
                    faces = varargin{1}.faces;
                elseif isfield(varargin{1}, 'face')
                    faces = varargin{1}.face;
                elseif isfield(varargin{1}, 'tri')
                    faces = varargin{1}.tri;
                end
                
            end
            
            % the next argument after this object
            nextArgIndex = 2;
            
        elseif size(varargin{1}, 1) == 3 || size(varargin{1}, 2) == 3
            % if vector/matrix where the size of one dimenion is 3
            
            % retrieve the first and second arguments
            vertices = varargin{1};    
            faces = varargin{2};
            
            % the next argument after these two arguments
            nextArgIndex = 3;

    end
    
    % the next argument are the faces that we wish to extract
    extractFaceIndices = varargin{nextArgIndex};

    % retrieve the faces that are relevant to the mesh that we want to extract
    extractFacesMatrix = faces(extractFaceIndices, :);

    % get the vertex indices of the vertices that should be extracted
    vertexConversion = unique(extractFacesMatrix(:));
    vertexConversion = sort(vertexConversion, 'asc');

    % extract the vertices based on the vertex indices
    % (which will effectively will create the new vertex matrix)
    vertexMatrix = vertices(vertexConversion, :);

    % generate new vertex id's (which, because of the sorting and extraction above, are linear ascending from 1)
    vertexConversion(1:length(vertexConversion), 2) = 1:1:length(vertexConversion);

    % convert the old vertex indices to the new ones
    % (using the counter iVertex, which corresponds with the second column of the vertexConversion matrix)
    facesMatrix = nan(size(extractFacesMatrix));
    for iVertex = 1:size(vertexConversion, 1)
        facesMatrix(extractFacesMatrix == vertexConversion(iVertex, 1)) = iVertex;
    end
    
end
