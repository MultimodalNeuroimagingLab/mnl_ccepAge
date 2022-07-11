% Display (gifti) surface
%
%   viewgii( surf1, [surf2], ..., pointMatrix1, lineMatrix1, ...)
%
%       surf#                       = surface(s) to display. Either a surface struct (e.g. gifti), handle or filename can be passed
%
%
%   Optional arguments:
%
%       <pointMatrix#>              = Nx3 matrix of points to plot. Each row in the
%                                     matrix (first dimension) represents a point; the
%                                     second dimension (colums) represent the x, y and z
%                                     coordinates of each point
%       <lineMatrix#>               = Nx6 matrix of lines to plot. Each row in the
%                                     matrix (first dimension) represents a line; the
%                                     second dimension (colums) represent the x1, y1, z1, x2, y2, z2 coordinates
%       <textCellArray>             = 
%
%       Trans or Transparent<alpha> = display the surface with transparency. Alpha value between 0 (transparent) and 1 (opaque)
%       Wireframe                   = display the wireframe of the surface
%       Vertices                    = display the vertices of the surface
%       FaceCenters                 = display a point at the center of each face
%       FaceNormals                 = display the normal of each face
%       FaceNormalsInv              = display the inverted normal of each face
%       VertexNormals               = display the normal of each vertex
%       VertexNormalsInv            = display the inverted normal of each vertex
%       FaceAreas                   = calculate the area of each face and display as text over each face
%       Merge                       = Merge the different surfaces together, each surface will be displayed in different color
%       Disks<radius>               = Show points as disks
%       PointProject                = Project the points onto surface
%       WireSpheres<radius>         = Show wire-spheres around each point
%       WireCylinders<radius>       = Show wire-cylinders around each line
%
%
%
%   Example 1:
%
%       gSurfPial       = gifti(giiPialFilepath);
%       viewgii(gSurfPial);
%
%
%   Example 2:
%
%       gSurfPial       = gifti(giiPialFilepath);
%       gSurfInflated   = gifti(giiInflatedFilepath);
%
%       viewgii(gSurfPial, gSurfInflated);
%
%
%   Example 3:
%
%       gSurfPial   = gifti(giiPialFilepath);
%       pointMatrix = [ 6.51, -47.34, 54.79; ...
%                       1.50, -46.12, 55.84];
%       lineMatrix  = [  6.51, -47.34, 54.79, 8.51, -69.34, 94.79; ...
%                        1.50, -46.12, 55.84, 1.50, -96.12, 95.84];
%       viewgii(gSurfPial, pointMatrix, lineMatrix);
%
%
%  Example 4:
%
%       gSurfPial   = gifti(giiPialFilepath);
%       pointMatrix = [ 6.51, -47.34, 54.79; ...
%                       1.50, -46.12, 55.84];
%       viewgii(gSurfPial, pointMatrix, 'WireSpheres7');    % radius of 7
%
%
%
%   Copyright (C) 2022, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function viewGii( varargin )
    giftis = {};
    pointMatrices = {};
    lineMatrices = {};
    textArrays = {};
    displayArgs = {};
    mergeGiftis = 0;
    showPointDisks = 0;
    showPointProject = 0;
    showPointWireSpheres = 0;
    showLineWireCylinders = 0;
    
    % todo, extend with more colors
    colors = {  [0 0 1], [1,0,0], [0 1 0], [1 1 0], [1,0,1], [0,1,1], [1, 0.65, 0], ...
                [0, 0.75, 0.75], [0.75, 0.75, 0], [0, 0.5, 0], [0.494, 0.184, 0.556], ...
                [0.850, 0.325, 0.098], [0.635, 0.078, 0.184], [0.929, 0.694, 0.125], ...
                [0, 0.447, 0.741], [0.301, 0.745, 0.933], [0.635, 0.078, 0.184]};
    
    % if an empty 3d object is given, then create a dummy object
    if nargin > 0 && isempty(varargin{1})
        dummy = [];
        dummy.vertices = [0 0 0; 0 0 1; 0 1 1];
        dummy.tri = [1,2,3];
        varargin{1} = gifti(dummy);
    end
    
    % loop through the input arguments
    for i = 1:nargin
        if isstruct(varargin{i}) || isa(varargin{i}, 'gifti') || (length(varargin{i}) == 1 && ishandle(varargin{i}))
            % surface, if struct or handle
            
            % add as gifti
            giftis{end + 1} = varargin{i};

        elseif isa(varargin{i}, 'char')
            % single string, check if the argument is a string (char array)
            
            % add as display argument
            displayArgs{end + 1} = varargin{i};
            
        elseif (size(varargin{i}, 2) == 3 || (size(varargin{i}, 1) == 3 && size(varargin{i}, 2) ~= 6)) && size(varargin{i}, 3) == 1
            % 2D point matrix, if 2D matrix where the size of the second dimenion is 3
            % or accidentally transposed (except for when the number of the second
            % dimension is six, because a second dimension of 6 signifies a
            % lineset)
            
            if size(varargin{i}, 2) ~= 3 && size(varargin{i}, 1) == 3
                varargin{i} = varargin{i}';
            end
            
            % add as point matrix
            pointMatrices{end + 1} = varargin{i};

        elseif (size(varargin{i}, 2) == 3) && size(varargin{i}, 3) > 1
            % multiple 2D point matrices, if vector/matrix where the size of either the 
            % first or second dimenion is 3 and the third dimension is more than 1
            
            if size(varargin{i}, 2) ~= 3 && size(varargin{i}, 1) == 3
                varargin{i} = permute(varargin{i}, [2 1 3]);
            end
            
            % add the point matrices
            for iPointMatrix = 1:size(varargin{i}, 3)
                pointMatrices{end + 1} = varargin{i}(:, :, iPointMatrix);
            end
            
        elseif size(varargin{i}, 1) == 6 || size(varargin{i}, 2) == 6
            % line matrix, if vector/matrix where the size of one dimenion is 6
            
            if size(varargin{i}, 2) ~= 6 && size(varargin{i}, 1) == 6
                varargin{i} = varargin{i}';
            end
            
            % add as line matrix
            lineMatrices{end + 1} = varargin{i};
            
        elseif iscell(varargin{i}) || isstring(varargin{i}) && (size(varargin{i}, 1) == 1 || size(varargin{i}, 2) == 1)
            % text, if vector is cell-array or a string-array and size of one dimenion is 1
            
            if size(varargin{i}, 2) ~= 1 && size(varargin{i}, 1) == 1
                varargin{i} = varargin{i}';
            end
            
            % add as text array
            textArrays{end + 1} = varargin{i};
            
        end
    end
    
    % retrieve the screen size and determine the display size
    screensize = get(0, 'Screensize');
    displaySize = (screensize(4) - 100) / 2;
	
    % loop through the giftis
    for i = 1:length(giftis)
        
        % set the config
        toolConfig = {};
        toolConfig.hideToolWindow           = 1;
        toolConfig.yokeCam                  = 1;
        
        % loop through the display arguments
        for j = 1:length(displayArgs)
            displayArgument = displayArgs{j};
            
            if contains(displayArgument, 'trans', 'IgnoreCase', true) || contains(displayArgument, 'transparent', 'IgnoreCase', true)
                
                
                try
                    if contains(displayArgument, 'transparent', 'IgnoreCase', true)
                        alpha = str2num(displayArgument(strfind(lower(displayArgument), 'transparent') + 11:end));
                    elseif contains(displayArgument, 'trans', 'IgnoreCase', true)
                        alpha = str2num(displayArgument(strfind(lower(displayArgument), 'trans') + 5:end));
                    end
                    if isempty(alpha),      alpha = .7;     end
                    if alpha < 0,           alpha = 0;      end
                    if alpha > 1,           alpha = 1;      end
                    
                    toolConfig.defaultBackgroundAlpha = alpha;
                    
                catch
                end

            end
            
            
            if strcmpi(displayArgument, 'showWireframe') || strcmpi(displayArgument, 'Wireframe')
                toolConfig.showWireframe = 1;
            end
            if strcmpi(displayArgument, 'showVertices') || strcmpi(displayArgument, 'Vertices')
                toolConfig.showVertices = 1;
            end
            if strcmpi(displayArgument, 'showFaceCenters') || strcmpi(displayArgument, 'FaceCenters')
                toolConfig.showFaceCenters = 1;
            end
            if strcmpi(displayArgument, 'showNormals') || strcmpi(displayArgument, 'normals') || strcmpi(displayArgument, 'showFaceNormals') || strcmpi(displayArgument, 'FaceNormals')
                toolConfig.showFaceNormals = 3;
            end
            if strcmpi(displayArgument, 'showNormalsInv') || strcmpi(displayArgument, 'showFaceNormalsInv') || strcmpi(displayArgument, 'FaceNormalsInv') || ...
                strcmpi(displayArgument, 'showInvNormals') || strcmpi(displayArgument, 'showFaceInvNormals') || strcmpi(displayArgument, 'FaceInvNormals')
                toolConfig.showFaceNormals = -3;
            end
            if strcmpi(displayArgument, 'showVertexNormals') || strcmpi(displayArgument, 'VertexNormals')
                toolConfig.showVertexNormals = 3;
            end
            if strcmpi(displayArgument, 'showVertexNormalsInv') || strcmpi(displayArgument, 'VertexNormalsInv') || ...
                strcmpi(displayArgument, 'showVertexInvNormals') || strcmpi(displayArgument, 'VertexInvNormals')
                toolConfig.showVertexNormals = -3;
            end
            if strcmpi(displayArgument, 'showFaceAreas') || strcmpi(displayArgument, 'showAreas') || strcmpi(displayArgument, 'FaceAreas')
                toolConfig.showFaceAreas = 1;
            end
            if strcmpi(displayArgument, 'merge') || strcmpi(displayArgument, 'mergegifti') || strcmpi(displayArgument, 'mergegiftis')
                mergeGiftis = 1;
            end
            
            if contains(displayArgument, 'showPointDisks', 'IgnoreCase', true) || contains(displayArgument, 'showDisks', 'IgnoreCase', true) || contains(displayArgument, 'Disks', 'IgnoreCase', true)
                try
                    radius = str2num(displayArgument(strfind(lower(displayArgument), lower('Disks')) + 5:end));
                    if isempty(radius) || radius < 1,   radius = 3; end
                    showPointDisks = radius;
                catch
                    showPointDisks = 3;
                end
            end
            if contains(displayArgument, 'showPointProject', 'IgnoreCase', true) || contains(displayArgument, 'PointProject', 'IgnoreCase', true) || contains(displayArgument, 'Project', 'IgnoreCase', true)
                showPointProject = 1;
            end
            if contains(displayArgument, 'showPointWireSpheres', 'IgnoreCase', true) || contains(displayArgument, 'WireSphere', 'IgnoreCase', true)
                try
                    radius = str2num(displayArgument(strfind(lower(displayArgument), lower('WireSpheres')) + 11:end));
                    if isempty(radius) || radius < 1,   radius = 3; end
                    showPointWireSpheres = radius;
                catch
                    showPointWireSpheres = 3;
                end
            end
            if contains(displayArgument, 'showLineWireCylinders', 'IgnoreCase', true) || contains(displayArgument, 'WireCylinder', 'IgnoreCase', true)
                try
                    radius = str2num(displayArgument(strfind(lower(displayArgument), lower('WireCylinders')) + 13:end));
                    if isempty(radius) || radius < 1,   radius = 3; end
                    showLineWireCylinders = radius;
                catch
                    showLineWireCylinders = 3;
                end
            end
            
        end
        
        % check if the giftis should be merged and there are multi giftis
        if mergeGiftis == 1 && length(giftis) > 1
            
            % variabl 
            numVertices = [];
            
            % loop through the gifti
            for iGifti = 1:length(giftis)
                
                % make sure the input is in gifti format
                giftis{iGifti} = gifti(giftis{iGifti});
                
                % store the number of vertices before merging
                numVertices(iGifti) = size(giftis{iGifti}.vertices, 1);
                
            end
            
            % calculate the total number of vertices in the final gifti
            totalVertices = sum(numVertices);
            
            % define the colors for each consecutive gifti
            colormaps = {'red', 'blue', 'green', 'violet', 'yellow', 'cyan', 'orange'};
            while (length(colormaps) < length(giftis) - 1)
                colormaps = repmat(colormaps, 1, 2);
            end
            
            % initiale a list of vertexlabels for each gifti (will result
            % in one overlay per gifti)
            vertexLabels = nan(totalVertices, length(giftis) - 1);

            % loop through the second and more giftis
            for iGifti = 2:length(giftis)


                % determine the offset of the indices for the vertices in the combined matrices.
                vertexIndexOffset = size(giftis{1}.vertices, 1);
                
                % increase the indices of the gifti's vertices with the offset
                newFaces = giftis{iGifti}.faces + vertexIndexOffset;

                % add the faces (with reindexed vertices) of the current gifti to the total
                giftis{1}.faces = [giftis{1}.faces; newFaces];

                % add the vertices of the current gifti to the total
                giftis{1}.vertices = [giftis{1}.vertices; giftis{iGifti}.vertices];
                
                % set an overlay the gifti
                vertexLabels(vertexIndexOffset + 1:vertexIndexOffset + numVertices(iGifti), iGifti - 1) = 1;
                toolConfig.(['overlay', num2str(iGifti - 1), 'PosEnabled']) = 1;
                toolConfig.(['overlay', num2str(iGifti - 1), 'PosColormap']) = colormaps(iGifti - 1);
                toolConfig.(['overlay', num2str(iGifti - 1), 'PosMin']) = 1;
                toolConfig.(['overlay', num2str(iGifti - 1), 'PosMax']) = 1;
                toolConfig.(['overlay', num2str(iGifti - 1), 'NegEnabled']) = 0;
                
            end
            for iGifti = 2:length(giftis)
                toolConfig.(['overlay', num2str(iGifti - 1)]) = vertexLabels(:, iGifti - 1);
            end
            
            % remove the rest of the gifti's (they are merged into the first)
            giftis(2:end) = [];
            
        end
        
        % initialize the color index
        colorIndex = 1;
        
        % loop through the point vectors/matrices
        for j = 1:length(pointMatrices)
            pointMarker = '*';
            pointSize = 8;
            pointColor = colors{colorIndex};
            colorIndex = colorIndex + 1;
            if colorIndex > length(colors), colorIndex = 1; end
            
            toolConfig.(['pointSet', num2str(j)]) = pointMatrices{j};
            toolConfig.(['pointSet', num2str(j), 'Marker']) = pointMarker;
            toolConfig.(['pointSet', num2str(j), 'Color']) = pointColor;
            if showPointDisks > 0
                toolConfig.(['pointSet', num2str(j), 'Type']) = 'Disks';
                toolConfig.(['pointSet', num2str(j), 'Size']) = showPointDisks;
            else
                toolConfig.(['pointSet', num2str(j), 'Size']) = pointSize;    
            end
            if showPointProject == 1
                toolConfig.(['pointSet', num2str(j), 'ProjectToHull']) = 1;    
            end
            if showPointWireSpheres > 0
                toolConfig.(['pointSet', num2str(j), 'WireSphereRad']) = showPointWireSpheres;
                
                % set the wirecolor and make the (center-)point color much brighter
                toolConfig.(['pointSet', num2str(j), 'WireSphereCol']) = pointColor;
                toolConfig.(['pointSet', num2str(j), 'Color']) = toolConfig.(['pointSet', num2str(j), 'Color']) + ((1 - toolConfig.(['pointSet', num2str(j), 'Color'])) * 0.5);
                
            end
            
        end
        
        % loop through the line vectors/matrices
        for j = 1:length(lineMatrices)
            lineSize = 1;
            lineColor = colors{colorIndex};
            colorIndex = colorIndex + 1;
            if colorIndex > length(colors), colorIndex = 1; end
            
            toolConfig.(['lineSet', num2str(j)]) = lineMatrices{j};
            toolConfig.(['lineSet', num2str(j), 'Color']) = lineColor;
            toolConfig.(['lineSet', num2str(j), 'Size']) = lineSize;
            if showLineWireCylinders > 0
                toolConfig.(['lineSet', num2str(j), 'WireCylinderRad']) = showLineWireCylinders;

                % set the wirecolor and make the (center-)line color much brighter
                toolConfig.(['lineSet', num2str(j), 'WireCylinderCol']) = lineColor;
                toolConfig.(['lineSet', num2str(j), 'Color']) = toolConfig.(['lineSet', num2str(j), 'Color']) + ((1 - toolConfig.(['lineSet', num2str(j), 'Color'])) * 0.5);
                
            end
            
        end
        
        % try to determine the number of faces and the number of vertices
        numFaces = 0;
        numVertices = 0;
        if isa(giftis{i}, 'gifti')
            % gifti object
            
            numFaces = size(giftis{i}.faces, 1);
            numVertices = size(giftis{i}.vertices, 1);
            
        elseif isstruct(giftis{i})
            % struct
            
            if isfield(giftis{i}, 'faces')
                numFaces = size(giftis{i}.faces, 1);
            elseif isfield(giftis{i}, 'face')
                numFaces = size(giftis{i}.face, 1);
            elseif isfield(giftis{i}, 'tri')
                numFaces = size(giftis{i}.tri, 1);
            end
            
            if isfield(giftis{i}, 'vertices')
                numVertices = size(giftis{i}.vertices, 1);
            elseif isfield(giftis{i}, 'vert')
                numVertices = size(giftis{i}.vert, 1);
            end
            
        end
        
        
        % loop through the text arrays
        for j = 1:length(textArrays)
            
            % check if the number of texts matches either the number of
            % faces or the number of vertices
            if length(textArrays{i}) == numFaces
                
                toolConfig.('faceText') = textArrays{j};
                toolConfig.('faceTextColor') = 'k';
                toolConfig.('faceTextSize') = 8;

            elseif length(textArrays{i}) == numVertices
            
                toolConfig.('vertexText') = textArrays{j};
                toolConfig.('vertexTextColor') = 'k';
                toolConfig.('vertexTextSize') = 8;
            
                % TODO: check if it matches any of the point matrices
                
                
            end
            
            for k = 1:length(pointMatrices)
                if length(textArrays{i}) == size(pointMatrices{k}, 1)
                    toolConfig.(['pointSet', num2str(k), 'Text']) = textArrays{j};
                    toolConfig.(['pointSet', num2str(k), 'TextColor']) = 'k';
                    toolConfig.(['pointSet', num2str(k), 'TextSize']) = 8;
                end
            end
            
        end        
        
        % display the gifti
        mx.three_dimensional.giftiTools(giftis{i}, toolConfig);
        
        % tile the position of the newly created windows
        current = gcf;
        y = mod(i - 1, 2) * (displaySize);
        y = screensize(4) - y - displaySize;
        x = floor((i - 1) / 2) * (displaySize + 10);
        set(current, 'OuterPosition', [x, y, displaySize, displaySize]);
        
        % weird behaviour in matlab, the for-loop will continue despite the
        % end condition being reached (in this case length(giftis))
        %
        % force the loop to end, do it here manually
        if i == length(giftis),  break; end
        
    end
    
end
