%   
%   Determine which points in 3D space are close to another set of points in the same space given a specific radius
%   [proximatePoints] = retrieveRadialProximatePoints(retrievalPoints, searchPoints, radius)
%
%   	retrievalPoints   = the points in 3D space (n-by-3 matrix) which are being searched in (to be retrieved)
%   	searchPoints      = the points in 3D space (n-by-3 matrix), where each point is a source in space from which to search
%   	radius            = the radius within each search point where a retrieval point will be considered close    
% 
%
%   Returns: 
%       proximatePoints   = cell array. Each cell represents an (input) search point, the values in the cell are the indices
%                           of the (input) retrieval points which are within the search point's radius
%   
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [proximatePoints] = retrieveRadialProximatePoints(retrievalPoints, searchPoints, radius)
    
    % allocate the return cell array
    n = size(searchPoints, 1);
    proximatePoints = cell(n, 1);
    
    % calculate the distance from each search point to each retrieval point
    dist = (searchPoints(:, 1)' - retrievalPoints(:, 1)) .^ 2 + ...
           (searchPoints(:, 2)' - retrievalPoints(:, 2)) .^ 2 + ...
           (searchPoints(:, 3)' - retrievalPoints(:, 3)) .^ 2;

    % determine which voxels are within radius distance
    dist = dist < (radius  .^ 2);
    
    % transfer the proximate point to a cell array
    [row, cols] = find(dist == 1);
    for i = 1:n
        proximatePoints{i} = row(cols == i);
    end
    
end
