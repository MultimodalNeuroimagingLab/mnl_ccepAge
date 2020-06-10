function el_add(els,elcol,msize)

%     Copyright (C) 2006  K.J. Miller
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
if ~exist('msize','var')
    msize=20; %marker size
end

if exist('elcol')==0, 
    elcol='r'; %default color if none input
end

hold on, plot3(els(:,1),els(:,2),els(:,3),'.','Color', elcol,'MarkerSize',msize)
