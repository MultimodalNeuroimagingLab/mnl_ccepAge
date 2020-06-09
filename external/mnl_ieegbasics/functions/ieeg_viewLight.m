function ecog_ViewLight(th, phi)
%function loc_view(theta, phi) 
%this function orients the brain and always puts the lighting behind you
%theta and phi are in degrees, not radians
%make sure the brain plot is your current axis
%this function rotates the brain and lighting to the spherical angle 
%inputted.   it is non-standard b/c of matlab.  so, phi is really "elevation" and not
%phi from standard physics coordinates.  (0,0) is at the back of the brain.  for example: 
%loc_view(180,90) views from the top with the front downward, and
%loc_view(180,90) has the front upward, loc_view(90,0) is from the right,
%
%     Copyright (C) 2009 K.J. Miller, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
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
%
% DH 2017, updated for GiftiRender
    

view(th,phi); % change viewing direction of object in angle

view_pt=[cosd(th-90)*cosd(phi) sind(th-90)*cosd(phi) sind(phi)]; % light position for viewing angle

%in order to change the direction of the light, get children of graphics:
a=get(gca,'Children');
for k=1:length(a) % run through to find the light 
    if strcmp(a(k).Type,'light') %find the correct child (the one that is the light)
        set(a(k),'Position',view_pt) % change light position
    end
end
