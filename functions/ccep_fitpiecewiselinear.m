function F = ccep_fitpiecewiselinear(pp,Y,x)

% fitting a piecewise linear function with 2 segments to ccep data
%     Copyright (C) 2020  D Hermes and D van Blooijs
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

F = Y - (pp(1) + pp(2)*min(pp(4),x) + pp(3)*max(x-pp(4),0));

% --> intercept = pp(1)
% --> tipping point = pp(4)
% --> slope before tipping point = pp(2)
% --> slope after tipping point = pp(3)

%% fit with:
% figure
% 
% x = my_output(~isnan(my_output(:,2)),1)';
% y = my_output(~isnan(my_output(:,2)),2)'*1000;
% plot(x,y,'o')
% 
% my_options = optimset('Display','off','Algorithm','trust-region-reflective');
% 
%         [pp] = lsqnonlin(@(pp) ccep_fitpiecewiselinear(pp,y,x),...
%             [40 -1 0 20],[0 -Inf -Inf 10],[40 0 Inf 30],my_options);
% 
% x_fit = [1:50];
% 
% y_fit = (pp(1) + pp(2)*min(pp(4),x_fit) + pp(3)*max(x_fit-pp(4),0));
% 
% hold on
% plot(x_fit,y_fit)
