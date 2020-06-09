function ieeg_elAdd_sizable(els,r2,varargin)

% function ieeg_elAdd_sizable(els,r2)
% 
% els: rows = electrodes, columns = xyz
% r2: for size
%
% ieeg_elAdd_sizable(els,r2,maxr2,max_elsize)  
%
% optional inputs maxr2 and max_elsize
%   maxr2 (varargin{1}) maximum for scale (if {''} absmax)
%       example: el_add_withr2size(els,r2,1)
%   max_elsize (varargin{2}) maximum electrode size: sometimes we want larger electrodes
%       default: 45 
%       example: el_add_withr2size(els,r2,1,45)
%
%     Copyright (C) 2006  K.J. Miller & D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
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

hold on

%gray2green
% cm1=[repmat([0 0 0],100,1)];
% cm1(1:10,2)=[0:0.7/9:0.7]';
% cm1(10:100,2)=[0.7:(1-0.7)/90:1]';
% cm1(1:10,1)=[0]';
% cm1(1:10,3)=[0]';
% cm1(20:100,1)=[0:1/80:1]';
% cm1=[repmat([0 1 0],100,1)];

cm1(1:100,2)=[0:1/99:1]';
cm1=[repmat([0 0 0],100,1)];
cm1(1:10,1)=[0:0.7/9:0.7]';
cm1(10:100,1)=[0.7:(1-0.7)/90:1]';
cm1(1:10,2)=[0]';
cm1(1:10,3)=[0]';
cm1(20:100,2)=[0:1/80:1]';


cm2=[repmat([0 0 0],100,1)];
cm2(1:10,3)=[0:0.7/9:0.7]';
cm2(10:100,3)=[0.7:(1-0.7)/90:1]';
cm2(1:10,2)=[0]';
cm2(1:10,1)=[0]';
cm2(20:100,2)=[0:1/80:1]';


maxr2=round(max(r2));
if abs(round(min(r2)))>maxr2
    maxr2=abs(round(min(r2)));
end

if isempty(varargin)
    r2 = r2/maxr2;% scaled to absmax
else
    r2 = r2/varargin{1};% scaled to varargin{1}
    
    % use varargin{2} to set maximum electrode size
    if isempty(varargin{2})
        max_elsize = 45;
    else
        max_elsize = varargin{2};
    end
end

elsize = [15:(max_elsize-15)/(100-1):max_elsize];

% electrode with r2:
for k=1:size(els,1)
    if ~isnan(r2(k))
        if abs(r2(k))>0.01
            ind_color=abs(round(100*r2(k)));
            if ind_color>100
                ind_color=100;
            end
            elsize_r2=elsize(ind_color);
            if r2(k)>0.01
                elcol_r2=cm1(ind_color,:); 
                plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize_r2)
                plot3(els(k,1),els(k,2),els(k,3),'.','Color',elcol_r2,'MarkerSize',elsize_r2-5)
            elseif r2(k)<0.01
                elcol_r2=cm2(ind_color,:);
                plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize_r2)
                plot3(els(k,1),els(k,2),els(k,3),'.','Color',elcol_r2,'MarkerSize',elsize_r2-5)
            end
        else
            plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize(1))
        end
    end
end


%% example code to plot color scale

 
% cm1(1:100,2)=[0:1/99:1]';
% cm1=[repmat([0 0 0],100,1)];
% cm1(1:10,1)=[0:0.7/9:0.7]';
% cm1(10:100,1)=[0.7:(1-0.7)/90:1]';
% cm1(1:10,2)=[0]';
% cm1(1:10,3)=[0]';
% cm1(20:100,2)=[0:1/80:1]';
% 
% cm2=[repmat([0 0 0],100,1)];
% cm2(1:10,3)=[0:0.7/9:0.7]';
% cm2(10:100,3)=[0.7:(1-0.7)/90:1]';
% cm2(1:10,2)=[0]';
% cm2(1:10,1)=[0]';
% cm2(20:100,2)=[0:1/80:1]';
% 
% % cm2=cm2(end:-1:1,:);
% % cm=[cm2; cm1];
% elsize=[15:(45-15)/(100-1):45];
% 
% figure('Color',[1 1 1],'Position',[30 50 50 300]),hold on
% for k=1:100
%     plot(1,k,'.','MarkerSize',elsize(k),'Color',cm1(k,:))
% end
%     
% for k=1:100
%     plot(1,-k,'.','MarkerSize',elsize(k),'Color',cm2(k,:))
% end
% 
% ylim([-120 120])
% set(gcf, 'PaperPositionMode', 'auto');
% print('-painters','-r300','-dpng',strcat(['./figures/colorscale_el_add_sizable']));
% print('-painters','-r300','-depsc',strcat(['./figures/colorscale_el_add_sizable']));

