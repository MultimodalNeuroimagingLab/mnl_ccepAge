function tH = ieeg_RenderGiftiLabels(g,vert_label,cmapInput,roiNames,varargin)
% function to render a gifti 
% 
% input:
%     g: gifti file with faces and vertices
%     vert_label
%     cmapInput
%     roiNames
%     varargin{1}: sulcal_map, can be loaded with read_curv('lh/rh.sulc')
%
% output:
%   th: returns trimesh handle so you can change it
%       for example, tH.FaceAlpha = 0.5 will make the rendering transparent 
%
% Viewing Angle: can be changed with ecog_ViewLight(90,0), changes both
% angle and light accordingly
%
% DH 2017

if ischar(cmapInput)
    eval(['cmap =' cmapInput '(max(vert_label));']);
elseif isnumeric(cmapInput)
    cmap = cmapInput;
end

% convert surface labels into colors for vertices in mesh (c)
if isempty(varargin)
    c = 0.7+zeros(size(vert_label,1),3);
elseif ~isempty(varargin)
    sulcal_labels = varargin{1};
    c = 0.5+zeros(size(vert_label,1),3);
    c(sulcal_labels<0,:) = 0.7;
end

for k = 1:max(vert_label)
    c(ceil(vert_label)==k,:) = repmat(cmap(k,:),length(find(ceil(vert_label)==k)),1);
end

subplot(1,5,1:4)
tH = trimesh(g.faces, g.vertices(:,1), g.vertices(:,2), g.vertices(:,3), c); axis equal; hold on
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
l1 = light;
lighting gouraud
material([.3 .9 .2 50 1]); 
axis off
set(gcf,'Renderer', 'zbuffer')
view(270, 0);
set(l1,'Position',[-1 0 1])

if exist('roiNames','var')
    if ~isempty(roiNames)
        subplot(1,5,5),hold on
        if iscell(roiNames)
            for k = 1:length(roiNames)
                plot(1,k,'.','Color',cmap(k,:),'MarkerSize',30)
                text(1.03,k,roiNames{k},'Color',cmap(k,:),'VerticalAlignment','middle')
            end
        elseif isnumeric(roiNames)
            for k = 1:length(roiNames)
                plot(1,k,'.','Color',cmap(k,:),'MarkerSize',20)
            end
            % if numbers text at bottom, middle and top
            text(1.03,1,int2str(roiNames(1)),'Color',[0 0 0],'VerticalAlignment','middle')
            text(1.03,round(length(roiNames)/2),int2str(roiNames(round(length(roiNames)/2))),'Color',[0 0 0],'VerticalAlignment','middle')
            text(1.03,length(roiNames),int2str(roiNames(end)),'Color',[0 0 0],'VerticalAlignment','middle')
        end
        xlim([0.8 1.2]),ylim([0 length(roiNames)+1])
        axis off
    end
end

subplot(1,5,1:4)
