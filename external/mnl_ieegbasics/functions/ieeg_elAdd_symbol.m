function el_add_symbol(els,elcol,msize,s_b)

if ~exist('msize','var')
    msize=20; %marker size
end

if exist('elcol')==0, 
    elcol='r'; %default color if none input
end

hold on, plot3(els(:,1),els(:,2),els(:,3),s_b,'Color', elcol,'MarkerSize',msize,'MarkerEdgeColor',elcol,'MarkerFaceColor',elcol)
%hold on, plot3(els(:,1),els(:,2),els(:,3),'s','MarkerSize',msize,'MarkerEdgeColor','r','MarkerFaceColor','r')
