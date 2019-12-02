function ea_slideslices(h,~,cmd)

% update toggle state
uibjs=getappdata(gcf,'uibjs');
set(uibjs.slide3dtog,'State','on');
set(uibjs.rotate3dtog,'State','off');
if strcmp(cmd,'off')
    return
end

% reset button down function
hfig = h.Parent.Parent;
set(hfig,'WindowButtonDownFcn', []);

% enable click actions on surfaces (image slices)
ax = findobj(hfig.Children,'Type','axes');
set(findobj(ax.Children,'Type','surface'),'HitTest','on'); 
