function ea_rotate(h,~,cmd)

% update toggle state
uibjs=getappdata(gcf,'uibjs');
set(uibjs.rotate3dtog,'State','on');
set(uibjs.slide3dtog,'State','off');
if strcmp(cmd,'off')
    return
end

% get figure and axes
hfig = h.Parent.Parent;
ax = findobj(hfig.Children,'Type','axes');

% disable click actions on surfaces (image slices)
set(findobj(ax.Children,'Type','surface'),'HitTest','off');

% set cam opts to the mouse
ea_mouse_camera(hfig);
