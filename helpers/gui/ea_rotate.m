function ea_rotate(h,~,cmd)

ea_distogslide;
% get figure and axes
hfig = h.Parent.Parent;
ax = findobj(hfig.Children,'Type','axes');
% disable click actions on surfaces (image slices)
set(findobj(ax.Children,'Type','surface'),'HitTest','off');
% set cam opts to the mouse
if strcmp(cmd,'on')
    ea_mouse_camera(hfig);
else
    set(hfig,'WindowButtonDownFcn', []);
end