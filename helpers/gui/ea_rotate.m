function ea_rotate(h,~,cmd)

ea_distogslide;
% get figure and axes
hfig = h.Parent.Parent;
ax = findobj(hfig.Children,'Type','axes');
% disable cliuck actions on surfaces (image slices)
set(findobj(ax.Children,'Type','surface'),'HitTest','off'); 
% set cam opts to the mouse
ea_mouse_camera(hfig);