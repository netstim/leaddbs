function ea_slideslices(h,~,cmd)

ea_distogrotate;
hfig = h.Parent.Parent;
set(hfig,'WindowButtonDownFcn', []); % reset button down function
ax = findobj(hfig.Children,'Type','axes');

set(findobj(ax.Children,'Type','surface'),'HitTest',cmd); % enable/disable click actions on surfaces (image slices)
