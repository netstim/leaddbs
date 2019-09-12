function ea_slideslices(h,~,cmd)

hfig = h.Parent.Parent;
set(hfig,'WindowButtonDownFcn', []); % reset button down function
ax = findobj(hfig.Children,'Type','axes');
set(findobj(ax.Children,'Type','surface'),'HitTest','on'); % enable click actions on surfaces (image slices)