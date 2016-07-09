function [] = get_quiver_data

ax = gca;
h = get(ax, 'children');
hq = findobj(h, 'flat', 'type', 'hggroup');

q = get(hq);

c = get(hq, 'children');

get(c(1) )
xb = get(c(1), 'XData');
yb = get(c(1), 'YData');
zb = get(c(1), 'ZData');

ax1 = newax;

plotmd(ax1, [xb; yb; zb] )

hold(ax1, 'on')

get(c(1) )
xh = get(c(2), 'XData');
yh = get(c(2), 'YData');
zh = get(c(2), 'ZData');

plotmd(ax1, [xh; yh; zh] )
axis equal
