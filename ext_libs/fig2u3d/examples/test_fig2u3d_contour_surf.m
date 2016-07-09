%function [] = test_fig2u3d_contour_surf

clf

res = [30, 31];
q = domain2vec([-2, 2, -2, 2], res);

x = q(1, :);
y = q(2, :);
z = -sin(x) .*cos(y);

ax = gca;
hold(ax, 'on')
%vcontourf(ax, q, z, res);
vcontour(ax, q, z, res);
vsurf(ax, q, z, res);
plot_scalings(ax, 0)

%{
h = get(gca, 'children');
hp = get(h(2), 'Children');
fc = get(hp(1), 'Faces');

x = v(:, 1);
y = v(:, 2);
x(isnan(x) ) = [];
y(isnan(y) ) = [];

dt = DelaunayTri(x, y);
triplot(dt, 'Parent', newax)
%}
fig2u3d(ax, 'test')

axis(ax, 'equal')
