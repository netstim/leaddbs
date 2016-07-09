function [] = heart3d
% File:      heart3d.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2011.02.07 - 2012.07.16
% Language:  MATLAB R2012a
% Purpose:   export 3D heart to u3d
% Copyright: Ioannis Filippidis, 2011-

dom = [-3, 3, -3, 3, -3, 3];
res = [100, 100, 100];

[X Y Z] = domain2meshgrid(dom, res);
F = -X.^2 .*Z.^3 -(9/80) .*Y.^2 .*Z.^3 ...
    +(X.^2 +(9/4) .*Y.^2 +Z.^2 -1).^3;

ax = newax;
p = patch(isosurface(X, Y, Z, F, 0), 'Parent', ax);
set(p, 'FaceColor', 'r', 'EdgeColor', 'r');
camlight
daspect(ax, [1 1 1] )
view(ax, 3)
axis(ax, 'tight')
axis(ax, 'equal')

fig2u3d(ax, 'heart3d', '-pdf')
