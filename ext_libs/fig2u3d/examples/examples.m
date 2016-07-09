function [] = examples
% File:      examples.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10 - 2012.07.24
% Language:  MATLAB R2012a
% Purpose:   examples for using fig2u3d
% Copyright: Ioannis Filippidis, 2012-

% depends
%   plotmd, quivermd

ax = gca;

%% plot
% uncomment each example to use it
fname = peaks_contour(ax);
% fname = peaks_facecolors(ax);
% fname = peaks_sphere(ax);
% fname = mixed_scene(ax);

%% export
fig2u3d(ax, fname)

% uncomment to use other variants,
% ensure LaTeX & required packages are installed
%fig2pdf3d(ax, fname, 'media9', 'xelatex')
%fig2pdf3d(ax, fname, 'media9', 'pdflatex')
%fig2pdf3d(ax, fname, 'movie15', 'pdflatex')

function [fname] = peaks_contour(ax)
[x, y, z] = peaks;
surf(ax, x, y, z+3)
hold(ax, 'on')
contour(ax, x, y, z)

fname = 'peaks_contour';

function [fname] = peaks_facecolors(ax)
[x, y, z] = peaks;
[n, m] = size(x);
c = 10 *rand(n-1, m-1);

surf(ax, x, y, z +5, z +5)

hold(ax, 'on')
surf(ax, x, y, z, c)

[x, y, z] = sphere(20);
[n, m] = size(x);
c = 10 *rand(n-1, m-1);
surf(ax, x, y, z, c)

fname = 'peaks_facecolors';

function [fname] = peaks_sphere(ax)
[x, y, z] = peaks;
surf(ax, x, y, z)

hold(ax, 'on')

[x, y, z] = sphere(20);
surf(ax, x, y, z)

fname = 'peaks_sphere';

function [fname] = mixed_scene(ax)
%% init
hold(ax, 'on')

%% surf
[x, y, z] = sphere(20);
surf(ax, x, y, z)

%% line
t = linspace(0, 10, 100);
x = [cos(t); sin(t); t];
plotmd(ax, x, 'b-*')
plotmd(ax, x+1, 'ro') % 2nd line

%% quiver3
n = 8;
x = 10 *rand(3, n);
v = rand(3, n);
quivermd(ax, x, v, 'g')

%% view
axis(ax, 'equal')
axis(ax, 'tight')
view(ax, 3)

fname = 'mixed_scene';
