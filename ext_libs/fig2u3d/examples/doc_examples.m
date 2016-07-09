function [] = doc_examples
% File:      doc_examples.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.14 - 2012.06.16
% Language:  MATLAB R2012a
% Purpose:   documentation examples of fig2u3d
% Copyright: Ioannis Filippidis, 2012-

%fname = two_variable_func;
%fname = three_variable_level_set;
%fname = molecule_structure;
%fname = robotic_manipulator;
fname = fluid_flow;
fig2u3d(gca, fname, '-pdf')

function [fname] = two_variable_func
%f(x,y) surface, gradient, contours
ax = newax;
hold(ax, 'on')

dom = 10 *[-1, 1, -1, 1];
res = [30, 40];

q = domain2vec(dom, res);
x = q(1, :);
y = q(2, :);
f = sin(sqrt(0.5 .*x.^2 +y.^2) ) +2;
v = bsxfun(@times, cos(sqrt(0.5 .*x.^2 +y.^2) ) .*(0.5 .*x.^2 +y.^2).^(-0.5), [0.5 .*x; y] );

vsurf(ax, q, f, res)
vcontour(ax, q, f, res)
quivermd(ax, q, v)
axis tight
view(ax, 3)

fname = 'two_var_func';

function [fname] = three_variable_level_set
% Elliptic supercylide level set surface, gradient and 2D field section
example_supercyclide
fname = 'example_supercyclide';

function [fname] = molecule_structure
draw_crystal_lattice
fname = 'diamond';

function [fname] = robotic_manipulator
% requires the Robotics Toolbox by Peter I. Corke
mdl_puma560
p560.plot(zeros(1,6) )
fname = 'puma560';

function [fname] = fluid_flow
Wind
fname = 'wind';
