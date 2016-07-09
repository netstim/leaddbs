function [] = draw_crystal_lattice
% File:      draw_crystal_lattice.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2011.01.30 - 
% Language:  MATLAB R2012a
% Purpose:   draw crystal lattices
% Copyright: Ioannis Filippidis, 2011-

selection = 'diamond'; % 'diamond', 'copper'

a = 1;
c = a /4;

r_diamond_cubic = c *[0, 0, 0;
        2, 0, 0;
        0, 2, 0;
        2, 2, 0;
        1, 1, 0; %
        1, 0, 1;
        0, 1, 1;
        1, 2, 1;
        2, 1, 1; %
        0, 0, 2;
        2, 0, 2;
        0, 2, 2;
        2, 2, 2;
        1, 1, 2; %
        0.5, 0.5, 0.5;
        1.5, 1.5, 0.5;
        0.5, 1.5, 1.5;
        1.5, 0.5, 1.5]';

r_fcc = c *[0, 0, 0;
        2, 0, 0;
        0, 2, 0;
        2, 2, 0;
        1, 1, 0; %
        1, 0, 1;
        0, 1, 1;
        1, 2, 1;
        2, 1, 1; %
        0, 0, 2;
        2, 0, 2;
        0, 2, 2;
        2, 2, 2;
        1, 1, 2]'; %

i_diamond_cubic =[15, 1;
     15, 5;
     15, 6;
     15, 7; %
     16, 5;
     16, 4;
     16, 8;
     16, 9; %
     17, 7;
     17, 12;
     17, 14;
     17, 8; %
     18, 6;
     18, 9;
     18, 14;
     18, 11];
%

i_fcc = [1, 2;
    2, 4;
    2, 3;
    1, 3;
    1, 4;
    3, 4; %
    1, 10;
    2, 10;
    2, 11;
    1, 11;
    10, 11; %
    3, 12;
    4, 13;
    12, 13;
    4, 12;
    3, 13; %
    10, 12;
    10, 13;
    13, 11;
    12, 11;
    1, 12;
    3, 10; %
    2, 13;
    4, 11];

switch selection
    case 'diamond'
        r = r_diamond_cubic;
        i1 = i_diamond_cubic;
    case 'copper'
        r = r_fcc;
        i1 = i_fcc;
    otherwise
        disp('Unknown material selection.')
end

ax = newax;
hold(ax, 'on')
for i=1:size(r, 2)
    h = drawSphere(r(:, i), 0.05);
    set(h, 'FaceColor', 'r')
end
grid(ax, 'on')
box(ax, 'on')
axis(ax, 'equal')
axis(ax, 'tight')
view(ax, 3)
%maximize(get(ax, 'Parent') )

for i=1:size(i1, 1)
    a = i1(i, 1);
    b = i1(i, 2);
    drawCylinder([r(:, a)', r(:, b)', 0.01] )
end

fname = selection;
fig2u3d(gca, fname, '-pdf');
