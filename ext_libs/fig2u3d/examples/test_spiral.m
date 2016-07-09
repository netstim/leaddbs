%function [] = test_fig2u3d_spiral

clf
t = linspace(0, 10, 100);
x = [t; cos(t); sin(t) ];

ax = gca;
plotmd(ax, x, 'h')

fig2u3d(ax, 'test')

axis(ax, 'equal')

copyfile('test.u3d', '..\tex\personal\3dheart\img\test.u3d')
copyfile('test.vws', '..\tex\personal\3dheart\img\test.vws')
