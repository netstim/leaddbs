function [] = test_hexagram_shape

r1 = 2;
r2 = 1;

n = 6;
t = linspace(0, 2*pi, n+1);
x1 = r1 *[cos(t); sin(t) ];

t2 = t +pi /n;
x2 = r2 *[cos(t2); sin(t2) ];

m = 2 *n +1;
x(:, 1:2:m) = x1;
x(:, 2:2:m) = x2(:, 1:n);

ax = gca;
plotmd(ax, x)
axis(ax, 'equal')