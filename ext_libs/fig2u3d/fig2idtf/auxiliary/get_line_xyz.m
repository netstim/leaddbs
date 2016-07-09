function [p] = get_line_xyz(h)
% File:      get_line_xyz.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.15
% Language:  MATLAB R2012a
% Purpose:   get line point coordinates and fix z if empty
% Copyright: Ioannis Filippidis, 2012-

% get defined data-points
x = get(h, 'XData');
y = get(h, 'YData');
z = get(h, 'ZData');

% 2d line ?
if isempty(z)
    z = zeros(size(x) );
end

p = [x; y; z];
