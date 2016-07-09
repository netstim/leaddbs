function [minC, maxC] = coor_extremals(C)
% usage
%   [minC, maxC] = coor_extremals(C)
%
% input
%   C = cell matrix of nested cell/numeric matrices, with numeric or empty
%       matrices at the bottom nesting level, of coordinates. These belong
%       to different graphics objects.
%
% output
%   minC = minimum coordinate over all graphics objects
%   maxC = maximum coordinate over all graphics objects
%
% See also PLOT_SCALINGS, MIN_CELL, MAX_CELL.

% depends
%   min_cell, max_cell

maxC = max_cell(C);
minC = min_cell(C);
