function [varargout] = quivermd(ax, x, v, varargin)
%QUIVERMD     multi-dimensional quiver
%   QUIVERMD(AX, X, V, VARARGIN) plots the column vectors in matrix V
%   at the points with coordinates the column vectors in matrix X
%   within axes object AX using plot formatting options in VARARGIN.
%
% usage:
%   H = QUIVERMD(AX, X, V, VARARGIN)
%
% input
%   ax = axes handle (e.g. ax = gca)
%   x = matrix of points where vectors are plotted
%     = [#dim x #points]
%   v = matrix of column vectors to plot at points x
%     = [#dim x #points]
%   varargin = plot formatting
%
% output
%   h = handle to plotted object(s)
%
% example
%   x = linspace(0, 10, 20);
%   y = linspace(0, 10, 20);
%   [X, Y] = meshgrid(x, y);
%   x = [X(:), Y(:) ].';
%   v = [sin(x(1, :) ); cos(x(2, :) ) ];
%   quivermd(gca, x, v)
%
% See also PLOTMD, QUIVER, QUIVER3.
%
% File:      quivermd.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.01.22 - 
% Language:  MATLAB R2011b
% Purpose:   multi-dimensional quiver
% Copyright: Ioannis Filippidis 2011-

ndim = size(x, 1);

if ndim > 3
    warning('quivermd:ndim', '#dimensions > 3, plotting only 3D component.')
end

if ndim == 2
    h = quiver(ax, x(1, :), x(2, :), v(1, :), v(2, :), varargin{:} );
elseif ndim >= 3
    h = quiver3(ax, x(1, :), x(2, :), x(3, :),...
                    v(1, :), v(2, :), v(3, :), varargin{:} );
end

if nargout == 1
    varargout{1, 1} = h;
end
