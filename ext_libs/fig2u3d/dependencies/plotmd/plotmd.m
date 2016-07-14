function [varargout] = plotmd(ax, x, varargin)
%PLOTMD     multi-dimensional plot of 2/3D column vectors
%   PLOTMD(AX, X, VARARGIN) plots the points in matrix X in the axes with
%   handle AX using the plot formatting options in VARARGIN. X must be
%   a matrix whose columns are the 2D or 3D vectors to plot.
%
% usage
%   h = plotmd(ax, x, varargin)
%
% input
%   ax = axes handle (e.g. ax = gca) | nan (to turn off plotting)
%   x = matrix of points to plot
%     = [#dim x #pnts]
%   varargin = plot formatting
%
% output
%   h = handle to plotted object(s)
%
% usage example: plot 10 random 3D points
%   ax = gca;
%   ndim = 3;
%   npoints = 10;
%   x = rand(ndim, npoints);
%   h = plotmd(ax, x, 'ro');
%
% See also PLOT, PLOT3.
%
% File:      plotmd.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2011.11.12 - 2012.07.02
% Language:  MATLAB R2011b
% Purpose:   plot or plot3 of matrix of column vector points
% Copyright: Ioannis Filippidis 2011-

% ax empty ?
if isempty(ax)
    warning('plotmd:ax', 'Axes handle ax is empty. No graphical output.')
    return
end

% no plotting output (silent)
% if isnan(ax)
%     %disp('Axes handle is nan. No graphical output.')
%     return
% end

% >3D ?
ndim = size(x, 1);
if ndim > 3
    warning('plotmd:ndim', '#dimensions > 3, plotting only 3D component.')
end

% select 2D or 3D
if ndim == 2
    h = plot(ax, x(1, :), x(2, :), varargin{:} );
elseif ndim >= 3
    h = plot3(ax, x(1, :), x(2, :), x(3, :), varargin{:} );
end

if nargout == 1
    varargout{1, 1} = h;
end
