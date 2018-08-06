function ea_zoomcenter(varargin)
% Zoom in and out on a specifeid point on a 2-D or 3-D plot without losing
% data in the original plot.
%
% EA_ZOOMCENTER(CENTER) zooms the current axis on the point CENTER by a factor
% of 2.
%
% EA_ZOOMCENTER(CENTER, FACTOR) zooms the current axis on the point CENTER by
% FACTOR.
%
% EA_ZOOMCENTER(AX, ...) zooms the specified axis
%
% Example:
% line
% ea_zoomcenter([0.5, 0.5], 3)

if ishandle(varargin{1})    % axis specified
    ax = varargin{1};
    center = varargin{2};
    if nargin == 3
        factor = varargin{3};
    else    % no factor specified, use 2
        factor = 2;
    end
else    % no axis specified, use gca
	ax = gca;
    center = varargin{1};
    if nargin == 2
        factor = varargin{2};
    else    % no factor specified, use 2
        factor = 2;
    end
end


    xrange = center(1) + [-1, 1] * max(abs(ax.XLim - center(1)));
    yrange = center(2) + [-1, 1] * max(abs(ax.YLim - center(2)));
    zrange = center(3) + [-1, 1] * max(abs(ax.ZLim - center(3)));
    axis(ax, [xrange, yrange,zrange]);

if nargin>2
zoom(factor);
end
